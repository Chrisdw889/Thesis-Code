clear all; close all;
%% Load in data

[filename, path] = uigetfile("*.mat");
data = load(strcat(path, filename));

%% Combine all groups into one dataset

Data = combineGroups(data.importedfiles);


%% Filter Data

% Input Parameters
Data.Zdelta = Data.Z(2) - Data.Z(1);

Data.Filters.ThresholdSamples=10; % Number of consecutive points to search with for limits
Data.Filters.Noise.Samples=10; % Number of points in the trace's tails to sample noise
Data.Filters.Noise.StDevMultiplier=1; % Number of standard deviations used to characterise the noises variation
Data.Filters.Noise.Threshold= -3.0; % Noise values larger than this fail the filter

Data.Filters.logG.StartLimit = -0.5; % Traces with the mean of samples just before break lower than this fail the filter

Data.Filters.logG.HiThreshold= -0.5; % Traces with a mean logG/G0 value greater than this fail the filter
Data.Filters.logG.LowThreshold= -5.0; % Traces with a mean logG/G0 value less than this fail the filter

Data.Filters.Z.LoLimit = 0E-3; % Traces with a plateau length shorter than this fail the filter (in um)
Data.Filters.Z.HiLimit = 30E-3; % Traces with a plateau length longer than this fail the filter (in um)

Data.ZeroMark = -0.1; % Current when plateau starts, logG

if Data.ZeroMark > Data.Filters.logG.HiThreshold
    Data.ZeroMark = Data.Filters.logG.HiThreshold;
end
Data.PlateauMark = -5.5; % Current when plateau ends, logG

if Data.PlateauMark < Data.Filters.Noise.Threshold
    Data.PlateauMark = Data.Filters.Noise.Threshold;
end

% Filter the Data

Data = filterTraces(Data);
Data = calculatePlateau(Data);

goodTraces = Data.logG(:, Data.passedTraces == 1);

%% Save the Data

savetxt = input("Please enter the filename for saving the data: ");

cd(path);

save(savetxt, "Data");



%% Functions

function Data = combineGroups(importedfiles)

Data = struct();

numGroups = length(importedfiles);
bias = importedfiles{1}.header.bias;
setpoint = importedfiles{1}.header.setpoint;
BJDist = importedfiles{1}.header.BJDist;
speed = importedfiles{1}.header.speed;
pointsPerTrace = importedfiles{1}.pointsPerTrace;

% Assert the file properties match

for i=1:numGroups
    groupStruct = importedfiles{i};
    if ~isequaln(round(groupStruct.header.bias, 3), round(bias, 3))
        error("ImportedFile does not have matching bias!!");
    end
    if ~isequaln(groupStruct.header.setpoint, setpoint)
%         error("ImportedFile does not have matching setpoint!!");
    end
    if ~isequaln(groupStruct.header.BJDist, BJDist)
        error("ImportedFile does not have matching BJDist!!");
    end
    if ~isequaln(groupStruct.header.speed, speed)
        error("ImportedFile does not have matching speed!!");
    end
    if ~isequaln(groupStruct.pointsPerTrace,pointsPerTrace)
        error("ImportedFile does not have matching pointsPerTrace!!");
    end
end

% Create Header Data

header = struct();
header.setpoint = setpoint;
header.bias = bias;
header.BJDist = BJDist;
header.speed = speed;
header.pointsPerTrace = pointsPerTrace;

% Merge Data

rawI = [];
numTraces = 0;
logI = [];
logG = [];

Z = linspace(0, header.BJDist, header.pointsPerTrace);

for i=1:numGroups
    rawI = [rawI, importedfiles{i}.data];
    numTraces = numTraces + importedfiles{i}.numTraces;
    logI = [logI, importedfiles{i}.logI];
    logG = [logG, importedfiles{i}.logG];
end

% Construct and return merged data structure

Data.rawI = rawI;
Data.logI = logI;
Data.logG = logG;
Data.Z = Z;
Data.numTraces = numTraces;
Data.header = header;


end

function rData = filterTraces(Data)

% Check the mean value of each trace is below a threshold value
Data.passedTraces = zeros(Data.numTraces, 1);

for i=1:Data.numTraces
    traceLogG = Data.logG(:, i);
    if mean(traceLogG) < Data.Filters.logG.HiThreshold && ...
            mean(traceLogG) > Data.Filters.logG.LowThreshold % Traces passes check
        Data.passedTraces(i) = 1;
    else
        % value stays as false / 0;
    end
end
Data.NumGoodTraces = sum(Data.passedTraces == 1);
fprintf("%i out of %i passed high conductance filter with threshold %1.2f log(G/G0). \r\n", ...
    Data.NumGoodTraces, Data.numTraces, Data.Filters.logG.HiThreshold);



% Check the noise level of the traces tails

for i=1:length(Data.passedTraces)
    passed = Data.passedTraces(i);
    
    if passed
        traceLogG = Data.logG(:, i);
        
        noisetail = Data.Filters.Noise.StDevMultiplier * ...
            std(traceLogG(end - Data.Filters.Noise.Samples : end));
        noiselevel = mean(traceLogG(end - Data.Filters.Noise.Samples : end));
        
        if noiselevel < Data.Filters.Noise.Threshold % Pass noise test
            % value stays as true / 1
        else
            Data.passedTraces(i) = 0;
        end
    else
        continue;
    end
end

Data.NumGoodTraces = sum(Data.passedTraces == 1);
fprintf("%i out of %i passed both conductance and now the noise filter with threshold %1.2f log(G/G0). \r\n", ...
    Data.NumGoodTraces, Data.numTraces, Data.Filters.Noise.Threshold);


% Return filter results
rData = Data;


end

function rData = calculatePlateau(Data)

% Fix Z offset based on breaking point

Data.ZOffsetIndices = zeros(Data.numTraces, 1);
Data.ZOffsetLengths = zeros(Data.numTraces, 1);

lessThanAll = Data.logG < Data.ZeroMark;
for i=1:Data.numTraces
    if Data.passedTraces(i) == 0
        continue;
    end
    lessThan = lessThanAll(:, i);
    for j=1:Data.header.pointsPerTrace - Data.Filters.ThresholdSamples
        test = sum(lessThan(j:j + Data.Filters.ThresholdSamples));
        if test == Data.Filters.ThresholdSamples
            if j == 1 || j > (Data.header.pointsPerTrace - Data.Filters.ThresholdSamples)
                Data.passedTraces(i) = 0;
                break;
            else
                Data.ZOffsetIndices(i) = j;
                Data.ZOffsetLengths(i) = j * Data.Zdelta;
                break;
            end
        end
    end
end

Data.NumGoodTraces = sum(Data.passedTraces == 1);
fprintf("%i out of %i have passed through the Z offset of %1.2f log(G/G0) in addition. \r\n", ...
    Data.NumGoodTraces, Data.numTraces, Data.ZeroMark);

 % Calculate Plateau Length
 Data.PlateauLengths = zeros(Data.numTraces, 1);
 Data.PlateauIndices = zeros(Data.numTraces, 1);
 
lessThanAll = Data.logG < Data.PlateauMark;
for i=1:Data.numTraces
    if Data.passedTraces(i) == 0
        continue;
    end
    lessThan = lessThanAll(:, i);
    for j=1:Data.header.pointsPerTrace - Data.Filters.ThresholdSamples
        test = sum(lessThan(j : j + Data.Filters.ThresholdSamples));
        
        if test == Data.Filters.ThresholdSamples
            if j == 1 || j == Data.header.pointsPerTrace - Data.Filters.ThresholdSamples
                Data.passedTraces(i) = 0;
                break;
            else
                Data.PlateauIndices(i) = j;
                Data.PlateauLengths(i) = (j - Data.ZOffsetIndices(i)) * Data.Zdelta;
                if Data.PlateauLengths(i) < Data.Filters.Z.LoLimit || ...
                        Data.PlateauLengths(i) > Data.Filters.Z.HiLimit
                    Data.passedTraces(i) = 0;
                    break;
                else
                    break;
                end
            end
        end
    end
    
end

Data.NumGoodTraces = sum(Data.passedTraces == 1);
fprintf("%i out of %i have passed plateau filetering with %1.2f log(G/G0) break distance and with %f and %f HiLo Length limits. \r\n", ...
    Data.NumGoodTraces, Data.numTraces, Data.PlateauMark, Data.Filters.Z.HiLimit, Data.Filters.Z.LoLimit);

% Check the first few points start at the right current

for i=1:length(Data.passedTraces)
    passed = Data.passedTraces(i);
    if passed
        traceLogG = Data.logG(:, i);
        traceZOffsetIndex = Data.ZOffsetIndices(i);
        if traceZOffsetIndex < 11
            start = traceLogG(1 : Data.Filters.Noise.Samples);
        else
            start = traceLogG(traceZOffsetIndex - Data.Filters.Noise.Samples: traceZOffsetIndex);
        end
        startlevel = mean(start);
        
        if startlevel > Data.Filters.logG.StartLimit
            % values stays as true / 1
        else
            Data.passedTraces(i) = 0;
        end
    else
        continue;
    end
end

Data.NumGoodTraces = sum(Data.passedTraces == 1);
fprintf("%i out of %i passed high conductance filter and start level filter with threshold %1.2f log(G/G0). \r\n", ...
    Data.NumGoodTraces, Data.numTraces, Data.Filters.logG.StartLimit);

rData = Data;

end



