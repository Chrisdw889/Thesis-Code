clear all; close all;

directory = uigetdir("Y:\CDW\Thesis Work by Chapter\Break Junctions\");
directory = strcat(directory, "\");

file_text = input("Please enter the pattern for filenames:");

BJDist = input("Please enter the Distance used in I(s) Spectroscopy in um:");
BJSpeed = input("Please enter the retraction speed in I(s) Spectroscopy in um/s:");

files = dir(strcat(directory, file_text));
numfiles = length(files);

importedfiles = cell(numfiles, 1);
for i=1:numfiles
    filename = files(i).name;
    importedfiles{i} = Import.import(directory, filename, 'I(s)');
end

%% Convert to real current and add in distance and Conductance 
load("CalibCurve.mat");

G0 = 7.74809173e-5; % in S

for i=1:numfiles
    for j=1:importedfiles{i}.numTraces
        
        realLogCurr = interp1(CalibCurve.VoltageOut,...
            CalibCurve.CurrentIn, importedfiles{i}.data(:, j));
        for k=1:length(realLogCurr)
            val = realLogCurr(k);
            if isnan(val)
                realLogCurr(k) = -12;
            end
        end
        
        importedfiles{i}.logI(:, j) = realLogCurr;
        
        importedfiles{i}.logG(:, j) = importedfiles{i}.logI(:, j) - log10(abs(importedfiles{i}.header.bias))...
            - log10(G0);
        
        importedfiles{i}.header.BJDist = BJDist;
        importedfiles{i}.header.speed = BJSpeed;
        importedfiles{i}.Z = linspace(0, importedfiles{i}.header.BJDist, 2000); % Distance in um
    end
end

%% Save data

savetxt = input("Please enter the filename for saving the data: ");

cd(directory);

save(savetxt, "importedfiles");


