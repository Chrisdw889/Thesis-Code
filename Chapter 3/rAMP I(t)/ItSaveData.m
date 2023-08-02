%%% Reads in either binary or ascii files produced by the STM
%%% Can import either single or multiple files
%%% Restructures imported data based on chunk sizes
%%% Recommended to only save data with the same operating params (setpoint)

clear all; close all;

directory = uigetdir("Z:\CDW\PhD\", "Please select the directory:");
directory = strcat(directory, "\");


file_text = input("Please enter the pattern for filenames:");

files = dir(strcat(directory, file_text));
numfiles = length(files);

importedfiles = cell(numfiles, 1);
for i=1:numfiles
    filename = files(i).name;
    importedfiles{i} = Import.import(directory, filename, 'I(t) Spec');
end


% Structure Data
for f = 1:length(importedfiles)
    data = importedfiles{f}.data;
    chunks = importedfiles{f}.header.chunks;
    
    i = 1;
    j = 1;
    tracenum = 1;
    while i <= length(chunks)
        m_size = chunks{i};
        newdata(1:m_size, tracenum) = data(j:j+m_size-1);
        i = i + 1;
        j = j + m_size - 1;
        tracenum = tracenum + 1;
    end
    importedfiles{f}.structdata = newdata;
end

% Filter Data

for i = 1:length(importedfiles)
    strucdata = importedfiles{i}.structdata;
    [rows, cols] = size(strucdata);
    
    filtereddata = zeros(rows, cols);
    
    for j = 1:cols
        trace = strucdata(:, j);
        filteredtrace = besselfilter(4, 2500, 25000, trace);
        filtereddata(:, j) = filteredtrace;
    end
    
    importedfiles{i}.filtereddata = filtereddata;
end

savetxt = input("Please enter the filename for saving the data: ");

cd(directory);

save(savetxt, "importedfiles");


%% Functions

function [filteredtrace, b, a] = besselfilter(order, high, sampling, data)
    [b, a] = besself(4, 2 * pi * high);
    [num, den] = bilinear(b, a, sampling);
    y = filter(num, den, data);
    filteredtrace = y;
end
