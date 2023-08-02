%%% Version - v1.00
%%% Author  - Christopher Weaver

%%% Summary:
%%% This script will extract and save all the relevent data out of a
%%% directory containing files that start with the "I(s)" header. Once
%%% extracted all data will be saved in a matlab structure .mat file.

%%% Instructions:
%%% Run and select options when prompted.

clear all; close all;
directory = uigetdir("Y:\CDW\Thesis Work by Chapter\I(t) Classification\Distance Calib", "Select the folder of study");
home = "Y:\CDW\Thesis Work by Chapter\I(t) Classification\Distance Calib\MATLAB Scripts";

study = IsGroup(directory, 0);
study.plotAllCurrentDistance(1);
betas = study.betaAnalysis();

struct.currents = study.current;
struct.time = study.time;
struct.distance = study.distance;
struct.directory = study.directory;
struct.traceNum = study.traceNum;
struct.pointsPerTrace = study.tracesPoints;
struct.betas = betas;

%%
cd(directory);

filename = input('Enter the name of the file to be saved:\r\n');

save(filename, 'struct');


