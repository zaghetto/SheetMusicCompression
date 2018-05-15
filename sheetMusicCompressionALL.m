%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alexandre Zaghetto                               %
% zaghetto@unb.br                                  %
% University of Brasília                           %
% Department of Computer Science                   %
% Laboratory of Images, Signals and Acoustics      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Encode all documents in dataset.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare environment
clear all;
close all;
clc;

% Select music folders
PATHNAME = uigetdir([], 'Select image folder');
docs = dir(PATHNAME);

% Select music subfolder
tic
for i = 3:length(docs)    
    disp(docs(i).name)    
    PATHMUSIC = [PATHNAME '\' docs(i).name];
    sheetMusicCALL(PATHMUSIC)   
end
toc


% Results folder
resdir = 'FINALRESULT';

% Current folder
currDir = pwd;

% Create results 
if ~(exist(resdir, 'dir') == 7)
    mkdir(resdir);
else
    delete([currDir '\' resdir '\*.mat'])
end

% Copy configuration files to results folder
copyfile([currDir '\*.cfg'],[currDir '\' resdir '\'])

% Move most recent file
dirMAT = dir('MAT')
for i = 3:length(dirMAT)        
    dirRes = dir([currDir '\MAT\' dirMAT(i).name '\*.mat']);
    if length(dirRes) == 5
        movefile([currDir '\MAT\' dirMAT(i).name '\' dirRes(5).name], [currDir '\' resdir '\' dirRes(5).name]);
    else
        break
    end
end

rmdir([currDir '\MAT'], 's')

cd FINALRESULT;

% Plot results
plot_result;

cd ..







