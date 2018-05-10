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
for i = 3:length(docs)
    disp(docs(i).name)
    tic
    PATHMUSIC = [PATHNAME '\' docs(i).name];
    sheetMusicCALL(PATHMUSIC)
    toc
end





