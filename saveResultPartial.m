function resultsave = saveResultPartial(FILENAME, MAT, DOC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alexandre Zaghetto                                %
% zaghetto@unb.br                                   %
% University of Brasília                            %
% Department of Computer Science                    %
% Laboratory of Images, Signals and Acoustics       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Save workspace with partial results. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = char(datetime('now','Format','y_M_d_HH_mm_ss'));

if ~(exist(MAT, 'dir') == 7)
    mkdir(MAT)
end

if ~(exist([MAT '/' DOC '/'], 'dir') == 7)    
    mkdir([MAT '/' DOC '/']);
end

resultsave = ['save ' MAT '/' DOC '/resultado_' FILENAME(1:end-2) '_' d  '.mat'];

end