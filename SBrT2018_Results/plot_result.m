%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alexandre Zaghetto                                          %
% zaghetto@unb.br                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intitute of Exact Sciences                                  %
% Department of Computer Science                              %
% University of Brasília                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version: 2018-04-24                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Plot results of document encoding.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare environment
clear all
close all

files = dir('*.mat');

% Interpolation interval
interi = [0.25:0.01:1];

pHEVCSum = zeros(1, length(interi));
pHEVCISum = zeros(1, length(interi));
pAVCSum = zeros(1, length(interi));
pAVCISum = zeros(1, length(interi));
pJP2Sum = zeros(1, length(interi));

for z = 1:9
    
    eval(['load ' files(z).name]);
    
    files(z).name
    
    pHEVCSum = pHEVCSum+pHEVC;
    pHEVCISum = pHEVCISum+pHEVCI;
    pAVCSum = pAVCSum+ pAVC;
    pAVCISum = pAVCISum+pAVCI;
    pJP2Sum = pJP2Sum+pJP2;
    
end

% Calculate average
pHEVCSum = pHEVCSum/6;
pHEVCISum = pHEVCISum/6;
pAVCSum = pAVCSum/6;
pAVCISum = pAVCISum/6;
pJP2Sum = pJP2Sum/6;

%%%%%%%%%%%%%%%%%
% Global plot   %
%%%%%%%%%%%%%%%%%
inte = [1:ceil(length(interi)/8):length(interi)-2*ceil(length(interi)/8)];
figure
plot(rHEVC(inte), pHEVCSum(inte),'r.-','Linewidth', 1,'MarkerSize',10)
set(gca,'fontsize', 12);
hold on
grid on
plot(rHEVC(inte), pHEVCISum(inte),'ks-','Linewidth', 1,'MarkerSize',8)
plot(rHEVC(inte), pAVCSum(inte),'-*','Linewidth', 1,'MarkerSize',8)
plot(rHEVC(inte), pAVCISum(inte),'-o','Linewidth', 1,'MarkerSize',8)
plot(rHEVC(inte), pJP2Sum(inte),'k->','Linewidth', 1,'MarkerSize',8)
legend('HEVC-INTER','HEVC-INTRA','AVC-INTER','AVC-INTRA','JPEG2000','Location','Best')
grid on
ylabel('Average PSNR (dB)')
xlabel('bitrate (bits/pixel)')
title('Comparison Between Encoders')
axis([0.2 0.8 34 64])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEVC-INTER vs. HEVC-INTRA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inte = 1:18:length(rHEVC)-18;

% Titles
tituloPSNR = 'Average PSNR difference';
tituloTaxa = 'Average bitrate difference';

% Curve 1
label1 = 'HEVC-INTER';  % Label for curve 1

% Curve 2
label2 = 'JPEG2000'; % Label for curve 2

% Ajusta label para eixo x e para eixo y
labelx = 'bitrate (bits/pixel)';
labely = 'Average PSNR (dB)';

% Method
method = 0
%   method = 0: Average PSNR difference;
%   method = 1; Average bitrate difference;

% Chama rdcurvediff (calcula a diferença média PSNR ou bitrate)
figure
pause(0.1)
[Coefs, diff] = rdcurvediff(method, rHEVC(inte), pHEVCSum(inte), rHEVC(inte), pHEVCISum(inte),label1,label2, labelx, labely, tituloPSNR, tituloTaxa);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEVC-INTER vs. JPEG2000 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
label2 = 'JPEG2000'; % Label for curve 2
pause(0.1)
[Coefs, diff] = rdcurvediff(method, rHEVC(inte), pHEVCSum(inte), rHEVC(inte), pJP2Sum(inte),label1,label2, labelx, labely, tituloPSNR, tituloTaxa);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEVC-INTER vs. JPEG2000 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
TableRes = zeros(5,4);
TableRes(1,:) = rHEVC(inte);
TableRes(2,:) = pJP2Sum(inte);
TableRes(3,:) = pAVCISum(inte);
TableRes(4,:) = pAVCSum(inte);
TableRes(5,:) = pHEVCISum(inte);
TableRes(6,:) = pHEVCSum(inte);

TableRes













