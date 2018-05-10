function [Coefs, diff] = rdcurvediff(method, bitrates1, PSNRs1, bitrates2, PSNRs2, label1, label2, labelx, labely, tituloPSNR, tituloTaxa)
% RDCURVEDIFF Realiza:
%
%   method = 0: Average PSNR difference;
%   method = 1; Average bitrate difference;
%
% 	According to:
%
%   [1] G. Bjontegaard,  “Calculation of average PSNR differences between RD-curves”, 
%   presented at the 13th VCEG-M33 Meeting, Austin, TX, Apr. 2001.

if nargin ~=11 
    error('Número de argumentos deve ser igual a 11.')
elseif (length(bitrates1) ~= 4) | (length(PSNRs1) ~= 4) | (length(bitrates2) ~= 4) | (length(PSNRs2) ~= 4)
    error('Os vetores bitrates e PSNRs devem conter exatamente 4 valores.')
elseif (method <0) | (method>1)
    error('Método deve ser 0: diferença PSNR em dB ou 1: diferença de bitrates em %')
else

% Monda a matriz de dados    
Dados = [PSNRs1; bitrates1; PSNRs2; bitrates2];

% Salva os dados no formato que o avsnr precisa
save snr.dat method -ASCII
save snr.dat Dados -ASCII -APPEND

% Chama o avsnr
[status,result_cod] = system('avsnr snr.dat');

% Carrega os resultados do avsnr
Coefs = double(load('log.dat'));
diff = double(load('diff.dat'));

% Faz a interpolação e imprime a curva 1
ini = min(bitrates1);
fim = max(bitrates1);
resolucao = (fim-ini)/100;
xx = ini:resolucao:fim;
yy1 = pchip(bitrates1,PSNRs1,xx);
plot(xx,yy1,'r-')
hold on

% Faz a interpolação e imprime a curva 2
ini = min(bitrates2);
fim = max(bitrates2);
resolucao = (fim-ini)/100;
xx = ini:resolucao:fim;
yy2 = pchip(bitrates2,PSNRs2,xx);
plot(xx,yy2,'k-')
grid on
%legend(label1, label2,' ',2)

Font_Size = 12;
set(gca,'fontsize', Font_Size);
% Aplica os labels
 xlabel(labelx);
 ylabel(labely);
 if method==0
      title([tituloPSNR ': ' num2str(diff) ' (dB)']);
 else
     title([tituloTaxa ': ' num2str(diff) ' (%)']);
 end
 
% Imprime os pontos sem a interpolação
plot(bitrates1,PSNRs1,'r.','MarkerSize',10);
plot(bitrates2,PSNRs2,'ks','MarkerSize',8);
p1 = plot(bitrates1(1),PSNRs1(1),'-r.','MarkerSize',10);
p2 = plot(bitrates2(1),PSNRs2(1),'-ks','MarkerSize',8);
legend([p1 p2], label1, label2,'Location', 'Best');
XMIN = min([bitrates1 bitrates2]) - (20/100)*min([bitrates1 bitrates2]);
XMAX = max([bitrates1 bitrates2]) + (20/100)*min([bitrates1 bitrates2]);
YMIN = min([PSNRs1 PSNRs2]);
YMAX = max([PSNRs1 PSNRs2]);
axis([XMIN XMAX YMIN YMAX+0.1])

delete *.dat;

end

return