function sheetMusicCALL(PATHNAME)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alexandre Zaghetto                               %
% zaghetto@unb.br                                  %
% University of Brasília                           %
% Department of Computer Science                   %
% Laboratory of Images, Signals and Acoustics      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Encode one document using JPEG2000,              %
% AVC, AVC-INTRA, HEVC and HEVC-INTRA              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of pages
RESDIR = dir([PATHNAME '\*.bmp']);
N = length(RESDIR);

% Remove file extension
FILENAME = RESDIR(1).name(1:length(RESDIR(1).name)-5);

% Get current folder
currentFolder = pwd;

% Get file size
i = 0;
F = imread([PATHNAME '\' FILENAME num2str(i) '.bmp']);
[h, w] = size(F);

% Get total number of bytes
OrigSize = N*RESDIR(1).bytes;
 
% Set econders on or off
jpeg2000  = 0;

h264intra = 0;

h264video = 1;

hevcintra = 0;

hevcvideo = 0;

% Turn h264video on if hevcvideo if on. Both will use the same YUV file.
if hevcvideo
   h264video = 1;
end

% Turn h264intra on if hevcintra if on. Both will use the same YUV file.
if hevcintra
    h264intra = 1;
end
    
% PSNR calculation mode
PSNR_MODE = 0;  % 0 - Global, 1 - Average

% Encoding parameters
QPI = 10;
QPF = 50;
QPStep = 10;
SearchRange = 64;

% Length of quantization paremeters vector
L = length(QPI:QPStep:QPF);

% Legend
leg = '';
q = char(39);

% Interpolation interval
interi = [0.25:0.01:1];

% Checkpoint nothing encoded
chkpt = 'NULL';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H.264 - x264  (AVC)        % Video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if h264video
    
    % Number of vertical (h) and horizontal (w) subpages: 0, 2 or 4.
    %prompt = {'Number of vertical (h) and horizontal (w) subpages'};
    %name = 'Number of subpages';
    %numlines = 1;
    %defaultanswer={'1'};
    
    %nSPh = str2double(cell2mat(inputdlg(prompt, name, numlines, defaultanswer)));
    nSPh = 2;
    nSPw = nSPh;
    
    % Total number of frames
    Np = N*nSPh*nSPw;
    
    % Generate qpfile for x264
    qpf = 'qpfile.txt';
    x264QPfile(Np, qpf)
             
    % Size of each frame
    H = h/nSPh;
    W = w/nSPw;
    
    % Initialize Y
    Y = zeros(H, W, Np);
    
    % YUV 4:2:0 format fos x264. U and V are constant
    U = 128*ones(H/2, W/2, Np);
    V = 128*ones(H/2, W/2, Np);
    
    k = 1;
    for frame = 1:N
        
        % Load each frame
        F = imread([PATHNAME '\' FILENAME num2str(frame-1) '.bmp']);
        
        % Frame limits
        seths = h/nSPh: h/nSPh :h;
        sethi = 1: h/nSPh: h;
        setws = w/nSPw: w/nSPw :w;
        setwi = 1: w/nSPw: w;
        
        for i = 0 : nSPh-1
            for j = 0 : nSPw-1
                Y(:,:,k) = F(sethi(i+1):seths(i+1),setwi(j+1):setws(j+1) );
                k = k + 1;
            end
        end
        
    end
      
    % Save YUV
    delete('*.yuv');         
    writeYUV400vid([FILENAME '_video'], Y, U, V, Np);        
   
    % Initialize PSNR and Bitrate vectors
    PSNRY_AVCTotal = zeros(1,L);
    Rate_AVCTotal  = zeros(1,L);
    
    % For different quantization parameters
    cont = 1;            
    
    for QPISlice = QPI:QPStep:QPF % 0 (Best quality) a 51 (Worst quality)

        % Encode
        [status, result] = system(['x264 --qp ' num2str(QPISlice) ' --qpfile ' qpf ' --merange ' num2str(SearchRange) ' --tune psnr --psnr --verbose --profile high --preset placebo -o ' ...
            FILENAME '_video.264 ' FILENAME '_video.yuv --input-res ' num2str(W) 'x' num2str(H) ' -r 4']);
                        
        % H.264 file size
        Info_Img_h264 = dir([FILENAME '_video.264']);
        Size_h264 = Info_Img_h264.bytes;
        
        % Calculate bitrate
        BpP = 8*Size_h264/(N*h*w);
        
        % Decode                
        [status,result_decod] = system(['ffmpeg -i ' FILENAME '_video.264 -c:v rawvideo -pix_fmt yuv420p '...
            FILENAME '_video_decod_x264.yuv']);
        
        % Read decoded yuv
        [Yrec, Urec, Vrec] = readyuv([FILENAME '_video_decod_x264.yuv'], W, H, Np, 0);
        
        % Reconstruct video
        k = 1;        
        for frame = 1:N
                                   
            % Reconstruct video
            seths = h/nSPh: h/nSPh :h;
            sethi = 1: h/nSPh: h;
            setws = w/nSPw: w/nSPw :w;
            setwi = 1: w/nSPw: w;
            
            for i = 0 : nSPh-1
                for j = 0 : nSPw-1
                    F(sethi(i+1):seths(i+1),setwi(j+1):setws(j+1), frame ) = Yrec(:,:,k);                    
                    k = k + 1;
                end
            end
            
        end
                         
        % Initialize decoded frames matrix
        Frec = zeros(h, w, N);
        
        % Calculate PSNR
        if PSNR_MODE == 0          
            PSNR = calcPSNR(Y(:), Yrec(:));
        else
            PSNR = 0;
            for frame = 1:N
                FOrig = imread([FILENAME num2str(frame-1) '.bmp']);
                FRec = F(:, :, frame);
                PSNR = PSNR + calcPSNR(FOrig(:), FRec(:));
            end                        
                PSNR = PSNR/N;
        end
        
        % Save bitreate and PSRN
        PSNRY_AVCTotal(cont) = PSNR
        Rate_AVCTotal(cont) = BpP      
        
        % Delete auxiliary files
        delete('*_video.264');
        delete('*_video_decod_x264.yuv');
        
        % Next quantization parameter
        cont = cont+1;
        
    end            
    
    rAVC = interi;
    pAVC = interp1(Rate_AVCTotal,PSNRY_AVCTotal,rAVC,'spline');
%     plot(Rate_AVCTotal,PSNRY_AVCTotal,'o',rAVC,pAVC,'--');
%     xlim([0 1])

    % Delete QP file
    delete(qpf);
            
else
    Rate_AVCTotal = [0.3000    0.7000    1.1000    1.5000];
end


% Checkpoint AVC-INTER
chkpt = 'AVC';

% Save partial results
eval(saveResultPartial(FILENAME, 'MAT', FILENAME(1:end-2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H.264 - JM (AVC-INTRA)     % Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if h264intra    
    
    % Load frames
    Page = zeros(h, w, N);
    for i = 1:N
        Page(:,:,i) = imread([PATHNAME '\' FILENAME num2str(i-1) '.bmp']);
    end
    
    % Generate qpfile
    qpfile = '  0 I';
    for frNum = 1:N-1
        if frNum < 10
          qpfile = [qpfile; '  ' num2str(frNum) ' I'];        
        elseif frNum < 100
          qpfile = [qpfile; ' ' num2str(frNum) ' I'];        
        else
          qpfile = [qpfile; num2str(frNum) ' I'];
        end
    end
    
    % Save qpfile    
    qpf = 'qpfile.txt';
    delete(qpf);
    fid = fopen('qpfile.txt','wt');
    for i = 1:N
        fprintf(fid,'%d %s\n', uint8(str2double(qpfile(i,1:3))), qpfile(i,4:end));            
    end
    fclose(fid);
                
    % Initialize U and V with 128
    U = 128*ones(h/2, w/2, N);
    V = 128*ones(h/2, w/2, N);
    
    % Save YUV
    writeYUV400vid([FILENAME '_image'], Page, U, V, N)
            
    % Encoding parameters    
%     QPI = QPI+5;
%     QPF = QPF+5;
        
    % Length of quantization paremeters vector
    L = length(QPI:QPStep:QPF);
    
    % Initialize PSNR and Bitrate vectors
    PSNRY_AVCITotal = zeros(1,L);
    Rate_AVCITotal  = zeros(1,L);
           
    % For different quantization parameters
    cont = 1;    
    for QPISlice = QPI:QPStep:QPF 
        
        QPISlice
                                          
        % Encode
        [status, result] = system(['x264 --qp ' num2str(QPISlice) ' --qpfile ' qpf ' --tune psnr --psnr --verbose --profile high --bframes 0 --preset placebo -o ' ...
            FILENAME '_image.264 ' FILENAME '_image.yuv --input-res ' num2str(w) 'x' num2str(h) '--keyint 1']);                        
        
        % Size in bytes
        Info_Img_h264 = dir([ FILENAME '_image.264']);
        Size_h264 = Info_Img_h264.bytes;
        
        % Calculate bitrate
        BpP = 8*Size_h264/(N*h*w);
        
        % Decode
        [status,result_decod] = system(['ffmpeg -i ' FILENAME '_image.264 -c:v rawvideo -pix_fmt yuv420p '...
            FILENAME '_image_decod_x264.yuv']);
        
        % Read decoded file
        [PageRecHEVCI, Urec, Vrec] = readyuv([FILENAME '_image_decod_x264.yuv'], w, h, N, 0);                             
        
        % Calculate PSNR
        if PSNR_MODE == 0
            PSNR = calcPSNR(Page(:), PageRecHEVCI(:));
        else
            PSNR = 0;
            for frame = 1:N
                FOrig = Page(:,:,frame);
                FRec = PageRecHEVCI(:,:, frame);
                PSNR = PSNR + calcPSNR(FOrig(:), FRec(:));
            end
            PSNR = PSNR/N;
        end   
                  
        % Save PSNR and bitrate
        PSNRY_AVCITotal(cont) = PSNR
        Rate_AVCITotal(cont) = BpP
        
        % Delete temporary files
        delete('*_image_decod_x264.yuv');
        delete('*.264');
                
        % Next quantization parameter
        cont = cont+1;        
    end   
    
    rAVCI = interi;
    pAVCI = interp1(Rate_AVCITotal,PSNRY_AVCITotal,rAVCI,'spline');  
    
end

% Checkpoint
chkpt = 'AVC-I';

% Save partial results
eval(saveResultPartial(FILENAME, 'MAT', FILENAME(1:end-2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kakadu (JPEG2000)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if jpeg2000
            
    % Load frames
    Page = zeros(h,w,N);
    for i = 1:N
        Page(:,:,i) = imread([PATHNAME '\' FILENAME num2str(i-1) '.bmp']);
    end

    % Bitrates
    Taxas = Rate_AVCTotal;
    
    % Initialize decoded frames
    PageRecJPG2000 = zeros(h,w,N);    
    
    % Initializ rates and PSNRs
    PSNRY_JP2Total = zeros(1, length(Taxas));
    Rate_JP2Total = zeros(1,length(Taxas));
    
    % For each rate
    cont  = 1;
    for TaxaKadu = Taxas
        
        Size_jp2 = 0;
        
        for i = 1:N           
                        
            % Save image as temporary pgm
            imwrite(uint8(Page(:,:,i)), 'temp.pgm' ,'pgm');
                                              
            % Encode
            [status,result] = system(['kdu_compress'...
                ' -i temp.pgm'...
                ' -o temp.j2c'...
                ' Cuse_sop=yes'...
                ' Cuse_eph=yes'...
                ' Creversible=no'...
                ' Cmodes=RESTART'...
                ' Qderived=yes'...
                ' -rate ' num2str(TaxaKadu)]);
                        
            % File size in bytes
            Info_Img_jp2 = dir(['temp.j2c']);            
            Size_jp2 = Size_jp2 + Info_Img_jp2.bytes; % Em bytes
            
            % Decode
            [status,result] = system(['kdu_expand'...
                ' -i temp.j2c'...
                ' -o temp_kakadu.pgm']);
            
            % Read decoded file
            PageRecJPG2000(:,:,i) = imread('temp_kakadu.pgm');
            
            % Delete temporary files
            delete('*.pgm');                                
            delete('*.j2c');
            
        end                
        
        % Calculate PSNR
        if PSNR_MODE == 0
            PSNR = calcPSNR(Page(:), PageRecJPG2000(:));
        else
            PSNR = 0;
            for frame = 1:N
                FOrig = Page(:,:,frame);
                FRec = PageRecJPG2000(:,:, frame);
                PSNR = PSNR + calcPSNR(FOrig(:), FRec(:));
            end
            PSNR = PSNR/N;
        end
            
        % Save bitrate and PSNR
        PSNRY_JP2Total(cont) = PSNR
        Rate_JP2Total(cont) = 8*Size_jp2/(N*h*w)
        
        % Next rate
        cont = cont+1;
                
    end
    
    rJP2 = interi;
    pJP2 = interp1(Rate_JP2Total,PSNRY_JP2Total,rJP2,'spline');
        
end

% Checkpoint JPEG2000
chkpt = 'JPEG2000';

% Save partial results
eval(saveResultPartial(FILENAME, 'MAT', FILENAME(1:end-2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H.265 (HEVC - HM)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hevcvideo
    
    % HEVC configuration files
    ARQCONFIGH265 = 'encoder_lowdelay_P_main.cfg';
       
    % Initialize PSNR and bitreate vectors
    PSNRY_TotalHEVC = zeros(1,L); % dB
    Rate_TotalHEVC  = zeros(1,L); % bpp
    
    % QPs
%     QPI = QPI-5;
%     QPF = 42-5;
        
    % Length of QP vector
    L = length(QPI:QPStep:QPF);
    
    % For each QP
    cont = 1;    
    for QPISlice = QPI:QPStep:QPF % de 0 (melhor qualidade) a 51 (pior qualidade)
        
        QPISlice                                      
                        
        % Encode
        [status, result] = system(['TAppEncoder -c ' ARQCONFIGH265...
            ' -i ' FILENAME '_video.yuv'...
            ' -o ' FILENAME '_recon_hevc.yuv'...
            ' -b ' FILENAME '_video.265'...
            ' -wdt ' num2str(W)...
            ' -hgt ' num2str(H)...
            ' -fr 1'...
            ' -f ' num2str(Np)...
            ' -sr ' num2str(SearchRange)...
            ' --SEIDecodedPictureHash'...
            ' --Level'...
            ' -q ' num2str(QPISlice)])                   
        
        % Size of H.265 file
        Info_Img_h265 = dir([FILENAME '_video.265']);
        Size_h265 = Info_Img_h265.bytes;
        
        % Calculate bitrate
        BpP = 8*Size_h265/(N*h*w);
                
        % Decode
        [status,result_decod] = system(['TAppDecoder -b ' FILENAME '_video.265 -d 8'...
            ' -o ' FILENAME '_video_decod_hevc.yuv']);
        
        % Read decoded file
        [Yrec, Urec, Vrec] = readyuv([FILENAME '_video_decod_hevc.yuv'], W, H, Np, 0);
        
        k = 1;
        for vista = 1:N
            
            % Constroi o video
            seths = h/nSPh: h/nSPh :h;
            sethi = 1: h/nSPh: h;
            setws = w/nSPw: w/nSPw :w;
            setwi = 1: w/nSPw: w;
            
            for i = 0 : nSPh-1
                for j = 0 : nSPw-1
                    F(sethi(i+1):seths(i+1),setwi(j+1):setws(j+1), vista ) = Yrec(:,:,k);
                    % imshow(uint8(Y(:,:,k)))
                    % pause
                    k = k + 1;
                end
            end
            
        end
        
        if PSNR_MODE == 0
            PSNR = calcPSNR(Y(:), Yrec(:));
        else
            PSNR = 0;
            for vista = 1:N
                FOrig = imread([FILENAME num2str(vista-1) '.bmp']);
                FRec = F(:,:, vista);
                PSNR = PSNR + calcPSNR(FOrig(:), FRec(:));
            end
            PSNR = PSNR/N;
        end
                
        % Save PSNR and bitrate
        PSNRY_TotalHEVC(cont) = PSNR
        Rate_TotalHEVC(cont) = BpP
        
         % Delete temporary files
        delete('*.265');
        delete([FILENAME '*_video_decod_hevc.yuv']);
        
        % Next QP
        cont = cont+1;
        
    end
    
    delete('*_video.yuv');
            
    rHEVC = interi;
    pHEVC = interp1(Rate_TotalHEVC, PSNRY_TotalHEVC ,rHEVC,'spline');
        
end

% Checkpoint HEVC-INTER
chkpt = 'HEVC';

% Save partial results
eval(saveResultPartial(FILENAME, 'MAT', FILENAME(1:end-2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H.265 - HM (HEVC-INTRA)    % Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if hevcintra
    
    % HEVC configuration files
    ARQCONFIGH265INTRA = 'encoder_intra_main.cfg';
    
    % Save decoded pages
    PageRecHEVCI = zeros(h,w,N);
         
    % Length of QP vector
    L = length(QPI:QPStep:QPF);
    
    % Initialize PSNR and bitreate vectors
    PSNRY_HEVCITotal = zeros(1,L);
    Rate_HEVCITotal = zeros(1, L);
    
    % For each QP
    cont = 1;    
    for QPISlice = QPI:QPStep:QPF % de QPI (melhor qualidade) a QPF (pior qualidade)
        
        QPISlice
        
        % Chama o JM via linha de comando
        [status, result] = system(['TAppEncoder -c ' ARQCONFIGH265INTRA...
            ' -i ' FILENAME '_image.yuv' ...
            ' -o ' FILENAME '_recon_hevc.yuv' ...
            ' -b ' FILENAME '_image.265'...
            ' -wdt ' num2str(w) ...
            ' -hgt ' num2str(h) ...
            ' -fr 1'...
            ' -f ' num2str(N) ...
            ' --SEIDecodedPictureHash'...
            ' --Level'...
            ' -q ' num2str(QPISlice)]);
        
        % Decodifica
        [status,result_decod] = system(['TAppDecoder -b ' FILENAME '_image.265 -d 8'...
            ' -o ' FILENAME '_image_decod_hevc.yuv']);
        
        % H.265 file size
        Info_Img_h265 = dir([ FILENAME '_image.265']);
        Size_h265 = Info_Img_h265.bytes; 
        BpP = 8*Size_h265/(N*h*w);
        
        % Read decoded file
        [PageRecHEVCI, Urec, Vrec] = readyuv([FILENAME '_image_decod_hevc.yuv'], w, h, N, 0);
                
        % Calculate PSNR
        if PSNR_MODE == 0
            PSNR = calcPSNR(Page(:), PageRecHEVCI(:));
        else
            PSNR = 0;
            for vista = 1:N
                FOrig = Page(:,:,vista);
                FRec = PageRecHEVCI(:,:, vista);
                PSNR = PSNR + calcPSNR(FOrig(:), FRec(:));
            end
            PSNR = PSNR/N;
        end
        
        % Save PSNR and bitrate
        PSNRY_HEVCITotal(cont) = PSNR
        Rate_HEVCITotal(cont) = BpP
        
        % Delete temporary files
        delete('*.265');
        delete([FILENAME '*_image_decod_hevc.yuv']);
        
        % Next QP
        cont = cont+1;
        
    end
    
    delete('*_image.yuv');
    
    rHEVCI = interi;
    pHEVCI = interp1(Rate_HEVCITotal, PSNRY_HEVCITotal,rHEVCI,'spline');
    
end

% Checkpoint HEVC-INTRA
chkpt = 'HEVC-I';

% Save partial results
eval(saveResultPartial(FILENAME, 'MAT', FILENAME(1:end-2)));

toc

return






