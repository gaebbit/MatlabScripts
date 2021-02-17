% Original script from Andreys computer

% Sample: 0.5uL of 1:1 (AngI 100uM in 0.1% TFA):(sDHB20 in TA30) on Stainless Steel, air dried
% where TA30 is 3:7 ACN:0.1% TFA in water
% Sample (UV 2-A): 16kV
% Extraction (UV 2-B): 11kV
function [tall,yall] = FCT_read_AndreysData(filename,logfilename,vari)


signalBounds= [];

%% data readout
% read logfile
hLogFile    = fopen(logfilename, 'r');
logdate     = fscanf(hLogFile, 'Date: %s %s\n'); % Zeile wird irgendwie gebraucht
logtriggersperframe = fscanf(hLogFile, 'Triggers per step: %u\n');
fclose(hLogFile);

% read signal file
hFile       = fopen(filename, 'r');
header      = AcqReadHeader(hFile);
t           = AcqGetTime(header);

figure
N           = floor(header.nFramesCount/logtriggersperframe);
yall        = nan*ones(N,40000); tall = yall;
% split the shots into groups to save memory.
if N>2000
    Nind    = vari;
else
    Nind    = 1:N;
end
for ii      = Nind%1:N
    frameBounds= [1+(ii-1)*logtriggersperframe, ii*logtriggersperframe];
    [~, yy] = AcqReadAverage(hFile, header, frameBounds, signalBounds, [], []);
    cla(gca),plot(t,yy), pause(.01)
    yall(ii,:)= yy;
    tall(ii,:)= t;
    
end

fclose(hFile);