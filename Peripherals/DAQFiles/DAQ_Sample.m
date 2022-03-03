core = DAQengine(); %   start enginge
%% Good old cup
mac = '00-11-1C-02-04-00'; % mac address of the low frequency daq
samples = 2030;
core.initDAQ( samples, 256, mac, 1:256, 0 ); % init daq with 2030 samples per channel
num_av = 1; % number of averages
num_frames = 64; % number of frames
clear sigMat; % clear sigMat if there is one
core.lockDAQ(1); % set up daq for recording
tic;
sigMat = core.DoAcquisition( num_av, num_frames ); toc; % do the acquisition with the defined number of averages and frames
% Close daq
core.lockDAQ(0);
core.closeDAQ();
% save('try1', 'sigMat');
sigMatCal = sigMat - repmat(mean(sigMat,1),samples,1); % DC shift per channel
% figure; imagesc(sigMatCal)

sigMean = mean(sigMat,3);
sigMatCal = sigMean - repmat(mean(sigMean,1),samples,1); % DC shift per channel
figure; imagesc(sigMatCal)

fs = 40e6;
dts=1/fs;
t = 0:dts:((samples-1)*dts);
figure();
plot(t*1e6,sigMatCal(:,124));
