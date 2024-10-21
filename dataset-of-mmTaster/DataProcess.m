%%%%%%%%%%%%%%%%
% generate data of different spectrum  
% input size: 7000 x 64 x 4
%       7000: frames at the speed of ~200 frames per second, collected while the radar rotates around the liquid target for about 540 degrees
%         64: samples collected in one chirp, from 57GHz to 63GHz
%          4: number of receiving antennas

% output: 7000 files of .mat
%       size: 8 x 2  x 10 x 4
%          8: number of spectrum channels 
%          2: magnitude and phase
%         10: range bins of 35-45cm
%          4: number of receiving antennas

%%%%%%%%%%%%%%%%
clear;clc;clear all;close all;
path= '.\rawdata\';
FrameNum = 7000; BW = 6e9; rxn = 4; %57GHz-63GHz
Rres = 3.0e8/(2*BW); samples = 64; Zeropad = samples; Rmax = Rres * samples;
range = linspace(0,1-1/Zeropad,Zeropad)*Rmax;
fs = 2e6; n2 = (0:(Zeropad-1))*fs/Zeropad;
subdirpath = fullfile( path, '*.mat' );  dat = dir( subdirpath );
%set(0,'DefaultFigureVisible','on');
%set(gca,'Box','off');
for j = 1  : size(dat,1)
    datpath = fullfile( path, dat( j ).name); %rawdata
    name = strsplit(dat( j ).name,'.');
    name = name{1}
    savedir = fullfile(path,name);
    if exist(savedir) == 0
        mkdir(savedir);
    end
    load(datpath);
    feature_list = permute(matData,[3,1,2]);%rxn * FrameNum * samples
    
    % 减均值去零频
    x1 = feature_list(1,:,:);%一根天线的数据
    meanx1 = mean(x1,3);
    meanx1 = reshape(meanx1,1,FrameNum,1);
    meanx1 = repmat(meanx1,1,1,samples);
    x1 = x1 - meanx1;
    
    x2 = feature_list(2,:,:);%一根天线的数据
    meanx1 = mean(x2,3);
    meanx1 = reshape(meanx1,1,FrameNum,1);
    meanx1 = repmat(meanx1,1,1,samples);
    x2 = x2 - meanx1;
    
    x3 = feature_list(3,:,:);%一根天线的数据
    meanx1 = mean(x3,3);
    meanx1 = reshape(meanx1,1,FrameNum,1);
    meanx1 = repmat(meanx1,1,1,samples);
    x3 = x3 - meanx1;
    
    x4 = feature_list(4,:,:);%一根天线的数据
    meanx1 = mean(x4,3);
    meanx1 = reshape(meanx1,1,FrameNum,1);
    meanx1 = repmat(meanx1,1,1,samples);
    x4 = x4 - meanx1;
    
    %加窗 做fft
    ends = 42;count = 1;
    Zeropad=ends;
    Rres = 3.0e8*samples/ends/(2*BW);
    Rmax = Rres * ends;
    range = linspace(0,1-1/Zeropad,Zeropad)*Rmax;
    
    Hann_window = hann(ends,'periodic');
    ScaleHannWin = 1/sum(Hann_window);
    Hann_window = repmat(Hann_window,1,FrameNum);
    feature_list8=zeros(8,4,7000,42);
    x1 = squeeze(x1);x2 = squeeze(x2);x3 = squeeze(x3);x4 = squeeze(x4);
    
    while(ends<=samples)
        feature_list8(count,1,:,:) = fft(x1(:,(ends-41):ends).*Hann_window',Zeropad,2)*1*ScaleHannWin;
        feature_list8(count,2,:,:) = fft(x2(:,(ends-41):ends).*Hann_window',Zeropad,2)*1*ScaleHannWin;
        feature_list8(count,3,:,:) = fft(x3(:,(ends-41):ends).*Hann_window',Zeropad,2)*1*ScaleHannWin;
        feature_list8(count,4,:,: ) = fft(x4(:,(ends-41):ends).*Hann_window',Zeropad,2)*1*ScaleHannWin;
%               xx1 = fft(x1(:,(ends-41):ends).*Hann_window',Zeropad,2)*1*ScaleHannWin;
%               plot(range,abs(squeeze(xx1(650,:))));hold on;
%              xlim([0.3,0.6])
        ends = ends + 3;
        count = count + 1;
    end
    
    % [0.304,0.342,0.380,0.419,0.45,0.49,0.53] 9:15 该实验距离为35-45cm
    k1 = 8;
    k2 = 17;
    % 每一帧数据进行处理 最后保存每一帧的数据
    for i = 1:FrameNum
        save_mat = zeros(8,2,10,4);
        for cha = 1:8
            Rx1 = squeeze(feature_list8(cha,1,i,:));
            Rx2 = squeeze(feature_list8(cha,2,i,:));
            Rx3 = squeeze(feature_list8(cha,3,i,:));
            Rx4 = squeeze(feature_list8(cha,4,i,:));
            %fft_all ( Rx1,Rx2,Rx3,Rx4,k1,k2,range);
            [save_mat(cha,:,:,:)] = fft_all2 ( Rx1,Rx2,Rx3,Rx4,k1,k2,range);
        end
        str = strcat(num2str(i),'.mat');
        savePath = fullfile(savedir, str);
        save(savePath,'save_mat');
    end
end
disp('******************************over******************************')