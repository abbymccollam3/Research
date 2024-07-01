%% DP burst chirp wave velocity
clc
clearvars
close all
datapath = '/Users/abbymccollam/Desktop/Research/230206 Baseball players/'; % path where to load data
addpath(datapath) % optional path where matlab scripts are located

fs = 100e3;
data1 = [];
data2 = [];
meantime_delay = [];
stdtime_delay = [];
time_delay = [];
sensor_distance = 0.005;
number_of_bursts = 30;
period = 0.1;

% single burst and envelope segmentation
for k = 1:3             % subjects
    for k1 = 1:2        % sides (2) left/right
        %for k5 = 1:2    % channels (2)

task{1} = 'sit'; %fill in the name of your file / task
% task{2} = 'task2';
% task{3} = 'task3';
% task{4} = 'task4';
% task{5} = 'task5';

side{1} = 'left';
side{2} = 'right';

normal = 1;

for k2 = 1:number_of_bursts
            
            filename = [datapath 'S' num2str(k) '/segmented/S' num2str(k) '_' side{k1} '_stat_sit_c1' '.mat']; %fill in the first part of the file name
            load(filename); %S1/segmented/S1_left_squat
            data1(k2,:) = data_f_seg(1+(k2-1)*period*fs:1+(k2-1)*period*fs + period*fs,1);          
            data2(k2,:) = data_f_seg(1+(k2-1)*period*fs:1+(k2-1)*period*fs + period*fs,2);  
            data1_n(k2,:) = normalize(data_f_seg(1+(k2-1)*period*fs:1+(k2-1)*period*fs + period*fs,1),'zscore');          
            data2_n(k2,:) = normalize(data_f_seg(1+(k2-1)*period*fs:1+(k2-1)*period*fs + period*fs,2),'zscore');
            cross1 = find(data1_n(k2,:) >= 0.5);
            cross2 = find(data2_n(k2,:) >= 0.5);
            [m,peak1] = min(data1_n(k2,1:150)); %optimize this
            [m,peak2] = min(data2_n(k2,1:150)); %optimize this
            diff(k2) = peak2 - peak1;
end

figure
plot(mean(data1',2))
hold on
plot(mean(data2',2))
legend
  
time_delay(:,k1) = (diff'.*(1/fs)).*1000; %ms
meantime_delay(k,k1) = (mean(diff)*(1/fs))*1000; %ms
stdtime_delay(k,k1) = (std(diff)*(1/fs))*1000; %ms

figure(111)
scatter(ones(length(diff')).*(k+k1/10),sensor_distance./(diff.*(1/fs)),'filled','SizeData',100) % wave velocity
hold on
figure(222)
viola_plot_JDR_QG(sensor_distance./(diff.*(1/fs)),k*3+k1,'right','#A2142F')
figure(333)
viola_plot_JDR_QG((diff.*(1/fs)).*1000,k*3+k1,'right','#A2142F')

% spectral analysis

 f = [500:1:2000];
% [pxx1,f1] = pwelch(data1',[],[],f,fs);
% [pxx2,f1] = pwelch(data2',[],[],f,fs);
% pxx1_m = mean(pxx1,2);
% pxx2_m = mean(pxx2,2);
% 
% figure(k*1000)
% plot(f,pxx2_m./pxx1_m,'LineWidth',3)
% hold on
% hold on
% semilogx(f,pxx2_m)

figure
modalfrf(mean(data1',2),mean(data2',2),fs,9000,9000-1);
xlim([0.5 2])
ylim([-3.5 3.5])

% [frf,f_v] = modalfrf(mean(data1',2),mean(data2',2),fs,9000,9000-1);
% figure(k*100)
% plot(f_v,abs(frf),'LineWidth',3)
% xlim([0.5e3 2e3])
% hold on
% 
% figure(k*100+1)
% plot(f_v,rad2deg(angle(frf)),'LineWidth',3)
% xlim([0.5e3 2e3])
% hold on


% [s1,w,t] = spectrogram(mean(data1',2),1500,1500-1,f,fs);
% [s2,w,t] = spectrogram(mean(data2',2),1500,1500-1,f,fs);
% figure
% imagesc(abs(s2)./abs(s1))

    end

p = signrank(time_delay(:,1),time_delay(:,2))

end
