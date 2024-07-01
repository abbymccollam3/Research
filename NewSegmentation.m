%% Achilles recording segmentation
% script to semi-automatically segment different tasks out of continuous
% recordings.

% Part 1 Load and visualize raw data
close all
clearvars

loadpath = '/Users/abbymccollam/Desktop/221104 burst chirp repeatability/';
addpath(loadpath);

filename = ['M1_sit_ext_stand_squat_right_chirp_r3.mat']; % fill in filename
load(filename)

%fs = 1e4;
dt = 1/fs;
t = [0:dt:dt*length(data)-dt]';

fc = 235;
f_sqr = 5;
w_sin = fc * (2*pi);
w_sqr = f_sqr * (2*pi);
a_sin = 3;
duty = 50;


data_filt = bandpass(data,[50 5e3], fs);

%burst = normalize(burst,'zscore'); % z score norm of signal
data_f = normalize(data_filt,'zscore'); % z score norm of signal

figure(1)
plot(t,data);
ylabel('Norm accel (g)','FontWeight','bold')
xlabel('Time (s)')

figure(1)
plot(t,data_f);
ylabel('Norm accel (g)','FontWeight','bold')
xlabel('Time (s)')

% [yupper,ylower] = envelope(data_filt,30,'peak');

% figure(2)
% subplot(2,1,1)
% plot(t,data(:,1));
% hold on
% plot(t,ylower(:,1),'Color',[0.5 0.1 1],'Linewidth',1);
% hold on
% plot(t,yupper(:,1),'Color',[0.5 0.1 1],'Linewidth',1);
% ylabel('acceleration envelope (-)','FontWeight','bold')
% xlabel('Time (s)','FontWeight','bold')
% title('Envelope1')
% 
% subplot(2,1,2)
% plot(t,data_f(:,2));
% hold on
% plot(t,ylower(:,2),'Color',[0.5 0.1 1],'Linewidth',1);
% hold on
% plot(t,yupper(:,2),'Color',[0.5 0.1 1],'Linewidth',1);
% ylabel('acceleration envelope (-)','FontWeight','bold')
% xlabel('Time (s)','FontWeight','bold')
% title('Envelope1')

% envelope_v = yupper-ylower;
% [envelope_envelope_v,dum] = envelope(envelope_v,1500,'peak');
% figure(3)
% p1 = plot(t,envelope_v);
% hold on
% plot(t,envelope_envelope_v,'Linewidth',2);
% ylim([0 inf]);
% ylabel('acceleration envelope (-)','FontWeight','bold')
% xlabel('Time (s)','FontWeight','bold')
% title('Envelope2')

%% Part 2 Semi-auto data segmentation
% Segments data automatically. Default 5s segmentation, input is
% approximate start time. Script will segment data at closest burst starting sample.

period = 1/f_sqr; % time of one burst period

clear segmentname data_seg data_f_seg startSample endSample
prompt = {'Start time segment (s)?'};
        dlg_title =  'Recording Setup';
        num_lines = 1;
        defaultans = {'00'};    
        N=30;
        answer1 = inputdlg(prompt,dlg_title,[num_lines,  length(dlg_title)+N],defaultans);
        
prompt = {'Task?'};
        dlg_title =  'Segment Name';
        num_lines = 1;
        defaultans = {'sit'};    
        N=30;
        answer2 = inputdlg(prompt,dlg_title,[num_lines,  length(dlg_title)+N],defaultans);
        
segmentname = answer2{1};

startSample= str2double(answer1{1}) * fs; % sets recording time

startSample = ceil(startSample/(fs*0.2))*fs*0.2;
endSample = startSample+5*fs;

data_seg = data(startSample:endSample,:);
data_f_seg = data_filt(startSample:endSample,:);

figure
plot(data_f_seg);
ylabel('Normalized acceleration (-)','FontWeight','bold')
xlabel('Time (s)')
title('Achilles Mic')


filename=['name.mat'];

% Ask if you'd like to save data: 
choice_to_save = questdlg('Save Data?','Save');
% Handle response
switch choice_to_save
    case 'Yes'
            prompt = {'Enter File Name:'};
            dlg_title = 'Input';
            num_lines = 1;
            dir=[loadpath 'segmented/'];
           fullFileName=[dir filename];
            defaultans = {fullFileName};
            answer = inputdlg(prompt,dlg_title,[num_lines, length(defaultans{1})+10], defaultans);
            fullFileName=answer{1}; %Input File Name to be saved.
           
           save(fullFileName,'data_seg','data_f_seg')
           case 'No'
end
