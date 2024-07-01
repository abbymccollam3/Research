data_filt = bandpass(data, [200 2e3], fs);
%% Achilles recording segmentation
% script to semi-automatically segment different tasks out of continuous
% recordings.

% Part 1 Load and visualize raw data
close all
clearvars

loadpath = '/Users/abbymccollam/Desktop/Research/230206 Baseball players/S1';
addpath(loadpath);

filename = 'S1_chirp_stat_r.mat'; % fill in filename
load(filename)

flow_filt = 500;
fhigh_filt = 2e3;

dt = 1/fs;
t = [0:dt:dt*length(data)-dt]';

data_filt = bandpass(data,[flow_filt fhigh_filt], fs);

data_f = normalize(data_filt,'zscore'); % z score norm of signal

figure(1)
plot(t,data_filt);
ylabel('Norm accel (g)','FontWeight','bold')
xlabel('Time (s)')

%% Part 2 Semi-auto data segmentation
% Segments data automatically. Default 5s segmentation, input is
% approximate start time. Script will segment data at closest burst starting sample.

%period = 1/f_sqr; % time of one burst period

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

filename=[segmentname '.mat'];

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
           
           save(fullFileName,'data_seg','data_f_seg','fs')
           case 'No'
end

