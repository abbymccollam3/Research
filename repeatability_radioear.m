%% Consider of medial and lateral signal seperately - envelope repeatabillity

%make sure to define file names

clc
clearvars
close all
addpath('/Users/abbymccollam/Desktop/221104 burst chirp repeatability/segmented/') % path where to load data
%addpath('fillinpath') % optional path where matlab scripts are located

% single burst and envelope segmentation
for k = 1:4             % tasks
    for k1 = 2:2        % sides (2) left/right
        for k5 = 1:2    % channels (2)

number_of_bursts = 20;
fs = 5e4;
period = 0.1;

task{1} = 'sit'; %fill in the name of you file / task
task{2} = 'extension';
task{3} = 'stand';
task{4} = 'squat';
%task{5} = 'task5';

side{1} = 'left';
side{2} = 'right';

normal = 1;

channel = k5;
            
for k3 = 1:3
            filename = ['Abby_' side{k1} '_' task{k} '_' num2str(k3) '.mat']; % fill in the first part of the file name
            load(filename);
            data(:,k3) = data_f_seg(1:(number_of_bursts+1)*period*fs,channel);        
            %data(:,k3) = data_f_seg(1:(number_of_bursts)*period*fs,channel);          
end 

[y_upper y_lower] = envelope(data,100,'peak');

figure
plot(data)
hold on
plot(y_lower)
hold on 
plot(y_upper)


data_env = y_upper%-y_lower; 
data_env = data_env(period*fs+1:end,:); % this line skips the first burts because of transient boundary issues. Check original data_env to see if its the case. This also why we consider one more burst than 'number_of_bursts' in line 31.


% figure
% plot(data_env)

% segment individual bursts and calc mean

for k4 = 1:number_of_bursts
    
    if normal == 1
    data1_env(:,k4) = normalize(data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1),"range");    % envelope rep 1
    data2_env(:,k4) = normalize(data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2),"range");    % envelope rep 2
    data3_env(:,k4) = normalize(data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3),"range");    % envelope rep 3
    else if normal == 0 
    data1_env(:,k4) = data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1);    % envelope rep 1
    data2_env(:,k4) = data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2);    % envelope rep 2
    data3_env(:,k4) = data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3);    % envelope rep 3
    end
    end

end

data_m_env = [mean(data1_env,2) mean(data2_env,2) mean(data3_env,2)]; % mean envelope per rep

% reproducibility analysis

data_mm_env = mean(data_m_env,2); % mean envelope over all reps

figure(k+1000) 
subplot(1,2,k5)
plot(data_m_env)
hold on
plot(data_mm_env,'LineWidth',3,'Color','k')

mae11_env(1) = mean(abs(data_mm_env-data_m_env(:,1)))./max(data_mm_env); 
mae11_env(2) = mean(abs(data_mm_env-data_m_env(:,2)))./max(data_mm_env); 
mae11_env(3) = mean(abs(data_mm_env-data_m_env(:,3)))./max(data_mm_env); 

if k5 == 1   
figure(111)
viola_plot_JDR_QG(mae11_env, k+(k1-1)*0.2, 'right', '#A2142F')
hold on 
xlim([0.5 5.0])
ylim([0 0.3])

elseif k5 == 2
figure(222)
viola_plot_JDR_QG(mae11_env, k+(k1-1)*0.2, 'right', '#0072BD') %where it is plotting on x-axis
hold on 
xlim([0.5 5.0])
ylim([0 0.3])
box on

%viola_plot_JDR_QG(mae_env(:,3), 2, 'right', '#A2142F')
    end
    end
    end
end
%% Consider mean of medial and lateral signal - envelope & PSD repeatabillity

clc
clearvars
close all
addpath('/Users/abbymccollam/Desktop/221104 burst chirp repeatability/segmented/') 
%addpath('/Users/abbymccollam/Desktop/data q/segmented/')

% single burst and envelope segmentation
for k = 1:4             % tasks
        for k1 = 2:2    % left or right

number_of_bursts = 20;
fs = 5e4;
period = 0.1;

task{1} = 'sit'; %fill in the name of you file / task
task{2} = 'extension';
task{3} = 'stand';
task{4} = 'squat';
%task{5} = 'task5';

side{1} = 'left';
side{2} = 'right';

normal = 1;
            
for k3 = 1:3
            filename = ['Abby_' side{k1} '_' task{k} '_' num2str(k3) '.mat']; 
            load(filename);
            data1(:,k3) = data_f_seg(1:(number_of_bursts+1)*period*fs,1);          
            filename = ['Abby_' side{k1} '_' task{k} '_' num2str(k3) '.mat']; 
            load(filename);
            data2(:,k3) = data_f_seg(1:(number_of_bursts+1)*period*fs,2);          
end 


[y_upper1 y_lower1] = envelope(data1,50,'peak');
[y_upper2 y_lower2] = envelope(data2,50,'peak');

figure
plot(data1)
hold on
plot(y_lower1)
hold on 
plot(y_upper1)

data_env1 = y_upper1;%-y_lower; 
data_env1 = data_env1(period*fs+1:end,:);
data_env2 = y_upper2;%-y_lower; 
data_env2 = data_env2(period*fs+1:end,:)

data_env = (data_env1+data_env2)./2; % mean envelope of ch1 & ch2

% figure
% plot(data_env)

% segment individ burst and calc mean

for k4 = 1:number_of_bursts
    
    if normal == 1
    data1_env(:,k4) = normalize(data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1),"range");    % individ envelopes rep 1
    data2_env(:,k4) = normalize(data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2),"range");    % individ envelopes rep 2
    data3_env(:,k4) = normalize(data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3),"range");    % individ envelopes rep 3

    data11(:,k4) =  normalize(data1(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1),1);    % individ burst rep 1 - ch1
    data12(:,k4) =  normalize(data1(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2),1);    % individ burst rep 2 - ch1
    data13(:,k4) =  normalize(data1(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3),1);    % individ burst rep 3 - ch1 

    data21(:,k4) =  normalize(data2(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1),1);    % individ burst rep 1 - ch2
    data22(:,k4) =  normalize(data2(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2),1);    % individ burst rep 2 - ch2
    data23(:,k4) =  normalize(data2(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3),1);    % individ burst rep 3 - ch2

    else if normal == 0 
    data1_env(:,k4) = data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1);    % individ envelopes rep 1
    data2_env(:,k4) = data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2);    % individ envelopes rep 2
    data3_env(:,k4) = data_env(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3);    % individ envelopes rep 3

    data11(:,k4) = data1(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1);    % individ burst rep 1 - ch1
    data12(:,k4) = data1(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2);    % individ burst rep 2 - ch1
    data13(:,k4) = data1(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3);    % individ burst rep 3 - ch1 

    data21(:,k4) = data2(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,1);    % individ burst rep 1 - ch2
    data22(:,k4) = data2(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,2);    % individ burst rep 2 - ch2
    data23(:,k4) = data2(period*fs*(k4-1)+1:period*fs*(k4-1)+period*fs,3);    % individ burst rep 3 - ch2 


    end
    end

end

f = [500:1:2000];

[pxx11,f1] = pwelch(data11,[],[],f,fs);
[pxx12,f1] = pwelch(data12,[],[],f,fs);
[pxx13,f1] = pwelch(data13,[],[],f,fs);

[pxx21,f1] = pwelch(data21,[],[],f,fs);
[pxx22,f1] = pwelch(data22,[],[],f,fs);
[pxx23,f1] = pwelch(data23,[],[],f,fs);

pxx11_m = mean(pxx11,2); % mean psd ch1 - rep 1
pxx12_m = mean(pxx12,2); % mean psd ch1 - rep 2
pxx13_m = mean(pxx13,2); % mean psd ch1 - rep 3

pxx21_m = mean(pxx21,2); % mean psd ch2 - rep 1
pxx22_m = mean(pxx22,2); % mean psd ch2 - rep 2
pxx23_m = mean(pxx23,2); % mean psd ch2 - rep 3

% averaging for ch1/2 is slightly different than for envelope. For single
% chirp spectrum here, while envelope is done for the entire chirp train.
% Would be possible do to the same for envlope, however envelope is prone
% to initialization artefacts.

pxx1_m = (pxx11_m + pxx21_m)./2; % mean psd rep 1 (ch1 & 2)
pxx2_m = (pxx12_m + pxx22_m)./2; % mean psd rep 2 (ch1 & 2)
pxx3_m = (pxx13_m + pxx22_m)./2; % mean psd rep 3 (ch1 & 2)


pxx_m = [pxx1_m pxx2_m pxx3_m]                                          % matrix with 3 mean psd's of 20 bursts per rep
pxx_mm = mean(pxx_m,2);                                                 % single mean psd over all reps

data_m_env = [mean(data1_env,2) mean(data2_env,2) mean(data3_env,2)];   % matrix with 3 mean envelopes of 20 bursts per rep
data_mm_env = mean(data_m_env,2);                                       % single mean envelope over all reps

figure(k+1000) 
%subplot(1,2,k1)
plot(data_m_env)
hold on
plot(data_mm_env,'LineWidth',3,'Color','k')
xlabel('Time Sample (-)')
ylabel('Normalized Amplitude (-)')

figure(k+2000) 
%subplot(1,2,k1)
semilogy(f,pxx_m)
hold on
semilogy(f,pxx_mm,'LineWidth',3,'Color','k')
xlabel('Frequency (Hz)')
ylabel('PSD magnitude (-)')

mae11_env(1) = mean(abs(data_mm_env-data_m_env(:,1)))./max(data_mm_env); 
mae11_env(2) = mean(abs(data_mm_env-data_m_env(:,2)))./max(data_mm_env); 
mae11_env(3) = mean(abs(data_mm_env-data_m_env(:,3)))./max(data_mm_env); 

mae11_pxx(1) = mean(abs(pxx_mm-pxx_m(:,1)))./max(pxx_mm); 
mae11_pxx(2) = mean(abs(pxx_mm-pxx_m(:,2)))./max(pxx_mm); 
mae11_pxx(3) = mean(abs(pxx_mm-pxx_m(:,3)))./max(pxx_mm); 

if k1 == 1   
figure(111)
viola_plot_JDR_QG(mae11_env, k, 'right', '#A2142F')
hold on 
xlim([0.5 5.75])
ylim([0 0.3])
xlabel('Task - Sit Extension Stand Squat')
ylabel('MAE (-)')

figure(1111)
viola_plot_JDR_QG(mae11_pxx, k, 'right', '#A2142F')
hold on 
xlim([0.5 5.75])
ylim([0 0.3])
xlabel('Task - Sit Extension Stand Squat')
ylabel('MAE PSD (-)')

elseif k1 == 2
figure(111)
viola_plot_JDR_QG(mae11_env, k+0.2, 'right', '#0072BD')
hold on 
xlim([0.5 5.0])
ylim([0 0.3])
box on
xlabel('Task - Sit Extension Stand Squat')
ylabel('MAE (-)')

figure(1111)
viola_plot_JDR_QG(mae11_pxx, k+0.2, 'right', '#0072BD')
hold on
xlim([0.5 5.0])
ylim([0 0.3])
xlabel('Task - Sit Extension Stand Squat')
ylabel('MAE PSD (-)')
box on

%viola_plot_JDR_QG(mae_env(:,3), 2, 'right', '#A2142F')
    end
        end
end