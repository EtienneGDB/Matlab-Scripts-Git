clear ;
close all ;
clc ;

% Task = {'TUG', '3min Free', '4min Free', '5min Free'};
Task = '4min Free';
Xunit = 'Frame';

% % Declare tables
% Percent_Time_Frozen = [NaN];
% Percent_Time_Tremor_Walking = [NaN];
% Percent_Time_Tremor_NOWalking = [NaN];
% Tremor_Level = [NaN];
% Brady_Level = [NaN];
% Dysk_Level = [NaN];
% 
% cd('F:\Projet RPQ\Extracted_Data')
% save([Task '_Percent_Time_Frozen.mat'],'Percent_Time_Frozen')
% save([Task '_Percent_Time_Tremor_Walking.mat'],'Percent_Time_Tremor_Walking')
% save([Task '_Percent_Time_Tremor_NOWalking.mat'],'Percent_Time_Tremor_NOWalking')
% save([Task '_Tremor_Level.mat'],'Tremor_Level')
% save([Task '_Brady_Level.mat'],'Brady_Level')
% save([Task '_Dysk_Level.mat'],'Dysk_Level')

% Participant to load
% Participant = 3;
for Participant = 5:14
    if Participant < 10
        Participant_Name = sprintf('P0%s',num2str(Participant));
    else
        Participant_Name = sprintf('P%s',num2str(Participant));
    end

    dataPath    = ['F:\Projet RPQ\RawData\',Participant_Name,'\',Participant_Name,'_Labo'];
    saveDir     = 'F:\Projet RPQ\ExportedData';
    addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
    addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
    addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\watch_extraction'))

    % %% Load FileNames
    cd([saveDir, '\FileNames'])
    load(['P', num2str(Participant), '_FileNames.mat'])
    for iFilNames = 1:length(keep_raw)
        if keep_raw{iFilNames,1}(1:3) == Task(1:3)
            FileName = keep_raw{iFilNames,1};
            break
        end
    end

    %% Load Watch Data
    cd([saveDir, '\Watch'])
    % load(['P', num2str(Participant), '_Watch_TUG_3 combinés.mat'])
    load(['P', num2str(Participant), '_Watch_', FileName, '.mat'])
    Freq = 50;

    %% Detect TUG parts
    % acc_vect = Module(Watch_Data.Ankle.Wrist_Raw_Accelerometer);
    % gyr_vect = Module(Watch_Data.Ankle.Wrist_Raw_Gyro);

    acc_ankle_vertical = Watch_Data.Ankle.Free_Acceleration(:,1);
    gyr_ankle_vertical = Watch_Data.Ankle.Wrist_Raw_Gyro(:,1);

    % figure; subplot(2,1,1); plot(acc_vect);
    % subplot(2,1,2); plot(gyr_vect);

    NpeaksToDetect = 2;
    RemoveThreePeaks = 0;
    showplot = 0;
    details = 1;
    if Participant == 2 & Task(1:3) == 'TUG'
        [TaskSeg, TurnSeg, SegNOWalking, Number_of_Steps] = RPQ_ThreePeaksDetection(acc_ankle_vertical, gyr_ankle_vertical, Freq, 9, 1, details); % P2
    else
        [TaskSeg, TurnSeg, SegNOWalking, Number_of_Steps] = RPQ_ThreePeaksDetection(acc_ankle_vertical, gyr_ankle_vertical, Freq, NpeaksToDetect, RemoveThreePeaks, details);  % others
    end
    % [TaskSeg] = RPQ_ThreePeaksDetection(acc_ankle_vertical, gyr_ankle_vertical, Freq, NpeaksToDetect, RemoveThreePeaks);  % others

% %% Detect walking periods
% config.acc.Fe = 50;
% 
% config.stepcounter.pks_min_dist = 50; % Inter peak distance (adj for walking)
% config.stepcounter.pks_min_height = 0.2; % Height of individual peaks (peaks = step)
% config.stepcounter.mov_avg = 30; % Window size for moving avg window
% config.stepcounter.nb_pts = 2048; % # PSD data points
% config.stepcounter.max_echelle = 20; % Maximal freq to keep at FFT output (Hz)
% config.stepcounter.win_length = 5*config.acc.Fe; % Length of window for FFT
% config.stepcounter.threshold_env = 0.05; % Threshold for the envelope (help find peaks)
% config.stepcounter.Fclow = 2.4; % Main filter lowpass freq (very important)
% config.stepcounter.filter_order = 4; % Order of butterworth filter
% config.stepcounter.freq_amp_coeff = 10; % Check the # of harmonics with high amp and high freq above this threshold
% config.stepcounter.fft_Hz_thresh = 3; % "high freq treshold to find harmonics for the freq_amp_coeff
% config.stepcounter.fft_Amp_thresh = 0.02; % "high amplitude threshold to find harmonics for the freq_amp_coeff
% 
% stepThresh = 1;
% 
% [totalStep, totalCad, numberOfStepPeriods, newLocs] = ...
%     step_detectAG_v2022(acc_ankle_vertical, config, stepThresh);
% 
% figure;
% plot(acc_vect(TaskSeg(1,1):TaskSeg(1,2)), 'color', [.5 .5 .5]);
% hold on,
% % plot(acc_vect_smooth, 'color', [.1 .7 .1 .7]);
% for i = 1:length(newLocs)
%     line([(newLocs(i)), (newLocs(i))],ylim,'Color',[.7 .1 .1 .5],'LineStyle','--')
% end
% title([num2str(totalStep) ' pas, cad:', num2str(totalCad), ', detected step periods: ' num2str(numberOfStepPeriods)]);

%% Data processing
Raw_acc_ankle = Watch_Data.Ankle.Wrist_Raw_Accelerometer;
acc_ankle = Watch_Data.Ankle.Free_Acceleration;
gyr_ankle = Watch_Data.Ankle.Wrist_Raw_Gyro;

Raw_acc_wrist = Watch_Data.Wrist.Wrist_Raw_Accelerometer;
acc_wrist = Watch_Data.Wrist.Free_Acceleration;
gyr_wrist = Watch_Data.Wrist.Wrist_Raw_Gyro;

% Filter Data
[b,a] = butter(2,2*[0.5 15]/Freq);
Filtre_acc_ankle = filtfilt(b,a,acc_ankle);
Filtre_gyr_ankle = filtfilt(b,a,gyr_ankle);

Filtre_acc_wrist = filtfilt(b,a,acc_wrist);
Filtre_gyr_wrist = filtfilt(b,a,gyr_wrist);

% for iaxe = 1:3
%     figure(1); subplot(3,1,iaxe); plot((acc_ankle(:,iaxe)))
%     figure(2); subplot(3,1,iaxe); plot((Filtre_acc_ankle(:,iaxe)))
%     figure(3); subplot(3,1,iaxe); plot((Raw_acc_ankle(:,iaxe)))
%     % figure(2); subplot(3,1,iaxe); plot(Filtre_gyr_ankle(:,iaxe))
% end
% spec_fft(acc_ankle(:,1),50,1)
% spec_fft(Raw_acc_ankle(:,1),50,1)

% Apply CWT to acc. data
axe = {'x','y','z'};
for iaxe = 1:length(axe)
    [coefcwt, freq] = cwt(acc_ankle(:,iaxe),'amor',Freq);
    assignin('base', ['PSD_' axe{iaxe} '_ankle'], abs(coefcwt).^2)
    assignin('base', ['f' axe{iaxe} '_ankle'], freq)

    [coefcwt, freq] = cwt(acc_wrist(:,iaxe),'amor',Freq);
    assignin('base', ['PSD_' axe{iaxe} '_wrist'], abs(coefcwt).^2)
    assignin('base', ['f' axe{iaxe} '_wrist'], freq)
end

% Calculate symptoms variables
for iaxe = 1:length(axe)
    [PowerBrady_ankle(:,iaxe), PowerMVT_ankle(:,iaxe), PowerFoG_ankle(:,iaxe), PowerTremor_ankle(:,iaxe), PowerTot_ankle(:,iaxe), ...
        Amplitude_ankle(:,iaxe), Amplitude_MVTband_ankle(:,iaxe), Amplitude_Tremorband_ankle(:,iaxe), ...
        PeakPower_Freq_ankle(:,iaxe), PeakPower_Freq_MVTband_ankle(:,iaxe), PeakPower_Freq_Tremorband_ankle(:,iaxe), ...
        PeakPower_ankle(:,iaxe), PeakPower_MVTband_ankle(:,iaxe), PeakPower_Tremorband_ankle(:,iaxe)] = ...
        RPQ_Compute_Power(eval(['PSD_' axe{iaxe} '_ankle']), eval(['f' axe{iaxe} '_ankle']), Filtre_acc_ankle(:,iaxe), false);

    [PowerBrady_wrist(:,iaxe), PowerMVT_wrist(:,iaxe), PowerFoG_wrist(:,iaxe), PowerTremor_wrist(:,iaxe), PowerTot_wrist(:,iaxe), ...
        Amplitude_wrist(:,iaxe), Amplitude_MVTband_wrist(:,iaxe), Amplitude_Tremorband_wrist(:,iaxe), ...
        PeakPower_Freq_wrist(:,iaxe), PeakPower_Freq_MVTband_wrist(:,iaxe), PeakPower_Freq_Tremorband_wrist(:,iaxe), ...
        PeakPower_wrist(:,iaxe), PeakPower_MVTband_wrist(:,iaxe), PeakPower_Tremorband_wrist(:,iaxe)] = ...
        RPQ_Compute_Power(eval(['PSD_' axe{iaxe} '_wrist']), eval(['f' axe{iaxe} '_wrist']), Filtre_acc_wrist(:,iaxe), false);
end

%% Detect Freezing
showplot = 0;
[PTF, Freezing_Index, Freezing_mask, FreezingSeg] = RPQ_FOG_Detection(acc_ankle_vertical, PowerFoG_ankle(:,1), PowerMVT_ankle(:,1), TaskSeg, showplot);

%% Detect Tremor
showplot =0;
[PTT_Walking, PTT_NOWalking, Tremor, Tremor_mask, TremorSeg] = RPQ_Tremor_Detection(acc_wrist, PowerTremor_wrist, sum(PowerMVT_wrist,2), PeakPower_Freq_wrist, TaskSeg, SegNOWalking, showplot);

%% Detect Bradykinesia
showplot = 0;
[Brady] = RPQ_Brady_Detection(acc_wrist, sum(PowerBrady_wrist,2), sum(PowerMVT_wrist,2), TaskSeg, SegNOWalking, showplot);

%% Detect Dyskinesia
showplot = 0;
[Dysk] = RPQ_Dysk_Detection(acc_wrist, sum(PowerMVT_wrist,2), Tremor_mask, TaskSeg, SegNOWalking, showplot);


%% Load and save sympt
% Load and save PTF
cd('F:\Projet RPQ\Extracted_Data')
load([Task '_Percent_Time_Frozen.mat'])
Percent_Time_Frozen(Participant,:) = mean(PTF);

cd('F:\Projet RPQ\Extracted_Data')
save([Task '_Percent_Time_Frozen.mat'],'Percent_Time_Frozen')


% Load and save PTF
cd('F:\Projet RPQ\Extracted_Data')
load([Task '_Percent_Time_Tremor_Walking.mat'])
load([Task '_Percent_Time_Tremor_NOWalking.mat'])
load([Task '_Tremor_Level.mat'])
Percent_Time_Tremor_Walking(Participant,:) = mean(PTT_Walking);
Percent_Time_Tremor_NOWalking(Participant,:) = mean(PTT_NOWalking);
Tremor_Level(Participant,:) = mean(Tremor);

cd('F:\Projet RPQ\Extracted_Data')
save([Task '_Percent_Time_Tremor_Walking.mat'],'Percent_Time_Tremor_Walking')
save([Task '_Percent_Time_Tremor_NOWalking.mat'],'Percent_Time_Tremor_NOWalking')
save([Task '_Tremor_Level.mat'],'Tremor_Level')


% Load and save Brady
cd('F:\Projet RPQ\Extracted_Data')
load([Task '_Brady_Level.mat'])
Brady_Level(Participant,:) = nanmean(Brady);

cd('F:\Projet RPQ\Extracted_Data')
save([Task '_Brady_Level.mat'],'Brady_Level')


% Load and save Dysk
cd('F:\Projet RPQ\Extracted_Data')
load([Task '_Dysk_Level.mat'])
Dysk_Level(Participant,:) = nanmean(Dysk);

cd('F:\Projet RPQ\Extracted_Data')
save([Task '_Dysk_Level.mat'],'Dysk_Level')




%% Plot Results
% Filter Data for plot
[b,a] = butter(2,2*[0.5 3]/Freq);
Filtre_acc_ankle_plot = filtfilt(b,a,acc_ankle);
[b,a] = butter(2,2*[0.5]/Freq,'low');
Filtre_gyr_ankle_plot = filtfilt(b,a,gyr_ankle);

Filtre_acc_ankle_plot = Filtre_acc_ankle_plot(:,1)./max(abs(Filtre_acc_ankle_plot(:,1)));
Filtre_gyr_ankle_plot = Filtre_gyr_ankle_plot(:,1)./max(abs(Filtre_gyr_ankle_plot(:,1)));

% Plot symptoms
Sympt = [Freezing_Index, Tremor, Brady, Dysk];
Symptoms = {'FoG', 'Tremor', 'Bradykinesia', 'Dyskinesia'};
figure;
for iSympt = 1:size(Sympt,2)
    subplot(size(Sympt,2),1,iSympt);
    plot(Sympt(:,iSympt))
    title(Symptoms{iSympt})
end


% Plot recap
figure;
subplot(3,1,1); plot(Dysk, '.', 'color', '#FFD700');
legend('Dysk.')
legend('boxoff')
title(['Participant ' num2str(Participant) ' - ' Task])

subplot(3,1,2); plot(Brady, '.', 'color', '#FF4500');
legend('Brady.')
legend('boxoff')
h = gca;
set( h, 'YDir', 'reverse' )

if Xunit(1) == 'T'
    dt = 1/Freq;
    Time = [0:dt:(length(Filtre_acc_ankle_plot)-1)/Freq];
    TaskSeg = TaskSeg/Freq;
    TurnSeg = TurnSeg/Freq;
    FreezingSeg = FreezingSeg/Freq;
    TremorSeg = TremorSeg/Freq;
else
    Time = 1:length(Filtre_acc_ankle_plot);
end


subplot(3,1,3); plot(Time, Filtre_acc_ankle_plot);
hold on
plot(Time, Filtre_gyr_ankle_plot);
TaskSeg_Legend = {'Walk. phase'};
% if length(TaskSeg) > 1
%     for i = 2:length(TaskSeg)
%         TaskSeg_Legend(i) = {''};
%     end
% end
%     ttt = [TaskSeg_Legend{1,1}, TaskSeg_Legend{1,2}, TaskSeg_Legend{1,3}];
% legend('Vertical acc.', 'Vertical gyr.', ttt, 'Turn. phase' )


% Show walking phases and symptomes
cSegActivity = '#00FF00';
cSegTurning  = '#696969';
cSegFreezing = '#000000';
cSegTremor   = '#1E90FF';
fill([TaskSeg(1,:) flip(TaskSeg(1,:),2)], [-1 -1 1 1], 'w', 'facecolor', cSegActivity, 'linestyle', 'none', 'facealpha', 0.4);
fill([TurnSeg(1,:) flip(TurnSeg(1,:),2)], [-1 -1 1 1], 'w', 'facecolor', cSegTurning, 'linestyle', 'none', 'facealpha', 0.4);

if size(FreezingSeg,1) >= 1
    fill([FreezingSeg(1,:) flip(FreezingSeg(1,:),2)], [-1 -1 1 1], 'black', 'linestyle', 'none', 'facealpha', 0.7);
    FOGLegend = 'FoG';
else
    FOGLegend = '';
end
if size(TremorSeg,1) >= 1
    fill([TremorSeg(1,:) flip(TremorSeg(1,:),2)], [-1 -1 1 1], 'w', 'facecolor', cSegTremor, 'linestyle', 'none', 'facealpha', 0.4);
    TremorLegend = 'Tremor';
    legend(TremorLegend, 'AutoUpdate', 'off')
else
    TremorLegend = '';
end

if size(FreezingSeg,1) < 1 & size(TremorSeg,1) >= 1
    legend('Vertical acc.', 'Vertical gyr.', 'Walk. phase', 'Turn. phase', TremorLegend , 'AutoUpdate', 'off')
else 
    legend('Vertical acc.', 'Vertical gyr.', 'Walk. phase', 'Turn. phase', FOGLegend, TremorLegend, 'AutoUpdate', 'off')
end
legend('boxoff')

if size(TaskSeg,1) > 1
    for iTaskSeg = 2:length(TaskSeg)
        fill([TaskSeg(iTaskSeg,:) flip(TaskSeg(iTaskSeg,:),2)], [-1 -1 1 1], 'w', 'facecolor', cSegActivity, 'linestyle', 'none', 'facealpha', 0.4);
    end
end
if size(TurnSeg,1) > 1
    for iTurnSeg = 2:length(TurnSeg)
        fill([TurnSeg(iTurnSeg,:) flip(TurnSeg(iTurnSeg,:),2)], [-1 -1 1 1], 'w', 'facecolor', cSegTurning, 'linestyle', 'none', 'facealpha', 0.4);
    end
end
if size(FreezingSeg,1) > 1
    for iFreezingSeg = 2:length(FreezingSeg)
        fill([FreezingSeg(iFreezingSeg,:) flip(FreezingSeg(iFreezingSeg,:),2)], [-1 -1 1 1], 'black', 'linestyle', 'none', 'facealpha', 0.7);
    end
end
if size(TremorSeg,1) > 1
    for iTremorSeg = 2:length(TremorSeg)
        fill([TremorSeg(iTremorSeg,:) flip(TremorSeg(iTremorSeg,:),2)], [-1 -1 1 1], 'w', 'facecolor', cSegTremor, 'linestyle', 'none', 'facealpha', 0.4);
    end
end

pause
clearvars -except Task Xunit

end





