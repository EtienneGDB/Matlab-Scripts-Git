clear ;
close all ;
clc ;

% Participant to load
Participant = 3;
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

% %% Load Excel File
% excelfile = dir(fullfile(dataPath, '*.xlsx'));
% if isempty(excelfile)
%     display('SVP vous assurez que la liste des essais est compilÃ©e dans un fichier excel!')
%     return
% elseif size(excelfile,1)>=1
%     for i=size(excelfile,1)
%         tst = excelfile(i).name;
%         if isempty(strfind(tst,'~$'))
%             excelfile = excelfile(i).name;
%             break
%         end
%     end
% end
% file_partinfo = fullfile(dataPath, excelfile);
% % [num, txt, raw] = xlsread(file_partinfo,'Feuille1');
% % [num, txt, raw] = xlsread(file_partinfo,'Sheet1');
% [num, txt, raw] = xlsread(file_partinfo);
% keep_raw = raw; keep_txt = txt; keep_num = num;
% for iL = 1:length(keep_raw)
%     if keep_raw{iL,1}(end) == '1' | keep_raw{iL,1}(end) == '2'
%         keep_raw(iL,:) = [];
%         keep_txt(iL,:) = [];
%         keep_num(iL-1,:) = [];
%     end
%     if iL >= length(keep_raw)
%         break
%     end
% end

%% Load Watch Data
cd([saveDir, '\Watch'])
load(['P', num2str(Participant), '_Watch_TUG_3 combinés.mat'])
Freq = 50;

%% Detect TUG parts
acc_vect = Module(Watch_Data.Ankle.Wrist_Raw_Accelerometer);
% [TaskSeg] = RPQ_ThreePeaksDetection(acc_vect, 50, 9, 1, 1, 1); % P2
[TaskSeg] = RPQ_ThreePeaksDetection(acc_vect, 50, 2, 0, 1, 1);  % others

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
%     step_detectAG_v2022(acc_vect(TaskSeg(1,1):TaskSeg(1,2)), config, stepThresh);
% 
% figure;
% plot(acc_vect(TaskSeg(1,1):TaskSeg(1,2)), 'color', [.5 .5 .5]);
% hold on,
% % plot(acc_vect_smooth, 'color', [.1 .7 .1 .7]);
% for i = 1:length(newLocs)
%     line([(newLocs(i)), (newLocs(i))],ylim,'Color',[.7 .1 .1 .5],'LineStyle','--')
% end
% title([num2str(totalStep) ' pas, cad:', num2str(totalCad), ', detected step periods: ' num2str(numberOfStepPeriods)]);

%% Detect Freezing
[b,a] = butter(2,2*[0.5 15]/50);
FiltreAcc_Ankle = filtfilt(b,a,Watch_Data.Ankle.Wrist_Raw_Accelerometer);
figure;
axe = {'x','y','z'};
for iaxe = 1:3
    subplot(3,1,iaxe); plot(FiltreAcc_Ankle(:,iaxe))
    title(['Filtered Ankle Accleration ', axe{iaxe}])
end

% Gyro
[b,a] = butter(2,2*[0.5 15]/50);
FiltreGyr_Ankle = filtfilt(b,a,Watch_Data.Ankle.Wrist_Raw_Gyro);
figure;
axe = {'x','y','z'};
for iaxe = 1:3
    subplot(3,1,iaxe); plot(FiltreGyr_Ankle(:,iaxe))
%     hold on
%     subplot(3,1,iaxe); plot(FiltreAcc_Ankle(:,iaxe))
    title(['Filtered Ankle Gyro ', axe{iaxe}])
end

% figure;
% for iWalking = 1:length(TaskSeg)
%     figure(1); subplot(3,1,iWalking); plot(FiltreAcc_Ankle(TaskSeg(iWalking,1):TaskSeg(iWalking,2),1))
%     title(['Filtered Ankle Accleration ', num2str(iWalking)])
%     spec_fft(FiltreAcc_Ankle(TaskSeg(iWalking,1):TaskSeg(iWalking,2),1),50,1)
% end
[coefAcc fAcc] = cwt(FiltreAcc_Ankle(:,1),'amor',50);
norm = abs(coefAcc).^2;

[PowerMVT, ...
PowerFoG, ...
PowerTremor, ...
PowerTot, ...
Amplitude, ...
Amplitude_MVTband, ...
Amplitude_Tremorband, ...
PeakPower_Freq, ...
PeakPower_Freq_MVTband, ...
PeakPower_Freq_Tremorband, ...
PeakPower, ...
PeakPower_MVTband, ...
PeakPower_Tremorband] = ...
RPQ_Compute_Power(norm, fAcc);

figure;
Thresh_Freezing = 3;
for iWalking = 1:length(TaskSeg)
    WalkingTime = length(TaskSeg(iWalking,1):TaskSeg(iWalking,2));
    Freezing_Index = PowerFoG(TaskSeg(iWalking,1):TaskSeg(iWalking,2))./PowerMVT(TaskSeg(iWalking,1):TaskSeg(iWalking,2));
    subplot(3,1,iWalking); plot(Freezing_Index)
    hold on;
    line([0, WalkingTime], [Thresh_Freezing, Thresh_Freezing], 'color', 'r')
    Freezing_State = find(Freezing_Index > Thresh_Freezing);
    if ~isempty(Freezing_State)
        PTF(:,iWalking) = (length(Freezing_State)/WalkingTime)*100;
    else
        PTF(:,iWalking) = 0;
    end
    text(length(Freezing_Index)-length(Freezing_Index)*0.2,max(Freezing_Index)-max(Freezing_Index)*0.2,['%timeFrozen:', num2str(PTF(iWalking))])
    title(['Freezing index TUG', num2str(iWalking)])
end

cd('F:\Projet RPQ\Extracted_Data')
load('TUG_Percent_Time_Frozen.mat')
Percent_Time_Frozen(Participant,:) = PTF;

% Save
% Percent_Time_Frozen = [NaN NaN NaN];
cd('F:\Projet RPQ\Extracted_Data')
save('TUG_Percent_Time_Frozen.mat','Percent_Time_Frozen')

%% Detect Tremor
[b,a] = butter(2,2*[0.5 15]/50);
FiltreAcc_Wrist = filtfilt(b,a,Watch_Data.Wrist.Wrist_Raw_Accelerometer);
figure;
axe = {'x','y','z'};
for iaxe = 1:3
    subplot(3,1,iaxe); plot(FiltreAcc_Wrist(:,iaxe))
    title(['Filtered Wrist Accleration ', axe{iaxe}])
end

% figure;
% for iWalking = 1:length(TaskSeg)
%     figure(1); subplot(3,1,iWalking); plot(FiltreAcc_Wrist(TaskSeg(iWalking,1):TaskSeg(iWalking,2),3))
%     title(['Filtered Wrist Accleration ', num2str(iWalking)])
%     spec_fft(FiltreAcc_Wrist(TaskSeg(iWalking,1):TaskSeg(iWalking,2),3),50,1)
% end

for iaxe = 1:3
    [coefAcc fAcc] = cwt(FiltreAcc_Wrist(:,iaxe),'amor',50);
    norm = abs(coefAcc).^2;

    [PowerMVT(:,iaxe),...
    PowerFoG(:,iaxe),...
    PowerTremor(:,iaxe),...
    PowerTot(:,iaxe),...
    Amplitude(:,iaxe),...
    Amplitude_MVTband(:,iaxe),...
    Amplitude_Tremorband(:,iaxe),...
    PeakPower_Freq(:,iaxe),...
    PeakPower_Freq_MVTband(:,iaxe),...
    PeakPower_Freq_Tremorband(:,iaxe),...
    PeakPower(:,iaxe),...
    PeakPower_MVTband(:,iaxe),...
    PeakPower_Tremorband(:,iaxe)] = RPQ_Compute_Power(norm, fAcc);

%     figure
%     subplot(2,1,1); plot(PowerTremor(:,iaxe))
%     ylim([0 1.5])
%     subplot(2,1,2); plot(PeakPower_Freq(:,iaxe))
%     hold on
%     line([TaskSeg(:,1), TaskSeg(:,1)], [0, 10], 'color', 'r')
%     line([TaskSeg(:,2), TaskSeg(:,2)], [0, 10], 'color', 'r')
end
% PeakPower_Tremorband_Sum = sum(PeakPower_Tremorband,2);
PowerTremor_Sum = sum(PowerTremor,2);
PowerMVT_Sum = sum(PowerMVT,2);
% R1 = PeakPower_Tremorband_Sum./PowerTremor_Sum;
% R2 = PowerTremor_Sum./PowerMVT_Sum;
% figure;plot(PeakPower_Freq(:,3))
% figure;plot(PeakPower_Tremorband_Sum)
% figure;plot(PowerTremor_Sum)
% figure;plot(PowerMVT_Sum)
% figure;plot(R1)
% figure;plot(R2)

figure;
for iWalking = 1:length(TaskSeg)
    WalkingTime = length(TaskSeg(iWalking,1):TaskSeg(iWalking,2));
    % Main Tremor Axis
    MedTremorAxis = median(PowerTremor(TaskSeg(iWalking,1):TaskSeg(iWalking,2),:));
    Main_Tremor_Axis = find(MedTremorAxis==max(MedTremorAxis));
    
    Tremor_Values = PowerTremor(TaskSeg(iWalking,1):TaskSeg(iWalking,2),Main_Tremor_Axis);
    NO_Tremor_Values = PowerMVT_Sum(TaskSeg(iWalking,1):TaskSeg(iWalking,2));
    
    % Percent time with Tremor
    Tremor_Freq = PeakPower_Freq(TaskSeg(iWalking,1):TaskSeg(iWalking,2),Main_Tremor_Axis);
    PTT(:,iWalking) = (length(Tremor_Freq(Tremor_Freq>=4 & Tremor_Freq<=7))/WalkingTime)*100;
%     figure; plot(Tremor_Freq)
    
    % Tremor values
    if ~isempty(Tremor_Values(Tremor_Freq>=4 & Tremor_Freq<=7))
        idx = find(Tremor_Values(Tremor_Freq>=4 & Tremor_Freq<=7));
%         Tremor(:,iWalking) = median(Tremor_Values(Tremor_Freq>=4 & Tremor_Freq<=7));
        Tremor(:,iWalking) = median(Tremor_Values(idx)./NO_Tremor_Values(idx));
    else
        Tremor(:,iWalking) = 0;
    end
    
    % NO Tremor values
    if ~isempty(NO_Tremor_Values)
        NO_Tremor(:,iWalking) = median(NO_Tremor_Values);
    else
        NO_Tremor(:,iWalking) = 0;
    end
    
    subplot(3,1,iWalking); plot(Tremor_Values./NO_Tremor_Values)
    text(length(Tremor_Values)-length(Tremor_Values)*0.2,max(Tremor_Values)-max(Tremor_Values)*0.2,['%timeTremor:', num2str(PTT(iWalking))])
    title(['Tremor TUG', num2str(iWalking)])

%     figure; plot(Tremor_Values)
%     figure; plot(NO_Tremor_Values)
end

cd('F:\Projet RPQ\Extracted_Data')
load('TUG_Percent_Time_Tremor.mat')
load('TUG_Tremor_Level.mat')
load('TUG_NO_Tremor_Level.mat')
Percent_Time_Tremor(Participant,:) = PTT;
Tremor_Level(Participant,:) = Tremor;
NO_Tremor_Level(Participant,:) = NO_Tremor;

% Save
% Percent_Time_Tremor = [NaN NaN NaN];
% Tremor_Level = [NaN NaN NaN];
% NO_Tremor_Level = [NaN NaN NaN];
cd('F:\Projet RPQ\Extracted_Data')
save('TUG_Percent_Time_Tremor.mat','Percent_Time_Tremor')
save('TUG_Tremor_Level.mat','Tremor_Level')
save('TUG_NO_Tremor_Level.mat','NO_Tremor_Level')




