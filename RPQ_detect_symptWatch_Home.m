clear ;
close all ;
clc ;

% Data = {'Watch', 'WatchXsens'};
Data = 'WatchXsens';
window = 10;

% Participant to load
% Participant = 3;
for Participant = 5%:9
    if Participant < 10
        Participant_Name = sprintf('P0%s',num2str(Participant));
    else
        Participant_Name = sprintf('P%s',num2str(Participant));
    end

    dataPath    = ['F:\Projet RPQ\ExportedData\Home\', Data];
    saveDir     = 'F:\Projet RPQ\ExtractedData\Home';
    addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
    addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
    addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\watch_extraction'))

    %% Load Watch Data
    if Data(end) == 's'
        FileName = '_WatchXsens_Home.mat';
        FileName_Time = '_TimeWatchXsens_Home.mat';
        FileName_rawTime = '_rawTimeWatchXsens_Home.mat';
    else
        FileName = '_Watch_Home.mat';
        FileName_Time = '_TimeWatch_Home.mat';
        FileName_rawTime = '_rawTimeWatch_Home.mat';
    end
    cd(dataPath)
    load(['P', num2str(Participant), FileName])
    load(['P', num2str(Participant), FileName_Time])
    load(['P', num2str(Participant), FileName_rawTime])
    Freq = 50;

    Location = fieldnames(Watch_allData);
    if Location{1}(1) == 'A'
        idx_Ankle = 1;
        idx_Wrist = 2;
    else
        idx_Ankle = 2;
        idx_Wrist = 1;
    end

    Sensors = {'Wrist_Raw_Accelerometer', 'Wrist_Raw_Gyro'};

    % Synchronisation signals
    for iLoc = 1:length(Location)
        for iSen = 1:length(Sensors)
            Temp_FirstMin(iLoc,iSen) = floor(Time_allData.(Location{iLoc}).(Sensors{iSen})(1))+1;
            Temp_LastMin(iLoc,iSen) = floor(Time_allData.(Location{iLoc}).(Sensors{iSen})(end))-1;
        end
    end

    for iLoc = 1:length(Location)
        for iSen = 1:length(Sensors)
            FirstMin = max(Temp_FirstMin);
            LastMin = min(Temp_LastMin);
            idx_FirstMin = find(floor(Time_allData.(Location{iLoc}).(Sensors{iSen})) == FirstMin(1));
            idx_LastMin = find(floor(Time_allData.(Location{iLoc}).(Sensors{iSen})) == LastMin(1));

            Watch_allData.(Location{iLoc}).(Sensors{iSen}) = Watch_allData.(Location{iLoc}).(Sensors{iSen})(idx_FirstMin(1):idx_LastMin(end),:);
            Time_allData.(Location{iLoc}).(Sensors{iSen}) = Time_allData.(Location{iLoc}).(Sensors{iSen})(idx_FirstMin(1):idx_LastMin(end),:);
            rawTime_allData.(Location{iLoc}).(Sensors{iSen}) = rawTime_allData.(Location{iLoc}).(Sensors{iSen})(idx_FirstMin(1):idx_LastMin(end),:);

            convTime.(Location{iLoc}).(Sensors{iSen}) = datetime(rawTime_allData.(Location{iLoc}).(Sensors{iSen}),'ConvertFrom','epochtime','Format','dd-MMM-yyyy HH:mm:ss.SSSSSSSSS');
        end
    end

    % Analyse
        acc_ankle = Watch_allData.(Location{idx_Ankle}).Wrist_Raw_Accelerometer;
        gyr_ankle = Watch_allData.(Location{idx_Ankle}).Wrist_Raw_Gyro;

        acc_wrist = Watch_allData.(Location{idx_Wrist}).Wrist_Raw_Accelerometer;
        gyr_wrist = Watch_allData.(Location{idx_Wrist}).Wrist_Raw_Gyro;


        %% Detect walking parts
        NpeaksToDetect = 2;
        RemoveThreePeaks = 0;
        showplot = 0;
        details = 0;
        [TaskSeg, TurnSeg, SegNOWalking, Number_of_Steps] = RPQ_ThreePeaksDetection(acc_ankle(:,1), gyr_ankle(:,1), Freq, NpeaksToDetect, RemoveThreePeaks, details);
        TurnSeg = round(TurnSeg);

        %% Data processing
        % Filter Data
        [b,a] = butter(2,2*[0.5 15]/Freq);
        Filtre_acc_ankle = filtfilt(b,a,acc_ankle);
        Filtre_gyr_ankle = filtfilt(b,a,gyr_ankle);

        Filtre_acc_wrist = filtfilt(b,a,acc_wrist);
        Filtre_gyr_wrist = filtfilt(b,a,gyr_wrist);

        % for iaxe = 1:3
        %     figure(1); subplot(3,1,iaxe); plot((acc_ankle(:,iaxe)))
        %     figure(2); subplot(3,1,iaxe); plot((Filtre_acc_ankle(:,iaxe)))
        %     % figure(2); subplot(3,1,iaxe); plot(Filtre_gyr_ankle(:,iaxe))
        % end

        % Apply CWT to acc. data
        axe = {'x','y','z'};
        for iaxe = 1:length(axe)
            [coefcwt, freq] = cwt(Filtre_acc_ankle(:,iaxe),'amor',Freq);
            assignin('base', ['PSD_' axe{iaxe} '_ankle'], abs(coefcwt).^2)
            assignin('base', ['f' axe{iaxe} '_ankle'], freq)

            [coefcwt, freq] = cwt(Filtre_acc_wrist(:,iaxe),'amor',Freq);
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
        [PTF, Freezing_Index, Freezing_mask, FreezingSeg] = RPQ_FOG_Detection(acc_ankle(:,1), PowerFoG_ankle(:,1), PowerMVT_ankle(:,1), TaskSeg, showplot);

        %% Detect Tremor
        showplot = 0;
        [PTT_Walking, PTT_NOWalking, Tremor, Tremor_mask, TremorSeg] = RPQ_Tremor_Detection(acc_wrist, PowerTremor_wrist, sum(PowerMVT_wrist,2), PeakPower_Freq_wrist, TaskSeg, SegNOWalking, showplot);

        %% Detect Bradykinesia
        showplot = 0;
        [Brady] = RPQ_Brady_Detection(acc_wrist, sum(PowerBrady_wrist,2), sum(PowerMVT_wrist,2), TaskSeg, SegNOWalking, showplot);

        %% Detect Dyskinesia
        showplot = 0;
        [Dysk] = RPQ_Dysk_Detection(acc_wrist, sum(PowerMVT_wrist,2), Tremor_mask, TaskSeg, SegNOWalking, showplot);




        %% Plot Results
        % Filter Data for plot
        [b,a] = butter(2,2*[0.5 3]/Freq);
        Filtre_acc_ankle_plot = filtfilt(b,a,acc_ankle);
        [b,a] = butter(2,2*[0.5]/Freq,'low');
        Filtre_gyr_ankle_plot = filtfilt(b,a,gyr_ankle);

        Filtre_acc_ankle_plot = Filtre_acc_ankle_plot(:,1)./max(abs(Filtre_acc_ankle_plot(:,1)));
        Filtre_gyr_ankle_plot = Filtre_gyr_ankle_plot(:,1)./max(abs(Filtre_gyr_ankle_plot(:,1)));

        % % Plot symptoms
        % Sympt = [Freezing_Index, Tremor, Brady, Dysk];
        % Symptoms = {'FoG', 'Tremor', 'Bradykinesia', 'Dyskinesia'};
        % figure;
        % for iSympt = 1:size(Sympt,2)
        %     subplot(size(Sympt,2),1,iSympt);
        %     plot(Sympt(:,iSympt))
        %     title(Symptoms{iSympt})
        % end


        % Plot recap
        figure;
        subplot(3,1,1); plot(convTime.WristHomeXsens.Wrist_Raw_Accelerometer,Dysk, '.', 'color', '#FFD700');
        legend('Dysk.')
        legend('boxoff')
        title(['Participant ' num2str(Participant) ' - ' Data])

        subplot(3,1,2); plot(convTime.WristHomeXsens.Wrist_Raw_Accelerometer,Brady, '.', 'color', '#FF4500');
        legend('Brady.')
        legend('boxoff')
        h = gca;
        set( h, 'YDir', 'reverse' )

        subplot(3,1,3); plot(convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer,Filtre_acc_ankle_plot);
        hold on
        plot(convTime.AnkleHomeXsens.Wrist_Raw_Gyro,Filtre_gyr_ankle_plot);
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
        fill([convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(TaskSeg(1,:))' flip(convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(TaskSeg(1,:))',2)], [-1 -1 1 1], 'w', 'facecolor', cSegActivity, 'linestyle', 'none', 'facealpha', 0.4);
        fill([convTime.AnkleHomeXsens.Wrist_Raw_Gyro(TurnSeg(1,:))' flip(convTime.AnkleHomeXsens.Wrist_Raw_Gyro(TurnSeg(1,:))',2)], [-1 -1 1 1], 'w', 'facecolor', cSegTurning, 'linestyle', 'none', 'facealpha', 0.4);

        if size(FreezingSeg,1) >= 1
            fill([convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(FreezingSeg(1,:))' flip(convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(FreezingSeg(1,:))',2)], [-1 -1 1 1], 'black', 'linestyle', 'none', 'facealpha', 0.7);
            FOGLegend = 'FoG';
        else
            FOGLegend = '';
        end
        if size(TremorSeg,1) >= 1
            fill([convTime.WristHomeXsens.Wrist_Raw_Accelerometer(TremorSeg(1,:))' flip(convTime.WristHomeXsens.Wrist_Raw_Accelerometer(TremorSeg(1,:))',2)], [-1 -1 1 1], 'w', 'facecolor', cSegTremor, 'linestyle', 'none', 'facealpha', 0.4);
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
                fill([convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(TaskSeg(iTaskSeg,:))' flip(convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(TaskSeg(iTaskSeg,:))',2)], [-1 -1 1 1], 'w', 'facecolor', cSegActivity, 'linestyle', 'none', 'facealpha', 0.4);
            end
        end
        if size(TurnSeg,1) > 1
            for iTurnSeg = 2:length(TurnSeg)
                fill([convTime.AnkleHomeXsens.Wrist_Raw_Gyro(TurnSeg(iTurnSeg,:))' flip(convTime.AnkleHomeXsens.Wrist_Raw_Gyro(TurnSeg(iTurnSeg,:))',2)], [-1 -1 1 1], 'w', 'facecolor', cSegTurning, 'linestyle', 'none', 'facealpha', 0.4);
            end
        end
        if size(FreezingSeg,1) > 1
            for iFreezingSeg = 2:length(FreezingSeg)
                fill([convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(FreezingSeg(iFreezingSeg,:))' flip(convTime.AnkleHomeXsens.Wrist_Raw_Accelerometer(FreezingSeg(iFreezingSeg,:))',2)], [-1 -1 1 1], 'black', 'linestyle', 'none', 'facealpha', 0.7);
            end
        end
        if size(TremorSeg,1) > 1
            for iTremorSeg = 2:length(TremorSeg)
                fill([convTime.WristHomeXsens.Wrist_Raw_Accelerometer(TremorSeg(iTremorSeg,:))' flip(convTime.WristHomeXsens.Wrist_Raw_Accelerometer(TremorSeg(iTremorSeg,:))',2)], [-1 -1 1 1], 'w', 'facecolor', cSegTremor, 'linestyle', 'none', 'facealpha', 0.4);
            end
        end

        pause
        % clearvars -except Task

    end




    % Load and save PTF
    % Percent_Time_Frozen = [NaN NaN NaN];
    % cd('F:\Projet RPQ\Extracted_Data')
    % load('TUG_Percent_Time_Frozen.mat')
    % Percent_Time_Frozen(Participant,:) = PTF;
    % 
    % cd('F:\Projet RPQ\Extracted_Data')
    % save('TUG_Percent_Time_Frozen.mat','Percent_Time_Frozen')
    % 
    % 
    % % Load and save PTF
    % Percent_Time_Tremor_Walking = [NaN NaN NaN];
    % Percent_Time_Tremor_NOWalking = [NaN NaN NaN];
    % Tremor_Level = [NaN NaN NaN];
    % 
    % cd('F:\Projet RPQ\Extracted_Data')
    % load('TUG_Percent_Time_Tremor_Walking.mat')
    % load('TUG_Percent_Time_Tremor_NOWalking.mat')
    % load('TUG_Tremor_Level.mat')
    % load('TUG_NO_Tremor_Level.mat')
    % Percent_Time_Tremor_Walking(Participant,:) = mean(PTT_Walking);
    % Percent_Time_Tremor_NOWalking(Participant,:) = mean(PTT_NOWalking);
    % Tremor_Level(Participant,:) = mean(Tremor);

    % cd('F:\Projet RPQ\Extracted_Data')
    % save('TUG_Percent_Time_Tremor_Walking.mat','Percent_Time_Tremor_Walking')
    % save('TUG_Percent_Time_Tremor_NOWalking.mat','Percent_Time_Tremor_NOWalking')
    % save('TUG_Tremor_Level.mat','Tremor_Level')

