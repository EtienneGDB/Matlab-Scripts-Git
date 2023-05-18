clear ;
close all ;
clc ;

% Participant to export
Participant = 9;
if Participant < 10
    Participant_Name = sprintf('P0%s',num2str(Participant));
else
    Participant_Name = sprintf('P%s',num2str(Participant));    
end

% Paths
projectPath = 'F:\Projet RPQ';
dataPath    = ['F:\Projet RPQ\RawData\',Participant_Name,'\',Participant_Name,'_Labo'];
saveDir     = 'F:\Projet RPQ\ExportedData';
libsPath = fullfile(projectPath, 'libs');
addpath(genpath(libsPath));
C3DserverPath   = 'C:\Users\Public\Documents\Motion Lab Systems\C3Dserver\';
addpath(genpath(C3DserverPath));
addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Codes_Analyse_Xsens_EMG
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))

XSENS = false;
OPTI = false;
WATCH = true;


%% Load Excel File
excelfile = dir(fullfile(dataPath, '*.xlsx'));
if isempty(excelfile)
    display('SVP vous assurez que la liste des essais est compilÃ©e dans un fichier excel!')
    return
elseif size(excelfile,1)>=1
    for i=size(excelfile,1)
        tst = excelfile(i).name;
        if isempty(strfind(tst,'~$'))
            excelfile = excelfile(i).name;
            break
        end
    end
end
file_partinfo = fullfile(dataPath, excelfile);
% [num, txt, raw] = xlsread(file_partinfo,'Feuille1');
% [num, txt, raw] = xlsread(file_partinfo,'Sheet1');
[num, txt, raw] = xlsread(file_partinfo);
keep_raw = raw; keep_txt = txt; keep_num = num;
for iL = 1:size(keep_raw,1)
    if keep_raw{iL,1}(end) == '1' | keep_raw{iL,1}(end) == '2'
        keep_raw(iL,:) = [];
        keep_txt(iL,:) = [];
        keep_num(iL-1,:) = [];
    end
    if iL >= length(keep_raw)
        break
    end
end
cd([saveDir '\FileNames'])
save(['P', num2str(Participant), '_FileNames', '.mat'], 'keep_raw')

%% Load and save OptiTrack Data into .mat
if OPTI
    optiFolder = fullfile(dataPath, 'Optitrack');
    fileList = dir([optiFolder, '\*', '.c3d']);
    Filenames = {fileList(1:end).name};
    while Filenames{1}(1) == '.'
        Filenames(1) = [];
    end
    
    % Options
    OPTITRACK_FORMAT = 1; %Tu devrais avoir toujours du C3D
    % 1 -> C3D
    % 2 -> CSV
    
    % CONSTANTS
    d2r    = pi/180;    %Degrees to radians
    r2d    = 180/pi;    %Radians to degrees
    
    % Check if C3D server is installed
    if OPTITRACK_FORMAT == 1
        if ~isempty(which('c3dserver'))
            opti = c3dserver;
        else
            %You must install C3D server from here before loading OptiTrack data: https://www.c3dserver.com/
            disp("Veuillez installer C3D server et ajoutez les fichiers dans le dossier libs");
            return;
        end
    end
    
    % Load optitrack data
    for iFile = 1:length(keep_num)
        switch OPTITRACK_FORMAT
            case 1 % Format C3D
                optiFile = [optiFolder, '\', Filenames{keep_num(iFile)}];
                if ~isempty(optiFile)
                    openc3d(opti, 0, optiFile);
                    rawOpti = get3dtargets(opti);
                end
                
            case 2 % Format CSV
                optiFile = [optiFolder, '\', Filenames{iFile}];
                optiFile = optiFile(1:end-4);
                optiFile = strcat(optiFile,'.csv');
                if ~isempty(optiFile)
                    rawOpti = parseOptiTrackCSV(optiFile);
                end
        end
        
%         figure
%         plot(rawOpti.LANK)
%         title('Marqueur Left Ankle... attention données erronées dans ce cas Droite/gauche')
        
        % Interpoler les donnees NaN de OptiTrack avec spline. Sinon, problÃ¨mes avec la fonction de filtrage filtfilt
        % Spline semble Ãªtre la meilleure methode generale pour les 3 axes. Validee avec une analyse visuelle.
        interp = @(data) fillmissing(data, 'spline'); %à évaluer pour nous...
        
        % Si format C3D, retirer l'unite puisque ce n'est pas un tableau, et n'est pas trÃ¨s utile de toute faÃ§on.
        if OPTITRACK_FORMAT == 1; rawOpti.units = []; end   %Tu devrais aovir des c3d
        
        % Interpoler les NaN OptiTrack
        rawOpti = structfun(interp, rawOpti, 'UniformOutput',false);
        
        clear interp;
        
%         figure
%         plot(rawOpti.LANK)
%         title('Marqueur Left Ankle... attention données erronées dans ce cas Droite/gauche')

        
        % Synchro des signaux
        %Voir fonction ***fctAlignTemporal*** si peut être utilisée...
        delay = 0; %Assume 0 juste parce que je voulais te donner le code après...
        
        %Aligner les signaux
        rawOpti_al = rawOpti;
        if delay < 0
            %Retirer toutes ls donnees de debut des signaux qui ne sont pas
            %utiles, pour que les donnees des 2 systÃ¨mes soient synchronisees
            rawOpti_al.RHEE = rawOpti.RHEE(-delay+1:end,:);
            rawOpti_al.RTOE = rawOpti.RTOE(-delay+1:end,:);
            rawOpti_al.RANK = rawOpti.RANK(-delay+1:end,:);
            rawOpti_al.RPSI = rawOpti.RPSI(-delay+1:end,:);
            
            aank_al = aank(-delay+1:end,:);
            aank_filt = aank_filt(-delay:end, :);
            
            
        elseif delay > 0
            acc_ankle = acc_ankle(delay+1:end,:);
            gyr_ankle = gyr_ankle(delay+1:end,:);
            mag_ankle = mag_ankle(delay+1:end,:);
            pressure_pp100 = pressure_pp100(delay+1:end,:);
            time = time(1:end-delay,:);
            qSensoria_I = qSensoria_I(delay+1:end,:);
            qIMU_ankle = qIMU_ankle(delay+1:end,:);
            
            acc_ankle_filt = acc_ankle_filt(delay+1:end,:);
        end
        cd([saveDir, '\Opti'])
        save(['P', num2str(Participant), '_Opti_', keep_raw{(iFile)+1,1}, '.mat'], 'rawOpti')
    end
end


%% Load and save XSens Data into .mat
if XSENS
    segmentData_of_Interest = {'orientation', 'angularVelocity', 'acceleration', 'position'};
    jointData_of_Interest = {'jointAngle', 'jointAngleXZY'};

    xsensFolder = fullfile(dataPath, 'Xsens');
    fileList = dir([xsensFolder, '\*', '.mvnx']);
    Filenames = {fileList(1:end).name};
    
    for iFile = 1:length(keep_num)
        mvnx_file = load_mvnx([xsensFolder, '\', (Filenames{keep_num(iFile)})]);
        nSeg = length(mvnx_file.segmentData);
        nJoint = length(mvnx_file.jointData);
        
        % Segment
        for iSegment = 1:nSeg
            tempLabel = mvnx_file.segmentData(iSegment).label;
            
            for isegmentData_of_Interest = 1:length(segmentData_of_Interest)
                DataXsens.(segmentData_of_Interest{isegmentData_of_Interest}).(tempLabel) = ...
                    mvnx_file.segmentData(iSegment).(segmentData_of_Interest{isegmentData_of_Interest});
            end
            
        end
        
        % JointAngle
        for iJoint = 1:nJoint
            tempLabel = mvnx_file.jointData(iJoint).label;
            
            for ijointData_of_Interest = 1:length(jointData_of_Interest)
                DataXsens.(jointData_of_Interest{ijointData_of_Interest}).(tempLabel) = ...
                    mvnx_file.jointData(iJoint).(jointData_of_Interest{ijointData_of_Interest});
            end
            
        end

        cd([saveDir, '\Xsens'])
        save(['P', num2str(Participant), '_Xsens_', keep_raw{(iFile)+1,1}, '.mat'], 'DataXsens')
    end
end

%% Load and save Watch Data into .mat
if WATCH
    Freq = 50;
    Type = {'Wrist', 'Ankle'};
    FileNames = {'Wrist_Raw_Accelerometer', 'Wrist_Raw_Gyro'};
    Watch_allData = [];
    Time_allData = [];
    for iType = 1:length(Type)
        for iFileNames = 1:length(FileNames)
            watchFolder = fullfile([saveDir, '\Watch\RPQProject\P', num2str(Participant), '\Recordsets\', Type{iType},'Labo']);
            Data = load([watchFolder, '\', FileNames{iFileNames}, '.mat']);
            Field_Data = fieldnames(Data);
            WatchData_of_Interest = Data.(Field_Data{1}).values(:,[2,4,6]);
            unixTime = Data.(Field_Data{1}).values(:,1);
            convTime = datetime(unixTime,'ConvertFrom','epochtime','Format','dd-MMM-yyyy HH:mm:ss.SSSSSSSSS');
            Hour = hour(convTime);
            Hour = (Hour-Hour(1))*60;
            Minutes = minute(convTime)+Hour;
            Minutes = ((Minutes-Minutes(1))*60);
            Time = second(convTime)+Minutes;
%             for iaxe = 1:3
%                 figure; plot(WatchData_of_Interest(:,iaxe))
%             end
            if iFileNames == 1
                % Filter gravity
                [b,a] = butter(2,2*0.5/50,'high');
                FiltreGravity = filtfilt(b,a,WatchData_of_Interest);
                Watch_allData.(Type{iType}).(FileNames{iFileNames}) = FiltreGravity;
                Time_allData.(Type{iType}).(FileNames{iFileNames}) = Time;
            else
                Watch_allData.(Type{iType}).(FileNames{iFileNames}) = WatchData_of_Interest;
                Time_allData.(Type{iType}).(FileNames{iFileNames}) = Time;
            end
            figure; plot(Time, sum(abs(Watch_allData.(Type{iType}).(FileNames{iFileNames})),2))
            title([Type{iType} ' ' FileNames{iFileNames}])
            
%             for iaxe = 1:3
%                 figure; plot(Time, Watch_allData.(Type{iType}).(FileNames{iFileNames})(:,iaxe))
%                 title([Type{iType} ' ' FileNames{iFileNames}])
%                 pause;
%                 close all;
%             end
            
        end
    end
    
    % Synchronize signaux
    BestSensor = 'Wrist_Raw_Accelerometer';
    for iType = 1:2
        Rectified_Sum = sum(abs(Watch_allData.(Type{iType}).(BestSensor)),2);
        [PKS,LOCS] = findpeaks(Rectified_Sum,'MinPeakDistance',250,'SortStr','descend','NPeaks',1);
        figure;
        plot(Rectified_Sum);
        hold on;
        plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
        
        LOCS_Peak(iType) = LOCS;
    end
    LOCS_Peak(1)=LOCS_Selec1(1);
    LOCS_Peak(2)=LOCS_Selec2(1);
    
    T1 = Time_allData.(Type{1}).Wrist_Raw_Accelerometer(LOCS_Peak(1),:);
    T2 = Time_allData.(Type{2}).Wrist_Raw_Accelerometer(LOCS_Peak(2),:);
    
    Watch_allData_Modif = Watch_allData;
    Time_allData_Modif = Time_allData;
    if LOCS_Peak(2) > LOCS_Peak(1)
%         Watch_allData_Modif.(Type{2}).Wrist_Raw_Accelerometer = Watch_allData.(Type{2}).Wrist_Raw_Accelerometer(diff(LOCS_Peak)+1:end,:);
%         Watch_allData_Modif.(Type{2}).Wrist_Raw_Gyro = Watch_allData.(Type{2}).Wrist_Raw_Gyro(diff(LOCS_Peak)+1:end,:);
        Time_allData_Modif.(Type{1}).Wrist_Raw_Accelerometer = ...
            Time_allData.(Type{1}).Wrist_Raw_Accelerometer + abs(T2-T1);
        Time_allData_Modif.(Type{1}).Wrist_Raw_Gyro = ...
            Time_allData.(Type{1}).Wrist_Raw_Gyro + abs(T2-T1);
    else
%         Watch_allData_Modif.(Type{1}).Wrist_Raw_Accelerometer = Watch_allData.(Type{1}).Wrist_Raw_Accelerometer(abs(diff(LOCS_Peak))+1:end,:);
%         Watch_allData_Modif.(Type{1}).Wrist_Raw_Gyro = Watch_allData.(Type{1}).Wrist_Raw_Gyro(abs(diff(LOCS_Peak))+1:end,:);
        Time_allData_Modif.(Type{2}).Wrist_Raw_Accelerometer = ...
            Time_allData.(Type{2}).Wrist_Raw_Accelerometer + abs(T2-T1);
        Time_allData_Modif.(Type{2}).Wrist_Raw_Gyro = ...
            Time_allData.(Type{2}).Wrist_Raw_Gyro + abs(T2-T1);
    end
    % check
    for iType = 1:length(Type)
        for iFileNames = 1:length(FileNames)
            figure;
            Rectified_Sum = sum(abs(Watch_allData_Modif.(Type{iType}).(FileNames{iFileNames})),2);
            Time = Time_allData_Modif.(Type{iType}).(FileNames{iFileNames});
            plot(Time, Rectified_Sum);
            title([Type{iType} ' ' FileNames{iFileNames}])
        end
    end

    % Find peaks feet tapping
    BestSensor = 'Wrist_Raw_Accelerometer';
    figure; plot(Watch_allData_Modif.(Type{2}).(BestSensor)(:,1)*9.81)
    Rectified_Sum = sum(abs(Watch_allData_Modif.(Type{2}).(BestSensor)),2);
%     Rectified_Sum = abs(Watch_allData.(Type{2}).(BestSensor)(:,3));
    [PKS,LOCS] = findpeaks(Rectified_Sum,'MinPeakDistance',250,'SortStr','descend','NPeaks',6);
    figure;
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
    LOCS = sort(LOCS);
%     %---Correct when necessary---
    % P2
%     PeakToDelete = [8];
%     LOCS(PeakToDelete) = [];
    LOCS(1) = 49200;
    LOCS(6) = 58390;
    LOCS(7) = 78490;
    LOCS(8) = 157500;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
    
    %P3
    PeakToDelete = [7; 8];
    LOCS(PeakToDelete) = [];
    LOCS(1) = 38210;
    LOCS(2) = 45380;
    LOCS(4) = 62920;
    LOCS(5) = 126100;
    LOCS(6) = 150600;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
    
    %P4
    PeakToDelete = [7; 8];
    LOCS(PeakToDelete) = [];
    LOCS(1) = 95350;
    LOCS(5) = 112500;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
    
    %P5
    PeakToDelete = [7; 8];
    LOCS(PeakToDelete) = [];
    LOCS(1) = 34560;
    LOCS(2) = 41910;
    LOCS(3) = 44010;
    LOCS(6) = 58740;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
    
    %P6
    LOCS(1) = 64190;
    LOCS(2) = 71730;
    LOCS(3) = 75490;
    LOCS(4) = 85310;
    LOCS(5) = 122000;
    LOCS(6) = 145600;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
    
    %P7
    LOCS(1) = 59370;
    LOCS(2) = 68100;
    LOCS(3) = 75950;
    LOCS(4) = 90770;
    LOCS(5) = 129100;
    LOCS(6) = 152400;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
    
    %P8
    LOCS(1) = 68700;
    LOCS(2) = 76100;
    LOCS(3) = 77190;
    LOCS(4) = 84560;
    LOCS(5) = 106900;
    LOCS(6) = 128800;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;

    %P9
    LOCS(1) = 57197;
    LOCS(2) = 60300;
    LOCS(3) = 70513;
    LOCS(4) = 75613;
    LOCS(5) = 94900;
    LOCS(6) = 117734;
    LOCS = sort(LOCS);
    plot(Rectified_Sum);
    hold on;
    plot(LOCS,PKS(1:length(LOCS)),'*','color','red')
    pause()
    close all;
%     %----------------------------

    % Define segments
    segments = [];
    inc = 1;
    % P2
    for iLOCS = 1:length(LOCS)-1
        if iLOCS == 1
            segments(1,inc) = LOCS(iLOCS)-200;
            segments(2,inc) = LOCS(iLOCS+2)+LOCS(iLOCS+1)-LOCS(iLOCS);
            PeakToDelete = [iLOCS+1; iLOCS+2];
            LOCS(PeakToDelete) = [];
        elseif iLOCS == 3
            segments(1,inc) = LOCS(iLOCS)-200;
            segments(2,inc) = LOCS(iLOCS)+LOCS(iLOCS)-LOCS(iLOCS-1);
        elseif iLOCS == 4
            segments(1,inc) = LOCS(iLOCS)-200;
            segments(2,inc) = LOCS(iLOCS)+3.2*60*Freq;
        elseif iLOCS == 5
            segments(1,inc) = LOCS(iLOCS)-200;
            segments(2,inc) = LOCS(iLOCS)+4.2*60*Freq;
        elseif iLOCS == length(LOCS)
            segments(1,inc) = LOCS(iLOCS)-200;
            segments(2,inc) = min(length(Watch_allData_Modif.Wrist.Wrist_Raw_Gyro),length(Watch_allData_Modif.Ankle.Wrist_Raw_Gyro));
            break
        else
            segments(1,inc) = LOCS(iLOCS)-200;
            segments(2,inc) = LOCS(iLOCS+1)-200;
        end
        inc = inc+1;
    end
    
    % P3, P4, P5, P6, P7, P8, P9
    segments = [];
    inc = 1;
    for iLOCS = 1:length(LOCS)
        if iLOCS == 1
            segments(1,inc) = LOCS(iLOCS)-100;
            segments(2,inc) = LOCS(iLOCS)+40*Freq;
        elseif iLOCS == 2
            segments(1,inc) = LOCS(iLOCS)-200;
            segments(2,inc) = LOCS(iLOCS)+20*Freq;
        elseif iLOCS == 3
            segments(1,inc) = LOCS(iLOCS)-100;
            segments(2,inc) = LOCS(iLOCS)+11*Freq;
        elseif iLOCS == 4
            segments(1,inc) = LOCS(iLOCS)-100;
            segments(2,inc) = LOCS(iLOCS)+3.5*60*Freq;
        elseif iLOCS == 5
            segments(1,inc) = LOCS(iLOCS)-100;
            segments(2,inc) = LOCS(iLOCS)+4.5*60*Freq;
        elseif iLOCS == length(LOCS)
            segments(1,inc) = LOCS(iLOCS)-100;
            segments(2,inc) = LOCS(iLOCS)+5.5*60*Freq;
            break
        end
        inc = inc+1;
    end
    
    
%     segments(1,1) = LOCS(1)-200;
%     segments(2,1) = LOCS(3)+LOCS(2)-LOCS(1);
%     segments(1,2) = LOCS(3)+LOCS(2)-LOCS(1);
%     segments(2,2) = LOCS(4)-200;
%     segments(1,3) = LOCS(4)-200;
%     segments(2,3) = LOCS(5)-200;
%     segments(1,4) = LOCS(5)-200;
%     segments(2,4) = LOCS(6)-200;
%     segments(1,5) = LOCS(6)-200;
%     segments(2,5) = LOCS(7)-200;
%     segments(1,6) = LOCS(7)-200;
%     segments(2,6) = LOCS(8)-200;
%     segments(1,7) = LOCS(8)-200;
%     segments(2,7) = min(length(Watch_allData.Wrist.Wrist_Raw_Gyro),length(Watch_allData.Ankle.Wrist_Raw_Gyro));
    
    % Segment Data
    Time_Ref = Time_allData_Modif.Ankle.Wrist_Raw_Accelerometer;
    for iFileSegmented = 1:6
        inc = 1;
        for iType = 1:length(Type)
            for iFileNames = 1:length(FileNames)
                Time = Time_allData_Modif.(Type{iType}).(FileNames{iFileNames});
                % Target Values of Time_Ref
                N(1) = Time_Ref(segments(1,iFileSegmented));
                N(2) = Time_Ref(segments(2,iFileSegmented));
                % idx of closest values of actual Time compared to Time_Ref
                for iN = 1:2
                    dist = abs(Time - N(iN));
                    minDist = min(dist);
                    idx(iN) = min(find(dist == minDist));
%                     Time(idx)
                end
                
                Watch_Data.(Type{iType}).(FileNames{iFileNames}) = ...
                    Watch_allData_Modif.(Type{iType}).(FileNames{iFileNames})...
                    (idx(1):idx(2),:);
                Time_Watch_Data.(Type{iType}).(FileNames{iFileNames}) = ...
                    Time(idx(1):idx(2),:);
                
%                 subplot(2,2,inc); 
%                 plot(Watch_allData_Modif.(Type{iType}).(FileNames{iFileNames})...
%                     (idx(1):idx(2),:))
%                 title([(Type{iType}) ' ' (FileNames{iFileNames})])
%                 pause()
%                 inc = inc +1;
            end
        end
        
        % Interpolate for same length
        Sig_Ref = Watch_Data.Wrist.Wrist_Raw_Accelerometer;
        inc = 1;
        for iType = 1:length(Type)
            for iFileNames = 1:length(FileNames)
                Sig = Watch_Data.(Type{iType}).(FileNames{iFileNames});
                Normalization = interp1(1:length(Sig),Sig,linspace(1,length(Sig),length(Sig_Ref)));
                
                Watch_Data.(Type{iType}).(FileNames{iFileNames}) = Normalization;

                subplot(2,2,inc); 
                plot(Normalization)
                title([(Type{iType}) ' ' (FileNames{iFileNames})])
                inc = inc +1;
            end
        end
        pause()
        close all;

        cd([saveDir, '\Watch'])
        % Before P6
%         save(['P', num2str(Participant), '_Watch_', keep_raw{iFileSegmented+2,1}, '.mat'], 'Watch_Data')
        % From P6
        save(['P', num2str(Participant), '_Watch_', keep_raw{iFileSegmented+3,1}, '.mat'], 'Watch_Data')
    end
        
end


% %% Extract InfoMontres Data
% watch_acc = load('F:\Projet RPQ\Infos montres\InfoMontres\TestMontresAntoine\Recordsets\StatiqueMontrePlat\Wrist_Raw_Accelerometer.mat');
% watch_gyr = load('F:\Projet RPQ\Infos montres\InfoMontres\TestMontresAntoine\Recordsets\StatiqueMontrePlat\Wrist_Raw_Gyro.mat');
% Acc = watch_acc.RawAccelerometer.values(:,[2, 4, 6]);
% Gyr = watch_gyr.RawGyro.values(:,[2, 4, 6]);
% Freq = 50;
% 
% % Filter gravity
% [b,a] = butter(2,2*0.5/Freq,'high');
% Acc_static = filtfilt(b,a,Acc);
% Gyr_static = Gyr;
% 
% figure; plot(Acc_static)
% figure; plot(Gyr_static)
% 
% cd('F:\Projet RPQ\Infos montres')
% save('Acc_static.mat','Acc_static')
% save('Gyr_static.mat','Gyr_static')

