clear ;
close all ;
clc ;

% Participant to export
for Participant = 4:6
if Participant < 10
    Participant_Name = sprintf('P0%s',num2str(Participant));
else
    Participant_Name = sprintf('P%s',num2str(Participant));    
end

% Paths
projectPath = 'F:\Projet RPQ';
dataPath    = ['F:\Projet RPQ\RawData\',Participant_Name,'\',Participant_Name,'_Domicile'];
saveDir     = 'F:\Projet RPQ\ExportedData';
libsPath = fullfile(projectPath, 'libs');
addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Codes_Analyse_Xsens_EMG
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))

XSENS = false;
WATCH = true;

%% Load and save XSens Data into .mat
if XSENS
    segmentData_of_Interest = {'orientation', 'angularVelocity', 'acceleration', 'position'};
    jointData_of_Interest = {'jointAngle', 'jointAngleXZY'};

    xsensFolder = fullfile(dataPath, 'Xsens');
    fileList = dir([xsensFolder, '\*', '.mvnx']);
    Filenames = {fileList(1:end).name};
    
    for iFile = 1:length(keep_num)
        mvnx_file = load_mvnx([xsensFolder, '\', (Filenames{iFile})]);
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
        save(['P', num2str(Participant), '_Xsens_Home', num2str(iFile), '.mat'], 'DataXsens')
    end
end

%% Load and save Watch Data into .mat
if WATCH
    Freq = 50;
    Type = {'Wrist', 'Ankle'};
    FileNames = {'Wrist_Raw_Accelerometer', 'Wrist_Raw_Gyro'};
    Watch_allData = [];
    Time_allData = [];
    rawTime_allData = [];
    watchFolder = dir(['F:\Projet RPQ\ExportedData\Watch\RPQProjectHomeXsens\P', num2str(Participant) ,'\Recordsets\*H*']);
    for iwatchFolder = 1:length(watchFolder)
        for iFileNames = 1:length(FileNames)
            Folder = fullfile([saveDir, '\Watch\RPQProjectHomeXsens\P', num2str(Participant), '\Recordsets\', watchFolder(iwatchFolder).name]);
            Data = load([Folder, '\', FileNames{iFileNames}, '.mat']);
            Field_Data = fieldnames(Data);
            WatchData_of_Interest = Data.(Field_Data{1}).values(:,[2,4,6]);
            unixTime = Data.(Field_Data{1}).values(:,1);
            convTime = datetime(unixTime,'ConvertFrom','epochtime','Format','dd-MMM-yyyy HH:mm:ss.SSSSSSSSS');
            Day = day(convTime);
            Day = (Day-Day(1))*24;
            Hour = hour(convTime)+Day;
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
                Watch_allData.(watchFolder(iwatchFolder).name).(FileNames{iFileNames}) = FiltreGravity;
                Time_allData.(watchFolder(iwatchFolder).name).(FileNames{iFileNames}) = Time;
                rawTime_allData.(watchFolder(iwatchFolder).name).(FileNames{iFileNames}) = unixTime;
            else
                Watch_allData.(watchFolder(iwatchFolder).name).(FileNames{iFileNames}) = WatchData_of_Interest;
                Time_allData.(watchFolder(iwatchFolder).name).(FileNames{iFileNames}) = Time;
                rawTime_allData.(watchFolder(iwatchFolder).name).(FileNames{iFileNames}) = unixTime;
            end
%             figure; plot(Time, sum(abs(Watch_allData.(watchFolder(iwatchFolder).name).(FileNames{iFileNames})),2))
%             title([(watchFolder(iwatchFolder).name) ' ' FileNames{iFileNames}])
            
%             for iaxe = 1:3
%                 figure; plot(Time, Watch_allData.(Type{iType}).(FileNames{iFileNames})(:,iaxe))
%                 title([Type{iType} ' ' FileNames{iFileNames}])
%                 pause;
%                 close all;
%             end
            
        end
    end
    
    cd([saveDir, '\Home\WatchXsens'])
    save(['P', num2str(Participant), '_WatchXsens_Home', '.mat'], 'Watch_allData')
    save(['P', num2str(Participant), '_TimeWatchXsens_Home', '.mat'], 'Time_allData')
    save(['P', num2str(Participant), '_rawTimeWatchXsens_Home', '.mat'], 'rawTime_allData')
        
end
end
