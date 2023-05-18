clear ;
close all ;
clc ;

addpath \\10.89.24.15\e\Bureau\Etienne\MATLAB\Functions
addpath(genpath('\\10.89.24.15\e\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('\\10.89.24.15\e\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

FoldersNames = dir('C:\Users\p1098713\Documents\test_velo_bras');
Subjects = {...
    FoldersNames(3:end).name...
    } ;

% Define muscles used
% Column 1 = muscle name, colonne 2 = location in EMG matrix
%warning('vérifier que le numéro correspond bien pour le pec')
Muscles = {...
    'biceps',    'Sensor_1_IM_EMG1';...
    'triceps',  'Sensor_2_IM_EMG2';...
    'deltant',  'Sensor_3_IM_EMG3';...
    'deltmed',  'Sensor_4_IM_EMG4';...
    'deltpost',    'Sensor_5_IM_EMG5';...
    'uptrap',  'Sensor_6_IM_EMG6';...
    'dorsal',  'Sensor_7_IM_EMG7';...
    } ;

FileNames = {...
    FoldersNames(3:length(FoldersNames)).name...
    } ;
for iFileNames = 1:length(FileNames)
    while strcmp(FileNames{iFileNames}(length(FileNames{iFileNames})-2:length(FileNames{iFileNames})), 'c3d') == 0
        FileNames(iFileNames) = [];
    end
    if iFileNames == length(FileNames)
        break
    end
end
cd('C:\Users\p1098713\Documents\test_velo_bras');
acq = btkReadAcquisition(FileNames{1,1}) ;
Data = btkGetAnalogs(acq) ;
Freq = btkGetAnalogFrequency(acq) ;

DataSelec = {};
for iM = 1:length(Muscles)
    DataSelec.(Muscles{iM,1}) = Data.(Muscles{iM,2});
end

cd('C:\Users\p1098713\Documents\test_velo_bras');
save(['RawEMG_Muscles_' (FileNames{2,1}) '.mat'],'DataSelec')
% save(['RawEMG_Muscles_Test_liveFabien_Biceps03.mat'],'DataSelec')
