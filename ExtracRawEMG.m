clear ;
close all ;
clc ;

addpath E:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

addpath H:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

FoldersNames = dir('F:\Data\IRSST\RAW');
Subjects = {...
    FoldersNames(4:36).name...
   } ;
Subjects(5) = [];...Supprime boiteg (essai à part)

% Define muscles used
% Column 1 = muscle name, colonne 2 = location in EMG matrix
%warning('vérifier que le numéro correspond bien pour le pec')
Muscles = {...
    'deltant',  'deltant_EMG1';...
    'deltmed',  'deltmed_EMG2';...
    'deltpost',    'deltpost_EMG3';...
    'biceps',    'biceps_EMG4';...
    'triceps',  'triceps_EMG5';...
    'uptrap',  'uptrap_EMG6';...
    'pect',  'pect_IM_EMG7';...
    'ssp',  'ssp_EMG8';... supra epineux
    'isp',  'isp_EMG9';... infra ep.
    'subs',       'subs_EMG10';... subscapulaire
    } ;

iSubjects=29;
for iSubjects = 28:length(Subjects)
    %Load MVC Data
    cd(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\mvc'])
    load(['CleanData_MVC_' (Subjects{iSubjects}) '.mat']);
% 
%     FolderContent = dir(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\fatigue']);
%     FileNames = {...
%     FolderContent(3:length(FolderContent)).name...
%    } ;

    % Les participants sur disque E
    FolderContent = dir(['X:\Projet_ExpertsNovices\data\raw\2017-09-28\samc\fatigue']);
    FileNames = {...
    FolderContent(3:length(FolderContent)).name...
   } ;
    for iFileNames = 1:length(FileNames)
        while strcmp(FileNames{iFileNames}(length(FileNames{iFileNames})-2:length(FileNames{iFileNames})), 'c3d') == 0
            FileNames(iFileNames) = [];
        end
        if iFileNames == length(FileNames)
            break
        end
    end

    while (str2num(FileNames{1}(1)) == 1) | (str2num(FileNames{1}(1)) == 2)
        FileNames(1) = [];
    end
    
    %Crée une deuxième ligne avec le nom du fichier sans .c3d
    FileNames{2,1}=[];
    for iFileNames = 1:length(FileNames)
        FileNames{2,iFileNames} = erase(FileNames{1,iFileNames},'.c3d');
    end
      
    %iFiles=1;
    % Load EMG data
    for iFiles = 1:length(FileNames)
%         cd(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\fatigue'])
        cd(['H:\Projet_ExpertsNovices\data\raw\2017-09-29\aled\fatigue'])
        acq = btkReadAcquisition(FileNames{1,iFiles}) ;
        Data = btkGetAnalogs(acq) ;
        Freq = btkGetAnalogFrequency(acq) ;
        
        DataSelec = {};
        for iM = 1:length(Muscles)
            DataSelec.(Muscles{iM,1}) = Data.(Muscles{iM,2});
        end

        cd(['H:\Bureau\Etienne\Extracted data\Fatigue\DataSelec'])
        save(['RawEMG_Muscles_' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'DataSelec')
    end
    clear functions;
end