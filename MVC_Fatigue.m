close all ;
clear all ;
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

Files = {...
    'biceps_1';'biceps_2';'deltant_1';'deltant_2';'deltmed_1';'deltmed_2';'deltpost_1';...
    'deltpost_2';'isp_1';'isp_2';'pect_1';'pect_2';'ssp_1';'ssp_2';'subs_1';'subs_2';'triceps_1';'triceps_2';'uptrap_1';'uptrap_2'} ; 

Muscles = {...
    'deltant',  'deltant_EMG1';...
    'deltmed',  'deltmed_EMG2';...
    'deltpost',    'deltpost_EMG3';...
    'biceps',    'biceps_EMG4';...
    'triceps',  'triceps_EMG5';...
    'uptrap',  'uptrap_EMG6';...
    'pect',  'pect_IM_EMG7';...
    'ssp',  'ssp_EMG8';...
    'isp',  'isp_EMG9';...
    'subs',  'subs_EMG10';...
    } ;

for iSubjects = 1:2
    cd(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\mvc'])
%     cd(['H:\Projet_ExpertsNovices\data\raw\2017-11-22\marc\mvc'])
%     iSubjects = 21;
    
    EMG_allMVC = [];
    
    %Loading EMG data
    for iMuscles = 1:length(Muscles)   
        EMG_MVC = [] ;
             
        for iFiles = 1:length(Files)
            FileName = ([Files{iFiles,1} '.c3d']) ;
            acq = btkReadAcquisition(FileName) ;
            DataMVC = btkGetAnalogs(acq) ;
            Freq = btkGetAnalogFrequency(acq) ; %Fréquence d'échantillonage de l'EMG
            EMG_MVC = [EMG_MVC ; DataMVC.(Muscles{iMuscles,2})] ;
        end
        
        EMG_allMVC = [EMG_allMVC EMG_MVC];
%         subplot(3,4,iMuscles) ; plot(EMG_allMVC(:,iMuscles)) ; title (Muscles{iMuscles}, 'interpreter', 'none') ;
    end
    
    % Preprocessing
    [b,a] = butter(2,2*[10 400]/Freq) ; % Parametre du filtre BP (Band Pass) 10-400 Hz
    EMGBP = filtfilt(b,a,EMG_allMVC) ;
    
    [b,a] = butter(2,2*[55 65]/Freq,'stop') ;  % Parametre du filtre BS (Band Stop) 55-65 Hz (prises électriques)
    EMGBS = filtfilt(b,a,EMGBP) ;
    EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)
    
%     for iMuscles = 1:length(Muscles)
%         subplot(3,4,iMuscles) ; plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%     end
    
    % Data analysis
    % EMG activation level
    [b,a] = butter(2,2*9/Freq) ; % Parametre du filtre Low pass 9 Hz
    EMG_MVCenvelop = filtfilt(b,a,abs(EMGBL)) ; %abs : on ne garde que les valeurs absolues de EMGBL donc signal redressé au-dessus de 0
    
%     for iMuscles = 1:length(Muscles)
%         subplot(3,4,iMuscles) ; plot(EMG_MVCenvelop(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%     end
    
    %Sort EMG values from highest to lowest (descending order)
    EMG_sort = sort(EMG_MVCenvelop,'descend') ;
    
    %Mean EMG over 2000 Hz for each muscle
    EMG_sort = EMG_sort(1:Freq,:) ;
    EMG_max = mean(EMG_sort,1);
    
    % Data saving
    save(['CleanData_MVC_' (Subjects{iSubjects}) '.mat'],'EMG_max')
    clear functions
end
