close all ;
clear all ;
clc ;

addpath X:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('X:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('X:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

Muscles = {'Dent1_IM_EMG1';'Dent2_IM_EMG2';'DeltA_IM_EMG5';'DeltM_IM_EMG6';'DeltP_IM_EMG7'; ...
    'Bi_IM_EMG11';'Tri_IM_EMG12';'TrapInf_IM_EMG10';'TrapSup_IM_EMG8';'TrapMed_IM_EMG9';};
 
FileName = ({'MVCpostfat.c3d', 'MVCprefat1.c3d', 'MVCprefat2.c3d'}) ;

for iSubjects = 1:24
    cd(['J:\IRSST_Fatigue\Pointage_repetitif\Data\P' num2str(iSubjects) '\Trial'])
             
        for iFile = 1:3
            MVC = [];
            acq = btkReadAcquisition(FileName{iFile}) ;
            DataMVC = btkGetAnalogs(acq) ;
            for iM = 1:length(Muscles)
                MVC = [MVC DataMVC.(Muscles{iM})];
            end
            
        end
        
        for iMuscles = 1:length(Muscles)   
        EMG_MVC = [] ;
        EMG_MVC = [EMG_MVC ; DataMVC.(Muscles{iMuscles,2})] ;

        EMG_allMVC = [EMG_allMVC EMG_MVC];
        subplot(3,4,iMuscles) ; plot(EMG_allMVC(:,iMuscles)) ; title (Muscles{iMuscles}, 'interpreter', 'none') ;
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
%     plot(EMG_max)
    
    % Data saving
    cd(['H:\Bureau\Etienne\Extracted data\Fatigue\MVC'])
    save(['CleanData_MVC_' (Subjects{iSubjects}) '.mat'],'EMG_max')
    clear functions
end
