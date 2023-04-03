clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Sample entropy Fouaz'))

% addpath H:\Bureau\Etienne\MATLAB\Functions
% addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\btk'))
% addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

FoldersNames = dir('F:\Data\IRSST\RAW');
Subjects = {...
    FoldersNames(4:36).name...
   } ;
Subjects(5) = [];
Subjects(5) = [];

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

%iSubjects=1;
for iSubjects = 2:length(Subjects)
    %Load MVC Data
    cd(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\mvc'])
    load(['CleanData_MVC_' (Subjects{iSubjects}) '.mat']);

    % Plus de place sur serveur F donc les données ne sont pas au même
    % endroit pour tout le monde
    if ismember(iSubjects,[2:16 18:19 20:24 26:length(Subjects)])
        FolderContent = dir(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\fatigue']);
        FileNames = {...
        FolderContent(3:length(FolderContent)).name...
        } ;
    end
    if ismember(iSubjects,[1 17 19 22 25 28])
        if iSubjects == 1
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-09-29\aled\fatigue'];
        end
        if iSubjects == 17
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-11-27\jono\fatigue'];
        end
        if iSubjects == 19
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-11-22\marc\fatigue'];
        end
        if iSubjects == 22
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-11-10\nicl\fatigue'];
        end        
        if iSubjects == 25
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-12-12\pasd\fatigue'];
        end
        if iSubjects == 28
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-09-28\samc\fatigue'];
        end
        FolderContent = dir(pathData);
        FileNames = {...
        FolderContent(3:length(FolderContent)).name...
        } ;
        for iFileNames = 1:length(FileNames)
            if iFileNames >= length(FileNames)
                while strcmp(FileNames{length(FileNames)}(length(FileNames{length(FileNames)})-2:length(FileNames{length(FileNames)})), 'c3d') == 0
                    FileNames(length(FileNames)) = [];
                end
                break
            end
            while strcmp(FileNames{iFileNames}(length(FileNames{iFileNames})-2:length(FileNames{iFileNames})), 'c3d') == 0
                FileNames(iFileNames) = [];
                if iFileNames >= length(FileNames)
                    break
                end
            end
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
        cd(['X:\Bureau\Etienne\Extracted data\Fatigue\DataSelec'])
        load(['RawEMG_Muscles_' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat']);
        Data = DataSelec ;
        Freq = 2000 ;
        
        EMG = [] ;
        % EMG
        for iMuscles = 1:size(Muscles,1)
            if contains(fieldnames(Data),'_EMG')
                EMG(:,iMuscles) = Data.(Muscles{iMuscles,2}) ;
            else
                EMG(:,iMuscles) = Data.(Muscles{iMuscles,1}) ;
            end
%             subplot(3,3,iMuscles) ; plot(EMG(:,iMuscles)) ; title(Muscles{iMuscles,1}, 'interpreter', 'none') ;
        end

        % Preprocessing
        % Freq : fréquence d'échantillonage de l'EMG
        [b,a] = butter(2,2*[10 400]/Freq) ; % Parametre du filtre BP 10-400 Hz
        EMGBP = filtfilt(b,a,EMG) ;
        EMGBS = EMGBP;
          
        % Filtre les données à 60Hz et 120Hz pour EMG_intra
        Peak_Surface = [60];
        Peak_Intra = [60 120];
        for iM = 1:7 ...les 7 premiers sont EMG surface
            for iStop = Peak_Surface
                [b,a] = butter(2,2*[iStop-1 iStop+1 ]/Freq,'stop') ;  % Parametre du filtre BS 55-65 Hz
                EMGBS(:,iM) = filtfilt(b,a,EMGBS(:,iM)) ;
            end
        end
        for iM = 8:10 ...les 3 derniers sont EMG intra
            for iStop = Peak_Intra
                [b,a] = butter(2,2*[iStop-1 iStop+1 ]/Freq,'stop') ;  % Parametre du filtre BS 55-65 Hz
                EMGBS(:,iM) = filtfilt(b,a,EMGBS(:,iM)) ;
            end
        end
        
        % Remove baseline
        EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)
                
%         for iMuscles = 1:size(Muscles,1)
%             subplot(3,4,iMuscles) ; plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%         end
        
        %% Data analysis
        
        % EMG normalization
        Normalization_toto = [];
        for iNormalization = 1:length(Muscles)
%             Normalization (:,iNormalization) = EMGBL(:,iNormalization)./EMG_max(1,iNormalization)*100;
%             Normalization (:,iNormalization) = EMGBL(:,iNormalization)./max(EMGBL(:,iNormalization))*100;
            Normalization_toto(:,iNormalization) = toto(:,iNormalization)./max(toto(:,iNormalization))*100;
        end
%         figure ; plot(Normalization(:,1))
        
        % Dection of activity
        cd(['X:\Bureau\Etienne\Extracted data\Fatigue\Signal Segments'])
        load(['Seg_' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'])
        
        % Compute EMG envelop
        [b,a] = butter(2,2*9/Freq) ; % Parametre du filtre Low pass 9 Hz
        EMG_envelop = filtfilt(b,a,abs(Normalization)) ; % On ne garde que les valeurs absolues de EMGBL donc signal redressé au-dessus de 0
        
        % Compute Activity (feature)
        Activity = [];
        for iM = 1:length(Muscles)
            Activity.Trial.(Muscles{iM}) = var(EMG_envelop(:,iM));
            for iSeg = 1:length(Seg)
                Activity.Seg.(Muscles{iM})(:,iSeg) = var(EMG_envelop(Seg(iSeg,1):Seg(iSeg,2),iM));
            end     
        end
        
        % Compute Mobility (feature)
        Mobility = [];
        for iM = 1:length(Muscles)
            % Derive Signal
            axT = [1:length(EMG_envelop(:,iM))]/Freq;
            Derivative = diff(EMG_envelop(:,iM))./diff(axT(1:2));

            Mobility.Trial.(Muscles{iM}) = sqrt(var(Derivative)/var(EMG_envelop(:,iM)));
            
            for iSeg = 1:length(Seg)
                Mobility.Seg.(Muscles{iM})(:,iSeg) = var(Derivative(Seg(iSeg,1):Seg(iSeg,2)))/var(EMG_envelop(Seg(iSeg,1):Seg(iSeg,2),iM));
            end     
        end
        
%         % EMG activation level
%         [b,a] = butter(2,2*9/Freq) ; % Parametre du filtre Low pass 9 Hz
%         EMG_envelop = filtfilt(b,a,abs(Normalization)) ; % On ne garde que les valeurs absolues de EMGBL donc signal redressé au-dessus de 0
%         Amplitude = [];
%         for iM = 1:length(Muscles)
%             for iSeg = 1:length(Seg)
%                 Amplitude.Seg.(Muscles{iM})(:,iSeg) = mean(EMG_envelop(Seg(iSeg,1):Seg(iSeg,2),iM));
%             end
%         end
        
        % Compute Sample Entropy (feature)
        SampleEntropy = [];
        for iM = 1:length(Muscles)
            for iSeg = 1:length(Seg)
                SampleEntropy.Seg.(Muscles{iM})(:,iSeg) = sampenc(Normalization(Seg(iSeg,1):Seg(iSeg,2),iM),2,0.2); % erreur pour r -> vérif article Kartick
                SampleEntropy.Seg.(Muscles{iM})(:,iSeg) = sampenc(Normalization_toto(Seg(iSeg,1):Seg(iSeg,2),iM),2,0.2);
            end     
        end
        plot(Normalization_toto(:,10))

        cd(['H:\Bureau\Etienne\Extracted data\Fatigue\Activity'])
        save(['Activity' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'Activity')

        cd(['H:\Bureau\Etienne\Extracted data\Fatigue\Mobility'])
        save(['Mobility' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'Mobility')
        
%         cd(['H:\Bureau\Etienne\Extracted data\Fatigue\Amplitude'])
%         save(['Amplitude_' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'Amplitude')
        
%         cd(['H:\Bureau\Etienne\Extracted data\Fatigue\SampleEntropy'])
%         save(['SampleEntropy_' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'SampleEntropy')

    end
end