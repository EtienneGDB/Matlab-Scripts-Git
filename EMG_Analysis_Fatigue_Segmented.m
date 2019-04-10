clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))

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
%warning('v�rifier que le num�ro correspond bien pour le pec')
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
for iSubjects = 21:length(Subjects)
    %Load MVC Data
    cd(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\mvc'])
    load(['CleanData_MVC_' (Subjects{iSubjects}) '.mat']);

    % Plus de place sur serveur F donc les donn�es ne sont pas au m�me
    % endroit pour tout le monde
    if ismember(iSubjects,[2:16 18:19 21:length(Subjects)])
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
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-12-12\jono\fatigue'];
        end
        if iSubjects == 28
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-09-28\jono\fatigue'];
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
    
    %Cr�e une deuxi�me ligne avec le nom du fichier sans .c3d
    FileNames{2,1}=[];
    for iFileNames = 1:length(FileNames)
        FileNames{2,iFileNames} = erase(FileNames{1,iFileNames},'.c3d');
    end
    
    %iFiles=1;
    % Load EMG data
    for iFiles = 21:length(FileNames)
        cd(['H:\Bureau\Etienne\Extracted data\Fatigue\DataSelec'])
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
        % Freq : fr�quence d'�chantillonage de l'EMG
        [b,a] = butter(2,2*[10 400]/Freq) ; % Parametre du filtre BP 10-400 Hz
        EMGBP = filtfilt(b,a,EMG) ;
        EMGBS = EMGBP;
          
        % Filtre les donn�es � 60Hz et 120Hz pour EMG_intra
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
        EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS � chacune de ses valeurs)
        
%         for iMuscles = 1:size(Muscles,1)
%             subplot(3,4,iMuscles) ; plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%         end
        
%         % Plot FFT
%         T = 1/Freq;             % Sampling period       
%         L = length(EMGBL(:,iMuscles));             % Length of signal
%         t = (0:L-1)*T;        % Time vector
% 
%         figure
%         for iM=1:length(Muscles)
%             Y = fft(EMGBL(:,iM));
% 
%             P2 = abs(Y/L);
%             P1 = P2(1:L/2+1);
%             P1(2:end-1) = 2*P1(2:end-1);
% 
%             f = Freq*(0:(L/2))/L;
%             subplot(3,4,iM); plot(f,P1) 
%             title(['Part' (num2str(iSubjects)) '_' (Muscles{iM})])
%             xlabel('f (Hz)')
%             ylabel('|P1(f)|')
%         end
        
        %% Data analysis
        %% EMG activation level
%         [b,a] = butter(2,2*9/Freq) ; % Parametre du filtre Low pass 9 Hz
%         EMG_envelop = filtfilt(b,a,abs(EMGBL)) ; % On ne garde que les valeurs absolues de EMGBL donc signal redress� au-dessus de 0
%         plot(EMG_envelop)
        
        % EMG normalization
        Normalization = [];
        for iNormalization = 1:length(Muscles)
            Normalization (:,iNormalization) = EMGBL(:,iNormalization)./EMG_max(1,iNormalization)*100;
        end
%         figure ; plot(Normalization(:,:))
        
        % Sum normalized signals
        for iSig = 1:length(Normalization)
            SumNorm(iSig) = sum(Normalization(iSig,4));... � modifier en fonction
        end
%         figure ; plot(SumNorm)

        % Dection of activity
        Seg = ActivityDetection2(SumNorm, Normalization, 1000, 2000, Muscles, FileNames{2,iFiles}(1), 1);
        pause
        
        % Cut movements
        EMG_envelop_trials = [] ; EMG_cycles = [] ;
        for iMovements = 1:length(Seg)
            Movement = Normalization(Seg(iMovements,1):Seg(iMovements,2),:) ;
            EMG_cycles{1,iMovements} = EMGBL(Seg(iMovements,1):(Seg(iMovements,2)),:) ;
            
            % Interpolation
            EMG_envelop_trials(:,:,iMovements) = interp1(1:length(Movement),Movement,1:length(Movement)/500:length(Movement)) ;
        end
%         figure ; plot(mean(EMG_envelop_trials,3)) ; legend(Muscles{:,1}) ;
%         plot(EMG_cycles{1}(:,1))
%         plot(EMGBL(:,1))
        
        %% EMG Median frequency
        % TFR
        FreqMin = 1 ;
        FreqMax = 400 ;
        Resolution = 1 ;... pas de 1Hz
        WaveNumber = 7 ;
        Args = WaveletParameters(FreqMin,FreqMax,Resolution,WaveNumber,Freq) ;
        FreqRange = length(FreqMin:Resolution:FreqMax) ;
        Nb_Interp_Pnts = 1 ;
                
        % Compute TFR
        % Cut movements
        TFR = [] ;
        for iMovements = 1:length(Seg)
            for iM = 1:length(Muscles)
                disp('.')
                if sum(EMGBL(:,iM))~=0
                    tic
                    [norm, Time, Wave_FreqS] = TimeFreqTransform(EMGBL(:,iM),Freq,Args,Nb_Interp_Pnts) ;
                    toc
                    TFR.TFR.(Muscles{iM}) = norm ;
                    TFR.MedianFreq.(Muscles{iM})(:,:,iMovements) = Compute_Median_Frequency(norm,Wave_FreqS) ;
                else
                    TFR.TFR.(Muscles{iM})(:,:,iMovements) = nan(400,500) ;
                    TFR.MedianFreq.(Muscles{iM})(:,:,iMovements) = nan(500,1)  ;
                end
            end
        end
        
        TFR.freq = Wave_FreqS ;
        TFR.time = 1:500 ;
        TFR.muscles = Muscles ;
        TFR.Seg = Seg ;
        MedianFreq.muscles = Muscles ;
        
%         figure
%         imagesc(TFR.TFR.(Muscles{iM})(:,:,1));

        % Data saving
        cd(['E:\Bureau\Etienne\Extracted data\Fatigue'])
        save(['CleanData_TFR_' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'TFR')
        save(['CleanData_EMGenvelop_' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'EMG_envelop')        

        close all

        TFR = [];
        EMG = [];
        EMG_cycles = [];
        EMG_envelop = [];
        EMG_envelop_trials = [];
        EMGBL = [];
        EMGBP = [];
        EMGBS = [];
        Movement = [];
        Normalization = [];
    end
    
end    
