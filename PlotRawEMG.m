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

    FolderContent = dir(['F:\Data\IRSST\RAW\' Subjects{iSubjects} '\fatigue']);
    FileNames = {...
    FolderContent(3:length(FolderContent)).name...
   } ;

    while (str2num(FileNames{1}(1)) == 1) | (str2num(FileNames{1}(1)) == 2)
        FileNames(1) = [];
    end
    
    %Crée une deuxième ligne avec le nom du fichier sans .c3d
    FileNames{2,1}=[];
    for iFileNames = 1:length(FileNames)
        FileNames{2,iFileNames} = erase(FileNames{1,iFileNames},'.c3d');
    end
    
    FiletoPeak = {};
    varInc = 1;
    for iFileNames = 1:length(FileNames)
        if str2num(FileNames{2,iFileNames}(5)) == 1
            FiletoPeak{varInc} = FileNames{2,iFileNames};
            varInc = varInc +1;
        end
        if str2num(FileNames{2,iFileNames}(5)) == 3
            FiletoPeak{varInc} = FileNames{2,iFileNames};
            varInc = varInc +1;
        end
    end
    
    %iFiles=1;
    % Load EMG data
    for iFiletoPeak = 1:length(FiletoPeak)
        cd(['H:\Bureau\Etienne\Extracted data\Fatigue\DataSelec'])
        load(['RawEMG_Muscles_' (FiletoPeak{1,iFiletoPeak}) '_' (Subjects{iSubjects}) '.mat']);
        Data = DataSelec ;
        Freq = 2000 ;
        
        EMG = [] ;
        % EMG
        for iMuscles = 1:size(Muscles,1)
            EMG(:,iMuscles) = Data.(Muscles{iMuscles,2}) ;
%             subplot(3,3,iMuscles) ; plot(EMG(:,iMuscles)) ; title(Muscles{iMuscles,1}, 'interpreter', 'none') ;
        end

        % Preprocessing
        % Freq : fréquence d'échantillonage de l'EMG
        [b,a] = butter(2,2*[10 400]/Freq) ; % Parametre du filtre BP 10-400 Hz
        EMGBP = filtfilt(b,a,EMG) ;
        EMGBS = EMGBP ;

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
%         figure
%         plot(f,P1); hold on
%         plot([0 1000], [mean(P1)+5*std(P1) mean(P1)+5*std(P1)], 'color','Red')
        
        EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)
        
        % Spectre de fréquence raw
        T = 1/Freq;             % Sampling period       
        L = length(EMGBL(:,iMuscles));             % Length of signal
        t = (0:L-1)*T;        % Time vector

        figure
        for iM=1:length(Muscles)
            Y = fft(EMGBL(:,iM));

            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);

            f = Freq*(0:(L/2))/L;
            subplot(3,4,iM); plot(f,P1) 
            title(['Part' (num2str(iSubjects)) '_' (Muscles{iM})])
            xlabel('f (Hz)')
            ylabel('|P1(f)|')
        end

        % Carte Temps-Fréquence
        % TFR
%         FreqMin = 1 ;
%         FreqMax = 400 ;
%         Resolution = 1 ;
%         WaveNumber = 7 ;
%         Args = WaveletParameters(FreqMin,FreqMax,Resolution,WaveNumber,Freq) ;
%         FreqRange = length(FreqMin:Resolution:FreqMax) ;
%         Nb_Interp_Pnts = length(EMGBL) ;
% 
%         TFR = [] ;
%            iM = 5;
%            [norm, Time, Wave_FreqS] = TimeFreqTransform(EMGBL(:,iM),Freq,Args,Nb_Interp_Pnts) ;
%            TFR.TFR.(Muscles{iM})(:,:,:) = norm ;
%            figure
%            subplot(2,1,1) ; plot(EMGBL(:,iM));
%            subplot(2,1,2) ; imagesc([1:length(EMGBL(:,iM))],Wave_FreqS,TFR.TFR.(Muscles{iM})(:,:,1));

        % Dection of activity
%         Seg = ActivityDetection2(EMGBL, 'biceps', 1000, 0, 2000, Muscles);
%         pause

        figure
        for iMuscles = 4:size(Muscles,1)
            subplot(3,4,iMuscles) ; plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             for iSeg = 1:length(Seg)
%                 hold on
%                 plot([Seg(iSeg,1) Seg(iSeg,1)], [min(EMGBL(:,iMuscles)) max(EMGBL(:,iMuscles))],'Color','Green')
%                 plot([Seg(iSeg,2) Seg(iSeg,2)], [min(EMGBL(:,iMuscles)) max(EMGBL(:,iMuscles))],'Color','yellow')
%             end
        end
        pause
    end
    close all
end

signal_EMG = EMGBL(1:14360,:);
plot(signal_EMG)
cd(['C:\Users\p1098713\Documents\2.Post-Doc\Cours KIN6839\Exemples'])
save(['EMGBL_all.mat'],'signal_EMG')