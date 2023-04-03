clear ;
close all ;
clc ;

addpath E:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\Matlab-Functions-Git'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\FFT'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\Entropy Fouaz'))

% addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))

FoldersNames = dir('V:\Violon_SAB\Results\EMG');
Filenames = {...
    FoldersNames(5:34).name...
   } ;
inc = 1;
for iFiles = 1:2:length(Filenames)
    Filenames_FirstSession(inc) = Filenames(iFiles);
    inc = inc +1;
end
inc = 1;
for iFiles = 2:2:length(Filenames)
    Filenames_SecondSession(inc) = Filenames(iFiles);
    inc = inc +1;
end
    
% Define muscles used
Muscles = {'InfEp_L';'SupEp_L';'SousScap_L';'SCM_L';'Biceps_L';'Delt_Med_L';'Trap_Sup_L';'Delt_Ant_R';  ...
    'Delt_Med_R';'Delt_Post_R';'Trap_Sup_R'} ;

for iSubjects = 15:15
    % Load EMG data
    Session = 'FirstSession'; % Change between FirstSession & SecondSession
    cd(['V:\Violon_SAB\Results\EMG\' eval(['Filenames_' Session '{iSubjects}']) '\Preprocessed']);
    
    if Session(1) == 'F'
        File = 'Filtered_downSampled_Gscale.mat';
    elseif Session(1) == 'S'
        File = 'Filtered_downSampled_GscaleSAB.mat';
    end
    load(File);
    Freq = 1000 ; % The filtered signal was down-sampled to 1000 Hz
        
    EMG = Norm_EMGBL ;
    EMG1 = EMG(1:length(EMG)/4,:);
    EMG2 = EMG(length(EMG)/4:length(EMG)/2,:);
    EMG3 = EMG(length(EMG)/2:3*length(EMG)/4,:);
    EMG4 = EMG(3*length(EMG)/4:length(EMG),:);
        
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
        for iM = 1:length(Muscles)
            disp('.')
            if sum(EMG(:,iM))~=0 & nansum(EMG(:,iM))~=0
%                 tic
                [norm, Time, Wave_FreqS] = TimeFreqTransform(EMG2(~isnan(EMG2(:,iM)),iM),Freq,Args,Nb_Interp_Pnts) ;
%                 toc
%                 TFR.TFR.(Muscles{iM}) = norm ;
                TFR.MedianFreq.(Muscles{iM})(:,:,:) = Compute_Median_Frequency(norm,Wave_FreqS) ;
                TFR.SpectralEntropy.(Muscles{iM})(:,:,:) = Compute_Spectral_Entropy(norm, Wave_FreqS) ;
            else
                TFR.TFR.(Muscles{iM})(:,:,:) = nan(400,500) ;
                TFR.MedianFreq.(Muscles{iM})(:,:,:) = nan(500,1)  ;
                TFR.SpectralEntropy.(Muscles{iM})(:,:,:) = nan(500,1)  ;
            end
        end
        
%         TFR.freq = Wave_FreqS ;
%         TFR.muscles = Muscles ;
%         MedianFreq.muscles = Muscles ;
        
        TFR_MedianFreq = {};
        TFR_MedianFreq.(['MedianFreq']) = TFR.MedianFreq ;
        
        TFR_SpectralEntropy = {};
        TFR_SpectralEntropy.(['SpectralEntropy']) = TFR.SpectralEntropy ;

        % Plot TFR/Median_Freq/FFT
%         for iM=1:length(Muscles)
%             figure(2);
%             subplot(2,1,1) ; imagesc(norm); title([HaFilenames{2,iSubjects} '__' Muscles{iM}])
% %             subplot(3,1,2) ; plot(TFR.MedianFreq.(Muscles{iM})(:,:),'.')
%             subplot(2,1,2) ; plot(Normalization(:,iM))
% %             pause
% %             close all
%         end

        % Data saving
%         cd(['J:\Piano_Fatigue\Data_Exported\TFR'])
%         save(['CleanData_TFR_' (HaFilenames{2,iSubjects}) '.mat'],'TFR')

%         cd(['V:\Piano_Fatigue\Data_Exported\TFR_MedianFreq'])
        save(['TFR_MedianFreq_' File],'TFR_MedianFreq')

%         cd(['V:\Piano_Fatigue\Data_Exported\TFR_SpectralEntropy'])
        save(['TFR_SpectralEntropy_' File],'TFR_SpectralEntropy')

        TFR = [];
        TFR_MedianFreq = [];
        TFR_SpectralEntropy = [];
        EMG = [];
        clear functions
end
    
