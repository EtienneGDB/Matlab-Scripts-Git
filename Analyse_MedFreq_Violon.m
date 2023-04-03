clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))

FoldersNames = dir('J:\Violon_SAB\Results\EMG');
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

FreqMed_Intervals_FirstSession = [];
FreqMed_Intervals_SecondSession = [];
for iM = 1:length(Muscles)
    FreqMed_Intervals_FirstSession.(Muscles{iM}) = [];
    FreqMed_Intervals_SecondSession.(Muscles{iM}) = [];
end
for iSubjects = 1:15
    if iSubjects == 6
    else
        % Load EMG data
        % FirstSession
        cd(['J:\Violon_SAB\Results\EMG\' eval(['Filenames_FirstSession' '{iSubjects}']) '\Preprocessed']);
        load('TFR_MedianFreq_Filtered_downSampled2_Gscale.mat')
        TFR_MedianFreq_FirstSession = TFR_MedianFreq.MedianFreq;
    end
    
    if iSubjects == 14
    else
        %SecondSession
        cd(['J:\Violon_SAB\Results\EMG\' eval(['Filenames_SecondSession' '{iSubjects}']) '\Preprocessed']);
        load('TFR_MedianFreq_Filtered_downSampled2_GscaleSAB.mat')
        TFR_MedianFreq_SecondSession = TFR_MedianFreq.MedianFreq;
    end
    clear TFR_MedianFreq
    
    % Normalize Second Session based on First Session
    for iM = 1:length(Muscles)
        if length(TFR_MedianFreq_FirstSession.(Muscles{iM})) > 500 & ...
                length(TFR_MedianFreq_SecondSession.(Muscles{iM})) > length(TFR_MedianFreq_FirstSession.(Muscles{iM}))
            TFR_MedianFreq_SecondSession.(Muscles{iM}) = ...
                TFR_MedianFreq_SecondSession.(Muscles{iM})(1:length(TFR_MedianFreq_FirstSession.(Muscles{iM})));
        end
    end

    Int = 1:100;
	for iM = 1:length(Muscles)
        FreqMed_Intervals_FirstSession.(Muscles{iM})(iSubjects,1) = iSubjects;
        FreqMed_Intervals_SecondSession.(Muscles{iM})(iSubjects,1) = iSubjects;
        if nansum(TFR_MedianFreq_FirstSession.(Muscles{iM}))~=0
            L = length(TFR_MedianFreq_FirstSession.(Muscles{iM}));
            for iInt = 1:99
                FreqMed_Intervals_FirstSession.(Muscles{iM})(iSubjects,iInt+1) =...
                    median(TFR_MedianFreq_FirstSession.(Muscles{iM})(Int(iInt)*L/100:Int(iInt+1)*L/100));
            end
        else
            FreqMed_Intervals_FirstSession.(Muscles{iM})(iSubjects,2:100) = NaN;
        end
        if nansum(TFR_MedianFreq_SecondSession.(Muscles{iM}))~=0
            L = length(TFR_MedianFreq_SecondSession.(Muscles{iM}));
            for iInt = 1:99
                FreqMed_Intervals_SecondSession.(Muscles{iM})(iSubjects,iInt+1) =...
                    median(TFR_MedianFreq_SecondSession.(Muscles{iM})(Int(iInt)*L/100:Int(iInt+1)*L/100));
            end
        else
            FreqMed_Intervals_SecondSession.(Muscles{iM})(iSubjects,2:100) = NaN;
        end
    end
end
% Discard some muscles weird in time/freq domains
for iM = 1:length(Muscles)
    if iM==10
        FreqMed_Intervals_FirstSession.(Muscles{iM})([3 7 9],2:100)=NaN;
        FreqMed_Intervals_SecondSession.(Muscles{iM})([4 10],2:100)=NaN;
    end
    if iM==4
        FreqMed_Intervals_FirstSession.(Muscles{iM})([2 10 11],2:100)=NaN;
        FreqMed_Intervals_SecondSession.(Muscles{iM})([5 9 13],2:100)=NaN;
    end
    if iM==7
        FreqMed_Intervals_FirstSession.(Muscles{iM})([11],2:100)=NaN;
        FreqMed_Intervals_SecondSession.(Muscles{iM})([11],2:100)=NaN;
    end
    if iM==3
        FreqMed_Intervals_SecondSession.(Muscles{iM})([3 15],2:100)=NaN;
    end
    if iM==8
        FreqMed_Intervals_SecondSession.(Muscles{iM})([7],2:100)=NaN;
    end
    if iM==11
        FreqMed_Intervals_SecondSession.(Muscles{iM})([7],2:100)=NaN;
    end
    if iM==9
        FreqMed_Intervals_SecondSession.(Muscles{iM})([9],2:100)=NaN;
    end
end

% Boxplot
iM=6;
for iM = 1:length(Muscles)
    subplot(2,1,1); boxplot(FreqMed_Intervals_FirstSession.(Muscles{iM})(:,2:end))
    title(['First Session ' Muscles{iM}])
    subplot(2,1,2); boxplot(FreqMed_Intervals_SecondSession.(Muscles{iM})(:,2:end))
    title(['Second Session ' Muscles{iM}])
    pause()
end
% Plot each participant
for iSubjects = 1:15
    for iM = 1:length(Muscles)
        figure(1); subplot(4,3,iM); plot(FreqMed_Intervals_FirstSession.(Muscles{iM})(iSubjects,2:end))
        title(['Subject ' num2str(iSubjects) ' ' Muscles{iM}])
        
        figure(2); subplot(4,3,iM); plot(FreqMed_Intervals_SecondSession.(Muscles{iM})(iSubjects,2:end))
        title(['Subject ' num2str(iSubjects) ' ' Muscles{iM}])
    end
    pause()
    clf(figure(1))
    clf(figure(2))
end

% Reduction to 10 intervals
FreqMed_10Intervals_FirstSession = [];
FreqMed_10Intervals_SecondSession = [];
Nparticipants = [];
for iM = 1:length(Muscles)
    FreqMed_10Intervals_FirstSession.(Muscles{iM}) = [];
    FreqMed_10Intervals_SecondSession.(Muscles{iM}) = [];
    Nparticipants.(Muscles{iM}) = [];
end
for iSubjects = 1:15
    Int = [0 10 20 30 40 50 60 70 80 90 100];
	for iM = 1:length(Muscles)
        FreqMed_10Intervals_FirstSession.(Muscles{iM})(iSubjects,1) = iSubjects;
        FreqMed_10Intervals_SecondSession.(Muscles{iM})(iSubjects,1) = iSubjects;
        if nansum(FreqMed_Intervals_FirstSession.(Muscles{iM}))~=0
            for iInt = 1:10
                FreqMed_10Intervals_FirstSession.(Muscles{iM})(iSubjects,iInt+1) =...
                    median(FreqMed_Intervals_FirstSession.(Muscles{iM})(iSubjects,Int(iInt)+1:Int(iInt+1)));
            end
        else
            FreqMed_10Intervals_FirstSession.(Muscles{iM})(iSubjects,2:length(Int)) = NaN;
        end
        if nansum(FreqMed_Intervals_SecondSession.(Muscles{iM}))~=0
            for iInt = 1:10
                FreqMed_10Intervals_SecondSession.(Muscles{iM})(iSubjects,iInt+1) =...
                    median(FreqMed_Intervals_SecondSession.(Muscles{iM})(iSubjects,Int(iInt)+1:Int(iInt+1)));
            end
        else
            FreqMed_10Intervals_SecondSession.(Muscles{iM})(iSubjects,2:length(Int)) = NaN;
        end
        
        % Si une des deux session est NaN -> l'autre est NaN
        if isnan(FreqMed_10Intervals_FirstSession.(Muscles{iM})(iSubjects,2))
            FreqMed_10Intervals_SecondSession.(Muscles{iM})(iSubjects,2:length(Int)) = NaN;
        elseif isnan(FreqMed_10Intervals_SecondSession.(Muscles{iM})(iSubjects,2))
            FreqMed_10Intervals_FirstSession.(Muscles{iM})(iSubjects,2:length(Int)) = NaN;
        end
    end
end

% Boxplot
iM=6;
for iM = 1:length(Muscles)
    subplot(2,1,1); boxplot(FreqMed_10Intervals_FirstSession.(Muscles{iM})(:,2:end))
    title(['First Session ' Muscles{iM}])
    subplot(2,1,2); boxplot(FreqMed_10Intervals_SecondSession.(Muscles{iM})(:,2:end))
    title(['Second Session ' Muscles{iM}])
    pause()
end

% Stats
% Test Variance equality
for iM = 4:length(Muscles)
    L = size(FreqMed_10Intervals_FirstSession.(Muscles{iM}),2);
    p = vartestn(FreqMed_10Intervals_FirstSession.(Muscles{iM})(:,[2 L]),'TestType','LeveneAbsolute');
    pValuesVar(iM,:) = p;
end

% Test Normality
for iM = 1:length(Muscles)
    FreqMed_10Intervals_SecondSession.(Muscles{iM})...
      (isnan(FreqMed_10Intervals_SecondSession.(Muscles{iM})(:,2)),:) = [];
    [h p] = kstest(FreqMed_10Intervals_FirstSession.(Muscles{iM})(:,2));
    pValuesNorm(iM,1) = p;

    [h p] = kstest(FreqMed_10Intervals_SecondSession.(Muscles{iM})(:,L));
    pValuesNorm(iM,2) = p;
end

variation = [];
% Wilcoxon Signed-rank test
for iM = 1:length(Muscles)
    if iM == 3
    else
        L = size(FreqMed_10Intervals_FirstSession.(Muscles{iM}),2);
        
        x1 = FreqMed_10Intervals_FirstSession.(Muscles{iM})(:,2);
        y1 = FreqMed_10Intervals_FirstSession.(Muscles{iM})(:,L);
        %     variation.(Muscles{iM})(:,1) = y1-x1;
        
        x2 = FreqMed_10Intervals_SecondSession.(Muscles{iM})(:,2);
        y2 = FreqMed_10Intervals_SecondSession.(Muscles{iM})(:,L);
        %     variation.(Muscles{iM})(:,2) = y2-x2;
        
        % Count N participants
        Nparticipants.(Muscles{iM}) = sum(~isnan(x2));
        
        [p,h,stats] = signrank(x1,y1);
        pValuesWSRT1(iM,1) = p;
        Colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
        figure(1); subplot(4,3,iM); UnivarScatter([x1 y1], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
        %     figure(1); subplot(4,3,iM); boxplot([x1 y1])
        if pValuesWSRT1(iM,1) < 0.05
            hold on
            plot(1.5, max([x1;y1]),'*', 'color', 'black', 'LineWidth', 1, 'markers', 30, 'MarkerSize', 10)
        end
        title(['First Session ' Muscles{iM}])
        
        [p,h,stats] = signrank(x2,y2);
        pValuesWSRT2(iM,1) = p;
        Colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
        figure(2); subplot(4,3,iM); UnivarScatter([x2 y2], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
        if pValuesWSRT2(iM,1) < 0.05
            hold on
            plot(1.5, max([x2;y2]),'*', 'color', 'black', 'LineWidth', 1, 'markers', 30, 'MarkerSize', 10)
        end
        title(['Second Session ' Muscles{iM}])
    end
end

% Test Variation
% Test Variance equality
for iM = 4:length(Muscles)
    p = vartestn(variation.(Muscles{iM}),'TestType','LeveneAbsolute');
    pValuesVar(iM,:) = p;
end

% Test Normality
for iM = 1:length(Muscles)
    variation.(Muscles{iM})...
      (isnan(variation.(Muscles{iM})(:,1)),:) = [];
    [h p] = kstest(variation.(Muscles{iM})(:,1));
    pValuesNorm(iM,1) = p;

    [h p] = kstest(variation.(Muscles{iM})(:,2));
    pValuesNorm(iM,2) = p;
end

% Wilcoxon Signed-rank test
for iM = 1:length(Muscles)
    [p,h,stats] = signrank(variation.(Muscles{iM})(:,1),variation.(Muscles{iM})(:,2));
    pValuesVar(iM,1) = p;
    figure(3); subplot(4,3,iM); boxplot([variation.(Muscles{iM})(:,1) variation.(Muscles{iM})(:,2)])
    if pValuesVar(iM,1) < 0.05
        hold on
        plot(1.5, max([variation.(Muscles{iM})(:,1);variation.(Muscles{iM})(:,1)]),'*', 'color', 'black', 'LineWidth', 1, 'markers', 30, 'MarkerSize', 10)
    end
    title(['Variation ' Muscles{iM}])
end