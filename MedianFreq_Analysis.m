clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

FoldersNames = dir('F:\Data\IRSST\RAW');
Subjects = {...
    FoldersNames(4:36).name...
   } ;
Subjects(5) = [];
Subjects(5) = [];

cd(['H:\Projet_ExpertsNovices\excel'])
[num, txt, raw] = xlsread('participants.xlsx');
partic = participants ;

% Cr�e une 2eme ligne avec l'expertise du participant
Subjects{2,1}=[];
Subjects{3,1}=[];
for iSubjects = 1:length(Subjects)
    Subjects{2,iSubjects} = txt{find(contains(partic,Subjects(1,iSubjects)))+3,10};
    if Subjects{2,iSubjects} == 'Novice'
        Subjects{3,iSubjects} = 0;
    else
        Subjects{3,iSubjects} = 1;
    end
end

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
    'ssp',  'ssp_EMG8';... supra epineux
    'isp',  'isp_EMG9';... infra ep.
    'subs',       'subs_EMG10';... subscapulaire
    } ;

Trials = {...
    'l_d_';...
    'l_u_';...
    'm_d_';...
    'm_u_';...
    's_d_';...
    's_u_';...
    };

TableData = [];
for iSubjects = 1:2
    for iTrials = 1:length(Trials)
        for iM = 1:length(Muscles)
%             MeanMedianFreqSeg.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM}) = [];
%             MeanMedianFreqTrial.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM}) = [];
% %             TableMeanMedianFreqSeg.(Trials{iTrials}).(Muscles{iM}) = [];
%             TableMeanMedianFreqSeg.(Muscles{iM}) = [];
%             TableMeanMedianFreqTrial.(Muscles{iM}) = [];
%             TableMeanMedianFreqPart.(Muscles{iM}) = [];
%             MedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM}) = [];
            TableData.(Muscles{iM}) = [];
        end
    end
end

tic
varInc = 1;
for iSubjects = 1:2
    for iTrials = 1:length(Trials)
        [strcTFR,strcSeg] = LoadTrials(iSubjects, Subjects{1,iSubjects}, Trials{iTrials});
        
        if Trials{iTrials}(1) == 'l'
            C1 = 1;
        elseif Trials{iTrials}(1) == 'm'
            C1 = 2;
        elseif Trials{iTrials}(1) == 's'
            C1 = 3;
        end
        
        if Trials{iTrials}(3) == 'd'
            C2 = 1;
        elseif Trials{iTrials}(3) == 'u'
            C2 = 2;
        end
        
            for iNumTrials = 1:length(fieldnames((eval(['strcTFR']))))
                NumberOfSeg = eval(['strcTFR.TFR_Trial' num2str(iNumTrials) '.NumberOfSeg']);
                for iSeg = 1:NumberOfSeg
                    Cycle = eval(['strcSeg.TFR_Trial' num2str(iNumTrials)]);
                    for iM = 1:length(Muscles)
                        
                    MedianFreq = eval(['strcTFR.TFR_Trial' num2str(iNumTrials) '.MedianFreq.' Muscles{iM}]);
                    MedianFreq = MedianFreq(Cycle(iSeg,1):Cycle(iSeg,2));
                    MedianFreq = MedianFreq(MedianFreq > (mean(MedianFreq)-3*std(MedianFreq)) & MedianFreq < (mean(MedianFreq)+3*std(MedianFreq)));...Supprime outliers
                    
                    TableData.(Muscles{iM})(varInc,1) = iSubjects;...Participant
                    TableData.(Muscles{iM})(varInc,2) = C1;...Condition essai (l,m,s)
                    TableData.(Muscles{iM})(varInc,3) = C2;...Condition essai (d,u)
                    TableData.(Muscles{iM})(varInc,4) = iNumTrials;...Quel essai
                    TableData.(Muscles{iM})(varInc,5) = iSeg;...Quel segment
                    TableData.(Muscles{iM})(varInc,6) = mean(MedianFreq);...Moyenne
                    TableData.(Muscles{iM})(varInc,7) = median(MedianFreq);...Mediane
                    end
                varInc = varInc + 1;
                end
            end
%         TFR.(Subjects{1,iSubjects}).(Trials{iTrials}) = strcTFR;
%         Seg.(Subjects{1,iSubjects}).(Trials{iTrials}) = strcSeg;
    end
end
toc

%% Moyenne de la fr�quence m�diane par segment (variable d�pendante)


iSubjects=1;
iM=1;
iTrials=1;
iNumTrials=1;
iSeg=1;

% Incr�mente les variables MedianFreq
for iM = 1:length(Muscles)
    varInc = 1;
    for iSubjects = 1:15
        for iTrials = 1:length(Trials)
            for iNumTrials = 1:length(fieldnames((eval(['TFR.' (Subjects{1,iSubjects}) '.' (Trials{iTrials})]))))
                NumberOfSeg = eval(['TFR.' (Subjects{1,iSubjects}) '.' (Trials{iTrials}) '.TFR_Trial' num2str(iNumTrials) '.NumberOfSeg']);
                MedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM}) = [MedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM}) permute(eval(['TFR.' (Subjects{1,iSubjects}) '.' (Trials{iTrials}) '.TFR_Trial' (num2str(iNumTrials)) '.MedianFreq.' (Muscles{iM}) '(:,:,1:' num2str(NumberOfSeg) ')']),[1 3 2])];
                for iSeg = 1:NumberOfSeg
                    
                    % Supprime les outliers
                    MedFreq = eval(['TFR.' (Subjects{1,iSubjects}) '.' (Trials{iTrials}) '.TFR_Trial' (num2str(iNumTrials)) '.MedianFreq.' (Muscles{iM}) '(:,:,' num2str(iSeg) ')']);
                    MedFreq = MedFreq(MedFreq > (mean(MedFreq)-3*std(MedFreq)) & MedFreq < (mean(MedFreq)+3*std(MedFreq)));
                    % Outliers supprim�s
                    % boxplot(MedFreq)
                    
                    TableData.(Muscles{iM})(varInc,1) = iSubjects;...Participant
                    TableData.(Muscles{iM})(varInc,2) = iTrials;...Condition essai (l_d)
                    TableData.(Muscles{iM})(varInc,3) = iNumTrials;...Quel essai
                    TableData.(Muscles{iM})(varInc,4) = iSeg;...Quel segment
                    TableData.(Muscles{iM})(varInc,5) = mean(MedFreq);...Moyenne
                    TableData.(Muscles{iM})(varInc,6) = median(MedFreq);...Mediane

                    varInc = varInc + 1;
                end
            end
        end
    end
end
VarNames = {'Participant','ConditionTrial','WhichTrial','WhichSeg','Mean','Median'};
save('H:\Bureau\Etienne\Extracted data\TableData_LMM.mat','TableData')
save('H:\Bureau\Etienne\Extracted data\VarNames_LMM.mat','VarNames')

%% PLOT DATA
% Clean Outliers
% CleanMedianFreqMuscles = MedianFreqMuscles;
% for iSubjects = 1:15
%     for iTrials = 1:length(Trials)
%        for iMuscles = 1:length(Muscles)
%            for iSeg = 1:size(MedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM}),2)
%                % Supprime les outliers
%                MedFreq = MedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM})(:,iSeg);
%                CleanMedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM})(MedFreq < mean(MedFreq)-3*std(MedFreq) | MedFreq > mean(MedFreq)+3*std(MedFreq),iSeg) = NaN;
%                % Outliers supprim�s
%            end
%        end
%     end
% end

% Boxlot des fr�quences m�dianes par muscles par participants
for iM = 1:size(Muscles,1)
    for iTrials = 1:length(Trials)
        for iSubjects = 1:1
            subplot(1,2,iSubjects) ; boxplot(MedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM})) ; title([Muscles{iM,1} '_' Trials{iTrials}], 'interpreter', 'none') ;
%             subplot(1,2,iSubjects+1) ; boxplot(CleanMedianFreqMuscles.(Subjects{1,iSubjects}).(Trials{iTrials}).(Muscles{iM})) ; title([Muscles{iM,1} '_' Trials{iTrials}], 'interpreter', 'none') ;
        end
        pause
    end
end
%% Linear Mixed Model
% Appelle le code R Ancova_EEG
cd('C:\Program Files\R\R-3.5.3\bin\x64\')
system('Rscript.exe H:\Bureau\Etienne\MATLAB\Functions\R_Function\LMM_Audrey.R') ;
% cd('C:\Users\p1038617\Desktop\RESULTS_EXERCISE_RestingState\scripts')

info = hdf5info('C:\Users\p1038617\Desktop\LMM_Audrey\tmpStats.h5') ;
for iVar=1:length(info.GroupHierarchy.Groups.Datasets)
    toto = strfind(info.GroupHierarchy.Groups.Datasets(1,iVar).Name,'/') ;
    Variable(iVar,1) = str2double(info.GroupHierarchy.Groups.Datasets(1,iVar).Name(1,toto(2)+1:end)) ;
end
for iVar = 1:length(Variable)
    line = Variable(iVar) ;
    dset = h5read(...
        info.GroupHierarchy.Groups.Datasets(iVar).Filename,...
        info.GroupHierarchy.Groups.Datasets(iVar).Name) ;
    
    P_values(iVar,1) = dset.Pr0x280x3EChisq0x29(2) ;
    %         clear dset
end
delete('C:\Users\p1038617\Desktop\LMM_Audrey\tmpStats.h5')... penser  supprimer la variable avant r�-ex�cution
delete('C:\Users\p1038617\Desktop\LMM_Audrey\tmpMatrix.mat')













