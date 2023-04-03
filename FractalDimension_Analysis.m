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

% Crée une 2eme ligne avec l'expertise du participant
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

Trials = {...
    'l_d_';...
    'l_u_';...
    'm_d_';...
    'm_u_';...
    's_d_';...
    's_u_';...
    };

TableData_FracDim = [];
for iSubjects = 1:length(Subjects)
    for iTrials = 1:length(Trials)
        for iM = 1:length(Muscles)
            TableData_FracDim.(Muscles{iM}) = [];
        end
    end
end

varInc = 1;
for iSubjects = 26:length(Subjects)
    for iTrials = 1:length(Trials)
        [strcTFR,strcSeg,strcAct,strcMob,strcSpecEnt,strcAmp,strcSampEnt,strcFracDim] = LoadTrials(iSubjects, Subjects{1,iSubjects}, Trials{iTrials});
        
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
                    
                    ActivitySeg = eval(['strcAct.Trial' num2str(iNumTrials) '.Seg.' Muscles{iM}]);
                    MobilitySeg = eval(['strcMob.Trial' num2str(iNumTrials) '.Seg.' Muscles{iM}]);

                    SpecEnt = eval(['strcSpecEnt.Trial' num2str(iNumTrials) '.SpectralEntropy.' Muscles{iM}]);
                    SampEnt = eval(['strcSampEnt.Trial' num2str(iNumTrials) '.Seg.' Muscles{iM}]);
                    
                    Amplitude = eval(['strcAmp.Trial' num2str(iNumTrials) '.Seg.' Muscles{iM}]);

                    FracDim = eval(['strcFracDim.Trial' num2str(iNumTrials) '.' Muscles{iM}]);

                    if isnan(MedianFreq) | MedianFreq == 0
                        TableData.(Muscles{iM})(varInc,:) = NaN;
                        
                    elseif (iSubjects==10 & C1==1 & C2==1 & iNumTrials==2 & iM==10) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==10 & C1==1 & C2==1 & iNumTrials==3 & iM==10) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==11 & C1==3 & C2==2 & iNumTrials==4 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==12 & C1==1 & C2==1 & iNumTrials==1 & iM==10) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==13 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==14 & iM==10) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==1 & C2==2 & iNumTrials==4 & iM==7) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==2 & C2==1 & iNumTrials==4 & iM==7) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==2 & C2==1 & iNumTrials==4 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==2 & C2==2 & iNumTrials==3 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==3 & C2==1 & iNumTrials==3 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==3 & C2==1 & iNumTrials==4 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==3 & C2==2 & iNumTrials==3 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==15 & C1==3 & C2==2 & iNumTrials==4 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==21 & C1==1 & C2==1 & iNumTrials==2) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==21 & C1==1 & C2==2 & iNumTrials==1 & iM==6) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==21 & C1==2 & C2==1 & iNumTrials==1 & iM==6) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==21 & C1==3 & C2==2 & iNumTrials==4 & iM==3) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==21 & C1==3 & C2==2 & iNumTrials==5 & iM==3) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==1 & C2==1 & iNumTrials==3) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==1 & C2==2 & iNumTrials==1 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==1 & C2==2 & iNumTrials==2 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==1 & C2==2 & iNumTrials==3 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==1 & C2==2 & iNumTrials==4 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==2 & C2==1 & iNumTrials==1 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==2 & C2==1 & iNumTrials==2 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==2 & C2==1 & iNumTrials==3 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==2 & C2==1 & iNumTrials==4 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==2 & C2==2 & iNumTrials==4 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==3 & C2==1 & iNumTrials==1 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==3 & C2==1 & iNumTrials==2 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==23 & C1==3 & C2==1 & iNumTrials==3 & iM==8) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==27 & C1==2 & C2==2 & iNumTrials==1 & iM==9) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==27 & C1==3 & C2==2 & iNumTrials==1 & iM==9) ; TableData.(Muscles{iM})(varInc,:) = NaN;
                    elseif (iSubjects==29 & iM==7) ; TableData.(Muscles{iM})(varInc,:) = NaN;

                    else
                        MedianFreq = MedianFreq(Cycle(iSeg,1):Cycle(iSeg,2));
                        MedianFreq = MedianFreq(MedianFreq > (mean(MedianFreq)-3*std(MedianFreq)) & MedianFreq < (mean(MedianFreq)+3*std(MedianFreq)));...Supprime outliers
            
                        SpecEnt = SpecEnt(Cycle(iSeg,1):Cycle(iSeg,2));
                        SpecEnt = SpecEnt(SpecEnt > (mean(SpecEnt)-3*std(SpecEnt)) & SpecEnt < (mean(SpecEnt)+3*std(SpecEnt)));...Supprime outliers
                            
                        TableData_FracDim.(Muscles{iM})(varInc,1) = iSubjects;...Participant
                        TableData_FracDim.(Muscles{iM})(varInc,2) = Subjects{3,iSubjects};...Expertise
                        TableData_FracDim.(Muscles{iM})(varInc,3) = C1;...Condition essai (l,m,s)
                        TableData_FracDim.(Muscles{iM})(varInc,4) = C2;...Condition essai (d,u)
                        TableData_FracDim.(Muscles{iM})(varInc,5) = iNumTrials;...Quel essai
                        TableData_FracDim.(Muscles{iM})(varInc,6) = iSeg;...Quel segment
                        
                        TableData_FracDim.(Muscles{iM})(varInc,7) = Amplitude(iSeg);...Amplitude
                        TableData_FracDim.(Muscles{iM})(varInc,8) = ActivitySeg(iSeg);...Activity
                        TableData_FracDim.(Muscles{iM})(varInc,9) = MobilitySeg(iSeg);...Mobility
                        
                        TableData_FracDim.(Muscles{iM})(varInc,10) = SampEnt(iSeg);...SampEn
                        TableData_FracDim.(Muscles{iM})(varInc,11) = mean(SpecEnt);...Moyenne SpecEnt
                        TableData_FracDim.(Muscles{iM})(varInc,12) = median(SpecEnt);...Median SpecEnt
                        TableData_FracDim.(Muscles{iM})(varInc,13) = mean(MedianFreq);...Moyenne MedianFreq
                        TableData_FracDim.(Muscles{iM})(varInc,14) = median(MedianFreq);...Mediane MedianFreq
                        
                        TableData_FracDim.(Muscles{iM})(varInc,15) = FracDim(iSeg);...FracDim
                    end
                end
            varInc = varInc + 1;
            end
        end
    end
    save('H:\Bureau\Etienne\Extracted data\TableData_FracDim.mat','TableData_FracDim')
end

TABLEDATA_FracDim = TableData_FracDim;

% Remove outliers
outDat = [];
for iM = 1:length(Muscles)
    outDat.(Muscles{iM})(1,1) = 0;
    for iC1 = 1:3
        for iC2 = 1:2
            if iC1 == 1
                NSeg = 6;
            end
            if iC1 == 2
                NSeg = 9;
            end
            if iC1 == 3
                NSeg = 12;
            end
            for iSeg = 1:NSeg
                idDat = find(TABLEDATA_FracDim.(Muscles{iM})(:,3)==iC1 & TABLEDATA_FracDim.(Muscles{iM})(:,4)==iC2 & TABLEDATA_FracDim.(Muscles{iM})(:,6)==iSeg);
                
                % Clean Mean MedianFreq
                Dat = TABLEDATA_FracDim.(Muscles{iM})(TABLEDATA_FracDim.(Muscles{iM})(:,3)==iC1 & TABLEDATA_FracDim.(Muscles{iM})(:,4)==iC2 & TABLEDATA_FracDim.(Muscles{iM})(:,6)==iSeg,7);
                Out = find(Dat > prctile(Dat,95));
                TABLEDATA_FracDim.(Muscles{iM})(idDat(Out),7) = NaN;
                outDat.(Muscles{iM})(1,1) = outDat.(Muscles{iM})(1,1)+length(Out);
                
                % Clean Median MedianFreq
                Dat = TABLEDATA_FracDim.(Muscles{iM})(TABLEDATA_FracDim.(Muscles{iM})(:,3)==iC1 & TABLEDATA_FracDim.(Muscles{iM})(:,4)==iC2 & TABLEDATA_FracDim.(Muscles{iM})(:,6)==iSeg,8);
                Out = find(Dat > prctile(Dat,95));
                TABLEDATA_FracDim.(Muscles{iM})(idDat(Out),8) = NaN;
                
                % Clean Activity
                Dat = TABLEDATA_FracDim.(Muscles{iM})(TABLEDATA_FracDim.(Muscles{iM})(:,3)==iC1 & TABLEDATA_FracDim.(Muscles{iM})(:,4)==iC2 & TABLEDATA_FracDim.(Muscles{iM})(:,6)==iSeg,9);
                Out = find(Dat > prctile(Dat,95));
                TABLEDATA_FracDim.(Muscles{iM})(idDat(Out),9) = NaN;

                % Clean Mobility
                Dat = TABLEDATA_FracDim.(Muscles{iM})(TABLEDATA_FracDim.(Muscles{iM})(:,3)==iC1 & TABLEDATA_FracDim.(Muscles{iM})(:,4)==iC2 & TABLEDATA_FracDim.(Muscles{iM})(:,6)==iSeg,10);
                Out = find(Dat > prctile(Dat,95));
                TABLEDATA_FracDim.(Muscles{iM})(idDat(Out),10) = NaN;
                
                % Clean Mean SpecEnt
                Dat = TABLEDATA_FracDim.(Muscles{iM})(TABLEDATA_FracDim.(Muscles{iM})(:,3)==iC1 & TABLEDATA_FracDim.(Muscles{iM})(:,4)==iC2 & TABLEDATA_FracDim.(Muscles{iM})(:,6)==iSeg,11);
                Out = find(Dat > prctile(Dat,95));
                TABLEDATA_FracDim.(Muscles{iM})(idDat(Out),10) = NaN;

                % Clean Median SpecEnt
                Dat = TABLEDATA_FracDim.(Muscles{iM})(TABLEDATA_FracDim.(Muscles{iM})(:,3)==iC1 & TABLEDATA_FracDim.(Muscles{iM})(:,4)==iC2 & TABLEDATA_FracDim.(Muscles{iM})(:,6)==iSeg,12);
                Out = find(Dat > prctile(Dat,95));
                TABLEDATA_FracDim.(Muscles{iM})(idDat(Out),10) = NaN;
            end
        end
    end
    TABLEDATA_FracDim.(Muscles{iM})(find(TABLEDATA_FracDim.(Muscles{iM})(:,15)==0),15)=NaN;
end
boxplot(Dat)
find(isnan(TABLEDATA_FracDim.(Muscles{iM})(:,10)))

% Boxplot
for iM = 1:length(Muscles)
    boxplot(TABLEDATA_FracDim.(Muscles{iM})(:,14),TABLEDATA_FracDim.(Muscles{iM})(:,[2 3 4 5]))
    title((Muscles{iM}))
    pause
end
for iM = 1:length(Muscles)
    boxplot(TABLEDATA_FracDim.(Muscles{iM})(:,11),TABLEDATA_FracDim.(Muscles{iM})(:,[6]))
    title((Muscles{iM}))
    pause
end

% VarNames = {'Participant','Expertise','BoxSize','UpDown','Trial','Seg','Amplitude','Activity','Mobility','SampEn','MeanSE','MedianSE','MeanMF','MedianMF'};
% save('H:\Bureau\Etienne\Extracted data\TableData_LMM.mat','TableData')
% save('H:\Bureau\Etienne\Extracted data\CleanTABLEDATA_LMM.mat','TABLEDATA')
% save('H:\Bureau\Etienne\Extracted data\VarNamesEN.mat','VarNames')
% 

% Mean Activation level
Act_lev = [];
for iM = 1:length(Muscles)
    for iP = 1:31
        temp1 = TABLEDATA.(Muscles{iM})(TABLEDATA.(Muscles{iM})(:,1)==iP,:);
        for iC = 1:3
            temp2 = temp1(temp1(:,3)==iC,:);
            Act_lev.(Muscles{iM})(iP,iC) = mean(temp2(:,7),'omitnan');
        end
    end
end

Mean_Act_lev = [];
SD_Act_lev = [];
for iM = 1:length(Muscles)
    Mean_Act_lev(iM,:) = mean(Act_lev.(Muscles{iM}),'omitnan');
    SD_Act_lev(,:) = std(Act_lev.(Muscles{iM}),'omitnan');
end










