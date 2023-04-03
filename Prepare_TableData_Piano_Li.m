clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

FoldersNames = dir('J:\Piano_Fatigue\Data_Exported\EMG_Max_MVC');
MVCFilenames = {...
    FoldersNames(3:52).name...
   } ;
    %Crée une deuxième ligne avec le nom du fichier sans .mat
    MVCFilenames{2,1}=[];
    for iMVCFilenames = 1:length(MVCFilenames)
        MVCFilenames{2,iMVCFilenames} = erase(MVCFilenames{1,iMVCFilenames},'.mat');
    end
    
FoldersNames = dir('J:\Piano_Fatigue\Data_Exported\EMG_Li');
HaFilenames = {...
    FoldersNames(3:52).name...
   } ;
    %Crée une deuxième ligne avec le nom du fichier sans .mat
    HaFilenames{2,1}=[];
    for iHaFilenames = 1:length(HaFilenames)
        HaFilenames{2,iHaFilenames} = erase(HaFilenames{1,iHaFilenames},'.mat');
    end

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

cd(['J:\Piano_Fatigue\Data_Exported'])
load(['cycles_Li.mat'])

cd(['J:\Piano_Fatigue\Matlab_matrix\Info_participants'])
load(['Info_participants_corrected'])
for iP = 1:50
    Demography(iP,1) = Info_participants(iP).Sexe;
    Demography(iP,2) = Info_participants(iP).Age;
    Demography(iP,3) = Info_participants(iP).Poids;
    Demography(iP,4) = Info_participants(iP).nb_annees;
    Demography(iP,5) = Info_participants(iP).nb_heures;
    Demography(iP,6) = mean(Info_participants(iP).MVC);

    Lat(iP,:) = Info_participants(iP).Lateralite;
    Sex(iP,:) = Info_participants(iP).Sexe;
end

TableData = [];
TABLEDATA_Int = [];
for iM = 1:length(Muscles)
    TableData.(Muscles{iM}) = [];
    TABLEDATA_Int.(Muscles{iM}) = [];
end

TimePlaying = [];
for iSubjects = 18:50
    TimePlaying(iSubjects,:) = participants(iSubjects).t_do(length(participants(iSubjects).t_do))--participants(iSubjects).t_do(1);
end
r = corr(Demography(part1,5),TimePlaying(part1,:),'Type','Pearson');
plot(Demography(part1,4),TimePlaying(part1,:),'o')
xlabel('Years of experience')
ylabel('Endurance (s)')
title(['Chord task : r² = ' num2str(r^2)])

plot(TimePlaying,'.')
idx = kmeans(TimePlaying,2);
E = evalclusters(TimePlaying ,'kmeans','silhouette', 'KList', [1:6]);
plot(E)
verif = unique(TABLEDATA_Int_LessSixMin.m1_2(:,1));
length(find(idx==1))
length(find(idx==2))
gscatter(TimePlaying(:),idx)
min(TimePlaying(idx==1))
std(TimePlaying(idx==1))
TimePlaying(17)=[];
idx(17)=[];
vartestn(TimePlaying(:,1),Group,'TestType','LeveneAbsolute')
[h,p,ci,stats]=ttest2(TimePlaying(idx==1),TimePlaying(idx==2),'Tail','both','Alpha',0.05,'Vartype','equal')

varInc = 1;
for iSubjects = 18:50
    tic
    [TFR_MedianFreq,TFR_SpectralEntropy,Activity,Mobility,Amplitude,SampleEntropy] = LoadTrials_Piano(iSubjects);
    toc
    
    RPE = Info_participants(iSubjects).Liszt;
    RPE = rmmissing(RPE);
    inc = 30;
    for iRPE = 1:length(RPE)
        TimeRPE = abs((participants(iSubjects).t_do-participants(iSubjects).t_do(1))-iRPE*inc);
        dist = min(TimeRPE);
        RPE(iRPE,2) = find(TimeRPE==dist);
    end
    participants(iSubjects).t_do = participants(iSubjects).t_do*(2048/2);

            for iSeg = 1:(length(participants(iSubjects).t_do)-1)
                for iM = 1:length(Muscles)
                        
                    MedianFreq = eval(['TFR_MedianFreq.MedianFreq.' Muscles{iM}]);
                    SpecEnt = eval(['TFR_SpectralEntropy.SpectralEntropy.' Muscles{iM}]);
                    MobilitySeg = eval(['Mobility.' Muscles{iM}]);
                    ActivitySeg = eval(['Activity.' Muscles{iM}]);
                    SampEntSeg = eval(['SampleEntropy.' Muscles{iM}]);
                    AmplitudeSeg = eval(['Amplitude.' Muscles{iM}]);

                        MedianFreq = MedianFreq(round(participants(iSubjects).t_do(iSeg)):round(participants(iSubjects).t_do(iSeg+1)));
                        MedianFreq = MedianFreq(MedianFreq > (mean(MedianFreq)-3*std(MedianFreq)) & MedianFreq < (mean(MedianFreq)+3*std(MedianFreq)));...Supprime outliers
            
                        SpecEnt = SpecEnt(round(participants(iSubjects).t_do(iSeg)):round(participants(iSubjects).t_do(iSeg+1)));
                        SpecEnt = SpecEnt(SpecEnt > (mean(SpecEnt)-3*std(SpecEnt)) & SpecEnt < (mean(SpecEnt)+3*std(SpecEnt)));...Supprime outliers
                            
                        TableData.(Muscles{iM})(varInc,1) = iSubjects;...Participant
                        TableData.(Muscles{iM})(varInc,2) = iSeg;...Quel segment
                        if any(iSeg==RPE(:,2))
                            TableData.(Muscles{iM})(varInc,3) = 1;...Quel segment
                            TableData.(Muscles{iM})(varInc,4) = RPE(RPE(:,2)==iSeg,1);...Quel segment
                        end
                    
                        TableData.(Muscles{iM})(varInc,7) = AmplitudeSeg(iSeg);...Amplitude
                        TableData.(Muscles{iM})(varInc,8) = ActivitySeg(iSeg);...Activity
                        TableData.(Muscles{iM})(varInc,9) = MobilitySeg(iSeg);...Mobility
                        
                        TableData.(Muscles{iM})(varInc,10) = SampEntSeg(iSeg);...SampEn
                        TableData.(Muscles{iM})(varInc,11) = mean(SpecEnt);...Moyenne SpecEnt
                        TableData.(Muscles{iM})(varInc,12) = mean(MedianFreq);...Moyenne MedianFreq
                end
            varInc = varInc + 1;
            end
end
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\TableData_Piano_Li.mat','TableData')

TABLEDATA = TableData;
% Suppression données abérentes
% Participant 22 wierdo
% for iM = 1:length(Muscles)
%     TABLEDATA.(Muscles{iM})(TABLEDATA.(Muscles{iM})(:,1)==22,:)=[];
% end

% Normalisation des données par rapport au temps
MaxPoint=0;
for iM = 1:length(Muscles)
    for iP = 1:50
        Temp =  TABLEDATA.(Muscles{iM})(TABLEDATA.(Muscles{iM})(:,1)==iP,:);
        if MaxPoint < length(Temp)
            MaxPoint = length(Temp);
        end
        Temp(:,5) = Temp(:,2)*100/length(Temp);
        TABLEDATA.(Muscles{iM})(TABLEDATA.(Muscles{iM})(:,1)==iP,5) = Temp(:,5);
    end
end
for iM = 1:length(Muscles)
    for iP = 1:50
        if iP ~= 17 %& iP ~= 22
            Temp =  TABLEDATA.(Muscles{iM})(TABLEDATA.(Muscles{iM})(:,1)==iP,:);
            Int = interp1(1:length(Temp),Temp,linspace(1,length(Temp),MaxPoint));
            Idx = find(Int(:,3));
            IdxS = Idx(1);
            vinc = 1;
            for L = 2:length(Idx)
                if Idx(L)-Idx(L-1) > 1
                    IdxS(vinc,2) = Idx(L-1);
                    vinc = vinc+1;
                    IdxS(vinc,1) = Idx(L);
                end
            end
            IdxS(length(IdxS),2) = Idx(length(Idx));
            for M = 1:length(IdxS)
                V = max(Int(IdxS(M,1):IdxS(M,2),3));
                Rp = find(Int(:,3)==V);
                Int(IdxS(M,1):IdxS(M,2),3) = 0;
                Int(Rp,3) = 1;
                Int(Rp,4) = ceil(Int(Rp,4));
                temp = Int(IdxS(M,1):IdxS(M,2),4);
                for j = 1:length(temp)
                    if temp(j) < max(temp)
                        temp(j) = 0;
                    end
                end
                Int(IdxS(M,1):IdxS(M,2),4) = temp;
                Int(:,6) = 1:MaxPoint;
            end
            TABLEDATA_Int.(Muscles{iM}) = [TABLEDATA_Int.(Muscles{iM}); Int];
        end
    end
end

% Gestion des outliers
Outlier = TABLEDATA_Int.('m6_7')(TABLEDATA_Int.('m6_7')(:,1)==5,:);
Outlier2 = TABLEDATA_Int.('m6_7')(TABLEDATA_Int.('m6_7')(:,1)==43,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m6_7')(TABLEDATA_Int.('m6_7')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m6_7')(TABLEDATA_Int.('m6_7')(:,1)==5,:)=Outlier;
TABLEDATA_Int.('m6_7')(TABLEDATA_Int.('m6_7')(:,1)==43,:)=Outlier2;

Outlier = TABLEDATA_Int.('m8_9')(TABLEDATA_Int.('m8_9')(:,1)==5,:);
Outlier2 = TABLEDATA_Int.('m8_9')(TABLEDATA_Int.('m8_9')(:,1)==43,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m8_9')(TABLEDATA_Int.('m8_9')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m8_9')(TABLEDATA_Int.('m8_9')(:,1)==5,:)=Outlier;
TABLEDATA_Int.('m8_9')(TABLEDATA_Int.('m8_9')(:,1)==43,:)=Outlier2;

Outlier = TABLEDATA_Int.('m1_2')(TABLEDATA_Int.('m1_2')(:,1)==21,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m1_2')(TABLEDATA_Int.('m1_2')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m1_2')(TABLEDATA_Int.('m1_2')(:,1)==21,:)=Outlier;

Outlier = TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==22,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==22,:)=Outlier;

Outlier = TABLEDATA_Int.('m31_32')(TABLEDATA_Int.('m31_32')(:,1)==18,:);
Outlier2 = TABLEDATA_Int.('m31_32')(TABLEDATA_Int.('m31_32')(:,1)==25,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m31_32')(TABLEDATA_Int.('m31_32')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m31_32')(TABLEDATA_Int.('m31_32')(:,1)==18,:)=Outlier;
TABLEDATA_Int.('m31_32')(TABLEDATA_Int.('m31_32')(:,1)==25,:)=Outlier2;

Outlier = TABLEDATA_Int.('m32_33')(TABLEDATA_Int.('m32_33')(:,1)==18,:);
Outlier2 = TABLEDATA_Int.('m32_33')(TABLEDATA_Int.('m32_33')(:,1)==25,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m32_33')(TABLEDATA_Int.('m32_33')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m32_33')(TABLEDATA_Int.('m32_33')(:,1)==18,:)=Outlier;
TABLEDATA_Int.('m32_33')(TABLEDATA_Int.('m32_33')(:,1)==25,:)=Outlier2;

Outlier = TABLEDATA_Int.('m15_16')(TABLEDATA_Int.('m15_16')(:,1)==25,:);
Outlier2 = TABLEDATA_Int.('m15_16')(TABLEDATA_Int.('m15_16')(:,1)==33,:);
Outlier3 = TABLEDATA_Int.('m15_16')(TABLEDATA_Int.('m15_16')(:,1)==40,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m15_16')(TABLEDATA_Int.('m15_16')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
    Outlier3(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m15_16')(TABLEDATA_Int.('m15_16')(:,1)==25,:)=Outlier;
TABLEDATA_Int.('m15_16')(TABLEDATA_Int.('m15_16')(:,1)==33,:)=Outlier2;
TABLEDATA_Int.('m15_16')(TABLEDATA_Int.('m15_16')(:,1)==40,:)=Outlier3;

Outlier = TABLEDATA_Int.('m16_17')(TABLEDATA_Int.('m16_17')(:,1)==25,:);
Outlier2 = TABLEDATA_Int.('m16_17')(TABLEDATA_Int.('m16_17')(:,1)==33,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m16_17')(TABLEDATA_Int.('m16_17')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m16_17')(TABLEDATA_Int.('m16_17')(:,1)==25,:)=Outlier;
TABLEDATA_Int.('m16_17')(TABLEDATA_Int.('m16_17')(:,1)==33,:)=Outlier2;

Outlier = TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==32,:);
Outlier2 = TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==33,:);
Outlier3 = TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==40,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
    Outlier3(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==32,:)=Outlier;
TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==33,:)=Outlier2;
TABLEDATA_Int.('m23_24')(TABLEDATA_Int.('m23_24')(:,1)==40,:)=Outlier3;

Outlier = TABLEDATA_Int.('m24_25')(TABLEDATA_Int.('m24_25')(:,1)==32,:);
Outlier2 = TABLEDATA_Int.('m24_25')(TABLEDATA_Int.('m24_25')(:,1)==33,:);
Outlier3 = TABLEDATA_Int.('m24_25')(TABLEDATA_Int.('m24_25')(:,1)==40,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m24_25')(TABLEDATA_Int.('m24_25')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
    Outlier2(i,7:12) = median(Temp(:,7:12));
    Outlier3(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m24_25')(TABLEDATA_Int.('m24_25')(:,1)==32,:)=Outlier;
TABLEDATA_Int.('m24_25')(TABLEDATA_Int.('m24_25')(:,1)==33,:)=Outlier2;
TABLEDATA_Int.('m24_25')(TABLEDATA_Int.('m24_25')(:,1)==40,:)=Outlier3;

Outlier = TABLEDATA_Int.('m40_41')(TABLEDATA_Int.('m40_41')(:,1)==48,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m40_41')(TABLEDATA_Int.('m40_41')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m40_41')(TABLEDATA_Int.('m40_41')(:,1)==48,:)=Outlier;

Outlier = TABLEDATA_Int.('m41_42')(TABLEDATA_Int.('m41_42')(:,1)==48,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('m41_42')(TABLEDATA_Int.('m41_42')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('m41_42')(TABLEDATA_Int.('m41_42')(:,1)==48,:)=Outlier;

Outlier = TABLEDATA_Int.('DeltAnt')(TABLEDATA_Int.('DeltAnt')(:,1)==50,:);
for i = 1:MaxPoint
    Temp = TABLEDATA_Int.('DeltAnt')(TABLEDATA_Int.('DeltAnt')(:,6)==i,:);
    Outlier(i,7:12) = median(Temp(:,7:12));
end
TABLEDATA_Int.('DeltAnt')(TABLEDATA_Int.('DeltAnt')(:,1)==50,:)=Outlier;

% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\TABLEDATA_Int_Li.mat','TABLEDATA_Int')

% Exclure en fonction des clusters
LowTime = [find(idx==1)];
HighTime = [find(idx==2)];
Group1 = TABLEDATA_Int;
Group2 = TABLEDATA_Int;
for iExcl = 1:length([find(idx==2)])
    for iM = 1:length(Muscles)
        Group1.(Muscles{iM})(Group1.(Muscles{iM})(:,1) == HighTime(iExcl,:),:)=[];
        Group1.(Muscles{iM})(:,13)=1;
    end
end
for iExcl = 1:length([find(idx==1)])
    for iM = 1:length(Muscles)
        Group2.(Muscles{iM})(Group2.(Muscles{iM})(:,1) == LowTime(iExcl,:),:)=[];
        Group2.(Muscles{iM})(:,13)=1;
    end
end
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\Group1_Li.mat','Group1')
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\Group2_Li.mat','Group2')


% % Exclure participant qui ont tenus plus de 6 minutes
% TABLEDATA_Int_LessSixMin = TABLEDATA_Int;
% TABLEDATA_Int_MoreSixMin = TABLEDATA_Int;
% for iExcl = 1:length(MoreSixMin)
%     for iM = 1:length(Muscles)
%         TABLEDATA_Int_LessSixMin.(Muscles{iM})(TABLEDATA_Int_LessSixMin.(Muscles{iM})(:,1) == MoreSixMin(:,iExcl),:)=[];
%         TABLEDATA_Int_LessSixMin.(Muscles{iM})(:,13)=1;
%     end
% end

% % Exclure participant qui ont tenus moins de 6 minutes
% for iExcl = 1:length(LessSixMin)
%     for iM = 1:length(Muscles)
%         TABLEDATA_Int_MoreSixMin.(Muscles{iM})(TABLEDATA_Int_MoreSixMin.(Muscles{iM})(:,1) == LessSixMin(:,iExcl),:)=[];
%         TABLEDATA_Int_MoreSixMin.(Muscles{iM})(:,13)=1;
%     end
% end
% % save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\TABLEDATA_Int_LessSixMin.mat','TABLEDATA_Int_LessSixMin')
% % save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\TABLEDATA_Int_MoreSixMin.mat','TABLEDATA_Int_MoreSixMin')


% % Lissage des données
% LissData=[];
% for iM = 1:length(Muscles)
%     LissData.(Muscles{iM}) = [];
% end
% for iM = 1:length(Muscles)
%     for iP = 1:50
%         Temp =  TABLEDATA_Int.(Muscles{iM})(TABLEDATA_Int.(Muscles{iM})(:,1)==iP,:);
%         V=[];
%         for id = 15:(length(Temp)-15)
%             V(id,:) = mean(Temp(id-14:id+14,:));
%         end
%         V(:,1:6) = Temp(1:length(V),1:6);
%         if length(V)~=0
%             V(V(:,12)==0,:)=NaN;
%         end
%         V=rmmissing(V);
%         LissData.(Muscles{iM}) = [LissData.(Muscles{iM}); V];
% %         Tri.(Muscles{iM})(Tri.(Muscles{iM})(:,1)==iParticipants,5) = 1:length(Tri.(Muscles{iM})(Tri.(Muscles{iM})(:,1)==iParticipants,5));
%     end
% end
% plot(TABLEDATA.('m31_32')(TABLEDATA.('m31_32')(:,1)==3,12))
% hold on
% plot(TABLEDATA_Int.('m31_32')(TABLEDATA_Int.('m31_32')(:,1)==3,12))
% plot(LissData.('m31_32')(LissData.('m31_32')(:,1)==1,12))

% VarNames = {'Participant','Segment','one','RPE','SegTri','Constante','Amplitude','Activity','Mobility','SampEn','MeanSE','MeanMF'};
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\VarNames_Piano.mat','VarNames')

% VarNames_TABLEDATA_Int_LessSixMin = {'Participant','Segment','PresenceRPE','RPE','TempsNorm','FrameNorm','Amplitude','Activity','Mobility','SampEn','MeanSE','MeanMF','Constante'};
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\VarNames_TABLEDATA_Int_LessSixMin.mat','VarNames_TABLEDATA_Int_LessSixMin')

% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\LissData_Piano_MoreSixMin.mat','LissData')
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\Tri_Piano_MoreSixMin.mat','Tri')

% Trie LissData pour garder les valeurs correspondant à l'instant de mesure du RPE
for iM = 1:length(Muscles)
    Tri_Group1.(Muscles{iM}) = Group1.(Muscles{iM})(Group1.(Muscles{iM})(:,3)==1,:);
    Tri_Group2.(Muscles{iM}) = Group2.(Muscles{iM})(Group2.(Muscles{iM})(:,3)==1,:);
end

for iParticipants = 1:50
    for iM = 1:length(Muscles)
        Temp1 =  Tri_Group1.(Muscles{iM})(Tri_Group1.(Muscles{iM})(:,1)==iParticipants,:);
        Temp2 =  Tri_Group2.(Muscles{iM})(Tri_Group2.(Muscles{iM})(:,1)==iParticipants,:);
        Tri_Group1.(Muscles{iM})(Tri_Group1.(Muscles{iM})(:,1)==iParticipants,6) = 1:length(Tri_Group1.(Muscles{iM})(Tri_Group1.(Muscles{iM})(:,1)==iParticipants,5));
        Tri_Group2.(Muscles{iM})(Tri_Group2.(Muscles{iM})(:,1)==iParticipants,6) = 1:length(Tri_Group2.(Muscles{iM})(Tri_Group2.(Muscles{iM})(:,1)==iParticipants,5));
    end
end
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\Tri_Group1_Li.mat','Tri_Group1')
% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\Tri_Group2_Li.mat','Tri_Group2')

