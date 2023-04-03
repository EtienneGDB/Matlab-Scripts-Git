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

% Load TableData_LMM
load('H:\Bureau\Etienne\Extracted data\TableData_LMM.mat')

% Grosses boîtes
DAT1 = TableData;
for iM = 1:length(Muscles)
   DAT1.(Muscles{iM}) = DAT1.(Muscles{iM})(DAT1.(Muscles{iM})(:,3)==1,:);...Grosses boîtes
%    DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,3)==2,:);...Moyennes boîtes
%    DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,3)==3,:);...Petites boîtes
    
    DAT1.(Muscles{iM}) = DAT1.(Muscles{iM})(DAT1.(Muscles{iM})(:,5)<5,:);...4 Premiers essais
end

% Order DAT per trial including Up/Down
DAT_Trial1 = DAT1;
for iM = 1:length(Muscles)
    varInc = 1;
    DAT_Trial1.(Muscles{iM}) = NaN;
    for iP = 1:31
        Vp = 1;
        for iT = 1:4
            for iUD = 2:-1:1
                V = DAT1.(Muscles{iM})(DAT1.(Muscles{iM})(:,1)==iP & DAT1.(Muscles{iM})(:,5)==iT & DAT1.(Muscles{iM})(:,4)==iUD,:);
                for iCol = 1:14
                    DAT_Trial1.(Muscles{iM})(varInc,iCol) = mean(V(:,iCol));
                    DAT_Trial1.(Muscles{iM})(varInc,6) = Vp;
                end
                Vp=Vp+1;
                varInc=varInc+1;
            end
        end
    end
end
for iM = 1:length(Muscles)
    DAT_Trial1.(Muscles{iM}) = DAT_Trial1.(Muscles{iM})(isfinite(DAT_Trial1.(Muscles{iM})(:,1)),:);
end

% Moyennes boîtes
DAT2 = TableData;
for iM = 1:length(Muscles)
%   DAT1.(Muscles{iM}) = DAT1.(Muscles{iM})(DAT1.(Muscles{iM})(:,3)==1,:);...Grosses boîtes
    DAT2.(Muscles{iM}) = DAT2.(Muscles{iM})(DAT2.(Muscles{iM})(:,3)==2,:);...Moyennes boîtes
%    DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,3)==3,:);...Petites boîtes
    
    DAT2.(Muscles{iM}) = DAT2.(Muscles{iM})(DAT2.(Muscles{iM})(:,5)<5,:);...4 Premiers essais
end

% Order DAT per trial including Up/Down
DAT_Trial2 = DAT2;
for iM = 1:length(Muscles)
    varInc = 1;
    DAT_Trial2.(Muscles{iM}) = NaN;
    for iP = 1:31
        Vp = 1;
        for iT = 1:4
            for iUD = 2:-1:1
                V = DAT2.(Muscles{iM})(DAT2.(Muscles{iM})(:,1)==iP & DAT2.(Muscles{iM})(:,5)==iT & DAT2.(Muscles{iM})(:,4)==iUD,:);
                for iCol = 1:14
                    DAT_Trial2.(Muscles{iM})(varInc,iCol) = mean(V(:,iCol));
                    DAT_Trial2.(Muscles{iM})(varInc,6) = Vp;
                end
                Vp=Vp+1;
                varInc=varInc+1;
            end
        end
    end
end
for iM = 1:length(Muscles)
    DAT_Trial2.(Muscles{iM}) = DAT_Trial2.(Muscles{iM})(isfinite(DAT_Trial2.(Muscles{iM})(:,1)),:);
end

% Petites boîtes
DAT3 = TableData;
for iM = 1:length(Muscles)
%   DAT1.(Muscles{iM}) = DAT1.(Muscles{iM})(DAT1.(Muscles{iM})(:,3)==1,:);...Grosses boîtes
%    DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,3)==2,:);...Moyennes boîtes
    DAT3.(Muscles{iM}) = DAT3.(Muscles{iM})(DAT3.(Muscles{iM})(:,3)==3,:);...Petites boîtes
    
    DAT3.(Muscles{iM}) = DAT3.(Muscles{iM})(DAT3.(Muscles{iM})(:,5)<5,:);...4 Premiers essais
end

% Order DAT per trial including Up/Down
DAT_Trial3 = DAT3;
for iM = 1:length(Muscles)
    varInc = 1;
    DAT_Trial3.(Muscles{iM}) = NaN;
    for iP = 1:31
        Vp = 1;
        for iT = 1:4
            for iUD = 2:-1:1
                V = DAT3.(Muscles{iM})(DAT3.(Muscles{iM})(:,1)==iP & DAT3.(Muscles{iM})(:,5)==iT & DAT3.(Muscles{iM})(:,4)==iUD,:);
                for iCol = 1:14
                    DAT_Trial3.(Muscles{iM})(varInc,iCol) = mean(V(:,iCol));
                    DAT_Trial3.(Muscles{iM})(varInc,6) = Vp;
                end
                Vp=Vp+1;
                varInc=varInc+1;
            end
        end
    end
end
for iM = 1:length(Muscles)
    DAT_Trial3.(Muscles{iM}) = DAT_Trial3.(Muscles{iM})(isfinite(DAT_Trial3.(Muscles{iM})(:,1)),:);
end

% Concatenate DAT1,2,3
DAT = DAT_Trial1;
for iM=1:length(Muscles)
    DAT.(Muscles{iM}) = [DAT_Trial1.(Muscles{iM}); DAT_Trial2.(Muscles{iM}); DAT_Trial3.(Muscles{iM})];
end

% Ordonne les essais
% Order DAT per trial including Up/Down
DAT_Trial = DAT;
for iM = 1:length(Muscles)
    W=[];
    DAT_Trial.(Muscles{iM})=NaN;
  for iP = 1:31
    for iT = 1:4
        V=DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,1)==iP & DAT.(Muscles{iM})(:,5)==iT,:);
        while size(V,1)<6
            V(size(V,1)+1,:)=NaN;
        end
        W=[W;V(3,:);V(5,:);V(6,:);V(4,:);V(2,:);V(1,:)];
    end
  end
  DAT_Trial.(Muscles{iM})=W;
end
for iM = 1:length(Muscles)
    DAT_Trial.(Muscles{iM}) = DAT_Trial.(Muscles{iM})(isfinite(DAT_Trial.(Muscles{iM})(:,1)),:);
end
DAT = DAT_Trial;
for iM = 1:length(Muscles)
  for iP = 1:31
    Nobs=size(DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,1)==iP,:));
    if Nobs~=0
        DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,1)==iP,6)=[1:Nobs];
    end
  end
end

% Réduire les données à 8 lignes
for iM = 1:length(Muscles)
    W=[];
    DAT_Trial.(Muscles{iM})=NaN;
    for iP = 1:31
        V=DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,1)==iP,:);
        for ig = 1:3:size(V,1)
            W = [W;V(ig,:)];
            if size(V,1)<22 & V(ig,:)==19
                W = [W;V(size(V,1),:)];
            end
        end
    end
    DAT_Trial.(Muscles{iM})=W;
end
DAT=DAT_Trial;

% Load RPE Data
cd(['H:\Projet_ExpertsNovices\excel'])
[num, txt, raw] = xlsread('rpe.xlsx');
RPE = num(2:59,6:14);
Sbj = transpose(txt(3:2:59,2));
Sbj{2,1} = [];
Sbj1 = [];
for i = 1:29
    for j = 1:31
        if contains(Subjects{1,j},Sbj{1,i})
            Sbj1(i) = j;
        end
    end
end
Sbj1 = Sbj1';
Sbj1 = repmat(Sbj1,1,2)';
Sbj1 = Sbj1(:)';

% Table RPE
TableRPE = [transpose(Sbj1) RPE];
TableRPE(TableRPE(:,1)==0,:) = [];
TableRPE = sortrows(TableRPE,1);
TableRPE(11:12,:) = [];...2x sujet 6
TableRPE(27:28,:) = [];...2x sujet 16

RPE_Shoulder=[];
varInc=1;
for i = 2:2:52 %1:2:51 pour shoulder & 2:2:52 pour overall
    for j = 2:9
        RPE_Shoulder(varInc,1) = TableRPE(i,1);
        RPE_Shoulder(varInc,2) = TableRPE(i,j);
        varInc=varInc+1;
    end
end
% RPE_Shoulder(RPE_Shoulder(:,1)==10,:)=[];

% Spearman correlations
% iMus = 7;
% DAT_Shoulder = DAT_Trial;
% Vemg = [7:11 13];
% 
% for iMus = 1:length(Muscles)
%     RPE_SH = RPE_Shoulder;
%     Del = [5 11 14 18 25];...missing RPE data
%     DAT_Shoulder.(Muscles{iMus})(ismember(DAT_Shoulder.(Muscles{iMus})(:,1),Del),:)=[];
% 
%     X1 = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
%     X2 = unique(RPE_SH(:,1));
%     N1 = numel(X1);
%     N2 = numel(X2);
%     count1 = zeros(N1,1);
%     count2 = zeros(N2,1);
%     for k = 1:N1
%         count1(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X1(k));
%     end
%     CT1 = [ X1(:) count1 ];
%     for k = 1:N2
%         count2(k) = sum(RPE_SH(:,1)==X2(k));
%     end
%     CT2 = [ X2(:) count2 ];
%     Vsync = setdiff(CT2,CT1);
%     if length(Vsync)~=0
%         for iVsync = 1:length(Vsync)
%             DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==Vsync(iVsync),:)=[];
%             RPE_SH(RPE_SH(:,1)==Vsync(iVsync),:)=[];
%         end
%     end
%     X = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
%     N = numel(X);
%     count = zeros(N,1);
%     for k = 1:N
%         count(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X(k));
%     end
%     CT = [ X(:) count ];
%     
%     for iP = 1:length(CT)
%         if CT(iP,2)<8
%             DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==CT(iP,1),:)=[];
%             RPE_SH(RPE_SH(:,1)==CT(iP,1),:)=[];
%         end
%     end
%     % Spearman correlations
%     for iVemg = 1:6
%         [rho,pval] = corr(RPE_SH(:,2),DAT_Shoulder.(Muscles{iMus})(:,Vemg(iVemg)),'Type','Spearman');
%         RHO(iMus,iVemg) = rho;
%         PVAL(iMus,iVemg) = pval;
%         plot(RPE_SH(:,2),DAT_Shoulder.(Muscles{iMus})(:,Vemg(iVemg)),'.')
%         pause
%     end
%     % Regression
%     [b,bint,r,rint,stats] = regress(RPE_SH(:,2),DAT_Shoulder.(Muscles{iMus})(:,Vemg));
%     aaa = DAT_Shoulder.(Muscles{iMus})(:,[7:11 13]);
%     for i = 1:length(aaa)
%         RPE_Shoulder(i,3) = aaa(i,1)*b(1)+aaa(i,2)*b(2)+aaa(i,3)*b(3)+aaa(i,4)*b(4)+aaa(i,5)*b(5)+aaa(i,6)*b(6);
%     end
%     plot(RPE_Shoulder(:,3),RPE_Shoulder(:,2),'.')
%     pause
% end

% [A,B,r,U,V,stats] = canoncorr(DAT_Shoulder.(Muscles{iMus})(:,[7:11 13]),RPE_Shoulder(:,2));

% Regression
DAT_Shoulder = DAT;
DAT_SH = [];
Vemg = [7:11 13];

% Armonisation DAT_Shoulder
% Étape 1-3
MinSize = 1000;
for iMus = 1:length(Muscles)
    Sz = size(DAT_Shoulder.(Muscles{iMus}));
    if Sz(1) < MinSize
        MinSize = Sz(1);
        IdRef = iMus;
    end
end
XRef = unique(DAT_Shoulder.(Muscles{IdRef})(:,1));
NRef = numel(XRef);
countRef = zeros(NRef,1);
for k = 1:NRef
    countRef(k) = sum(DAT_Shoulder.(Muscles{IdRef})(:,1)==XRef(k));
end
CTRef = [ XRef(:) countRef ];
    
for iMus = 1:length(Muscles)
    X2 = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
    N2 = numel(X2);
    count2 = zeros(N2,1);
    for k = 1:N2
        count2(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X2(k));
    end
    CT2 = [ X2(:) count2 ];
    Vsync = setdiff(CT2,CTRef);
    if length(Vsync)~=0
        for iVsync = 1:length(Vsync)
            DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==Vsync(iVsync),:)=[];
        end
    end
end

% Étape 2-4
for iMus = 1:length(Muscles)
    RPE_SH = RPE_Shoulder;
    Del = [5 11 14 18 25];...missing RPE data
    DAT_Shoulder.(Muscles{iMus})(ismember(DAT_Shoulder.(Muscles{iMus})(:,1),Del),:)=[];

    X1 = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
    X2 = unique(RPE_SH(:,1));
    N1 = numel(X1);
    N2 = numel(X2);
    count1 = zeros(N1,1);
    count2 = zeros(N2,1);
    for k = 1:N1
        count1(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X1(k));
    end
    CT1 = [ X1(:) count1 ];
    for k = 1:N2
        count2(k) = sum(RPE_SH(:,1)==X2(k));
    end
    CT2 = [ X2(:) count2 ];
    Vsync = setdiff(CT2,CT1);
    if length(Vsync)~=0
        for iVsync = 1:length(Vsync)
            DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==Vsync(iVsync),:)=[];
            RPE_SH(RPE_SH(:,1)==Vsync(iVsync),:)=[];
        end
    end
    X = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
    N = numel(X);
    count = zeros(N,1);
    for k = 1:N
        count(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X(k));
    end
    CT = [ X(:) count ];
    
    for iP = 1:length(CT)
        if CT(iP,2)<8
            DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==CT(iP,1),:)=[];
            RPE_SH(RPE_SH(:,1)==CT(iP,1),:)=[];
        end
    end
end

% Suite
DAT_SH = [];
for iMus = 1:length(Muscles)
    DAT_SH = [DAT_SH DAT_Shoulder.(Muscles{iMus})(:,[11])];
end

[b,bint,r,rint,stats] = regress(RPE_SH(:,2),DAT_SH(:,:));
stats
[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(DAT_SH(:,:),RPE_SH(:,2),'penter',0.05,'premove',0.1);

mdl = fitlm(DAT_SH(:,:),RPE_SH(:,2))
mdl = stepwiselm(DAT_SH(:,:),RPE_SH(:,2),'constant','Upper','linear','PEnter',0.05)

% % Essai
% RPE_SH_Modif = [];
% DAT_Shoulder_Modif = [];
% for iMus = 1:10
%     RPE_SH_Modif = [RPE_SH_Modif ; RPE_SH(:,2)];
%     DAT_Shoulder_Modif = [DAT_Shoulder_Modif ; DAT_Shoulder.(Muscles{iMus})(:,[7 8 9 10 11 13])];
% end
% 
% [b,bint,r,rint,stats] = regress(RPE_SH_Modif(:),DAT_Shoulder_Modif(:,5));
% stats
% [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(DAT_Shoulder_Modif(:,5),RPE_SH_Modif(:),'penter',0.05,'premove',0.1);

% Names table
RepVar = repmat({'Amp','Act','Mob','SampEn','SpecEn','MedFreq'},1,10);
k=1;
Names=[];
for i = 1:10
    for j = 1:6
        Names{1,k} = strcat(RepVar{j},'_',Muscles{i});
        k=k+1;
    end
end
DAT_SHN = array2table(DAT_SH,'variableNames',Names);


aaa = DAT_Shoulder.(Muscles{iMus})(:,[7:11 13]);
for i = 1:length(aaa)
    RPE_Shoulder(i,3) = aaa(i,1)*b(1)+aaa(i,2)*b(2)+aaa(i,3)*b(3)+aaa(i,4)*b(4)+aaa(i,5)*b(5)+aaa(i,6)*b(6);
end
plot(RPE_Shoulder(:,3),RPE_Shoulder(:,2),'.')



