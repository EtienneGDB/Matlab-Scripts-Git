clear all;
close all;
clc;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\UnivarScatter-master'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\waterloo\Waterloo_MATLAB_Library'))

Task = 'SP';

cd('\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted')
RMSE_Effort = load(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted\RMSE_CrossVal_LOO_' Task '_effort.mat']);
RMSE_Pain = load(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted\RMSE_CrossVal_LOO_' Task '_pain.mat']);
VarExpl_Effort = load(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted\VarExplained_LOO_' Task '_effort.mat']);
VarExpl_Pain = load(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted\VarExplained_LOO_' Task '_pain.mat']);

Participants_effort = [1,2,7,8,9,10,11,14,15,20,22,23,26];
Participants_NOeffort = [3,4,5,6,12,13,16,17,18,19,21,24,25,27,28,29,30,31];

Participants_pain = [1,2,5,7,8,9,11,14,15,16,17,22,26];
Participants_NOpain = [3,4,6,10,12,13,18,19,20,21,23,24,25,27,28,29,30,31];

nanmean(RMSE_Effort.RMSE_CrossVal(Participants_effort,:))
nanstd(RMSE_Effort.RMSE_CrossVal(Participants_effort,:))

nanmean(RMSE_Pain.RMSE_CrossVal(Participants_pain,:))
nanstd(RMSE_Pain.RMSE_CrossVal(Participants_pain,:))

nanmean(VarExpl_Effort.VarExplained_LOO(Participants_effort,:))
nanstd(VarExpl_Effort.VarExplained_LOO(Participants_effort,:))

nanmean(VarExpl_Pain.VarExplained_LOO(Participants_pain,:))
nanstd(VarExpl_Pain.VarExplained_LOO(Participants_pain,:))

% if Task == 'Ha'
%     Colors = repmat([0 0.4470 0.7410],20,1);
% else
%     Colors = repmat([0.8500 0.3250 0.0980],100,1);
% end
% [p,tbl,stats] = anova1(RMSE_CrossVal(Participants_effort,:));
% multcomp = multcompare(stats);
figure; UnivarScatter(RMSE_CrossVal(Participants_effort,:), 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
% for iNcomp = 1:size(fullModel_VAREXPLAINED_IT,2)
%     iMod = find(multcomp(:,1)==iNcomp & multcomp(:,2)==iNcomp+1);
%     if multcomp(iMod,6) < 0.05
%         hold on;
%         plot(iNcomp+1, 105, '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%     end
% end
% xlabel('Number of latent variables')
% ylabel('RPE explained (%)')
% axis([0 21 0 105])




