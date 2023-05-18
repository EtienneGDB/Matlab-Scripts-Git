clear ;
close all ;
clc ;

% Task = {'TUG', '3min Free', '4min Free', '5min Free'};
Task = 'TUG';

cd('F:\Projet RPQ\Extracted_Data')
load('TUG_Percent_Time_Frozen.mat')
load([Task '_Percent_Time_Tremor_Walking.mat'])
load([Task '_Percent_Time_Tremor_NOWalking.mat'])
load([Task '_Tremor_Level.mat'])
load([Task '_Brady_Level.mat'])
load([Task '_Dysk_Level.mat'])

Sympt = [Percent_Time_Frozen Tremor_Level Brady_Level Dysk_Level];

