clear all;
close all;
clc;

Subjects = {'P1' 'P2' 'P3' 'P4' 'P5' 'P6' 'P7' 'P8' 'P9' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21'...
    'P22' 'P23' 'P24' 'P25' 'P26' 'P27' 'P28' 'P29' 'P30' 'P31'};

Variables = {'Effort', 'Arm pain', 'Shoulder pain', 'Neck pain', 'Cycle/s'} ;
Effort = [];
Pain = [];
inc = 1;
cd('\\10.89.24.15\q\IRSST_DavidsTea')
load('Longueur_15_cycles.mat')

for iS = 1:31
    CR100Table = [];
    cd('\\10.89.24.15\q\IRSST_DavidsTea\CR100')
    load(['CR100_' Subjects{iS} '.mat'])
    % Y_MatKernel = [Y_MatKernel; CR100Table(:,:)];
    Effort(iS,:) = CR100Table(:,1)';
    Arm_Pain(iS,:) = CR100Table(:,2)';
    Shoulder_Pain(iS,:) = CR100Table(:,3)';
    Neck_Pain(iS,:) = CR100Table(:,4)';
end
Effort_blocks = [];
ArmPain_blocks = [];
ShoulderPain_blocks = [];
NeckPain_blocks = [];
Perf_blocks = [];
MeanPain_blocks = [];
seq = 1:2:10;
for iseq = 1:length(seq)
    Effort_blocks(:,iseq) = nanmean(Effort(:,seq(iseq):seq(iseq)+1),2);
    ArmPain_blocks(:,iseq) = nanmean(Arm_Pain(:,seq(iseq):seq(iseq)+1),2);
    ShoulderPain_blocks(:,iseq) = nanmean(Shoulder_Pain(:,seq(iseq):seq(iseq)+1),2);
    NeckPain_blocks(:,iseq) = nanmean(Neck_Pain(:,seq(iseq):seq(iseq)+1),2);
    Perf_blocks(:,iseq) = nanmean(MatLengthWT(:,seq(iseq):seq(iseq)+1),2);
    MeanPain_blocks(:,iseq) = mean([ArmPain_blocks(:,iseq) ShoulderPain_blocks(:,iseq) NeckPain_blocks(:,iseq)],2);
end


% Save
cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted'])
save(['Effort_blocks.mat'],'Effort_blocks');
save(['ArmPain_blocks.mat'],'ArmPain_blocks');
save(['ShoulderPain_blocks.mat'],'ShoulderPain_blocks');
save(['NeckPain_blocks.mat'],'NeckPain_blocks');
save(['Perf_blocks.mat'],'Perf_blocks');
save(['MeanPain_blocks.mat'],'MeanPain_blocks');

