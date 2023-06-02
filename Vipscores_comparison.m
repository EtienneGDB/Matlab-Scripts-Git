clear all;
close all;
clc;

cd('\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted')
vip_xyz = load('VIPSCORE_xyz.mat');
vip_Meanxyz = load('VIPSCORE_Meanxyz.mat');
vip_Modulexyz = load('VIPSCORE_Modulexyz.mat');

Participants_effort = [1,2,7,8,9,10,11,14,15,20,22,23,26];
Participants_pain = [1,2,5,7,8,9,11,14,15,16,17,22,26];
Participants = Participants_effort;
Participants = Participants_pain;

Data_prep = [];
Part = [];
Var = [repmat(1,numel(vip_xyz.VIPSCORE),1); repmat(2,numel(vip_Meanxyz.VIPSCORE),1); repmat(3,numel(vip_Modulexyz.VIPSCORE),1)];
Var = {};
inc = 1;
for i = 1:numel(vip_xyz.VIPSCORE)
    Var(inc,1) = {'xyz'};
    inc = inc +1;
end
for i = 1:numel(vip_Meanxyz.VIPSCORE)
    Var(inc,1) = {'Mean xyz'};
    inc = inc +1;
end
for i = 1:numel(vip_Modulexyz.VIPSCORE)
    Var(inc,1) = {'Magnitude xyz'};
    inc = inc +1;
end

for i = 1:size(vip_xyz.VIPSCORE,2)
    Data_prep = [Data_prep; vip_xyz.VIPSCORE(:,i)];
    Part = [Part; repmat(Participants(i),length(vip_xyz.VIPSCORE),1)];
end
for i = 1:size(vip_Meanxyz.VIPSCORE,2)
    Data_prep = [Data_prep; vip_Meanxyz.VIPSCORE(:,i)];
    Part = [Part; repmat(Participants(i),length(vip_Meanxyz.VIPSCORE),1)];
end
for i = 1:size(vip_Modulexyz.VIPSCORE,2)
    Data_prep = [Data_prep; vip_Modulexyz.VIPSCORE(:,i)];
    Part = [Part; repmat(Participants(i),length(vip_Modulexyz.VIPSCORE),1)];
end

Data_Stats = [Part, Data_prep];

Var(Data_Stats(:,2) < 1,:) = [];
Data_Stats(Data_Stats(:,2) < 1,:) = [];









% -----All VipScores-----
clear all;
close all;
clc;

Task = 'SP';
Variable = 'Effort';

cd('\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted')
vip = load(['VIPSCORE_All_' Variable '_' Task '.mat']);

Participants_effort = [1,2,7,8,9,10,11,14,15,20,22,23,26];
% Participants_pain = [1,2,5,7,8,9,11,14,15,16,17,22,26];
Participants = Participants_effort;
% Participants = Participants_pain;

parametre  = {'Acceleration','Angular_Velocity','Jerk','Module_Acceleration',...
    'Module_Angular_Velocity','Module_Jerk'};

indicateur = {'MedianFreq','SpectralEntropy','PowerLF','PowerHF','PowerTot', ...
    'PeakPower','PeakPower_Freq'};

indicateur_time = {'Peak', 'Mean', 'STD', 'Median', 'IQR', 'Percentile10',...
    'Percentile25', 'Percentile75', 'Percentile90'};

Seg_of_interest = {'Pelvis', 'T8', 'Head', 'RightShoulder', 'RightUpperArm',...
    'RightForeArm', 'RightHand', 'LeftShoulder', 'LeftUpperArm', 'LeftForeArm', 'LeftHand'};

Var = {};
inc = 1;
for ipar = 1:length(parametre)
    for iind = 1:length(indicateur)
        for iSeg = 1:length(Seg_of_interest)

            if ipar == 4 | ipar == 5 | ipar == 6
                Var(inc,1) = {'Magnitude xyz'};
                inc = inc + 1;

            else

                for iaxe = 1:3
                    Var(inc,1) = {'xyz'};
                    inc = inc + 1;
                end
                Var(inc,1) = {'Mean xyz'};
                inc = inc + 1;

            end
        end
    end
end
for ipar = 1:length(parametre)
    for iind = 1:length(indicateur_time)
        for iSeg = 1:length(Seg_of_interest)

            if ipar == 4 | ipar == 5 | ipar == 6
                Var(inc,1) = {'Magnitude xyz'};
                inc = inc + 1;

            else

                for iaxe = 1:3
                    Var(inc,1) = {'xyz'};
                    inc = inc + 1;
                end
                Var(inc,1) = {'Mean xyz'};
                inc = inc + 1;

            end
        end
    end
end
Var = repmat(Var,length(Participants),1);

Data_prep = [];
Part = [];
inc = 1;
for i = 1:size(vip.VIPSCORE,2)
    Data_prep = [Data_prep; vip.VIPSCORE(:,i)];
    Part = [Part; repmat(Participants(i),length(vip.VIPSCORE),1)];
end

Data_Stats = [Part, Data_prep];

Var(Data_Stats(:,2) < 1,:) = [];
Data_Stats(Data_Stats(:,2) < 1,:) = [];
