clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

cd(['X:\Bureau\Etienne\Extracted data\Piano_Fatigue\Data\Piano'])
load(['Group1_Ha.mat'])
load(['Group2_Ha.mat'])

Measurment = [];
for iM = 1:length(Muscles)
    Measurment.(Muscles{iM}) = [];
    iM
end

graph = 0;

%% Extract data
for iGroup = 1:2
    Group_Selected = eval(['Group' num2str(iGroup) '_Ha']);... '_Ha']);
        
    part = unique(Group_Selected.(Muscles{1})(:,1));

    % Define id%
    for iM = 1:length(Muscles)
        id10 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==10,6);
        id20 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==20,6);
        id30 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==30,6);
        id40 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==40,6);
        id50 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==50,6);
        id60 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==60,6);
        id70 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==70,6);
        id80 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==80,6);
        id90 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==90,6);
        id100 = Group_Selected.(Muscles{iM})(floor(Group_Selected.(Muscles{iM})(:,5))==100,6);
        for j = length(id10):-1:2
            if id10(j)==id10(j-1)+1
                id10(j)=NaN;
            end
        end
        for j = length(id10):-1:2
            if id10(j)==id10(j-1)+1
                id10(j)=NaN;
            end
        end
        for j = length(id20):-1:2
            if id20(j)==id20(j-1)+1
                id20(j)=NaN;
            end
        end
        for j = length(id30):-1:2
            if id30(j)==id30(j-1)+1
                id30(j)=NaN;
            end
        end
        for j = length(id40):-1:2
            if id40(j)==id40(j-1)+1
                id40(j)=NaN;
            end
        end
        for j = length(id50):-1:2
            if id50(j)==id50(j-1)+1
                id50(j)=NaN;
            end
        end
        for j = length(id60):-1:2
            if id60(j)==id60(j-1)+1
                id60(j)=NaN;
            end
        end
        for j = length(id70):-1:2
            if id70(j)==id70(j-1)+1
                id70(j)=NaN;
            end
        end
        for j = length(id80):-1:2
            if id80(j)==id80(j-1)+1
                id80(j)=NaN;
            end
        end
        for j = length(id90):-1:2
            if id90(j)==id90(j-1)+1
                id90(j)=NaN;
            end
        end    
        id10 = round(mean(id10(~isnan(id10))));
        id20 = round(mean(id20(~isnan(id20))));
        id30 = round(mean(id30(~isnan(id30))));
        id40 = round(mean(id40(~isnan(id40))));
        id50 = round(mean(id50(~isnan(id50))));
        id60 = round(mean(id60(~isnan(id60))));
        id70 = round(mean(id70(~isnan(id70))));
        id80 = round(mean(id80(~isnan(id80))));
        id90 = round(mean(id90(~isnan(id90))));
        id100 = round(mean(id100(~isnan(id100))));    
    end

    % Ref_0-10 ; 10-20 ; 30-40 ; 50-60 ; 70-80 ; 90-100
    % Color Map 0-10%
    Dat = [];
    meas = [];
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = [];
        meas.(Muscles{iM}) = [];
    end
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = Group_Selected.(Muscles{iM})(Group_Selected.(Muscles{iM})(:,6) < id10,:);
        M1(iM) = mean(Dat.(Muscles{iM})(:,12));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,2) = mean(temp(:,12));
        end
    end
    seq = 1:6:length(M1);
    M1_MedFreq = [];
    for i = 1:7
        M1_MedFreq(:,i) = M1(seq(i):seq(i)+5)';
    end
    for i = 8
        M1_MedFreq_Arm = M1(seq(i):seq(i)+4)';
        M1_MedFreq_Arm = flip(M1_MedFreq_Arm);
    end

    % Color Map 10-20%
    Dat = [];
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = [];
    end
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = Group_Selected.(Muscles{iM})(Group_Selected.(Muscles{iM})(:,6)>id10 & Group_Selected.(Muscles{iM})(:,6)<id20,:);
        M2(iM) = mean(Dat.(Muscles{iM})(:,12));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,3) = mean(temp(:,12));
        end
    end
    seq = 1:6:length(M2);
    M2_MedFreq = [];
    for i = 1:7
        M2_MedFreq(:,i) = M2(seq(i):seq(i)+5)';
    end
    for i = 8
        M2_MedFreq_Arm = M2(seq(i):seq(i)+4)';
        M2_MedFreq_Arm = flip(M2_MedFreq_Arm);
    end

    % Color Map 30-40%
    Dat = [];
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = [];
    end
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = Group_Selected.(Muscles{iM})(Group_Selected.(Muscles{iM})(:,6)>id30 & Group_Selected.(Muscles{iM})(:,6)<id40,:);
        M3(iM) = mean(Dat.(Muscles{iM})(:,12));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,4) = mean(temp(:,12));
        end
    end
    seq = 1:6:length(M3);
    M3_MedFreq = [];
    for i = 1:7
        M3_MedFreq(:,i) = M3(seq(i):seq(i)+5)';
    end
    for i = 8
        M3_MedFreq_Arm = M3(seq(i):seq(i)+4)';
        M3_MedFreq_Arm = flip(M3_MedFreq_Arm);
    end

    % Color Map 50-60%
    Dat = [];
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = [];
    end
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = Group_Selected.(Muscles{iM})(Group_Selected.(Muscles{iM})(:,6)>id50 & Group_Selected.(Muscles{iM})(:,6)<id60,:);
        M4(iM) = mean(Dat.(Muscles{iM})(:,12));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,5) = mean(temp(:,12));
        end
    end
    seq = 1:6:length(M4);
    M4_MedFreq = [];
    for i = 1:7
        M4_MedFreq(:,i) = M4(seq(i):seq(i)+5)';
    end
    for i = 8
        M4_MedFreq_Arm = M4(seq(i):seq(i)+4)';
        M4_MedFreq_Arm = flip(M4_MedFreq_Arm);
    end

    % Color Map 70-80%
    Dat = [];
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = [];
    end
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = Group_Selected.(Muscles{iM})(Group_Selected.(Muscles{iM})(:,6)>id70 & Group_Selected.(Muscles{iM})(:,6)<id80,:);
        M5(iM) = mean(Dat.(Muscles{iM})(:,12));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,6) = mean(temp(:,12));
        end
    end
    seq = 1:6:length(M5);
    M5_MedFreq = [];
    for i = 1:7
        M5_MedFreq(:,i) = M5(seq(i):seq(i)+5)';
    end
    for i = 8
        M5_MedFreq_Arm = M5(seq(i):seq(i)+4)';
        M5_MedFreq_Arm = flip(M5_MedFreq_Arm);
    end

    % Color Map 90-100%
    Dat = [];
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = [];
    end
    for iM = 1:length(Muscles)
        Dat.(Muscles{iM}) = Group_Selected.(Muscles{iM})(Group_Selected.(Muscles{iM})(:,6)>id90,:);
        M6(iM) = mean(Dat.(Muscles{iM})(:,12));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,7) = mean(temp(:,12));
        end
        meas.(Muscles{iM})(:,8) = iGroup;
    end
    seq = 1:6:length(M6);
    M6_MedFreq = [];
    for i = 1:7
        M6_MedFreq(:,i) = M6(seq(i):seq(i)+5)';
    end
    for i = 8
        M6_MedFreq_Arm = M6(seq(i):seq(i)+4)';
        M6_MedFreq_Arm = flip(M6_MedFreq_Arm);
    end

    % Concatenate meas data
    for iM = 1:length(Muscles)
        Measurment.(Muscles{iM}) = [Measurment.(Muscles{iM}) ; meas.(Muscles{iM})];
    end
end

%% ANOVA repeated-measures
RepeatedMeasuresResults = [];
GroupResults = [];
PostHocResults = [];
PostHocResults2 = [];
for iM = 1:length(Muscles)
    t = table(Measurment.(Muscles{iM})(:,1),Measurment.(Muscles{iM})(:,2),Measurment.(Muscles{iM})(:,3),...
        Measurment.(Muscles{iM})(:,4),Measurment.(Muscles{iM})(:,5),Measurment.(Muscles{iM})(:,6),...
        Measurment.(Muscles{iM})(:,7),Measurment.(Muscles{iM})(:,8),'VariableNames',{'Part','MFRef','MF10_20',...
        'MF30_40','MF50_60','MF70_80','MF90_100','Group'});
    t.Group = categorical(t.Group);
    Meas = table([1 2 3 4 5 6]','VariableNames',{'Measurement'});
%     Meas2 = Meas;
%     Meas2.Measurement = categorical(Meas2.Measurement);
%     rm = fitrm(t,'MFRef-MF90_100~Part','WithinDesign',Meas);
    rm = fitrm(t,'MFRef-MF90_100~Group','WithinDesign',Meas);
%     ranovatbl = ranova(rm);
    anovatbl = anova(rm);
    posthoc = multcompare(rm,'Measurement');
    RepeatedMeasuresResults.(Muscles{iM}) = ranovatbl;
    GroupResults.(Muscles{iM}) = anovatbl;
    PostHocResults.(Muscles{iM}) = posthoc([1 7 13 19 25],:);
    PostHocResults2.(Muscles{iM}) = posthoc([2 8 14 20],:);
end
computeCohen_d(table2array(t(1:end,2)),table2array(t(1:end,7)),'paired')

% Anova group pour toutes conditions
% DetailGroupResults = [];
% for iM = 1:length(Muscles)
%     for iTime = 1:6
%         t = table(Measurment.(Muscles{iM})(:,1),Measurment.(Muscles{iM})(:,2),Measurment.(Muscles{iM})(:,3),...
%             Measurment.(Muscles{iM})(:,4),Measurment.(Muscles{iM})(:,5),Measurment.(Muscles{iM})(:,6),...
%             Measurment.(Muscles{iM})(:,7),Measurment.(Muscles{iM})(:,8),'VariableNames',{'Part','MFRef','MF10_20',...
%             'MF30_40','MF50_60','MF70_80','MF90_100','Group'});
%         group = {'G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2'};
%         p=anova1(table2array(t(:,iTime+1))',group,'off');
%         DetailGroupResults.(Muscles{iM})(iTime,1) = p;
%     end
% end

Correction = [];
Correction2 = [];
Correction_interaction = [];
Correction_Group = [];
% Correction_DetailGroupResults = [];
for iM = 1:length(Muscles)
    Correction = [Correction; table2array(PostHocResults.(Muscles{iM})(:,5))];
    Correction2 = [Correction2; table2array(PostHocResults2.(Muscles{iM})(:,5))];
    Correction_interaction = [Correction_interaction; table2array(RepeatedMeasuresResults.(Muscles{iM})(1,5))]; %(1,5)_Time; (2,5)_Interaction
    Correction_Group = [Correction_Group; table2array(GroupResults.(Muscles{iM})(2,7))];
%     Correction_DetailGroupResults = [Correction_DetailGroupResults; DetailGroupResults.(Muscles{iM})];
end
[FDR,Qtime] = mafdr(Correction);
PCorrected = [];
seq = 1:5:length(Qtime);
for iM = 1:length(Muscles)
    PCorrected.(Muscles{iM}) = Qtime(seq(iM):seq(iM)+4);
end
[FDR,Qtime2] = mafdr(Correction2);
PCorrected2 = [];
seq = 1:4:length(Qtime2);
for iM = 1:length(Muscles)
    PCorrected2.(Muscles{iM}) = Qtime2(seq(iM):seq(iM)+3);
end
[Qinter] = mafdr(Correction_interaction,'BHFDR',true);
Qinter;
[Qgroup] = mafdr(Correction_Group,'BHFDR',true);
Qgroup;
% [Qdetgroup] = mafdr(Correction_DetailGroupResults);
% Qdetgroup;

% Matrice Q-values et plot
% interaction
seq = 1:6:length(M1);
Qplot = [];
for i = 1:7
    Qplot(:,i) = Qinter(seq(i):seq(i)+5)';
end
Qplot = rot90(Qplot,2);
figure
imagesc(Qplot)
set(gca,'ColorScale','log')
set(gca,'clim',[10^-4 0.2])
set(colorbar,'Ticks',[10^-4 10^-3 10^-2 0.05 0.2])
set(gca,'XTick',[1.5 2.5 3.5 4.5 5.5 6.5])
set(gca,'YTick',[1.5 2.5 3.5 4.5 5.5 6.5])
grid on
set(gca,'GridColor', [1 1 1])
colormap('hot')

% % group
seq = 1:6:length(M1);
Qplot = [];
for i = 1:7
    Qplot(:,i) = Qgroup(seq(i):seq(i)+5)';
end
Qplot = rot90(Qplot,2);
figure
imagesc(Qplot)
set(gca,'ColorScale','log')
set(gca,'clim',[3E-52 0.05])
set(colorbar,'Ticks',[logspace(-60,-1,6)])
colormap('hot')

% time
seq = 1:6:length(M1);
min(Qtime)
for itime = 1:5
    Qplot = [];
    temp = Qtime(itime:5:length(Qtime));
    for i = 1:7
        Qplot(:,i) = temp(seq(i):seq(i)+5)';
    end
    Qplot = rot90(Qplot,2);
    subplot(4,5,itime)
    imagesc(Qplot)
    set(gca,'ColorScale','log')
    set(gca,'clim',[10^-6 0.1])
    set(colorbar,'Ticks',[logspace(-8,-1,8)])
    colormap('hot')
    set(gca,'XTick',[1.5 2.5 3.5 4.5 5.5 6.5])
    set(gca,'YTick',[1.5 2.5 3.5 4.5 5.5 6.5])
    grid on
    set(gca,'GridColor', [0 0 0])
%     pause()
end


    % Plot
    if graph == 1       
        %% Plot forearm map
        schem = [1:6;7:12;13:18;19:24;25:30;31:36;37:42]';

        figure
        imagesc(M2_MedFreq-M1_MedFreq,[-20 0])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(1) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
        imagesc(M3_MedFreq-M1_MedFreq,[-20 0])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(2) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(1) < 0.05 & PCorrected.(Muscles{iM})(1) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
        imagesc(M4_MedFreq-M1_MedFreq,[-20 0])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(3) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(2) < 0.05 & PCorrected.(Muscles{iM})(2) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
        imagesc(M5_MedFreq-M1_MedFreq,[-20 0])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(4) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(3) < 0.05 & PCorrected.(Muscles{iM})(3) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
%         imagesc(rot90(rot90(M6_MedFreq))-rot90(rot90(M1_MedFreq)),[-20 0])
        subplot(2,1,1);
        imagesc(M6_MedFreq-M1_MedFreq,[-20 0])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(5) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(4) < 0.05 & PCorrected.(Muscles{iM})(4) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end
        view([-180 90])



        %% Plot arm map
        schem = [47:-1:43]';

        figure
        imagesc(M2_MedFreq_Arm-M1_MedFreq_Arm,[-20 0])
        hold on
        for iM = 43:length(Muscles)
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(1) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
        imagesc(M3_MedFreq_Arm-M1_MedFreq_Arm,[-20 0])
        hold on
        for iM = 43:length(Muscles)
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(2) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(1) < 0.05 & PCorrected.(Muscles{iM})(1) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
        imagesc(M4_MedFreq_Arm-M1_MedFreq_Arm,[-20 0])
        hold on
        for iM = 43:length(Muscles)
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(3) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(2) < 0.05 & PCorrected.(Muscles{iM})(2) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
        imagesc(M5_MedFreq_Arm-M1_MedFreq_Arm,[-20 0])
        hold on
        for iM = 43:length(Muscles)
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(4) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(3) < 0.05 & PCorrected.(Muscles{iM})(3) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end

        figure
        imagesc(M6_MedFreq_Arm-M1_MedFreq_Arm,[-20 0])
        hold on
        for iM = 43:length(Muscles)
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(5) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(4) < 0.05 & PCorrected.(Muscles{iM})(4) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end
        pause()
    end
