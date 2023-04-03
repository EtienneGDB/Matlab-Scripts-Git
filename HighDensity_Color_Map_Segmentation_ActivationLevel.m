clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\emgpartition'))

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

cd(['C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano'])
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
        M1(iM) = mean(Dat.(Muscles{iM})(:,7));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,2) = mean(temp(:,7));
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
        M2(iM) = mean(Dat.(Muscles{iM})(:,7));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,3) = mean(temp(:,7));
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
        M3(iM) = mean(Dat.(Muscles{iM})(:,7));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,4) = mean(temp(:,7));
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
        M4(iM) = mean(Dat.(Muscles{iM})(:,7));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,5) = mean(temp(:,7));
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
        M5(iM) = mean(Dat.(Muscles{iM})(:,7));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,6) = mean(temp(:,7));
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
        M6(iM) = mean(Dat.(Muscles{iM})(:,7));
        for iP = 1:length(part)
            temp = Dat.(Muscles{iM})(Dat.(Muscles{iM})(:,1)==part(iP),:);
            meas.(Muscles{iM})(iP,1) = part(iP);
            meas.(Muscles{iM})(iP,7) = mean(temp(:,7));
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

%% ANOVA effet electrode
t=[];
for iM = 42:-1:1%length(Muscles)
    t = [t; repmat(iM,length(Measurment.(Muscles{iM})(:,7)),1),Measurment.(Muscles{iM})(:,7),Measurment.(Muscles{iM})(:,8)];%,...
end
[p,tbl,stats] = anova1(t(:,2),t(:,1));
[posthoc MEANS] = multcompare(stats);
Toto(1,1) = 0;
for iM = 1:42
    temp = posthoc(posthoc(:,1)==iM,:);
    Toto(iM,(iM+1):42) = temp(:,6)';
end
imagesc(Toto)
set(gca,'ColorScale','log')
set(gca,'clim',[10^-4 0.2])
set(colorbar,'Ticks',[10^-4 10^-3 10^-2 0.05 0.2])
set(gca,'XTick',[1.5:42.5])
set(gca,'YTick',[1.5:42.5])
grid on
set(gca,'GridColor', [1 1 1])
colormap('hot')

Total = 0;
for iE = 1:18
    for iF = 19:42
        if Toto(iE,iF)<0.05 & (MEANS(iE,1)-MEANS(iF,1))<0
            Total = Total + 1;
        end
    end
end

    % Plot
    if graph == 1       
        %% Plot forearm map
        schem = [1:6;7:12;13:18;19:24;25:30;31:36;37:42]';

        figure
        XXX = rot90(rot90(M6_MedFreq))-rot90(rot90(M1_MedFreq));
        imagesc(XXX)
        imagesc(rot90(rot90(M6_MedFreq))-rot90(rot90(M1_MedFreq)),[-20 0])
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        colormap('gray')
        cd(['C:\Users\p1098713\Documents\2.Post-Doc\Articles\Article Piano Fatigue\Figures'])
        rgb = imread('Seg_Colormap_Ha.png');
        
        toto = emgpartition(XXX);
        imshow(toto)
imi = XXX
        I = rgb2gray(rgb);
        imshow(I)

        % Use the Gradient Magnitude as the Segmentation Function
        figure
        gmag = imgradient(I);
        imshow(gmag,[])
        title('Gradient Magnitude')
        
        % Mark the Foreground Objects
        % Opening
        figure
        se = strel('disk',20);
        Io = imopen(I,se);
        imshow(Io)
        title('Opening')
        
%         % opening-by-reconstruction
%         figure
%         Ie = imerode(I,se);
%         Iobr = imreconstruct(Ie,I);
%         imshow(Iobr)
%         title('Opening-by-Reconstruction')

        % Closing
        figure
        Ioc = imclose(Io,se);
        imshow(Ioc)
        title('Opening-Closing')
        
%         % opening-closing-by-reconstruction
%         figure
%         Iobrd = imdilate(Iobr,se);
%         Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
%         Iobrcbr = imcomplement(Iobrcbr);
%         imshow(Iobrcbr)
%         title('Opening-Closing by Reconstruction')

        % regional maxima
        figure
%         fgm = imregionalmax(Iobrcbr);
        fgm = imregionalmax(Ioc);
        imshow(fgm)
        title('Regional Maxima of Opening-Closing by Reconstruction')

        I2 = labeloverlay(I,fgm);
        imshow(I2)
        title('Regional Maxima Superimposed on Original Image')
        
        figure
        se2 = strel(ones(5,5));
        fgm2 = imclose(fgm,se2);
        fgm3 = imerode(fgm2,se2);
        fgm4 = bwareaopen(fgm3,20);
        I3 = labeloverlay(I,fgm4);
        imshow(I3)
        title('Modified Regional Maxima Superimposed on Original Image')
        
        % Compute Background Markers
        figure
        bw = imbinarize(Ioc);
        imshow(bw)
        title('Thresholded Opening-Closing by Reconstruction')
        
        D = bwdist(~bw);
        imshow(D,[])
title('Distance Transform of Binary Image')

D = -D;
imshow(D,[])
title('Complement of Distance Transform')

L = watershed(D);
L(~bw) = 0;

rgb = label2rgb(L,'jet',[.5 .5 .5]);
imshow(rgb)
title('Watershed Transform')

        D = bwdist(bw);
        DL = watershed(D);
        bgm = DL == 0;
        imshow(bgm)
        title('Watershed Ridge Lines')

        gmag2 = imimposemin(gmag, bgm | fgm4);
        L = watershed(gmag2);

        labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm4;
        I4 = labeloverlay(I,labels);
        imshow(I4)
        title('Markers and Object Boundaries Superimposed on Original Image')

        Lrgb = label2rgb(L,'jet','w','shuffle');
        imshow(Lrgb)
        title('Colored Watershed Label Matrix')
        
        figure
        imshow(I)
        hold on
        himage = imshow(Lrgb);
        himage.AlphaData = 0.3;
        title('Colored Labels Superimposed Transparently on Original Image')


pause()
    end
