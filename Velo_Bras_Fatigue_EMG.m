clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))

Muscles = {'biceps','triceps','deltant','deltmed','deltpost','uptrap','dorsal'} ;

cd('C:\Users\p1098713\Documents\test_velo_bras');
load('RawEMG_Muscles_Etienne_30_endu.mat')
for iM = 1:length(Muscles)
    EMG(:,iM) = DataSelec.(Muscles{iM});
end
Freq = 2000;

% Preprocessing
% Freq : fréquence d'échantillonage de l'EMG
[b,a] = butter(2,2*[10 450]/Freq) ; % Parametre du filtre BP 10-400 Hz
EMGBP = filtfilt(b,a,EMG) ;
EMGBS = EMGBP;
EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)
Normalization2 = interp1(1:length(EMGBL),EMGBL,linspace(1,length(EMGBL),length(EMGBL)/2)) ;

figure;
for iM = 1:length(Muscles)
    subplot(2,5,iM);
    plot(Normalization2(:,iM));
    title(Muscles{iM})
end

% TFR
Freq = 1000;
FreqMin = 1 ;
FreqMax = 400 ;
Resolution = 1 ;... pas de 1Hz
WaveNumber = 7 ;... nombre de pics dans l'ondelette de Morlet
Args = WaveletParameters(FreqMin,FreqMax,Resolution,WaveNumber,Freq) ;
FreqRange = length(FreqMin:Resolution:FreqMax) ;
Nb_Interp_Pnts = 1 ;

% Compute TFR And MedFreq
MF_TFR = [] ;
x1 = [26193 36193];
x2 = [1240654 1250654];
Seg = [x1 ; x2];
figure;
for iM = 1:length(Muscles)
    for iSeg = 1:2
        [norm, Time, Wave_FreqS] = TimeFreqTransform(Normalization2(Seg(iSeg,1):Seg(iSeg,2),iM),Freq,Args,Nb_Interp_Pnts);
        PSD = norm.^2;
        MF_TFR.(Muscles{iM})(:,iSeg) = Compute_Median_Frequency(PSD,Wave_FreqS);
    end
    subplot(2,4,iM);
    boxplot(MF_TFR.(Muscles{iM}))
    title(Muscles{iM})
end

% Average MF_TFR
figure;
for iM = 4%1:10
    inc=1;
    for iSig = 1:1000:length(MF_TFR(:,1))-1000
        AvMF(inc,iM) = mean(MF_TFR(iSig:iSig+1000,iM));
        inc=inc+1;
    end
    % subplot(2,5,iM)
    plot(AvMF(:,iM))
end

% Compute MedFreq using fft
figure;
MF_FFT = [];
Baseline_MF = [];
for iM = 4%1:10
    inc = 1;
    for iSig = 1:1000:length(Normalization2(:,1))-1000
        [FFT, f, Mean_Freq, Med_Freq] = spec_fft(Normalization2(iSig:iSig+1000,iM),1000,0);
        MF_FFT(inc,iM) = Med_Freq;
        inc=inc+1;
    end
%     subplot(2,5,iM)
    plot(MF_FFT(:,iM))
end

figure;
for iM = 4%1:10
    inc = 1;
    for iSig = 1:5:length(MF_FFT(:,1))-5
        Av_MF_FFT(inc,iM) = mean(MF_FFT(iSig:iSig+5,iM));
        inc=inc+1;
    end
    subplot(2,5,iM)
    plot(Av_MF_FFT(:,iM))
end

[FFT, f, Mean_Freq, Med_Freq] = spec_fft(EMGBL(:,1),2000,1);

%% test threshold
envelop = abs(EMGBL);
[b,a] = butter(2,2*3/Freq,'low') ; % Parametre du filtre BP 10-400 Hz
envelop_filtered = filtfilt(b,a,envelop);
figure;
plot(envelop_filtered(:,4))
plot(sum(envelop_filtered,2))

% thresh by muscle
for iM = 1:length(Muscles)
    SigForThresh.(Muscles{iM}) = envelop_filtered(1:10000,iM);
    MedmaxVal = median(maxk(SigForThresh.(Muscles{iM}),5000));
    IQRmaxVal = iqr(maxk(SigForThresh.(Muscles{iM}),5000));
    thresh(:,iM) = MedmaxVal+1.5*IQRmaxVal;
    subplot(2,5,iM)
    plot(SigForThresh.(Muscles{iM}),'color','b')
    title(Muscles{iM})
end
hold on
window = 1:1000:length(envelop_filtered);
for iM = 1:length(Muscles)
    subplot(2,5,iM);
    line([0 length(envelop_filtered(:,iM))],[thresh(:,iM) thresh(:,iM)],'color','r')
    hold on;
    for i = 11:length(window)-10
        if sum(envelop_filtered(window(i):window(i)+9999,iM) > thresh(:,iM)) > 500
            a = double(envelop_filtered(window(i):window(i)+9999,iM) > thresh(:,iM));
            b = diff(a);
            c = find(b);
            d = diff([0;c]);
            if any(d(2:2:length(d)) > 500)
                plot(window(i):window(i)+9999,envelop_filtered(window(i):window(i)+9999,iM),'color','y');
            else
                plot(window(i):window(i)+9999,envelop_filtered(window(i):window(i)+9999,iM),'color','b');
                SigForThresh.(Muscles{iM}) = [SigForThresh.(Muscles{iM}); envelop_filtered(window(i):window(i)+999,iM)];
                MedmaxVal = median(maxk(SigForThresh.(Muscles{iM}),5000));
                IQRmaxVal = iqr(maxk(SigForThresh.(Muscles{iM}),5000));
                thresh(:,iM) = MedmaxVal+1.5*IQRmaxVal;
            end
        else
            plot(window(i):window(i)+9999,envelop_filtered(window(i):window(i)+9999,iM),'color','b');
            SigForThresh.(Muscles{iM}) = [SigForThresh.(Muscles{iM}); envelop_filtered(window(i):window(i)+999,iM)];
            MedmaxVal = median(maxk(SigForThresh.(Muscles{iM}),5000));
            IQRmaxVal = iqr(maxk(SigForThresh.(Muscles{iM}),5000));
            thresh(:,iM) = MedmaxVal+1.5*IQRmaxVal;
        end
        line([0 length(envelop_filtered(:,iM))],[thresh(:,iM) thresh(:,iM)],'color','r')
%         pause()
    end
    line([0 length(envelop_filtered(:,iM))],[thresh(:,iM) thresh(:,iM)],'color','r')
end

% one thresh for everymuscle
SigForThresh = sum(envelop_filtered(:,4),2);
MedmaxVal = median(maxk(SigForThresh,5000));
IQRmaxVal = iqr(maxk(SigForThresh,5000));
thresh = MedmaxVal+1.5*IQRmaxVal;
plot(SigForThresh,'color','b')
title(Muscles{iM})
hold on
window = 1:1000:length(envelop_filtered);

line([0 length(envelop_filtered(:,iM))],[thresh thresh],'color','r')
hold on;
for i = 11:length(window)-10
    if sum(sum(envelop_filtered(window(i):window(i)+9999,:),2) > thresh) > 500
        a = double(sum(envelop_filtered(window(i):window(i)+9999,:),2) > thresh);
        b = diff(a);
        c = find(b);
        d = diff([0;c]);
        if any(d(2:2:length(d)) > 500)
            plot(window(i):window(i)+9999,sum(envelop_filtered(window(i):window(i)+9999,:),2),'color','y');
        else
            plot(window(i):window(i)+9999,sum(envelop_filtered(window(i):window(i)+9999,:),2),'color','b');
            SigForThresh = [SigForThresh; sum(envelop_filtered(window(i):window(i)+999,:),2)];
            MedmaxVal = median(maxk(SigForThresh,5000));
            IQRmaxVal = iqr(maxk(SigForThresh,5000));
            thresh = MedmaxVal+1.5*IQRmaxVal;
        end
    else
        plot(window(i):window(i)+9999,sum(envelop_filtered(window(i):window(i)+9999,:),2),'color','b');
        SigForThresh = [SigForThresh; sum(envelop_filtered(window(i):window(i)+999,:),2)];
        MedmaxVal = median(maxk(SigForThresh,5000));
        IQRmaxVal = iqr(maxk(SigForThresh,5000));
        thresh = MedmaxVal+1.5*IQRmaxVal;
    end
    line([0 length(envelop_filtered(:,iM))],[thresh thresh],'color','r')
%         pause()
end

