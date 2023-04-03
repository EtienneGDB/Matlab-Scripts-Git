clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))

Freq = 40;
for iS = 17:24
    cd(['J:\IRSST_Fatigue\Pointage_repetitif\Compiled_data\3-Filtre'])
    Filenames(iS) = {['P' num2str(iS) '_Pointage.mat']};
    load([Filenames{iS}])
    AccPelvis = Data.lAcc.Pelvis(:,:);
    
    cd(['J:\IRSST_Fatigue\Pointage_repetitif\Compiled_data\Decoupe'])
    FilenamesFrame(iS) = {['FramesFab_P' num2str(iS) '.mat']};
    load([FilenamesFrame{iS}])
    
    % Filtre
    [b,a] = butter(2,2*[0.5 15]/Freq) ; % Parametre du filtre BP 0.5-15 Hz
    AccPelvis_Filtered = filtfilt(b,a,AccPelvis) ;
    
    % Magnitude
    AccMag = Module(AccPelvis_Filtered(:,1:3));
    
    % RMS début/fin
    PostSway(iS,1) = rms(AccMag(round(Dist(1)*40/2000):round(Dist(11)*40/2000)));
    PostSway(iS,2) = rms(AccMag(round(Dist(end-11)*40/2000):round(Dist(end)*40/2000)));
    
%     xxx = sqrt([Data.lAcc.RHand(:,1)].^2 + [Data.lAcc.RHand(:,2)].^2);
%     figure; plot(xxx)
%     for i = 1:length(Dist)
%         line([Dist(i)*40/2000 Dist(i)*40/2000], [0 10],'color','green')
%         line([Prox(i)*40/2000 Prox(i)*40/2000], [0 10],'color','red')
%     end
%     

end




    
figure; plot(Data.lAcc.Pelvis(:,1:3))
figure; plot(Data.lAcc.Pelvis(:,4))

rawAccGR = GlobalRef(Data.rawAcc.RHand(:,1:3),[Data.rawQuat.RHand(:,2:4) Data.rawQuat.RHand(:,1)]);
figure; plot(rawAccGR(:,1:3))
lAccLR = LocalRef([Data.rawAcc.RHand(:,1:2) Data.rawAcc.RHand(:,3)+9.81],[Data.rawQuat.RHand(:,2:4) Data.rawQuat.RHand(:,1)]);
plot(lAccLR)

figure; plot(Data.lAcc.RHand(:,1:3))
lAccLR = LocalRef(Data.lAcc.RHand(:,1:3),[Data.rawQuat.RHand(:,2:4) Data.rawQuat.RHand(:,1)]);
figure; plot(lAccLR)
Acc = Module(Data.rawAcc.RHand(:,1:3))-9.81;
hold on ; plot(Acc)

plot(Acc(end-2500:end-1000,1))
hold on
plot(Vel(1:600,1))
plot((Pos(1:600,1)-mean(Pos(10:200,1))))
t=(0:length(Acc)-1)*(1/40) ;

spec_fft(Data.lAcc.Pelvis(:,4),40,1)
axis([0 20 0 5])
spec_fft(Acc(end-2500:end-1000,1),40,1)
axis([0 20 0 5])

for i = 2:length(t)
	xxx(i)=trapz(t(1:i),Acc(1:i,1));
end
plot(xxx)
spec_fft(xxx,40,1)


rms(Data.lAcc.RHand(500:1000,4))
rms(Data.lAcc.RHand(11500:12000,4))





