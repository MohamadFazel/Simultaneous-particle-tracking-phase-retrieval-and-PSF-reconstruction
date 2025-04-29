%This script showcases an example of simultaneous singl eparticle tracking,
%phase retrieval and PSF reconstruction using synthetic data. The data is 
%generated assuming a single particle with aberrated PSF and a 2 planes 
%setup, 100 frames from each plane and frmes of 32x32 pixels. 
%
%INPUTS
% Data: input data is a 4D array with the 1st and 2nd dimensions are frame
%       sizes, the 3rd dimension is number of frames from each plane, and
%       the fourth dimension is number of planes.
% PSFstruct: structure containing experimental parameters
%    Na: numerical aperture
%    N: refractive index
%    Lambda: fluorescence emission wavelength (nm)
%    PixelSize: pixel size (nm)
%    NPix: length of frames assumed to be square (pixel)
%    DelX: spatial (X,Y,Z) shifts across planes (nm)
%    ZImage: Zernike polynomials
% BNP: structure containing math framework parameters
%    NJumps: number of samples from the posterior (chain length)
%    T_A: GP prior parameter on pupil magnitude (Default: 1)
%    L_A: GP prior parameter on pupil magnitude related to magnitude correlation
%    T_Phi: GP prior parameter on pupil phase (Default: 1)
%    L_Phi: GP prior parameter on pupil phase related to phase correlation
%    Dt: frame exposure time (ms)
%    SubPixel: subpixel resolution used in model calculations (Default: 1)
% Xstart: initial X-trajectory
% Ystart: initial Y-trajectory
%
%OUTPUTS
% Chain: chain of samples containing the following parameters
%    Mag: pupil magnitude
%    Phase: pupil phase
%    PSFstack: PSF model reconstructed from the last state of the chain
%    Bg: offset background (photons)
%    I: eitter intensity (photons)
%    X: X-trajectory (nm)
%    Y: Y-trajectory (nm)
%    Z: Z-trajectory (nm)
%    D: diffusion constant (nm^2/ms)
%    LogLike: log-likelihood
%    LogPost: log-posterior
%    SigConv: Sigma of the Gaussian convolved with the PSF
%
% Author:
%   Mohamadreza Fazel, Presse lab, 2025
%
%% add path to the codes

clear 
close all
addpath('Codes')

%% load data (4D array of simulated 100 frames of 32x32 pixels assuming 2 planes)

load('Data/Data_Simulated.mat')

%% all the experimental parameters are fields in "PSFstruct" structure

PSFstruct.Na = 1.30; %Numerical apperture
PSFstruct.N = 1.4; %Refractive index
PSFstruct.Lambda = 690; %Emission wavelength (nm)
PSFstruct.PixelSize = 100; %Spatial pixel size (nm)
PSFstruct.NPix = 32; %frame size assumed to be sqaure (pixel)
PSFstruct.DelX = [0 0 400]; %[X, Y, Z] shifts from first (refernce) to the second plane

%generating the Zernike polynomials associated to the phase due to particle location
NZernike = 5; %the first 4 Zernike polynomial (offset, X-tilt, Y-tilt, defocus)
PSFstruct = genZernikeIm(PSFstruct,NZernike);

%% all the math framework parameters are fields in "BNP" structure

BNP.NJump = 25000; %number of samples to be taken from posterior
BNP.T_A = 1; %GP prior parameter associated to the magnitude (see the SI)
%Larger values result in smoother magnitudes
BNP.L_A = 4/(PSFstruct.PixelSize*PSFstruct.NPix); 
BNP.T_Phi = 1; %GP prior parameters associated to the phase (see the SI)
%Larger values result in smoother phases
BNP.L_Phi = 4/(PSFstruct.PixelSize*PSFstruct.NPix);
BNP.Dt = 20; %exposure time (ms)
BNP.SubPixel = 4; %subpixel resolution used in constructing the frame models

%Initializing the X and Y trajectories by finding the brightest pixel in
%each frame.
NFrame = size(Data,3);
Xstart = zeros(1,NFrame);
Ystart = zeros(1,NFrame);
for ii = 1:NFrame
    tData = Data(:,:,ii,1);
    [Ystart(ii),Xstart(ii)] = find(tData==max(tData(:)),1);
end
Xstart = (Xstart - PSFstruct.NPix/2)*PSFstruct.PixelSize;
Ystart = (Ystart - PSFstruct.NPix/2)*PSFstruct.PixelSize;

%This takes about 1 CPU hours to run on an AMD Ryzen 3.8 GHz 12-core processor.
tic();
[Chain,PSFstack] = runPhaseRetrieval(Data,PSFstruct,BNP,Xstart,Ystart);
T=toc();
sprintf('It took %0.5ds per iteration.',T/BNP.NJump)

%% generating the results using the last 5000 samples samples

Nsamples = 5000;
Mag = 0;
Phase = 0;
FoundX = zeros(size(Chain(1).Z));
FoundY = zeros(size(Chain(1).Z));
FoundZ = zeros(size(Chain(1).Z));
for ii = BNP.NJump-(Nsamples-1):BNP.NJump
    Mag = Mag + Chain(ii).Mag;
    Phase = Phase + Chain(ii).Phase;
    FoundX = FoundX + Chain(ii).X;
    FoundY = FoundY + Chain(ii).Y;
    FoundZ = FoundZ + Chain(ii).Z;
end
Phase = Phase/Nsamples;
Mag = Mag/Nsamples;
FoundX = FoundX/Nsamples;
FoundY = FoundY/Nsamples;
FoundZ = FoundZ/Nsamples;

%subtracting the contribution of particle location from the obtained phase
Z0 = PSFstruct.ZImages(:,:,1);
Zx = PSFstruct.ZImages(:,:,2);
Zy = PSFstruct.ZImages(:,:,3);
Zz = PSFstruct.ZImages(:,:,4);
KPixelSize = 1/(PSFstruct.PixelSize*PSFstruct.NPix); 
PupilRadius = PSFstruct.Na/(PSFstruct.Lambda*KPixelSize); 
ADisk = pi*PupilRadius^2;
Mask = PSFstruct.Mask;
tmpPhase = Phase;
for ii = 1:20
    XShift = sum(sum(tmpPhase.*Zx.*Mask))/ADisk;
    YShift = sum(sum(tmpPhase.*Zy.*Mask))/ADisk;
    ZShift = sum(sum(tmpPhase.*Zz.*Mask))/ADisk;
    Offset = sum(sum(tmpPhase.*Z0.*Mask))/ADisk;
    tmpPhase = (tmpPhase - XShift*Zx ...
        - YShift*Zy - ZShift*Zz - Offset*Z0).*Mask;
end

%plotting the learned phase
figure;pcolor(tmpPhase(8:27,8:27).*Mask(8:27,8:27));A2=colorbar();
title('Learned-Phase');set(gca,'Xtick',[]);set(gca,'Ytick',[])
xlabel('K_x(1/\mum)','FontSize',20);ylabel('K_y(1/\mum)','FontSize',20)
A2.Ticks = [-2.5 0 1.5];A2.FontSize = 15;

%plotting the ground truth phase
figure;pcolor(TruePhase(8:27,8:27).*Mask(8:27,8:27));A2=colorbar();
title('True-Phase');set(gca,'Xtick',[]);set(gca,'Ytick',[])
xlabel('K_x(1/\mum)','FontSize',20);ylabel('K_y(1/\mum)','FontSize',20)
A2.Ticks = [-2.5 0 1.5];A2.FontSize = 15;

%plotting the groundtruth magnitude
tMag=(TrueMag);
tMag = tMag./(max(tMag(:).*Mask(:))); %scaling the magnitude between 0 and 1
figure;pcolor(tMag(8:27,8:27).*PSFstruct.Mask(8:27,8:27));colorbar();
A3=colorbar();title('True-Magnitude')
set(gca,'Xtick',[]);set(gca,'Ytick',[])
xlabel('K_x(1/\mum)','FontSize',20);ylabel('K_y(1/\mum)','FontSize',20)
A3.Ticks = [0 1];A3.FontSize = 15;

%plotting the learned magnitude
tmpMag = Mag;
tmpMag = tmpMag./(max(tmpMag(Mask==1)));
figure;pcolor(tmpMag(8:27,8:27).*PSFstruct.Mask(8:27,8:27));A4=colorbar();
title('Learned-Magnitude');A4.Ticks = [0 1];A4.FontSize = 15;
set(gca,'Xtick',[]);set(gca,'Ytick',[])
xlabel('K_x(1/\mum)','FontSize',20);ylabel('K_y(1/\mum)','FontSize',20)

%plotting trajectories using 2500 last samples
Xchain = zeros(Nsamples,100,'single');
Ychain = zeros(Nsamples,100,'single');
Zchain = zeros(Nsamples,100,'single');
Ind = 0;
for ii = BNP.NJump-(Nsamples-1):BNP.NJump
    Ind = Ind + 1;
    Xchain(Ind,:) = Chain(ii).X;
    Ychain(Ind,:) = Chain(ii).Y;
    Zchain(Ind,:) = Chain(ii).Z;
end
Xfill = [(1:100),fliplr(1:100)];
XBetween = [quantile(Xchain,0.0025),fliplr(quantile(Xchain,0.9975))]+1600;
YBetween = [quantile(Ychain,0.0025),fliplr(quantile(Ychain,0.9975))]+1600;
ZBetween = [quantile(Zchain,0.0025),fliplr(quantile(Zchain,0.9975))];

%plotting X-trajectory
figure;hold;
fill(Xfill,XBetween,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
plot(Xtrue+1600,'--m');ylim([1500 2500]);xlabel('time(frame)')
set(gca,'Xtick',[0 100]);set(gca,'Ytick',[1500 2500],'FontSize',15);ylabel('X(nm)','FontSize',20)
legend({'95% confidence interval','ground truth'},'FontSize',12,'Location','northwest')

%plotting Y-trajectory
figure;hold;
fill(Xfill,YBetween,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
plot(Ytrue+1600,'--m');ylim([1700 2500]);xlabel('time(frame)')
set(gca,'Xtick',[0 100]);set(gca,'Ytick',[1700 2500],'FontSize',15);ylabel('Y(nm)','FontSize',20)
legend({'95% confidence interval','ground truth'},'FontSize',12)

%plotting Z-trajectory
figure;hold;
fill(Xfill,ZBetween,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
plot(Ztrue,'--m');ylim([-700 500])
set(gca,'Xtick',[0 100],'FontSize',15);set(gca,'Ytick',[-700 500],'FontSize',15)
ylabel('Z(nm)','FontSize',20);xlabel('time(frame)','FontSize',20)
legend({'95% confidence interval','ground truth'},'FontSize',12)

%Diffusion coefficient, intensity and background
TotPhoton = zeros(length(Chain),2);
TotDD = zeros(length(Chain),1);
TotBg = zeros(length(Chain),2);
for ii = 1:length(Chain)
    TotPhoton(ii,:) = Chain(ii).I';
    TotDD(ii) = Chain(ii).D;
    TotBg(ii,:) = Chain(ii).Bg';
end

%Histogram of 5000 last sampled diffuion constant 
figure;histogram(TotDD(end-Nsamples:end)/1000,20,'normalization','pdf') %divided by zero to convert unit
hold;plot((TrueD)*ones(size(0:250)),0:250,'--r','linewidth',2)
xlabel('diffusion constant(\mum^2/s)','FontSize',20);ylabel('probabiliy','FontSize',20)
legend('sampled values','ground truth');set(gca,'FontSize',15)

%plotting the entire sample chain of diffuion constant 
figure;plot((1:BNP.NJump)',TotDD/1000,'linewidth',1.2)
hold;plot(1:BNP.NJump,TrueD*ones(size(TotPhoton)),'r--','linewidth',1.5)
xlabel('samples');ylabel('Diffusion constant(\mum^2/s)')
legend('smapled values','ground truth')
set(gca,'Fontsize',12)

%Histogram of 5000 last sampled intensity
tPhotons = TotPhoton(end-Nsamples:end,:);
figure;histogram(tPhotons(:),20,'normalization','pdf')
hold;plot(TrueIntensity*ones(size(0:0.002:0.02)),0:0.002:0.02,'--r','linewidth',2)
xlabel('intensity(photons)','FontSize',20);ylabel('probability','FontSize',20)
legend('sampled values','ground truth');set(gca,'FontSize',15)

%plotting the entire sample chain of intensities
figure;plot((1:BNP.NJump)',TotPhoton,'linewidth',1.2)
hold;plot(1:BNP.NJump,TrueIntensity*ones(size(TotPhoton)),'r--','linewidth',1.5)
xlabel('samples');ylabel('Photons')
legend('intensity plane#1','intensity plane#2','ground truth')
set(gca,'Fontsize',12)

%Histogram of last 5000 sampled background
tBg = TotBg(end-Nsamples:end,:);
figure;histogram(tBg(:),20,'normalization','pdf')
hold;plot(TrueBg*ones(size(0:0.2:10)),0:0.2:10,'--r','linewidth',2)
xlabel('background(photons)','FontSize',20);ylabel('probability','FontSize',20)
legend('sampled values','ground truth');set(gca,'FontSize',15)

%plotting the entire sample chain of background
figure;plot((1:BNP.NJump)',TotBg,'linewidth',1.2)
hold;plot((1:BNP.NJump)',TrueBg*ones(length(TotPhoton),1),'r--','linewidth',1.5)
xlabel('samples');ylabel('background(photons)')
legend('background plane#1','background plane#2','ground truth')
set(gca,'Fontsize',12)
