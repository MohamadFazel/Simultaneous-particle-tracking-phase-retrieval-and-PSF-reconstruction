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
% TempFlag: 1 implement temporal sampling, 0 regular sampling (Default: 0)
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

%Data is already corrected for gain and offset
load('Data/Data_Experimental.mat')

%% load affine transfrom used to register the two planes

load('PlaneRegistration_ExpDataOnly/tform.mat')

%% setting the experimental parameters

PSFstruct.Na = 1.35; %Numerical apperture
PSFstruct.N = 1.449; %Refractive index
PSFstruct.Lambda = 680; % Wavelength (nm)
PSFstruct.PixelSize = 120; %Spatial pixel size (nm)
PSFstruct.NPix = size(Data,1); %Length of image
PSFstruct.SubPixel = 2; %Subpixel resolution 
PSFstruct.DelX = [0 0 -370]; %[X Y Z]-shifts across the two planes
PSFstruct.Tform(1) = affine2d(eye(3)); %affine transforn of first plane with respect to the first frame
PSFstruct.Tform(2) = tform; %affine transform of the second plane with respect to the first plane

%Generating first 5 Zernike polynomial
NZernike = 5;
PSFstruct = genZernikeIm(PSFstruct,NZernike);

%% Analyze Data

BNP.NJump = 40000; %Number of samples
BNP.T_A = 1; %Magnitude GP prior parmeters (larger values impose smoothness)
BNP.L_A = 6/(PSFstruct.PixelSize*PSFstruct.NPix);
BNP.T_Phi = 1; %Phase GP prior parmeters (larger values impose smoothness)
BNP.L_Phi = 4/(PSFstruct.PixelSize*PSFstruct.NPix);
Dt = 20; %frame exposure time (ms)
BNP.Dt = Dt;
BNP.SubPixel = 2; %number of PSF samples inside each data pixel

%initializing the X and Y chains to the brightest pixels in each frame
Xstart = zeros(1,size(Data,3));
Ystart = zeros(1,size(Data,3));
for ii = 1:size(Data,3)
    tData = Data(:,:,ii,1);
    [Ystart(ii),Xstart(ii)] = find(tData==max(tData(:)),1);
end
Xstart = (Xstart - PSFstruct.NPix/2)*PSFstruct.PixelSize;
Ystart = (Ystart - PSFstruct.NPix/2)*PSFstruct.PixelSize;
TempFlag = 1;

%This takes about 5.5 CPU hours to run on an AMD Ryzen 3.8 GHz 12-core processor.
tic();
[Chain,PSFstack] = runPhaseRetrieval(Data,PSFstruct,BNP,Xstart,Ystart,TempFlag);
T = toc();

%% generating results using the last 10K samples

Mag = 0;
Phase = 0;
Photons = [];
Diff = [];
Xchain = zeros(10000,150,'single');
Ychain = zeros(10000,150,'single');
Zchain = zeros(10000,150,'single');
for ii = 30001:40000
    Mag = Mag + Chain(ii).Mag;
    Phase = Phase + Chain(ii).Phase;
    Diff = cat(1,Diff,Chain(ii).D);
    Photons = cat(1,Photons,Chain(ii).I);
    Xchain(ii-30000,:) = Chain(ii).X;
    Ychain(ii-30000,:) = Chain(ii).Y;
    Zchain(ii-30000,:) = Chain(ii).Z;
end
Phase = Phase/10000;
Mag = Mag/10000;

%pupil phase nd magnitude
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
    XShift = sum(sum(tmpPhase.*Zx.*Mask))/sum(sum(Zx.*Zx.*Mask));
    YShift = sum(sum(tmpPhase.*Zy.*Mask))/sum(sum(Zy.*Zy.*Mask));
    ZShift = sum(sum(tmpPhase.*Zz.*Mask))/sum(sum(Zz.*Zz.*Mask));
    Offset = sum(sum(tmpPhase.*Z0.*Mask))/sum(sum(Z0.*Z0.*Mask));
    tmpPhase = (tmpPhase - XShift*Zx ...
        - YShift*Zy - ZShift*Zz - Offset*Z0).*Mask;
end

figure;Ap=pcolor(-tmpPhase(14:51,14:51).*Mask(14:51,14:51));set(Ap,'EdgeColor','none')
A=colorbar();clim([-3.3 2]);A.Ticks = [-3 0 2];set(gca,'FontSize',15)
set(gca,'Xtick',[]);set(gca,'Ytick',[])

MaxMag = max(Mag(:).*Mask(:));
figure;Ap=pcolor(Mag(14:51,14:51).*Mask(14:51,14:51)/MaxMag);set(Ap,'EdgeColor','none')
A=colorbar();A.Ticks = [0 1];set(gca,'FontSize',15)
set(gca,'Xtick',[]);set(gca,'Ytick',[])

%trajectories
Xfill = [(1:150),fliplr(1:150)];
XBetween = [quantile(Xchain,0.0025),fliplr(quantile(Xchain,0.9975))]+3100;
YBetween = [quantile(Ychain,0.0025),fliplr(quantile(Ychain,0.9975))]+3100;
ZBetween = [quantile(Zchain,0.0025),fliplr(quantile(Zchain,0.9975))];

figure;hold;
fill(Xfill,XBetween,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
plot(median(Xchain)+3100,'m--','linewidth',1.10);ylim([2200 3800])
set(gca,'Xtick',[0 50 100 150],'FontSize',15);set(gca,'Ytick',[2200 3800],'FontSize',15)
ylabel('X(nm)','FontSize',15);xlabel('time(frame)','FontSize',15)
legend({'95% confidence interval','median'},'FontSize',15)
figure;hold;
fill(Xfill,YBetween,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
plot(median(Ychain)+3100,'m--','linewidth',1.10);ylim([1800 2400])
set(gca,'Xtick',[0 50 100 150],'FontSize',15);set(gca,'Ytick',[1800 2400],'FontSize',15)
ylabel('Y(nm)','FontSize',15);xlabel('time(frame)','FontSize',15)
figure;hold;
fill(Xfill,ZBetween,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
plot(median(Zchain),'m--','linewidth',1.10);ylim([-1000 500])
set(gca,'Xtick',[0 50 100 150],'FontSize',15);set(gca,'Ytick',[-1000 500],'FontSize',15)
ylabel('Z(nm)','FontSize',15);xlabel('time(frame)','FontSize',15)

%diffusion constant, particle intensities and background
figure;histogram(Diff/1000,20,'normalization','pdf')
xlabel('Diffusion constant (\mum^2/s)');ylabel('probability')
set(gca,'FontSize',15)

figure;histogram(Photons,50,'normalization','pdf')
xlabel('Photons');ylabel('Probabilioty')
set(gca,'FontSize',15)

TotPhoton = zeros(length(Chain),2);
TotDD = zeros(length(Chain),1);
TotBg = zeros(length(Chain),2);
for ii = 1:BNP.NJump
    TotPhoton(ii,:) = Chain(ii).I;
    TotDD(ii) = Chain(ii).D;
    TotBg(ii,:) = Chain(ii).Bg';
end

figure;plot(TotPhoton)
xlabel('samples');ylabel('Photons')
legend({'plane 1','plane 2'},'Location','South')
set(gca,'FontSize',15)

figure;plot(TotDD/1000)
xlabel('samples');ylabel('D(\mu m^2/s)')
set(gca,'FontSize',15)

figure;plot(TotBg,'linewidth',1.4)
xlabel('samples');ylabel('Background')
legend({'plane 1','plane 2'})
set(gca,'FontSize',15)

