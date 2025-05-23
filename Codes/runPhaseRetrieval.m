function [Chain,PSFstack] = runPhaseRetrieval(Data,Struct,BNP,Xstart,Ystart,TempFlag,CMOS_Noise)
%This is the main function that implements the the Gibbs sampling as 
% explained in the paper. Here, the chain is first initialized and then
% each samoling function is called to iteratively sample each parameter.
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
% TempFlag: 0 no temperature, 1 temporal sampling (DEfault: 0)
% CMOS_Noise: pixel-by-pixel map of CMOS noise
%
%OUTPUTS
% Chain: chain of samples including
%           Mag: sampled magnitudes
%           Phase: sampled phases
%           PSFstack: current model
%           Bg: sampled background
%           I: sampled photon counts
%           X: sampled particle X-locations
%           Y: sampled particle Y-locations
%           Z: sampled particle Z-locations
%           D: sampled diffuion coefficients
%           LogLike: log-likelihood
%           LogPost: log-posterior
%           SigConv: sigma of the Gaussian used in smoothening the PSF.
% PSFstack: model associated to the last state of the chain.
% Author:
%   Mohamadreza Fazel, Presse lab, 2024.
%

%setting up the chain
Chain(BNP.NJump).Mag = [];
Chain(BNP.NJump).Phase = [];
Chain(BNP.NJump).PSFstack = [];
Chain(BNP.NJump).Bg = [];
Chain(BNP.NJump).I = [];
Chain(BNP.NJump).X = [];
Chain(BNP.NJump).Y = [];
Chain(BNP.NJump).Z = [];
Chain(BNP.NJump).D = [];
Chain(BNP.NJump).LogLike = [];
Chain(BNP.NJump).LogPost = [];
Chain(BNP.NJump).SigConv = zeros(BNP.NJump,1,'single');
Chain(1).SigConv = 1;

if ~isfield(Struct,'SuPixel')
    Struct.SubPixel = 1;
end
if ~isfield(Struct,'Tform')
    for nn = 1:size(Struct.DelX,1)+1
        Struct.Tform(nn) = affine2d(eye(3));
    end
end

%size of frame model with subpixel resolution
SZ = Struct.SubPixel*Struct.NPix;
SubPixelZeros = complex(zeros(SZ)); 
StartInd = SZ/2-Struct.NPix/2+1;
EndInd = SZ/2+Struct.NPix/2;

%Sampling distance in frequency domain
KPixelSize = 1/(Struct.PixelSize*Struct.NPix); 
%Largest frequency available via pupil function in unit of KPixelSize
PupilRadius = Struct.Na/(Struct.Lambda*KPixelSize);
[Kx,Ky]=meshgrid((-Struct.NPix/2:Struct.NPix/2-1),(-Struct.NPix/2:Struct.NPix/2-1));
%Covariance Matrix
Coordinates = KPixelSize*[Kx(:),Ky(:)];
K_A = calCOV(Coordinates,Coordinates,BNP.T_A,BNP.L_A);
K_Phi = calCOV(Coordinates,Coordinates,BNP.T_Phi,BNP.L_Phi);
%Cholesky decomposition for covariance matrix. Here, a very small value is
%added to diagonal elements to avoid zero determinant.
Chol_A = cholcov(K_A+1000*eps*eye(size(K_A)));
Chol_Phi = cholcov(K_Phi+1000*eps*eye(size(K_Phi)));
if isempty(Chol_A) || isempty(Chol_Phi)
    Chol_A = cholcov(K_A+(10*size(K_A,1)*eps+eps)*eye(size(K_A)));
    Chol_Phi = cholcov(K_Phi+(10*size(K_Phi,1)*eps+eps)*eye(size(K_Phi)));
end

Rho = sqrt(Kx.^2+Ky.^2);
%Pupil function does not pass frequencies larger than PupilRadius
Mask = single(Rho <= PupilRadius);
%R^2 up to the radius of pupil
Kr_Image = Rho.*Mask*KPixelSize; 
DefocusK = 2*pi*sqrt((Struct.N/Struct.Lambda)^2-Kr_Image.^2);

XOffsetPhase = 2*pi*Kx*KPixelSize;
YOffsetPhase = 2*pi*Ky*KPixelSize;

NStack = size(Data,3);
DelX = cat(1,[0 0 0],Struct.DelX);
%initializing the chain
Chain(1).Mag = ones(size(Mask));
Chain(1).Phase = zeros(size(Mask));%Struct.Pupil_Phase;
Chain(1).Bg = gamrnd(1,20,size(DelX,1),1);
Chain(1).I = gamrnd(5,2000,size(DelX,1),1);
Chain(1).D = 100*rand();
if nargin < 4
    Chain(1).X = zeros(1,NStack);
else
    Chain(1).X = Xstart;
end
if nargin < 5
    Chain(1).Y = zeros(1,NStack);
else
    Chain(1).Y = Ystart;
end
if nargin < 6
   TempFlag = 0; 
end
if nargin < 7
   SZData = size(Data);
   CMOS_Noise = zeros([SZData(1:2),2]); 
end
Chain(1).Z = zeros(1,NStack);

%CMOS noise is added to the data for noise modeling in CMOS cameras
for nn = 1:NStack
    for ii = 1:size(DelX,1)
        Data(:,:,nn,ii) = Data(:,:,nn,ii) + CMOS_Noise(:,:,ii);
    end
end

%calculating the model using the initial values
PSFstack = [];
for ii = 1:size(DelX,1)
    PSFstack = cat(4,PSFstack,findPSF(Chain(1).Mag,Chain(1).Phase, ...
    Chain(1).Bg(ii),Chain(1).I(ii),DefocusK,Chain(1).Z+DelX(ii,3),Mask,...
    SubPixelZeros,StartInd,EndInd,Struct.SubPixel,Chain(1).X+DelX(ii,1),...
    Chain(1).Y+DelX(ii,2),XOffsetPhase,YOffsetPhase,Struct.Tform(ii),Chain(1).SigConv,CMOS_Noise(:,:,ii)));
end

Z0 = Struct.ZImages(:,:,1); %offset Zernike
Zx = Struct.ZImages(:,:,2); %X-tilt Zernike
Zy = Struct.ZImages(:,:,3); %Y-tilt Zernike
Zz = Struct.ZImages(:,:,4); %Z-tilt Zernike

AcceptPup = 0;
AcceptTraj = 0;
AcceptOffset = 0;
AcceptI = 0;
AcceptBg = 0;
AcceptSig = 0;
tmpPhase = Chain(1).Phase;
%Gibbs sampling iteratively samples the posterior, one parameter at a time
for jj = 2:BNP.NJump
    if jj/1000 == floor(jj/1000)
        fprintf('Jumps: %d out of %d\n',jj,BNP.NJump)
        fprintf('Accepted GP jumps: %d\n',AcceptPup)
        fprintf('Accepted X jumps: %d\n',AcceptTraj)
        fprintf('Accepted I jumps: %d\n',AcceptI)
        fprintf('Accepted Bg jumps: %d\n',AcceptBg)
        fprintf('Accepted SigCov jumps: %d\n',AcceptSig)
    end
    %setting temperature
    if TempFlag ~= 0
        if jj < 2000
            Temp = 1000;
        elseif jj < 10000
            Temp = 100;
        elseif jj < 20000
            Temp = 10;
        else
            Temp = 1;
        end
    else
       Temp = 1; 
    end
    tChain = Chain(jj-1);
    
    %sampling the standard deviation of the Gaussian convolved with PSF
    [tChain.SigConv,PSFstack,AcceptSig] = sampleSigConv(Data,tChain,PSFstack,...
    DefocusK,Mask,DelX,AcceptSig,XOffsetPhase,YOffsetPhase,tmpPhase,...
    SubPixelZeros,StartInd,EndInd,Struct.SubPixel,Struct.Tform,tChain.SigConv,CMOS_Noise);
    
    %sampling pupil function (phase and magnitude)
    [tChain.Mag,tChain.Phase,tmpPhase,PSFstack,AcceptPup] = samplePupil(Data,tChain,PSFstack,Chol_A,Chol_Phi,...
        DefocusK,Mask,DelX,AcceptPup,XOffsetPhase,YOffsetPhase,tmpPhase,Z0,Zx,Zy,Zz,...
        SubPixelZeros,StartInd,EndInd,Struct.SubPixel,Struct.Tform,tChain.SigConv,Temp,CMOS_Noise);
    
    %sampling tajectories
    [tChain.X,tChain.Y,tChain.Z,PSFstack,AcceptTraj] = sampleTraj(Data,tChain,...
              PSFstack,BNP,Struct,DefocusK,Mask,DelX,AcceptTraj,XOffsetPhase,...
              YOffsetPhase,tmpPhase,SubPixelZeros,StartInd,EndInd,...
              Struct.SubPixel,Struct.Tform,tChain.SigConv,Temp,CMOS_Noise);
    
    %sampling trajectories
    [tChain.X,tChain.Y,tChain.Z,PSFstack,Norm1_PSFstack,AcceptOffset] = ...
              sampleTrajOffset(Data,tChain,PSFstack,BNP,Struct,DefocusK,...
              Mask,DelX,AcceptOffset,XOffsetPhase,YOffsetPhase,...
              tmpPhase,SubPixelZeros,StartInd,EndInd,Struct.SubPixel,...
              Struct.Tform,tChain.SigConv,Temp,CMOS_Noise);      
    
    %sampling diffusion constant
    tChain.D = sampleDiff(tChain,DelX,BNP.Dt);      
          
    %sampling emitter intensity
    [tChain.I,PSFstack,AcceptI] = sampleIntensity(Data,tChain,PSFstack,Norm1_PSFstack,DelX,AcceptI,CMOS_Noise); 
    
    %sampling offset background
    [tChain.Bg,PSFstack,AcceptBg] = sampleBg(Data,tChain,PSFstack,Norm1_PSFstack,DelX,AcceptBg,CMOS_Noise);

    %finding log-posterior
    [tChain.LogPost,tChain.LogLike] = calLogPost(Data,Chain,PSFstack,Struct,Chol_A,Chol_Phi,BNP);      
          
    Chain(jj) = tChain;
end

end
