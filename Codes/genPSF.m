function [PSF,Struct] = genPSF(Struct,PSFtype,Scale)
%genPSF generates different PSFs for given pupil or from a list
%
%INPUTS:
%   Struct:    Structure containing refractive index(N), numerical
%              apperture (Na), wavelength (Lambda), Pixelsize, length of
%              the image (BPix), sampled planes along axial direction (Z)
%   PSFtype:   A string indicating PSF typr
%              AI: airry PSF
%              AS: astigmatic PSF
%              DH: double helix PSF
%              TP: tetrapod PSF
%              X0: X-offset from the center (Default = 0) (nm)
%              Y0: Y-offset from the center (Default = 0) (nm)
%
%OUTPUTS:
%   PSF:       PSF stack at different axial locations
%

if nargin < 3
    Scale = 1;
end

if isfield(Struct,'X')
    X0 = Struct.X;
else
    X0 = zeros(size(Struct.Z)); 
end
if isfield(Struct,'Y')
    Y0 = Struct.Y;
else
    Y0 = zeros(size(Struct.Z)); 
end
%Sampling distance in frequency domain
KPixelSize = 1/(Struct.PixelSize*Struct.NPix); 
%Largest frequency available via pupil function in unit of KPixelSize
PupilRadius = Struct.Na/(Struct.Lambda*KPixelSize);

[Kx,Ky]=meshgrid((-Struct.NPix/2:Struct.NPix/2-1),(-Struct.NPix/2:Struct.NPix/2-1));
Rho = sqrt(Kx.^2+Ky.^2);
%Pupil function does not pass frequencies larger than PupilRadius
Mask = Rho <= PupilRadius;

if nargin>1
    Theta = atan2(Ky,Kx).*Mask;
    NRho = Rho/PupilRadius;
    if strcmp(PSFtype,'AS')
        Struct.Pupil_Mag = (1-0.2*(2*NRho.^2-1)).*Mask;
        %Z_4^-2
        Struct.Pupil_Phase = -2*Mask.*(4*NRho.^4-3*NRho.^2).*cos(2*Theta);
    elseif strcmp(PSFtype,'DH')
        L = 3;
        %Cannot be expressed in term of Zernike functions
        for ll=1:L
            M=(NRho>=sqrt((ll-1)/L))&(NRho<sqrt(ll/L));
            Struct.Pupil_Phase(M)=mod((2*(ll-1)+1)*Theta(M),2*pi); %make hole
        end
        Struct.Pupil_Mag = Mask;
    elseif strcmp(PSFtype,'TP')
        Struct.Pupil_Mag = Mask;
        % Z_2^-2, Z_4^-2, Z_6^-2, Z_6^-6
        Struct.Pupil_Phase = (4.5*(NRho.^2).*sin(2*Theta)-4.5*(4*NRho.^4-3*NRho.^2)...
            .*sin(2*Theta)-0.5*NRho.^6.*sin(6*Theta)-0.1*(15*NRho.^6-20*NRho.^4+6*NRho.^2)...
            .*sin(2*Theta)).*Mask;
    elseif strcmp(PSFtype,'AI')
        Struct.Pupil_Phase = zeros(size(Mask));
        Struct.Pupil_Mag = Mask;
    elseif strcmp(PSFtype,'RND') 
        [Xg,Yg] = meshgrid(Struct.PixelSize*(0.5:size(Mask,2)),Struct.PixelSize*(0.5:size(Mask,1)));
        T = 1;L = Struct.PixelSize*size(Mask,1)/9;
        K = calCOV([Xg(:),Yg(:)],[Xg(:),Yg(:)],T,L);
        Struct.Pupil_Phase = Scale*reshape(mvnrnd(zeros(length(Mask(:)),1),K),size(Mask));
        Expon = reshape(mvnrnd(zeros(length(Mask(:)),1),K/40),size(Mask));
        Struct.Pupil_Mag = Mask.*exp(Expon);
    elseif strcmp(PSFtype,'GivenTraj')
        [Xg,Yg] = meshgrid(Struct.PixelSize*(0.5:size(Mask,2)),Struct.PixelSize*(0.5:size(Mask,1)));
        T = 1;L = Struct.PixelSize*size(Mask,1)/9;
        K = calCOV([Xg(:),Yg(:)],[Xg(:),Yg(:)],T,L);
        Struct.Pupil_Phase = Scale*reshape(mvnrnd(zeros(length(Mask(:)),1),K),size(Mask));
        Expon = reshape(mvnrnd(zeros(length(Mask(:)),1),K/40),size(Mask));
        Struct.Pupil_Mag = Mask.*exp(Expon);
    elseif strcmp(PSFtype,'Given') && isfield(Struct,'Pupil_Phase')
       %Do nothing
    else
        error('This PSF does not exist on the list.\n You must pick one of these options: AS, DH, TP, AI, RND');
    end
end

%R^2 up to the radius of pupil
Kr_Image = Rho.*Mask*KPixelSize; 
DefocusK = 2*pi*sqrt((Struct.N/Struct.Lambda)^2-Kr_Image.^2);
XOffsetPhase = 2*pi*Kx*KPixelSize;
YOffsetPhase = 2*pi*Ky*KPixelSize;

Z0 = Struct.ZImages(:,:,1);
Zx = Struct.ZImages(:,:,2);
Zy = Struct.ZImages(:,:,3);
Zz = Struct.ZImages(:,:,4);
Ax = sum(sum(Zx.*Zx.*Mask));
Ay = sum(sum(Zy.*Zy.*Mask));
Az = sum(sum(Zz.*Zz.*Mask));
A0 = sum(sum(Z0.*Z0.*Mask));

for ii = 1:5
    XShift = sum(sum(Struct.Pupil_Phase.*Zx.*Mask))/Ax;
    YShift = sum(sum(Struct.Pupil_Phase.*Zy.*Mask))/Ay;
    ZShift = sum(sum(Struct.Pupil_Phase.*Zz.*Mask))/Az;
    Offset = sum(sum(Struct.Pupil_Phase.*Z0.*Mask))/A0;
    Struct.Pupil_Phase = (Struct.Pupil_Phase - XShift*Zx ...
        - YShift*Zy - ZShift*Zz - Offset*Z0).*Mask;
end

SZ = Struct.SubPixel*Struct.NPix;
SubPixelZeros = complex(zeros(SZ)); 
StartInd = SZ/2-Struct.NPix/2+1;
EndInd = SZ/2+Struct.NPix/2;
PSF = zeros(Struct.NPix,Struct.NPix,length(Struct.Z),'single');
for zz = 1:length(Struct.Z)
    DefocusPhase = DefocusK*Struct.Z(zz);
    Phase = Mask.*(DefocusPhase+Struct.Pupil_Phase ...
        +XOffsetPhase*X0(zz)+YOffsetPhase*Y0(zz));
    %Optical transfer function (propagator)
    OTF = Mask.*Struct.Pupil_Mag.*exp(1i*Phase);
    %Parseval Normalizaton
    Norm = sqrt(sum(sum(abs(OTF))))*Struct.NPix;
    %Padding zeros for subpixel resolution
    SubPixelZeros(StartInd:EndInd,StartInd:EndInd) = OTF;
    tmpPSF = single((abs(fftshift(fft2(SubPixelZeros/Norm)))).^2);
    PSF(:,:,zz) = blockresample(tmpPSF,1/Struct.SubPixel)/Struct.SubPixel;
    PSF(:,:,zz) = PSF(:,:,zz)/sum(sum(PSF(:,:,zz)));
end
Struct.Mask = single(Mask);
Struct.Pupil_Phase = single(Struct.Pupil_Phase);
Struct.Pupil_Mag = single(Struct.Pupil_Mag);
end