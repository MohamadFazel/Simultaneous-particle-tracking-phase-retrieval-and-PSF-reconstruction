function Struct = genZernikeIm(Struct,NollIndex)
%This function generates Zernike polynomials up to the given number
%(NollIndex) using the parameters provided inside "Struct".
%
%INPUTS
% Struct: structure containing parameters required to generate Zernike
%         polynomials including:
%    Na: numerical aperture
%    N: refractive index
%    Lambda: fluorescence emission wavelength (nm)
%    PixelSize: pixel size (nm)
%    NPix: length of frames assumed to be square (pixel)
% NollIndex: number of Zernike polynomials.
%OUTPUTS
% Struct: the output structure contains an additional field "ZImages" which
%         are Zernike polynomials.
%
% Author:
%   Mohamadreza Fazel, Presse lab, 2024.
%

if nargin < 2
   NollIndex = 3; 
end

[Kx,Ky]=meshgrid((-Struct.NPix/2:Struct.NPix/2-1),(-Struct.NPix/2:Struct.NPix/2-1));
Rho=sqrt(Kx.^2+Ky.^2);
%Sampling distance in frequency domain
KPixelSize = 1/(Struct.PixelSize*Struct.NPix); 
%Largest frequency available via pupil function in unit of KPixelSize
PupilRadius = Struct.Na/(Struct.Lambda*KPixelSize);

Mask = Rho <= PupilRadius;
Theta=(atan2(Ky,Kx)); 
NRho = Rho/PupilRadius;

Struct.ZImages = zeros(Struct.NPix,Struct.NPix,NollIndex,'single');

for ii = 1:NollIndex+1

    %convert to NM
    N = ceil((-3 + sqrt(1 + 8*ii)) / 2);
    M = ii - N * (N + 1) / 2 - 1;
    if mod(N, 2) ~= mod(M, 2)
        M = M + 1;
    end
    if mod(ii, 2) == 1
        M = -M;
    end

    if M<0
        Image=Mask.*sin(M*Theta);
    else
        Image=Mask.*cos(M*Theta);
    end

    N=abs(N);
    M=abs(M);

    %calculate the radial polynomial
    tmp=0;
    for kk=0:(N-M)/2
        Poly=NRho.^(N-2*kk);
        Coef = (-1)^kk * prod((N - M)/2 - kk + 1 : N - kk) / ...
                         (factorial(kk) * factorial((N + M)/2 - kk));
        tmp=tmp+Coef*Poly;
    end

    Struct.ZImages(:,:,ii) =Mask.*Image.*tmp;
    Struct.Mask = Mask;
end
end