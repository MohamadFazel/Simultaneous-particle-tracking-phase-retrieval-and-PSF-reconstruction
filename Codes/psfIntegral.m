function PSFout = psfIntegral(PSFin,NSamPerPix,SZ)

[a,b,c]=size(PSFin);
Filter = ones(NSamPerPix);
PSF = zeros(a+NSamPerPix-1,b+NSamPerPix-1,c);
for ii=1:c
    PSF(:,:,ii) = (1/NSamPerPix)^2*conv2(PSFin(:,:,ii),Filter);
end
PSF = PSF(NSamPerPix:end-(NSamPerPix-1),NSamPerPix:end-(NSamPerPix-1),:);

Elems = floor(size(PSF,1)/NSamPerPix);
PSFout = zeros(SZ,SZ,c,NSamPerPix,NSamPerPix);
Indnn = 0;
for nn = 1:NSamPerPix
    Indmm = 0;
    Indnn = Indnn + 1;
    for mm = NSamPerPix:-1:1
        Indmm = Indmm + 1;
        PSFout(:,:,:,Indmm,Indnn) = PSF(nn:NSamPerPix:end,mm:NSamPerPix:end,:); 
    end
end
Cent = round(size(PSF,1)/2) + ceil(NSamPerPix/2);
IndArray = [];
for nn = 0:NSamPerPix-1
%     Inds = (Cent-nn) - NSamPerPix*floor((Cent-nn)/NSamPerPix):NSamPerPix:size(PSF,1);
    Inds = (Cent-nn)-floor(SZ/2)*NSamPerPix:NSamPerPix:(Cent-nn)+(floor(SZ/2)-1)*NSamPerPix;
    if Inds(1) <= 0; Inds(1) = []; end
    IndArray = cat(1,IndArray,Inds);
end

Indnn = 0;
for nn = 1:NSamPerPix
    Indmm = 0;
    Indnn = Indnn + 1;
    for mm = 1:NSamPerPix
        Indmm = Indmm + 1;
        PSFout(:,:,:,mm,nn) = PSF(IndArray(nn,:),IndArray(mm,:),:); 
    end
end

end