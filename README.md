# Simultaneous single particle tracking, phase retrieval and PSF reconstruction
Mohamadreza Fazel, Reza Hoseini, Maryam Mahmoodi, Kevin Lynn Scrudders, Lance W. Q. Xu, Zeliha Kilic, Julian Antolin, Shalini Low-Nam, Fang Huang, Steve Pressé

# Overview of the framework
3D localization and tracking of particles, typically fluorescently labeled biomolecules, provides a direct means to monitor details within nano-scale environments. However, sample induced aberrations of point spread functions (PSFs) are a challenge for 3D particle tracking often relying on pre-calibrated PSFs. For instance, errors in PSF calibration can lead to inaccurate mislocalizations, which may ultimately prevent surpassing the diffraction limit. In order to retain superresolution in tracking, it is often required to simultaneously retrieve the pupil function and, thus, the resulting PSF while tracking. In practice, sample-induced aberrations are often corrected using adaptive optics requiring additional hardware and must be repeated for each sample. Yet, assuming point emitters, the information on the pupil function is already encoded from data collected on a multiplane setup. Here, we provide a software package to track and simultaneously learn the pupil phase and amplitude from the input data. In order to achieve this, we use a Bayesian method placing continuous 2D priors on all possible phase and amplitudes warranted by the data without limiting ourselves to a Zernike set. Moreover, our method is data efficient by taking into account all existing sources of uncertainty in the problem, including detector noise, Poisson shot noise, pixelization, and other sources. 

# How to use the software


# Input parameters
Data: input data is a 4D array with the 1st and 2nd dimensions are frame
       sizes, the 3rd dimension is number of frames from each plane, and
       the fourth dimension is number of planes.
Struct: structure containing experimental parameters
    Na: numerical aperture
    N: refractive index
    Lambda: fluorescence emission wavelength (nm)
    PixelSize: pixel size (nm)
    NPix: length of frames assumed to be square (pixel)
    DelX: spatial (X,Y,Z) shifts across planes (nm)
BNP: structure containing math framework parameters
    NJumps: number of samples from the posterior (chain length)
    T_A: GP prior parameter on pupil magnitude (Default: 1)
    L_A: GP prior parameter on pupil magnitude related to magnitude correlation
    T_Phi: GP prior parameter on pupil phase (Default: 1)
    L_Phi: GP prior parameter on pupil phase related to phase correlation
    Dt: frame exposure time (ms)
    SubPixel: subpixel resolution used in model calculations (Default: 1)
Xstart: initial X-trajectory
Ystart: initial Y-trajectory
TempFlag: indicating to use temporal sampling (1) or not (0) (Default: 0)
CMOS_Noise: CMOS noise variance used in modeling data from CMOS cameras (Default: 0)

# Outputs
Chain: chain of samples containing the following parameters
    Mag: pupil magnitude
    Phase: pupil phase
    PSFstack: PSF model reconstructed from the last state of the chain
    Bg: offset background (photons)
    I: eitter intensity (photons)
    X: X-trajectory (nm)
    Y: Y-trajectory (nm)
    Z: Z-trajectory (nm)
    D: diffusion constant (nm^2/ms)
    LogLike: log-likelihood
    LogPost: log-posterior
    SigConv: Sigma of the Gaussian convolved with the PSF
