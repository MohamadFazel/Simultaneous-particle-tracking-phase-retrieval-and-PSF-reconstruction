# Simultaneous particle tracking, phase retrieval and PSF reconstruction
Mohamadreza Fazel, Reza Hoseini, Maryam Mahmoodi, Lance W. Q. Xu, Ayush Saurabh, Zeliha Kilic,
Julian Antolin, Kevin Lynn Scrudders, Douglas Shepherd, shalini Low-Nam, Fang Huang, Steve Presse

# Overview of the framework
3D localization and tracking of particles, typically fluorescently labeled biomolecules, provides a direct means to monitor details within nano-scale environments. However, sample induced aberrations of point spread functions (PSFs) are a challenge for 3D particle tracking often relying on pre-calibrated PSFs. For instance, errors in PSF calibration can lead to inaccurate mislocalizations, which may ultimately prevent surpassing the diffraction limit. In order to retain superresolution in tracking, it is often required to simultaneously retrieve the pupil function and, thus, the resulting PSF while tracking. In practice, sample-induced aberrations are often corrected using adaptive optics requiring additional hardware and must be repeated for each sample. Yet, assuming point emitters, the information on the pupil function is already encoded from data collected on a multiplane setup. Here, we provide a software package to track and simultaneously learn the pupil phase and amplitude from the input data. In order to achieve this, we use a Bayesian method placing continuous 2D priors on all possible phase and amplitudes warranted by the data without limiting ourselves to a Zernike set. Moreover, our method is data efficient by taking into account all existing sources of uncertainty in the problem, including detector noise, Poisson shot noise, pixelization, and other sources. 

# Software package overview and how to use it
The software package contains two example matlab scripts ("example_experimentalData.m" and "example_syntheticData.m") to analyze the experimental data and the synthetic data. The package also contains 4 folders: 
"Codes" containing algorithm functions and methods; 
"Data" containing a synthetic data set simulated assuming a bi-focal setup, 100 frames of data from each focal plane; and frames of 32x32 pixels. This is the data set used in Fig. S6. The experimental data set contains an in vitro data set obtained using a bi-focal setup, 150 frames from each focal plane, and 62x62 pixels per frame. This is the data set used to make Fig. 3. 
"PlaneRegistration_ExpDataOnly" contains the tForm matrix characterizing relative shift and rotation between the two focal planes whcih was obtained using calibration data from very bright beads. This is only used for the experimental data.
"ExpectedResults" contains two sub-folder each containing the expected outputs of the experimental and synthetic data analyses.

To use this software package you need matlab 2023 or higher versions, matlab image processing toolbox and matlab statistics and machine learning toolbox. To run the software package, just open the example scripts and run them. It took almost one hour to finish analyzing the synthetic data and approximately 5 hours to analyze the experimental data on an AMD Ryzen 3.8 GHz 12-core processor. When done analyzing data the example scripts also plot the results. 

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
    ZImage: Zernike polynomials
    Tform: The model for each plane is first generated with no XY-shift and rotation with 
           respect to the reference (first) plane. Tform is an array with the same size as the 
           number of planes where each element contains the affine transform (a 3x3 matrix) that
           transforms the generated model to match the input data. As such, its first element is 
           always a diagonal matrix corresponding to no XY-shift and rotation. Note that the 
           Z-differences across planes are given in the DelX parameter.
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

# OUTPUT PLOTS
 Background_Samples: plot of the entire chain of sampled backgrounds 
 Background_Hist: histogram of the sampled backgrounds after chain convergence
 DiffusionConst_Hist: histogram of the sampled diffusion constants after chain convergence
 DiffusionConst_Samples: plot of the entire chain of sampled diffusion constants
 Learned_Magnitude: median of the sampled pupil amplitudes after chain convergence
 True_Magnitude: ground truth pupil amplitude
 Learned_Phase: median of the sampled pupil phases after chain convergence
 True_Phase: ground truth phase
 ParticleIntensity_Hist: histogram of the sampled particle intensities after chain convergence
 ParticleIntensity_Samples: plot of the entire chain of sampled particle intensities
 X-trajectory: plot of the median and 95% confidence interval of the sampled X-trajectories
 Y-trajectory: plot of the median and 95% confidence interval of the sampled Y-trajectories
 Z-trajectory: plot of the median and 95% confidence interval of the sampled Z-trajectories

# How to process your data
To process your data, we would suggest modifying the script "example_experimentalData.m" and adjusting the experimental parameters (fields of PSFstruct and Dt in BNP) for your data. The other parameters should be mostly fine as they are. To reduce the run-time pick a small ROI around the target particle so the PSF does not move outside the ROI for both planes. The number of frames should be selected so that the particle ideally explores a few hunderd nano-meters below and above the focal plane. We also emphasize that the affine transform should be calibrated so that it transforms the reference plane into the second plane NOT the opposite.

# Support
If you happen to hve any questions please shoot us an email at:

Mohamadreza Fazel, fazel.mohamadreza@gmail.com

Steve Presse, spresse@asu.edu
