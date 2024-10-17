# Simultaneous single particle tracking, phase retrieval and PSF reconstruction
Mohamadreza Fazel, Reza Hoseini, Maryam Mahmoodi, Kevin Lynn Scrudders, Lance W. Q. Xu, Zeliha Kilic, Julian Antolin, Shalini Low-Nam, Fang Huang, Steve Pressé

# Overview of the framework
3D localization and tracking of particles, typically fluorescently labeled biomolecules, provides a direct means to monitor details within nano-scale environments. However, sample induced aberrations of point spread functions (PSFs) are a challenge for 3D particle tracking often relying on pre-calibrated PSFs. For instance, errors in PSF calibration can lead to inaccurate mislocalizations, which may ultimately prevent surpassing the diffraction limit. In order to retain superresolution in tracking, it is often required to simultaneously retrieve the pupil function and, thus, the resulting PSF while tracking. In practice, sample-induced aberrations are often corrected using adaptive optics requiring additional hardware and must be repeated for each sample. Yet, assuming point emitters, the information on the pupil function is already encoded from data collected on a multiplane setup. Here, we provide a software package to track and simultaneously learn the pupil phase and amplitude from the input data. In order to achieve this, we use a Bayesian method placing continuous 2D priors on all possible phase and amplitudes warranted by the data without limiting ourselves to a Zernike set. Moreover, our method is data efficient by taking into account all existing sources of uncertainty in the problem, including detector noise, Poisson shot noise, pixelization, and other sources. 

# How to use the software


# Input parameters


# Outputs

