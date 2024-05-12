# PSFToolbox-Matlab
Includes various ways to calculate point spread functions (PSFs) with a focus on high NA fluorescence microscopy. It also includes the ability to estimate aberrations caused by refractive index mismatch.


## Installation

To run this toolbox, the user needs to install [DipImage](https://diplib.org/).  

The MATLAB m files for generating the PSF are stored in the folder "Code". To run the code, change the current directory to the directory of the toolbox. This can be done using the code below:
```bash
cd <PSFToolbox-Matlab directory> 
```
If necessary, add the path of each subfolder in the Code folder into the MATLAB current directory. 

## Toolbox description

Use the test.m file to run a test of the computation. 

The main function for computing the PSF is called GenericPSFSim. It takes the following as input and output argument.

INPUT
    ImageParam: a structure containing the following fields:
        Sampling := voxel size represented in a vector format of three elements along x, y and, z
        Size := window size along x, y and, z  
    PSFParam: a structure containing the following minimum fields:
        polarization = 'circular';  % default if not mentionned in the structure
        polarization = 'linearX' or 'linearY' or 'radial' or 'azimuthal'; % are other options
        Aplanatic = -1;  % -1 means emission, 0: no aplanatic factor, 1: excitation. 
            Note: Some function can also support the PSFParam: 'Aplanatic', 2: % precomensate twice for excitation when projecting onto a sphere in 3D (SincR method), 3: % this means that the intensity will be modified with cosalpha^2 
        NA := numerical aperture of the objective lens
        n := refractive index of the immersion medium
        lambdaEm := emission wavelength
    Method: a string for selecting the method to use in the simulation. Choices are:
        'RichardsWolff' (default): Using the method from the Richards & Wolff paper [1, 2].
        'RichardsWolffInt' : A fast method, if only the intensity is needed [1, 2].
        'RWFast': Reworked Richards & Wolf bast on an chirp-Z transform along kz [3, 4].
        'SlicePropagation' : Simulates the in-plane ATF based on FFTs and the jinc function and propagates it to all planes [4].
        'CZT' chirped-Z transform-based (== Zoomed FFT) method avoiding the out-of-focus wrap-around problems of SlicePropagation [3, 4].
        'VolumeShell' : Creates a part of the McCutchen Pupil directly in Fourier space using pre-computed Fourier-space interpolators. Pretty good in the central region of the PSF but (intentionally) degrading near the edge of the volume [4].
        'SincR' :  Creates the Fourier shell by a 3D FFT of a sinc(abs(R))
        function and modifies it from there on [4].
   AddParams: structure having as fields
                        e.g. AddParams=struct('ns',1.33,'ng',1.518,'ni',1.516,'ni0',1.518,'ng0',1.518,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5); 
                          'ns': ri of sample, can be in a vector format if the sample is in a stratified medium.  First element of ns (if a vector) is the ri of where the emitter is at, it is the ri of the medium which the furthest from the coverslip. The last element of the vector is the ri of the medium closest to the coverslip
                          'ng': ri of coverslip in real condition
                          'ni': ri of immersion medium in real condition
                          'ni0': ri of immersion medium in design condition
                          'ts': thickness of the sample or position where the emitter is at, this should be with the same length as the vector of ns and within the same order as ns.
                          'tg': thickness of coverslip in real condition
                          'tg0': thickness of coverslip in design condition
                          'wd': working distance which is supposed to be the thickness of the immersion medium in design condition
          * if AddParams == [] or not added to the function, the code will run as there is no interfaces and no aberrations due to that, n = PSFParam.n will be used troughout the calculation in this case.
   AddPhase : additional phase that can be added into the function. In Method = 'ZoomedPupilProp', it is prefereable to have this AddPhase variable as a cell for computing the ZernikePoly instead while it can only be in 2D image for the other methods. 
  
  The Sampling, thicknesses and wavelengths must be in the same units.

  OUTPUT
    h : resulting intensity PSF
    amp4d : 4D complex amplitude distribution (if supported).
  
  EXAMPLE OF USE:
  ImageParam=struct('Sampling',[40 40 50],'Size',[128 128 32]);
  PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
  AddParams=struct('ns',1.33,'ng',1.518,'ni',1.518,'ng0',1.518,'ni0',1.518,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
  [h,amp4d]=GenericPSFSim(ImageParam,PSFParam,'VolumeShell',AddParams)

  REFERENCES
  [1] Wolf E. Electromagnetic diffraction in optical systems-I. An integral representation of the image field. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences. 1959 Dec 15;253(1274):349-57.
  [2] Richards B, Wolf E. Electromagnetic diffraction in optical systems, II. Structure of the image field in an aplanatic system. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences. 1959 Dec 15;253(1274):358-79.
  [3] Rabiner L, Schafer RW, Rader C. The chirp z-transform algorithm. IEEE transactions on audio and electroacoustics. 1969 Jun;17(2):86-92.
  [4] Miora RH, Rohwer E, Kielhorn M, Sheppard CJ, Bosman G, Heintzmann R. Calculating point spread functions: Methods, pitfalls and solutions. arXiv preprint arXiv:2301.13515. 2023 Jan 31.