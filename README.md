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

### INPUT <br>
> <font color='blue'>ImageParam</font>: a structure containing the following fields: <br>
* Sampling : voxel size represented in a vector format of three elements along x, y and, z <br>
* Size : window size along x, y and, z  <br>

> <font color='blue'>PSFParam</font>: a structure containing the following minimum fields: <br>
* polarization : polarization state of the field. It can take as an input integers 'linearX' or 'linearY' or 'radial' or 'azimuthal' or 'circular'. The default polarization stated is assumed to be 'circular' if not mentionned <br>
* Aplanatic : value can be $`\beta = -1, 0`$, or 1 <br>
This corresponds to the aplanatic correction applied in the pupil plane. For the correction, the pupil plane is multiplied by $`cos^{\beta/2}(\theta)`$, $`\beta`$ being the input argument in the Aplanatic and $`\theta`$ is the opening or elevation angle of the beam. A value of $`\beta = 0`$ means there is no aplanatic correction applied onto the field. By default, $`\beta = -1`$ <br>
* NA : numerical aperture of the objective lens <br>
* n : refractive index of the immersion medium <br>
* lambdaEm : emission wavelength <br>

> <font color='blue'>Method</font>: a string for selecting the method to use in the simulation. Choices are: <br>
- 'RichardsWolff' (default): Using the method from the Richards & Wolff paper to generate the electric field as well as the image intensity [1, 2]. <br>
- 'RichardsWolffInt' : A fast method, if only the intensity is needed [1, 2]. <br>
- 'RWFast': Reworked Richards & Wolf bast on an chirp-Z transform along kz [3, 4]. This method is still under revision. <br>
- 'SlicePropagation' : uses the angular spectrum method to propagate the field from the in-plane amplitude transfer function (ATF) using a fast Fourier transform (FFT). A jinc function is used to simulate a smooth aperture to reduce Fourier artefacts [4]. <br>
- 'CZT' : a chirped-Z transform-based (== Zoomed FFT) method to avoid the out-of-focus wrap-around problems of SlicePropagation [3, 4]. <br>
- 'VolumeShell' : this creates a part of the McCutchen Pupil directly in Fourier space using pre-computed Fourier-space interpolators which is pretty good in the central region of the PSF but (intentionally) degrading near the edge of the volume [4]. <br>
- 'SincR' :  this creates the Fourier shell by a 3D FFT of a sinc(abs(R)) <br>

> <font color='blue'>AddParams</font>: a structure with the details on the refractive indices and thickeness of each layer in the system <br>
* 'ns': refractive index (ri) of the sample. If the sample is a stratified medium, value attributed to this variable can be in a vector format listing the refractive index of each layer such that the first element is the ri of where the emitter is at i.e. the the ri of the layer furtherst away from the coverslip. The last element of the vector is the ri of the medium closest to the coverslip. <br>
* 'ng': ri of the coverslip in real condition. <br>
* 'ni': ri of the immersion medium in real condition. <br>
* 'ni0': ri of the immersion medium in design condition. <br>
* 'ts': thickness of the sample or position where the emitter is at, this should be with the same length as the vector of ns and within the same order as ns. <br>
* 'tg': thickness of the coverslip in real condition. <br>
* 'tg0': thickness of the coverslip in design condition. <br>
* 'wd': working distance which is supposed to be the thickness of the immersion medium in design condition. <br>

If AddParams is empty i.e. AddParams == [] or not added to the function, the code will run as with the assumption that there is no interfaces with different refractive index mismatches yielding to aberrations into the system. In this case, the refractive index that is input in the dictionary PSFParam will be used troughout the calculation as the ri of the medium in which light is propagating. <br>

Example of usage: AddParams = struct('ns', 1.33, 'ng', 1.518, 'ni', 1.516, 'ni0', 1.518, 'ng0', 1.518, 'ts', 2e3, 'tg', 1.7e5, 'tg0', 1.7e5, 'wd', 1.5e5); <br>

> <font color='blue'>AddPhase</font> : additional phase that can be added into the function. In Method = 'ZoomedPupilProp', it is prefereable to have this AddPhase variable as a cell for computing the ZernikePoly instead while it can only be in 2D image for the other methods. 
  
The Sampling, thicknesses and wavelengths must be in the same units.

### OUTPUT <br>
<font color='blue'>h</font>: resulting intensity PSF <br>
<font color='blue'>amp4d</font> : 4D complex amplitude distribution (if supported). <br>
  
### EXAMPLE OF USE: <br>
```bash
ImageParam = struct('Sampling', [40 40 50], 'Size', [128 128 32]); 
PSFParam = struct('NA', 1.2, 'n', 1.33, 'lambdaEm', 520); 
AddParams = struct('ns', 1.33, 'ng', 1.518, 'ni', 1.518, 'ng0', 1.518, 'ni0', 1.518, 'ts', 2e3, 'tg', 1.7e5, 'tg0', 1.7e5, 'wd', 1.5e5);
[h,amp4d] = GenericPSFSim(ImageParam, PSFParam, 'VolumeShell', AddParams);
```

### REFERENCES

[1] Wolf E. Electromagnetic diffraction in optical systems-I. An integral representation of the image field. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences. 1959 Dec 15;253(1274):349-57.

[2] Richards B, Wolf E. Electromagnetic diffraction in optical systems, II. Structure of the image field in an aplanatic system. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences. 1959 Dec 15;253(1274):358-79.

[3] Rabiner L, Schafer RW, Rader C. The chirp z-transform algorithm. IEEE transactions on audio and electroacoustics. 1969 Jun;17(2):86-92.
  
[4] Miora RH, Rohwer E, Kielhorn M, Sheppard CJ, Bosman G, Heintzmann R. Calculating point spread functions: Methods, pitfalls and solutions. arXiv preprint arXiv:2301.13515. 2023 Jan 31.