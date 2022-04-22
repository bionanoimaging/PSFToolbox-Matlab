Draft after the dashed line
--------------------
R Holinirina Dina Miora, G Bosman, E Rohwer, M Kielhorn, R Heintzmann (#YearOfPublication)
Calculating Point Spread Functions - Methods, Pitfalls and Solutions. #Journal and number
=================================================
Model version: 1.0
-------------------------------------------------
Abstract: 

The knowledge of the exact structure of the Point Spread Function (PSF) for a given optical system
is of great interest in fluorescence microscopy to improve the image resolution. However,
the exact structure of the PSF is generally unknown and particularly difficult to construct due to
the ever-present optical aberrations in real systems. Moreover, the vector feature of the PSF is
important to consider for the correct prediction of intensity-based measurements. This includes
the polarization state and the direction of the energy flow in the image plane. The scalar model
for computing a PSF is common due to the ease of implementation but limited, whereas the vectorial
model is computationally more expensive but more accurate, especially for high numerical
apertures. In this work, we present new techniques for computing vector PSF using the Fourier
transform. The fast Fourier-transformation is a very handy tool to speed up PSF calculations, but
its pitfalls have to be avoided. We show that even without assuming radial symmetry in the field,
that these Fourier-based methods are computationally inexpensive. The accuracy and quality of
each model is compared with the well-known PSF modelled by Richards andWolf for the ability
of this last one to represent an ideal field with a very high accuracy.

-------------------------------------------------
Usage of the proposed method:

* Install DipImage under Matlab http://www.diplib.org/download
* Follow the instruction in the test.m file

-------------------------------------------------

Contact: Prof. Dr. R Heintzmann (heintzmann@gmail.com), University of Jena, Germany
-------------------------------------------------
Copyright (c) Stellenbosch University and Friedrich Schiller University Jena, 2022