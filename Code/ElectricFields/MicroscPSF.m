function PSF = MicroscPSF(params)
%MICROPSF Compute the 3D PSF model described by Gibson and Lanni (JOSA 1992).
%   PSF = MICROPSF(params) return a 3D PSF given parameters
%
%   Parameters include:
%   (1) image properties
%           'size'  -  the size of the 3D PSF, e.g. params.size = [256 256 128];
%   (2) precision control
%           'numBasis'      - the number of approximation basis, default '100'
%           'numSamp'       - the number of sampling to determine the basis
%                             coefficients, default '1000'
%           'overSampling'  - the oversampling ratio, default 2
%   (3) microscope parameters
%        'NA'        - numerical aperture of the microscope, default 1.4
%        'lambda'    - Emission wavelength in vacuum, default 610nm
%        'M'         - magnification factor, default 100
%        'ns'        - specimen refractive index (RI), default 1.33
%        'ng0'       - coverslip RI, design value, default 1.5
%        'ng'        - coverslip RI, experimental, default 1.5
%        'ni0'       - immersion RI, design value, default 1.5
%        'ni'        - immersion RI, experimental, defualt 1.5
%        'ti0'       - working distance, design, default 150um
%        'tg0'       - coverslip thickness, design value, default 170um
%        'tg'        - coverslip thickness, experimental, default 170um
%        'resLateral' - lateral pixel size, default 100nm
%        'resAxial'  - axial pixel size, default 250nm
%        'pZ'        - position of particle, default 2000nm
%
%   Reference:
%       [1] Gibson, S.F. & Lanni, F., 1992.
%           Experimental test of an analytical model of aberration in an
%           oil-immersion objective lens used in three-dimensional light
%           microscopy. JOSA A, 9(1), pp.154-166.
%       [2] Li, J., Xue, F. and Blu, T. Fast and accurate 3D PSF
%           computation for fluorescence microscopy. JOSA A. Accepted.
%
%   See also AUX_BESSEL, AUX_SHOWPSF
%
%   Acknowledgement: PSFgenerator (http://bigwww.epfl.ch/algorithms/psfgenerator/)
%
%   Copyright © Jizhou Li, Feng Xue and Thierry Blu, 2017
%   Update date: 4 May, 2017

warning off;
if ~isfield(params,'size')
    error('Please set the size of PSF model');
end

sizeData = params.size;
params.nx = sizeData(1);
params.ny = sizeData(2);
params.nz = sizeData(3);

%% default parameters
if ~isfield(params,'numBasis')
    params.numBasis = 100;
end
if ~isfield(params,'numSamp')
    params.numSamp = 1000;
end
if ~isfield(params,'fastcom')
    params.fastcom = 0;
end
if ~isfield(params,'overSampling')
    params.overSampling = 2;
end
if ~isfield(params,'NA')
    params.NA = 1.4;
end
if ~isfield(params,'lambda')
    params.lambda = 610e-9;
end
if ~isfield(params,'M')
    params.M = 100;
end
if ~isfield(params,'ns')
    params.ns = 1.33;
end
if ~isfield(params,'ng0')
    params.ng0 = 1.5;
end
if ~isfield(params,'ng')
    params.ng = 1.5;
end
if ~isfield(params,'ni0')
    params.ni0 = 1.5;
end
if ~isfield(params,'ni')
    params.ni = 1.5;
end
if ~isfield(params,'ti0')
    params.ti0 = 150e-6;
end
if ~isfield(params,'tg0')
    params.tg0 = 170e-6;
end
if ~isfield(params,'tg')
    params.tg = 170e-6;
end
if ~isfield(params,'resLateral')
    params.resLateral = 100e-9;
end
if ~isfield(params,'resAxial')
    params.resAxial = 250e-9;
end
if ~isfield(params,'pZ')
    params.pZ = 2000e-9;
end

x0 = (params.nx-1)/2;
y0 = (params.ny-1)/2;
xp = x0; yp=y0;
maxRadius = round(sqrt((params.nx - x0).^2 + (params.ny - y0).^2)) + 1; % maximum radius along the diagonal centered at the zero central of the window // unit in pixel
R = [0:params.overSampling*maxRadius-1]./params.overSampling; % sampling here is 1/overSampling and length is overSampling*maxRadius // unit in pixel

Ti = params.ti0 + params.resAxial*([0:params.nz-1] - ((params.nz - 1.0) / 2.0)); % this is related to thickness, ri, z, and OPD

a = 0; % lower integration limit 
b = min([1, params.ns/params.NA, params.ni/params.NA, params.ni0/params.NA,...
    params.ng0/params.NA,params.ng/params.NA]); % upper integration limit (borne d'integration)
% b = min([1,params.NA/ params.ns, params.NA/params.ni, params.NA/params.ni0,...
%     params.NA/params.ng0,params.NA/params.ng]); % upper integration limit (borne d'integration)

L = params.numSamp; % this defines the number of coefficient parameters denoted by c_m in the manuscript 

Rho = linspace(a,b,L)'; % uniform samples points of rho_l 

%% 1. approximate function exp(iW) as Bessel series
NN = params.numBasis;    % number of Bessel series basis
k0 = 2*pi/params.lambda; % wavenumber

r = R*params.resLateral;  % radial position in nm // origin point is zero center of image and end point is edge of the diagonal // length : overSampling*maxRadius*resLateral

A = k0*params.NA*r; % argument of the Bessel function in Eq. 2 of the manuscript // denoted by beta in the Eq. of R_m in Eq. 5
A2 = A.^2; % beta^2

Ab = A.*b; % beta*a in the manuscript

%%%%%%%
% min wavelength
k00 = 2*pi/(545e-9);
factor1 = k0./k00;
% max numerical aperture
NA0 = 1.4;
factor2 = params.NA./NA0;
%%%%%%%
an = (3*[1:NN]-2); % empirical parameter that was defined as 3m-2 in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
an = an.*(factor1)*(factor2); % empirical parameters as sigma_m

anRho = bsxfun(@times,an,Rho); % sigma_m * rho // argument of the Bessel function in Eq. 4
J = besselj(0, anRho); % Bessel function J_0(sigma_m rho_l)

J0A = besselj(0, Ab);
J1A = A.*besselj(1, Ab);

anJ0A = bsxfun(@times,J0A,an');
anb = an.*b; % sigma_m * a (a being the upper integration limit in the notation in the manuscript)
an2 = an.^2; % sigma_m^2
B1anb = besselj(1, anb);
B0anb = besselj(0, anb);

Ele = bsxfun(@times,anJ0A,B1anb') - bsxfun(@times,J1A,B0anb'); % numerator of R_m(r,p) in Eq. 5
domin = bsxfun(@minus, an2', A2); % denominator 
Ele = Ele.*b./domin; % R_m(r,p)

C1 = params.ns*params.pZ; % ns: sample ri and params.pZ is sample thickness from the coverslip
C2 = params.ni*(Ti - params.ti0);
C3 = params.ng*(params.tg - params.tg0);

OPDs = C1*sqrt(1-(params.NA*Rho/params.ns).^2);
OPDi = bsxfun(@times, C2,sqrt(1-(params.NA*Rho/params.ni).^2));
OPDg = C3*sqrt(1-(params.NA*Rho/params.ng).^2);

OPD = bsxfun(@plus, OPDi, OPDs+OPDg); % optical path difference

% determine the coefficients
W = k0*OPD;
Ffun = cos(W) + 1i*sin(W); % this is exp(iW) // size L (sample points of rho) x z_range
Ci = J\Ffun; % coefficient c_m , \ is a matrix left division // size of J is NN x L and size of Ci is NN x z_range

%% 2. get PSF in each slice
ciEle = Ele'*Ci;
PSF0 = abs(ciEle).^2;

%% 3. apply axial asymmetry
%
% The 2D component is resampled to a Cartesian grid using
% piecewise-linear interpolation
[X,Y] = meshgrid(0:params.nx-1,0:params.ny-1);
rPixel = sqrt((X-xp).^2 + (Y-yp).^2);
index = floor(rPixel*params.overSampling);
disR = (rPixel - R(index+1))*params.overSampling;
index1 = index+1;
index2 = index+2;
disR1 = 1-disR;

for zi = 1:params.nz
    h = PSF0(:,zi);
    slice = h(index2).*disR +  h(index1).*disR1;
    PSF(:,:,zi) = slice;
end

PSF = PSF./max(PSF(:));

end
