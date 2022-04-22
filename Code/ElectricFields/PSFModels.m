% function h = PSFModels(ImageParam,PSFParam, scal_vec,AddParams)
% Generate a perfect PSF scalar or vector based on the algorithm of Aguet 
% Values input in the arguments are in nanometer but they will be converted in meter in this function. 
% If there is a need of change of the microscope parameters or any other parameters, it is better to directly use scalarPSF (scalar) or vectorialPSF(vector)
% Source: Aguet, François. "Super-ResolutionFluorescence MicroscopyBasedonPhysical Models." (2009).
% link: http://www.francoisaguet.net/software.html
% Example:
% ImageParam=struct('Sampling',[80 80 150],'Size',[256 256 32]);
% PSFParam=struct('NA',1.4,'n',1.518,'MinOtf',1.2e-3,'lambdaEm',580,'Aplanatic',-1);
% h_scal=PSFModels(ImageParam,PSFParam, 'ScalarPSF');
% h_vec=PSFModels(ImageParam,PSFParam, 'VectorPSF');

function h = PSFModels(ImageParam,PSFParam, scal_vec,AddParams)
% scales=ImageParam.Sampling.*1e-9;
xp = [0 0 0]; % Source position, 3-element vector [xp yp zp]
if nargin<4 || isempty(AddParams)
    n = PSFParam.n;
    p.ti0 = 0.15e-3;    % working distance of the objective 
    p.ni0 = n;             % immersion medium refractive index, design value
    p.ni = n;               % immersion medium refractive index, experimental value
    p.tg0 = 0.17e-3;   % coverslip thickness, design value
    p.tg = 0.17e-3;     % coverslip thickness, experimental value
    p.ng0 = n;             % coverslip refractive index, design value
    p.ng = n;               % coverslip refractive index, experimental value
    p.ns = n;               % sample refractive index
else
    p = AddParams;
    p.ti0 = p.wd.*1e-9;  % working distance of the objective
    p.tg0 = p.tg0*1e-9; % coverslip thickness, design value in meters
    p.tg = p.tg*1e-9;     % coverslip thickness, experimental value in meters
    xp(3)= p.ts*1e-9;    % Source position, 3-element vector [xp yp zp] // I assume this is how it should be 
%     p.pZ = p.ts*1e-9;
end

p.lambda = PSFParam.lambdaEm.*10^(-9);%/p.ni;      % emission wavelength in the real medium
p.M = 1;                % magnification
p.NA = PSFParam.NA;%/p.ni; %min(PSFParam.NA/p.ni0, PSFParam.NA/p.ni);%PSFParam.NA;     % numerical aperture
p.pixelSize = ImageParam.Sampling(1).*10^(-9);  % physical size (width) of the camera pixels in metres
p.sf = 3;               % (optional, default: 3) oversampling factor to approximate pixel integration
p.mode = 1;             % (optional, default: 1) if 0, returns oversampled PSF

sz = ImageParam.Size;
nsz = sz + mod(sz+1,2);
zcale=1;%0.95;%sz(3)/nsz(3);
z = double(xx(sz(3)+mod(sz(3)+1,2)).*ImageParam.Sampling(3).*zcale.*10^(-9));% Vector of z-plane positions
nx =nsz(1);              % xy-window which size is forced it to be odd

% p.pixelSize = (p.pixelSize).*ImageParam.Size(1)/nx;

switch scal_vec
    case 'ScalarPSF'
        [h, dxp, dyp, dzp] = scalarPSF(xp, z, nx, p);
    case 'VectorPSF'
        h = vectorialPSF(xp, z, nx, p);
    otherwise
        error('Choose ScalarPSF or VectorPSF as last ')
end

h = extract(dip_image(h), ImageParam.Size);
h = h./sum(h);