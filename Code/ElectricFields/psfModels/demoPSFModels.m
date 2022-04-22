% This script contains an example illustrating how to use the 
% scalar and vectorial PSF models provided as part of the psfModels package.
% For more information, see the documentation of scalarPSF.m and vectorialPSF.m
% as well as refs. 1 and 2
%
% [1] F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
% [2] F. Aguet, Ph.D Thesis, Swiss Federal Institute of Technology, Lausanne (EPFL), 2009

% Francois Aguet 2013


% Objective parameters:
p.ti0 = 0.19e-3;        % working distance of the objective
p.ni0 = 1.518;          % immersion medium refractive index, design value
p.ni = 1.518;           % immersion medium refractive index, experimental value
p.tg0 = 0.17e-3;        % coverslip thickness, design value
p.tg = 0.17e-3;         % coverslip thickness, experimental value
p.ng0 = 1.515;          % coverslip refractive index, design value
p.ng = 1.515;           % coverslip refractive index, experimental value
p.ns = 1.33;            % sample refractive index
p.lambda = 550e-9;      % emission wavelength
p.M = 100;              % magnification
p.NA = 1.45;            % numerical aperture
p.pixelSize = 6.45e-6;  % physical size (width) of the camera pixels
p.sf = 3;               % (optional, default: 3) oversampling factor to approximate pixel integration
p.mode = 1;             % (optional, default: 1) if 0, returns oversampled PSF


zp = 1e-6; % z-position of the source
zv = zp + (-2000:100:4000)*1e-9; % position of the z-planes
nz = numel(zv);
nx = 41; % window width; the volume size is nx x nx x nz

% calculate models
spsf = scalarPSFFun([0 0 zp], zv, nx, p);
vpsf = vectorialPSF([0 0 zp], zv, nx, p);

nxi = size(spsf,1);
zi = floor((nz+1)/2);
xi = (nxi+1)/2; % center pixel coordinate for plotting

% lateral coordinate vector
xa = (-xi+1:xi-1)*p.pixelSize/p.M*1e6;
if ~p.mode
    xa = xa/p.sf;
end

figure;
colormap(gray(256));
subplot(1,2,1);
imagesc(xa, zv*1e6, sqrt((squeeze(spsf(:,xi,:))))'); axis equal tight;
title('Scalar model');
xlabel('x [um]');
ylabel('z [um]');

subplot(1,2,2);
imagesc(xa, zv*1e6, sqrt((squeeze(vpsf(:,xi,:))))'); axis equal tight;
title('Vectorial model');
xlabel('x [um]');
ylabel('z [um]');

