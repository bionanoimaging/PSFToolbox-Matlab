% [h,amp4d]=GenericPSFSim(ImageParam,PSFParam,Method) : This function simulates a PSF using one of several possible approaches.
% 
% INPUT
%   ImageParam: a structure containing the following fields:
%       Sampling := voxel size represented in a vector format of three elements along x, y and, z
%       Size := window size along x, y and, z  
%   PSFParam: a structure containing the following minimum fields:
%       polarization = 'circular';  % default if not mentionned in the structure
%       polarization = 'linearX' or 'linearY' or 'radial' or 'azimuthal'; % are other options
%       Aplanatic = -1;  % -1 means emission, 0: no aplanatic factor, 1: excitation. 
%           Note: Some function can also support the PSFParam: 'Aplanatic', 2: % precomensate twice for excitation when projecting onto a sphere in 3D (SincR method), 3: % this means that the intensity will be modified with cosalpha^2 
%       NA := numerical aperture of the objective lens
%       n := refractive index of the immersion medium
%       lambdaEm := emission wavelength
%   Method: a string for selecting the method to use in the simulation. Choices are:
%       'RichardsWolff' (default): Using the method from the Richards & Wolff paper [1, 2].
%       'RichardsWolffInt' : A fast method, if only the intensity is needed [1, 2].
%       'RWFast': Reworked Richards & Wolf bast on an chirp-Z transform along kz [3, 4].
%       'SlicePropagation' : Simulates the in-plane ATF based on FFTs and the jinc function and propagates it to all planes [4].
%       'CZT' chirped-Z transform-based (== Zoomed FFT) method avoiding the out-of-focus wrap-around problems of SlicePropagation [3, 4].
%       'VolumeShell' : Creates a part of the McCutchen Pupil directly in Fourier space using pre-computed Fourier-space interpolators. Pretty good in the central region of the PSF but (intentionally) degrading near the edge of the volume [4].
%       'SincR' :  Creates the Fourier shell by a 3D FFT of a sinc(abs(R))
%       function and modifies it from there on [4].
%  AddParams: structure having as fields
%                       e.g. AddParams=struct('ns',1.33,'ng',1.518,'ni',1.516,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5); 
%                         'ns': ri of sample, can be in a vector format if the sample is in a stratified medium.  First element of ns (if a vector) is the ri of where the emitter is at, it is the ri of the medium which the furthest from the coverslip. The last element of the vector is the ri of the medium closest to the coverslip
%                         'ng': ri of coverslip in real condition
%                         'ni': ri of immersion medium in real condition
%                         'ni0': ri of immersion medium in design condition
%                         'ts': thickness of the sample or position where the emitter is at, this should be with the same length as the vector of ns and within the same order as ns.
%                         'tg': thickness of coverslip in real condition
%                         'tg0': thickness of coverslip in design condition
%                         'wd': working distance which is supposed to be the thickness of the immersion medium in design condition
%         * if AddParams == [] or not added to the function, the code will run as there is no interfaces and no aberrations due to that, n = PSFParam.n will be used troughout the calculation in this case.
%  AddPhase : additional phase that can be added into the function. In Method = 'ZoomedPupilProp', it is prefereable to have this AddPhase variable as a cell for computing the ZernikePoly instead while it can only be in 2D image for the other methods. 
% 
% OUTPUT
%   h : resulting intensity PSF
%   amp4d : 4D complex amplitude distribution (if supported).
% 
% REFERENCES
% [1] Wolf E. Electromagnetic diffraction in optical systems-I. An integral representation of the image field. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences. 1959 Dec 15;253(1274):349-57.
% [2] Richards B, Wolf E. Electromagnetic diffraction in optical systems, II. Structure of the image field in an aplanatic system. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences. 1959 Dec 15;253(1274):358-79.
% [3] Rabiner L, Schafer RW, Rader C. The chirp z-transform algorithm. IEEE transactions on audio and electroacoustics. 1969 Jun;17(2):86-92.
% [4] Calculating Point Spread Functions - Methods, Pitfalls and Solutions
% 
% See README.pdf for license agreement and disclaimer.
% 
% Last modified: 23/08/2021 by R. Dina
% Copyright (c) 2021
% Rainer Heintzmann
% Nano-imaging Research group, 
% Friedrich Schiller University and  Leibniz Institute of Photonic Technology, Jena, Germany
% Stellenbosch University, South Africa
% 
% EXAMPLES:
% ImageParam=struct('Sampling',[40 40 50],'Size',[128 128 32]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
% AddParams=struct('ns',1.33,'ng',1.518,'ni',1.518,'ng0',1.518,'ni0',1.518,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
% [h,amp4d]=GenericPSFSim(ImageParam,PSFParam,'VolumeShell',AddParams)
%_________________________________________________________________________
%
function [h,amp3d]=GenericPSFSim(ImageParam,PSFParam,Method,AddParams,AddPhase,noCheck)
if nargin<4 || isempty(AddParams) % parameters containing the ri and thicknesses parameters
    AddParams=[];
end

if nargin<5 || isempty(AddPhase) %|| all(AddPhase==0)
    PlusPhase=0;
elseif iscell(AddPhase)
%     switch Method
%         case 'CZT'
%             PlusPhase = AddPhase; %
%         otherwise
            PlusPhase = ZernikePoly(AddPhase,ImageParam,PSFParam); 
%     end
else
    PlusPhase = AddPhase; 
%     switch Method
%         case 'CZT'
%             PlusPhase = []; %
%             fprintf('WARNING: Additional phase is not added. Input the Zernike polynomial instead for CZT method \n');
%         otherwise
%             PlusPhase = AddPhase;
%     end
end

if nargin < 6 
    noCheck=0;
end
h=[];
amp3d=[];
% if ~CheckPSFParams(ImageParam,PSFParam,noCheck) % make some sanity checks
%     return;
% end
% noCheck = 1; % just ignore all other checks.

if isfield(PSFParam,'Mode')
    myMode = PSFParam.Mode;
else
    myMode = 'Widefield';
end

if strcmp(myMode,'Confocal')
        h=PSF_Confocal(ImageParam,PSFParam,PSFParam.PinholeSize,Method,1);
        amp3d =[];
        return;
end

%% First calculate the wiedefield emission PSF

% add conditions for AddParams here if any field in the structure is missing

if ~isfield(PSFParam, 'polarization')
    PSFParam.polarization='circular';  % default
end

if ~isfield(PSFParam, 'circPol')
    PSFParam.circPol = 1;
end

if ~isfield(PSFParam, 'computeASF')
    PSFParam.computeASF = 1;
end

if ~isfield(PSFParam, 'Aplanatic')
    PSFParam.Aplanatic = -1;  % default is emission PSF
end

if PSFParam.NA > PSFParam.n
    error('Impossible PSF: PSFParam.NA has to be below PSFParam.n\n');
end

if nargin < 3
    Method='RichardsWolff';
end
if nargin < 2 || isempty(PSFParam)
    PSFParam=struct('NA',1.2,'n',1.33,'MinOtf',1.2e-3,'lambdaEm',520);
end
if nargin < 1 || isempty(ImageParam)
	ImageParam=struct('Sampling',[40 40 100],'Size',[256 256 16]);
end

if strcmp(Method,'RichardsWolff')  % to remain compatible with old (wrong spelling)
    Method = 'RichardsWolf';
end
if strcmp(Method,'RichardsWolffInt')  % to remain compatible with old (wrong spelling)
    Method = 'RichardsWolfInt';
end

tic
switch Method
    case 'RichardsWolfInt'
        switch PSFParam.polarization
            case 'circular'
                PSFParam.circPol=1;
%             case 'linearX' % DOES NOT WORK for Intensity-only method
%                 PSFParam.circPol=0;
             otherwise
                error('For RichardsWolfInt polarization, only cicular and linearX are supported.')
        end
        % Note below that lambdaEx and lambdaEm are set to the emission wavelength to get the correct result
        if ~isfield(PSFParam,'Aplanatic') || PSFParam.Aplanatic==-1 % (emission)
            h=kSimPSF({'scaleX',ImageParam.Sampling(1);'scaleY',ImageParam.Sampling(2);'scaleZ',ImageParam.Sampling(3);  ...
            'sX',ImageParam.Size(1);'sY',ImageParam.Size(2);'sZ',ImageParam.Size(3);'na',PSFParam.NA;'ri',PSFParam.n; ...
            'lambdaEm',PSFParam.lambdaEm;'lambdaEx',PSFParam.lambdaEm;'circPol',PSFParam.circPol});
        elseif PSFParam.Aplanatic==1 % excitation:  Calucate a confocal PSF with open pinhole!
            error('In the RickardsWolfInt method only the PSFParam.Aplanatic values -1 (emission) is supported.');  % The confocal does not calculate properly
            h=kSimPSF({'scaleX',ImageParam.Sampling(1);'scaleY',ImageParam.Sampling(2);'scaleZ',ImageParam.Sampling(3);  ...
            'sX',ImageParam.Size(1);'sY',ImageParam.Size(2);'sZ',ImageParam.Size(3);'na',PSFParam.NA;'ri',PSFParam.n; ...
            'lambdaEm',PSFParam.lambdaEm;'lambdaEx',PSFParam.lambdaEm;'circPol',PSFParam.circPol;'confocal',1;'pinhole',1e6});  % There is a BUG in circpol and amplitudes!!
        else
            error('In the RichardsWolfInt method only the PSFParam.Aplanatic values -1 (emission) and 1 (excitation) are supported.')
        end
    case 'RichardsWolf'
        if isempty(AddParams)
            facNA=1; facLambda=1;
        else
            facNA=PSFParam.n/(max([AddParams.ni0 AddParams.ni])); % The maximum allowed NA is related to the refractive index of the immersion medium in real and design condition
            facLambda=PSFParam.n/AddParams.ni; % The wavelength depends on the medium in which the field is propagating
        end
%         if ~isempty(AddParams)
%             PSFParam.n=AddParams.ni;
%             PSFParam.NA=PSFParam.NA*PSFParam.n/(max([AddParams.ni0 AddParams.ni]));
% %             PSFParam.lambdaEm=PSFParam.lambdaEm*PSFParam.n/AddParams.ni;
%         end
        % Note below that lambdaEx and lambdaEm are set to the emission wavelength to get the correct result
        if ~isfield(PSFParam,'Aplanatic') || PSFParam.Aplanatic==-1 % (emission)
            amp3dX=kSimPSF({'scaleX',ImageParam.Sampling(1);'scaleY',ImageParam.Sampling(2);'scaleZ',ImageParam.Sampling(3);  ...
            'sX',ImageParam.Size(1);'sY',ImageParam.Size(2);'sZ',ImageParam.Size(3);'na',PSFParam.NA*facNA;'ri',PSFParam.n; ...
            'lambdaEm',PSFParam.lambdaEm*facLambda;'lambdaEx',PSFParam.lambdaEm*facLambda;'circPol',0;'computeASF',1});  % There is a BUG in circpol==1 and amplitudes!!
            amp3dX=reshape(amp3dX,[size(amp3dX,1) size(amp3dX,2) size(amp3dX,3) size(amp3dX,5)]);
        else
            error('In the Richards Wolf method only the PSFParam.Aplanatic values -1 (emission) is supported.')
        end
        % amp3dY=extract(flipud(permute(amp3dX,[2 1 3 4 5])),size(amp3dX),[floor(size(amp3dX,2)/2) floor(size(amp3dX,1)/2)-mod(size(amp3dX,1)+1,2) floor(size(amp3dX,3)/2) floor(size(amp3dX,4)/2) floor(size(amp3dX,5)/2)]);
        if size(amp3dX,1)==size(amp3dX,2)
            amp3dY=permute(amp3dX,[2 1 3 4]);
        else
            amp3dY=kSimPSF({'scaleX',ImageParam.Sampling(2);'scaleY',ImageParam.Sampling(1);'scaleZ',ImageParam.Sampling(3);  ...
            'sX',ImageParam.Size(2);'sY',ImageParam.Size(1);'sZ',ImageParam.Size(3);'na',PSFParam.NA;'ri',PSFParam.n; ...
            'lambdaEm',PSFParam.lambdaEm;'lambdaEx',PSFParam.lambdaEm;'circPol',0;'computeASF',1});  % There is a BUG in circpol==1 and amplitudes!!
            amp3dY=reshape(amp3dY,[size(amp3dY,1) size(amp3dY,2) size(amp3dY,3) size(amp3dY,5)]);
            amp3dY=permute(amp3dY,[2 1 3 4]);
        end
        switch PSFParam.polarization
            case 'circular'
                amp3d=amp3dX+1i*amp3dY(:,:,:,[1 0 2]);
            case 'linearX'
                amp3d=amp3dX;
            case 'linearY'
                amp3d=amp3dY(:,:,:,[1 0 2]);
            otherwise
                error('for RichardsWolff PSF calculations, only circular, linearX and linearY are supported.')
        end
        % Why is the Y-pol different than the X-pol??
        % add phase and Fresnel coefficients 23.08.21
        if ~isempty(AddParams) 
            ktotal = 2*pi/(PSFParam.lambdaEm/AddParams.ni);
            TheoPhase =  ktotal.*OPDTsp(AddParams,PSFParam,ImageParam); 
%             TheoPhase = ZernPhase(AddParams,ImageParam,PSFParam,6) ;%
        else
            AddParams.ni=PSFParam.n;
            TheoPhase = 0; 
        end
%         flipdim(PlusPhase,1)
        TotalPhase = -PlusPhase + TheoPhase; % Plusphase shoud also be an array not an integer value. Otherwise it will be considered in the same way as if it is 0. 
        if any(size(TotalPhase)>1)
            RWPupil = ft3d(amp3d);%.*exp(-1i.*TotalPhase);
            amp3d = ift3d(RWPupil.*exp(1i*TotalPhase));
        end
        
        % end add
        h=squeeze(sum(abssqr(amp3d),[],4));
    case 'SlicePropagation'
        [h,amp3d]=SlicePSFSim(ImageParam,PSFParam,AddParams,PlusPhase);
    case 'SPIterative'
        [h,amp3d]=PSFsIterative(ImageParam,PSFParam,AddParams,PlusPhase,1);  % last argument 1 means
    case 'RWFast' % uses the facto that you can write eq. 3.2 in teh RicharsWolf paper as an FFT.        
        [h,amp3d]=RWFastProp(PSFParam,ImageParam);  % centered at zero        
    case 'CZT'
        [amp3d,h,ZoomedPupil]=ZoomedPupilProp(PSFParam,ImageParam,0,[],[],AddParams,PlusPhase);  % centered at zero
    case 'VolumeShell'  % Use the code with the interpolation in Fourier space
        % error('Not implemented yet.')
        [h,amp3d]=VolumePSFSim(ImageParam,PSFParam,[],1,AddParams,PlusPhase);  % calculates an extended region and cuts it
    case 'SincR'
        PSFParam.Aplanatic=PSFParam.Aplanatic+2; % the additional aplanatic factor (+2) is to compensate for the projection of the propagator (sinc function) onto a sphere in 3D
%         PSFParam.Aplanatic=10;
        [h,amp3d]=mPSFSimSincR(ImageParam,PSFParam,1,AddParams,PlusPhase);  %  calculates an extended region and cuts it out // older version is PSFSimSincR
    case 'GibsonLanni' % based on MicroscopePSF from other group.
 
        if isempty(AddParams)
            n = PSFParam.n;
            params.ti0 = 0.15e-3;    % working distance of the objective 
            params.ni0 = n;             % immersion medium refractive index, design value
            params.ni = n;               % immersion medium refractive index, experimental value
            params.tg0 = 0.17e-3;   % coverslip thickness, design value
            params.tg = 0.17e-3;     % coverslip thickness, experimental value
            params.ng0 = n;            % coverslip refractive index, design value
            params.ng = n;              % coverslip refractive index, experimental value
            params.ns = n;              % sample refractive index
            params.pZ = 0;              % source z-position
        else
            params = AddParams;
            params.ti0 = params.wd.*1e-9;  % working distance of the objective % express values in meters
            params.tg0 = params.tg0*1e-9; % coverslip thickness, design value in meters
            params.tg = params.tg*1e-9;     % coverslip thickness, experimental value in meters
            params.pZ = params.ts*1e-9;    % source z-position
        end

        params.size = ImageParam.Size+mod(ImageParam.Size+1,2); % to avoid even window size, algorithm has problem with even window size but works good with odd numbers
        params.resLateral = ImageParam.Sampling(1)*10^(-9);
        params.resAxial = ImageParam.Sampling(3)*10^(-9);
        params.lambda = PSFParam.lambdaEm.*10^(-9);
        params.NA = PSFParam.NA;
        params.M = 1;
        
        h = dip_image(MicroscPSF(params));
        h = extract(h,ImageParam.Size); % back to the wanted size
        h = h./sum(h);
        
    case 'ScalarPSF' % based on PSFModels from other group % PSFmodels-py https://pypi.org/project/psfmodels/
        h = PSFModels(ImageParam,PSFParam, 'ScalarPSF',AddParams);
        
    case 'VectorPSF' % from other group PSFmodels-py https://pypi.org/project/psfmodels/
        h=PSFModels(ImageParam,PSFParam, 'VectorPSF',AddParams);
        
    otherwise
        error('Unknown PSF generation method. Try: RichardsWolffInt, RichardsWolffInt, SlicePropagation, CZT, VolumeShell, SincR')
end
h=squeeze(h);
% IS=sum(MidPos(h,3));  % normalization factor as the integrated intensity in focal plane 
IS=sum(h); % normalize by the integrated intensity
h=h./IS; % normalize to the integrated intensity in the focal plane
if exist('amp3d','var') % this variable does not exist for some models
    amp3d=amp3d./sqrt(IS);    % normalize afterwards to avoid indidual slice normalization
end
toc

switch myMode
    case 'Widefield'
        return
    case '2-Photon'
        h = h .^ 2.0;
    case 'Lightsheet'
        psf_det = h;

        LS_PSFParam=PSFParam; % struct('NA',NA_illu,'n',ri,'lambdaEm',lambdaEx);
        LS_PSFParam.LambdaEm = LS_PSFParam.LambdaEx;
        % change Dina 16.12.21 as the next line LS_ImageParam leads to an error
        %         LS_ImageParam=struct('Sampling',LS_PSFParam.Sampling([1,3,2]),'Size',LS_PSFParam.Size([1,3,2]));  % reaches border at 262 µm
        LS_ImageParam=struct('Sampling',ImageParam.Sampling([1,3,2]),'Size',ImageParam.Size([1,3,2]));  % reaches border at 262 µm
        CalcMethod=Method;%
        LS_PSFParam.Mode='Widefield';
        % end change
        psf_illu = GenericPSFSim(LS_ImageParam,LS_PSFParam,CalcMethod);
        h = permute(sum(psf_illu,[],3),[1,3,2]) .* psf_det;
    otherwise
        error('unknown mode')
end

end
