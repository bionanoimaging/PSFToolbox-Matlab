% [FPlane,Aperture]=SimLens(BFPAperture,lambda,scales,NA,aplanar,smoothAperture)  : Simulates the Electric field vectors corresponding to Ex, Ey and Ez spatial frequency distribution in the focal plane of a lens
% BFPAperture : Either XY size or an XY-polarisation map, e.g. {polimX,polimY}. Default is linear polarization along X
% lambda : wavelength in the medium (supply lambda/n to account for the refractive index in the medium).
% scales : pixelsizes
% NA: numerical aperture
% aplanar: flag for the aplanar factor. 1 means multiply intensity with
% cos(alpha) as in illumination, -1 means divide intensity by cos(alpha) as in detection, 0 means no aplanatic factor (default), 2 precomensate twice, i.e. intensity by cos^2(alpha), for excitation when projecting onto a sphere in 3D (SincR method)
% smoothAperture : optional flag to define whether the aperture is a hard discretised circle or a smooth version obtained by Fourier-transforming the jinc function. default: 0 (hard aperture)
%                  if smoothAperture > 1, a Damping is applied the the jinc function
% AddParams : contains additional parameters about refractive indexes (ri) and thicknesses of the sample, coverlsip and immersion medium in real and design conditions. 
%                       The emission wavelength is function of the medium in which the light is propagating (i.e. ri of immersion medium in real condition) whereas the 
%                       opening aperture is limited by how the optical system was designed (i.e. function of the ri of immersion medium in design condition).
%                       e.g. AddParams=struct('ns',1.33,'ng',1.518,'ni',1.516,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5); 
%                         'ns': ri of sample, can be in a vector format if the sample is in a stratified medium.  First element of ns (if a vector) is the ri of where the emitter is at, it is the ri of the medium which the furthest from the coverslip. The last element of the vector is the ri of the medium closest to the coverslip
%                         'ng': ri of coverslip in real condition
%                         'ni': ri of immersion medium in real condition
%                         'ni0': ri of immersion medium in design condition
%                         'ts': thickness of the sample or position where the emitter is at, this should be with the same length as the vector of ns and within the same order as ns.
%                         'tg': thickness of coverslip in real condition
%                         'tg0': thickness of coverslip in design condition
%                         'wd': working distance which is supposed to be the thickness of the immersion medium in design condition
%
% Alternative usage:
% FPlane=SimLens(BFPAperture,PSFParam,ImageParam,aplanar,smoothAperture)  
% Here PSFParam.lambdaEm refers to the wavelength in vaccum and it will be automatically divided by PSFParam.n to obtain lambda in the medium
% BFPAperture can be empty [] here as well
%
% Example1:
% NA=0.6; lambda=500; scales=[100 100];
% FPlane=SimLens([256 256],[lambda NA],scales,-1); % the last argument here is defining the aplanatic factor
% 
% Example2:
% ImageParam=struct('Sampling',[100 100 100],'Size',[256 256 16]);
% PSFParam=struct('NA',0.6,'n',1.33,'lambdaEm',500);
% % AddParams=struct('RIReal',[1.33 1.518 1.516],'RIDesign',[1.518 PSFParam.n],'ThicknessReal',[2e3 1.7e5 1.5e5],'ThicknessDesign',[1.7e5 1.5e5]);
% AddParams=struct('ns',1.33,'ng',1.518,'ni',1.516,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
% orderlist={[3,1];[1 -1 1 60]}; % order of the Zernike phase: coma and tilt aberration
% [TiltComaPhase,rho,phi] = ZernikePoly(orderlist,ImageParam,PSFParam);
% FPlane=mSimLens([],PSFParam,ImageParam,1,1,AddParams,TiltComaPhase)  % simulates with smooth aperture and aplanatic factor 1 which means multiply the intensity at higher angle by cos(theta)
% [RPlane,I,FPlane]=mPropagateFieldBlock([],0.0,FPlane,PSFParam,ImageParam) % constructs the vectorial PSF and its intensity distribution
%
% Copyright (C) 2021, Rainer Heintzmann,  All rights reserved.
% This is a modified version of SimLens [edit: 23.08.21 - R. Dina]
% Aberration due to refractive index mismatch is included
%_________________________________________________________________________
function [FPlane,OpticalPathDiff]=mSimLens(BFPAperture,lambdaNA,scales,aplanar,smoothAperture,AddParams,AddPhase)
if nargin < 7 || isempty(AddPhase)
    AddPhase=0;
end

if nargin < 6
    smoothAperture=0;
end
if isstruct(lambdaNA)  % to enable the alternate usage
    PSFParam=lambdaNA;
    ImageParam=scales;
    if nargin<5 || isempty(smoothAperture)     
        smoothAperture=0;
    end
    if nargin < 4 || ~isfield(PSFParam,'Aplanatic')
        aplanar=[];
    else 
        aplanar=PSFParam.Aplanatic;
    end
    if isempty(BFPAperture)
        BFPAperture=ImageParam.Size(1:2);
    elseif isvector(BFPAperture) && numel(BFPAperture) == 1
            BFPAperture=newim(ImageParam.Size(1:2))+ BFPAperture;  % BFP Polarization is given by complex number
    end
    scales=ImageParam.Sampling;
elseif isvector(lambdaNA) % lamba is a vector containg as first element the emission wavelength and second element the NA
    PSFParam=struct('lambdaEm',lambdaNA(1), 'NA', lambdaNA(2),'n',1); % input lambda and NA are already in the corresponding medium where it should be
end

if nargin<4 || isempty(aplanar)
    aplanar=0;  % do not use the aplanatic factor
end
if isvector(BFPAperture) % User means size instead of BFP image
    BFPAperture=newim(BFPAperture)+1;  % assumes linear X polarisation only
end
BFPAperture=expanddim(BFPAperture,4);
if size(BFPAperture,3) ~= 1
    if size(BFPAperture,3) == 2
        BFPAperture=permute(BFPAperture,[1 2 4 3]);
    else
        error('BFP Aperture needs to have two elements along dimension 3 or 4');
    end
end

if size(BFPAperture,4) < 2
    BFPAperture=cat(4,BFPAperture,BFPAperture(:,:,:,0)*0);
end

% RI and thicknesses data
% -- start --
if nargin<6 || isempty(AddParams)  % i.e. if the structure is not given
    tp=1; ts=1; % transmission coefficients (Fresnel)
    TFactor=1;
    OpticalPathDiff=0;
    AddParams=[]; % for our next reference
    ndes=PSFParam.n; % ri in the immersion medium in a design or ideal condition
    nreal=ndes; % ri in the immersion medium in real condition is the same as in design condition
    J=1;
   nlist=[ndes, nreal];
else
    ndes=AddParams.ni0;
    nreal=AddParams.ni; % refractive index in real condition
    PSFParam.n=nreal; % the field is propagating in the immersion medium whose refractive index is given in real condition
%     PSFParam.n=ndes; % 
    nlist=[AddParams.ni, AddParams.ni0, AddParams.ng, AddParams.ng0, AddParams.ns];
end

% clear lambda NA
PSFParam.NA=min([PSFParam.NA, nlist]); % the minimum aperture here is to avoid any total internatl
NA=min(PSFParam.NA/ndes, PSFParam.NA/nreal); % the minimum aperture is defined account for what can be achieved by the optical system

lambda=PSFParam.lambdaEm / PSFParam.n;  % to obtain lambda in the medium

sz=size2d(BFPAperture);
ktotal=2*pi/lambda; 
ktotalsqr=ktotal*ktotal;
kxysqr=(abssqr(2*pi*xx(sz,'freq')/scales(1))+abssqr(2*pi*yy(sz,'freq')/scales(2)));
kzsqr=(ktotalsqr-kxysqr);
kzsqr(kzsqr < 0)=0;
kz=sqrt(kzsqr);

if ~isempty(AddParams) && isstruct(AddParams)
    if isfield(AddParams,'tilt')
        alpha = AddParams.tilt(1)*pi/180; cosa = cos(alpha); sina=sin(alpha); % elevation angle of the coverslip from y-axis about x-axis in counterclockwise direction and in degrees
        beta = AddParams.tilt(2)*pi/180;   cosb = cos(beta);  sinb=sin(beta);
        kx0=2*pi*xx(sz,'freq')/scales(1); 
        ky0=2*pi*yy(sz,'freq')/scales(2);
        kx=cosb*kx0 +sinb*ky0;
        ky=-sinb*cosa*kx0 + cosb*cosa*ky0 + sina*kz;
        kxysqr=abssqr(kx) + abssqr(ky);
        kzsqr=(ktotalsqr-kxysqr);
        kzsqr(kzsqr < 0)=0;
        kz=sqrt(kzsqr);
        J=abs(cosa - cosb*sina*div(ky,kz) + sinb*sina*div(kx,kz)); % Jacobian
    else 
        J=1;
    end
else
    J=1;
end

cosalpha=kz./ktotal; % angle in the immersion medium in real condition
sinalpha=sqrt(kxysqr)./ktotal;

HardAperture= kxysqr/ktotalsqr < NA*NA;
if (smoothAperture)
%         PSFParamAp=PSFParam; PSFParamAp.NA=NA; % define properly the amount of light going through the system
        jpsf=jincPSF(ImageParam,PSFParam,ndes); % n defines in PSFParam corresponds to the ri of the medium through which light propagates in reality and ndes defines the ri of the medium through which the medium is supposed to propagate as described by the manufacturer
    if smoothAperture>1  % means: Do Damping
        jpsf=DampEdge(jpsf,0.05,2,1,3);
    end
    Aperture= real(ft(jpsf));  % force it to be real
else
    Aperture= HardAperture;
end
    
% calculate the transmission coefficients and optical path difference if there are interfaces
if ~isempty(AddParams) && isstruct(AddParams)
        [OpticalPathDiff, tp, ts]=OPDTsp(AddParams,PSFParam,ImageParam); % 
%     else
%         [OpticalPathDiff, tp, ts]=OPDDirect(AddParams,ImageParam,PSFParam); % faster than the above // I prefer having two different functions here as it facilitates the job somehow than combining it
%     end
    % TFactor : ( n_t cos(theta_t) )./ ( n_i cos(theta_i) ) . This factor is
    % to account for normal incidence where n_t/n_i is the reciprocal of
    % the media wave impedance and cos(theta_t)/cos(theta_i) is needed from
    % expressing the power in the direction normal to the interface
end


if (smoothAperture)
    AMask=cosalpha>1e-6;  % to avoid division by zero
else
    AMask=HardAperture>0;
end
Aperture= dip_image(Aperture,'sfloat') .* AMask; % Aperture+0.0;

if any(-2:1:3==aplanar)
    pow = aplanar/2;
    Aperture(AMask) = Aperture(AMask) .* (cosalpha(AMask)).^pow;  
else
    error('Such an aplanatic factor is not allowed.')
end

FPlane=newim([sz 1 3],'scomplex');

myx=ramp(sz,1)./(sz(1).*scales(1));
myy=ramp(sz,2)./(sz(2).*scales(2));

myphiphi = atan2(myy,myx);
cosphi=cos(myphiphi); % phiphi(sz,'freq')
sinphi=sin(myphiphi);

% tp=1; ts=1;

BFPRadial=tp.*((BFPAperture(:,:,:,0).*cosphi+ BFPAperture(:,:,:,1).*sinphi)) *Aperture;
FPlane(:,:,:,2)=  sinalpha * BFPRadial; % Z-component of the vector 

BFPTangential=ts.*((BFPAperture(:,:,:,1).*cosphi - BFPAperture(:,:,:,0).*sinphi)) *Aperture;
RadialNew=cosalpha * BFPRadial;  % This is what is left over for the XY radial component

FPlane(:,:,:,0)=  cosphi * RadialNew - sinphi * BFPTangential; % X-component of the vector 
FPlane(:,:,:,1)=  sinphi * RadialNew + cosphi * BFPTangential; % Y-component of the vector 

% Add aberrations phase
FPlane=J.*FPlane.*exp(1i.*((ktotal.*OpticalPathDiff) - AddPhase)); % the reason why it is not needed to include ktotal anymore is .... and the negative sign is because .... (it works like this for now but the reason I am still figuring it out :D)

end

function res = div(num,denom) % division to avoid dividing by zero
    res=newim(denom);
    msk=denom~=0;
    if length(num)<length(denom) % is there a better way to check whether the numerator is within the same size as the denominator (matrix) or just a number
        res(msk)=num./denom(msk);
    else
        res(msk)=num(msk)./denom(msk);
    end
end
