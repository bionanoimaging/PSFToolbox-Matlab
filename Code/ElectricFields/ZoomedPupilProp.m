% [Propagated,I,ZoomedPupil]=ZoomedPupilProp(OrigPupil,PSFParam,ImageParam,ZDist,ZoomedPupil,WantedZoom,AddParams,AddPhase) : propagtes a (4D) pupil by a distance using chirped Z-transforms
% INPUT
%     OrigPupil : original Pupil (not resampled)
%     PSFParam : List of PSF parameters
%     ImageParam: List of image parameters
%     ZDist: Z-Position from the focus to propagate to. If ImageParam.Size(3) is larger than one, a stack is generated with ZDist being the position of the middle (floor(size/2)) slice.
%     AddParams: this includes the ri and thicknesses of samples, coverslip and immersion medium. See help GenericPSFSim.m to learn more
%     AddPhase: can be a direct additional phase or a cell containing the Zernike order for phase aberration which formati is {[n m aberration_coeff azimuth_rotation_angle]}.
%                      Since this parameter goes into the loop where there might be different zooming factor, a care must be taken when inputing it. The best might be to just input the Zernike order in cells and let the code manages it. 
% 
% OUTPUT
%     Propagated: Amplitude distribution
%     I: Intensity distribution
%     ZoomedPupil : resampled version of the pupil (if not empty [])
% 
% EXAMPLE:
% 
% ImageParam=struct('Sampling',[40 40 100],'Size',[256 256 64]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
% AberrationList={[4 0 1];[1 -1 1 60]}; % spherical aberration of power 1 and tilt by 60degrees 
% AddParams=struct('ns',1.518,'ng',1.518,'ni',1.518,'ng0',1.518,'ni0',1.518,'ts',0,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
% tic;[Propagated,I,ZoomedPupil]=ZoomedPupilProp(PSFParam,ImageParam,0.0,[], [],AddParams,AberrationList);t1=toc
% tic;[PropagatedNonAberr,INonAberrated,ZoomedPupilNonAberr]=ZoomedPupilProp(PSFParam,ImageParam,0.0);t2=toc
% cat(4,I,INonAberrated)
%
% Edit 30.06.21 -- R. Dina: wavelength emission is in the ri of the immersion medium in real condition whereas the opening aperture (NA) is defined as how the 
%                                        optics were manufactured (ri of immersion medium in design condition). Functions which use these are AbbeMaxRadiusFromPSF.m ;
%                                        jincPSF ; ...
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [Propagated,I,ZoomedPupil]=ZoomedPupilProp(PSFParam,ImageParam,ZDist,ZoomedPupil, WantedZoom,AddParams,AddPhase)

if nargin<7
    AddPhase=0;
end

if nargin<6
    AddParams=[];
end

if nargin < 5 || isempty(WantedZoom)
    WantedZoom=[];
end

if ~isfield(PSFParam,'Aplanatic')
    PSFParam.Aplanatic = 1;  % 1 means emission and multiply intensity by cos(theta), 0: no aplanatic factor
end

% add effect of ri mismatch -----
if nargin<6 || isempty(AddParams)  % i.e. if the structure is not given
    AddParams=[]; % for our next reference
    ndes=PSFParam.n; % ri in the immersion medium in a design or ideal condition
    nreal=ndes; % ri in the immersion medium in real condition is the same as in design condition
else
    ndes=AddParams.ni0;
    nreal=AddParams.ni;
    PSFParam.n=AddParams.ni; % the field is propagating in the immersion medium whose refractive index is given in real condition
end
% end -----

if numel(ImageParam.Size)>2 && ImageParam.Size(3) > 1
    Propagated=newim([ImageParam.Size,3],'scomplex');
    I=newim(ImageParam.Size);
    ZoomedPupil=[];  % only for the first call.
    ImageParam2D=ImageParam;
    ImageParam2D.Size(3)=1;
    MidP=floor(ImageParam.Size(3)/2);
    for n=0:ImageParam.Size(3)-1
        MyDist=ZDist+(n-MidP)*ImageParam.Sampling(3);
        [myP,myI,ZoomedPupil]=ZoomedPupilProp(PSFParam,ImageParam2D,MyDist,ZoomedPupil,[],AddParams,AddPhase);
        Propagated(:,:,n,:)=myP;
        I(:,:,n)=myI;
    end
    return
end

%*******
lambda=PSFParam.lambdaEm ./ PSFParam.n; % emission wavelength in the immersion medium in real condition

[IncoherentRadius,PupilRadius]=AbbeMaxRadiusFromPSF(PSFParam,ImageParam,0,ndes);

OrigSize=ImageParam.Size(1:2);
HalfImgSize=floor(OrigSize/2);
AllowedZoom = HalfImgSize ./ PupilRadius;

if nargin < 5 || isempty(WantedZoom)  % The wanted zoom is calculated with a heuristics, based on the hour-glass model of the PSF to avoid wrap-around
    kxymax=PupilRadius./OrigSize ./ ImageParam.Sampling(1:2);  % pupil radius in reciprocal space units
    kz=sqrt(1/(lambda*lambda)-kxymax.*kxymax);        % calculate kz in reciprocal space units (using the lambda modified by n)
    tanTheta=kxymax/kz;
    MaxTravel = abs(ZDist) * tanTheta ./ ImageParam.Sampling(1:2) + HalfImgSize; % How far do the wave present in the in-focus slice maximally travel out to the side (in destination pixels)
    WantedZoom = MaxTravel ./ HalfImgSize .*1.3;   % This zoom should be obtained, to avoid wrap-around. The factor 1.5 was multiplied as a heiristic margin
    
    % WantedZoom(WantedZoom < 1.0)=1.0;
    WantedZoom(WantedZoom < AllowedZoom)=AllowedZoom(WantedZoom < AllowedZoom);  % if we can get the zoom "for free" (not making the currently allowed image bigger, we just zoom in as much as possible
    WantedZoom=floor(WantedZoom);  % round the Wanted Zoom down, to not clash with border pixels and to be integer
end    
% WantedZoom

if any(AllowedZoom - WantedZoom > 0)  % TargetSize ist the size in pixels with which we have to calculate now
    TargetSize=OrigSize;
else
    TargetSize=ceil(PupilRadius .* WantedZoom * 2.0);  % calculate the TargeSize such that the pupil just about fits inside
    TargetSize=TargetSize+mod(TargetSize,2);  % force to be even, due to a problem of czt for uneven sizes.
end
ZoomedImageParam=ImageParam;
ZoomedImageParam.Sampling(1:2)=ZoomedImageParam.Sampling(1:2);  % OrigSize./TargetSize
ZoomedImageParam.Size(1:2)=TargetSize;
[IncoherentRadius,PupilRadiusBig]=AbbeMaxRadiusFromPSF(PSFParam,ZoomedImageParam,0,ndes);
WantedZoom=floor(floor(TargetSize/2.0)./PupilRadiusBig);  % No matter what WantedZoom was calculated, always make it fit just inside.

ZoomToApply = WantedZoom; % floor(TargetSize/2) ./ PupilRadius;

% The speed reusability improvement is anyway marginal (2 sec out of 14)
if nargin < 4 || isempty(ZoomedPupil) || norm(size2d(ZoomedPupil)-TargetSize)~=0
    ZoomedImageParam=ImageParam;
    ZoomedImageParam.Sampling(1:2)=ZoomedImageParam.Sampling(1:2).*ZoomToApply;  % OrigSize./TargetSize
    ZoomedImageParam.Size(1:2)=TargetSize;
    InputAmp=Make2DPolPhase(ZoomedImageParam,PSFParam);  % generates the polarization dependent on the ImageParam and PSFParam. Specifically using PSFParam.polarization

    ZoomedPupil =mSimLens(InputAmp,PSFParam,ZoomedImageParam,PSFParam.Aplanatic,2,AddParams);  % jinc-based aperture for the last parameter==2
%    fprintf('recalculated pupil, zoom [%d,%d], TargetSize [%d,%d]\n',ZoomToApply,TargetSize);
    % add phase aberration using the Zernike Polynomial--------------------

    if ~iscell(AddPhase) 
        if all(AddPhase==0 | isempty(AddPhase))
            AddPhase=0;
        end
    else %all(iscell(PlusPhase) | (size(PlusPhase,1)==1 & size(PlusPhase,2)~=1) | (size(PlusPhase,1)~=1 & size(PlusPhase,2)==1))  % i.e. Zernike radial order and azimuth order
        AddPhase=ZernikePoly(AddPhase,ZoomedImageParam,PSFParam); 
    end
    % AddPhase end---------------------------------------------------------
    ZoomedPupil=ZoomedPupil.*exp(-1i*AddPhase);
end

% apply the propagator
scales=ImageParam.Sampling(1:2)  .* ZoomToApply;
sz=size2d(ZoomedPupil);
kxysqr=(abssqr(xx(sz,'freq')/scales(1))+abssqr(yy(sz,'freq')/scales(2)));
kz=2*pi*sqrt(1/(lambda*lambda)-kxysqr);
kz(1/(lambda*lambda)-kxysqr < 0)=0;

% change Dina 08.06.20 @17:53----------------------------------------------
Propagator = exp(-1i*kz*ZDist);
% shiftpos=norm([2 2 2].*ImageParam.Sampling); k0=2*pi/lambda;
% Propagator = exp(-1i*kz*ZDist).*exp(-1i*k0*shiftpos);
% end editing--------------------------------------------------------------

ZoomedPupilPropagated = ZoomedPupil .* Propagator;

% Propagated =  czt2d(ZoomedPupilProp,ZoomToApply(1),ZoomToApply(2),0,'inverse').*prod(TargetSize); 
Propagated =  czt2d(ZoomedPupilPropagated,ZoomToApply(1),ZoomToApply(2),0,'inverse'); % no scaling needed any more

Propagated=extract(Propagated,OrigSize);
% Propagated=reshape(Propagated,[size2d(Propagated) 1 size(Propagated,3)*size(Propagated,4)]);
I=squeeze(sum(abssqr(Propagated),[],4));

suma=sum(I,[],[1 2]); 
I=I/suma; % scaling is needed to avoid discontinuties in the zooming region

% single slice cannot be normalized as it would skew the 3D PSF
