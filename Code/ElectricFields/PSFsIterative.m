% function [h,amp3d]=PSFsIterative(ImageParam,PSFParam,PupilPhase,Options)
% % Propagate the field in an iterative way instead of a PropagateFieldBlock
% Options: an optional flag. 
%         Options = 0: compute psfs with the given window but do not dim the window
%         Options = 1: dim the window with a pre-calculated factor
%         Options = 2: expand the window to accomodate the wrap-around and dim accordingly
%
%
% Example:
% ImageParam=struct('Sampling',[100 100 100],'Size',[256 256 16]);
% PSFParam=struct('NA',0.6,'n',1.33,'lambdaEm',500);
% AddParams=struct('ns',1.518,'ng',1.518,'ni',1.518,'ng0',1.518,'ni0',1.518,'ts',0,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
% orderlist={[3,1];[1 -1 1 60]}; % order of the Zernike phase: coma and tilt aberration
% [TiltComaPhase,rho,phi] = ZernikePoly(orderlist,ImageParam,PSFParam);
% [h,amp3d]=PSFsIterative(ImageParam,PSFParam,AddParams,TiltComaPhase);

function [h,amp3d]=PSFsIterative(ImageParam,PSFParam,AddParams,AddPhase,Options)
if nargin<5
    Options=1;
end

if nargin<3 || isempty(AddParams)  % i.e. if the structure is not given
    AddParams=[]; % for our next reference
    ndes=PSFParam.n; % ri in the immersion medium in a design or ideal condition
    nreal=ndes; % ri in the immersion medium in real condition 
else
    ndes=AddParams.ni0;
    nreal=AddParams.ni;
    PSFParam.n=nreal; % the field is propagating in the immersion medium whose refractive index is given in real condition
end

if nargin<4 || isempty(AddPhase)
    AddPhase=0;
end

if ~isfield(PSFParam,'Aplanatic')
    PSFParam.Aplanatic = 1;  % 1 means emission the intensity is scaled by cos(theta) and the electric field is scaled by sqrt(cos(theta)), 0: no aplanatic factor, updated today 29.06.21
end
if ~isfield(PSFParam, 'PosMidZ')
    PSFParam.PosMidZ = 0.0;
end
% mSimLens(BFPAperture,lambda,scales,NA,aplanar,smoothAperture,AddParams)

InputPol=Make2DPolPhase(ImageParam,PSFParam);  % generates the polarization dependent on the ImageParam and PSFParam. Specifically using PSFParam.polarization
FPlane=mSimLens(InputPol,PSFParam,ImageParam,PSFParam.Aplanatic,1,AddParams,AddPhase);  % Circular input polarisation, Emission aplanatic factor, jinc-based-Aperture

sz=ImageParam.Size;
amp3d=newim([sz 3],'scomplex');

zfoc=(sz(3)-mod(sz(3),2))/2; % focal position // -mod(x,2) is added to account for even and odd length in z

% calculate parameters to calculate kz then prop_pupil
scales=ImageParam.Sampling;
NA=min(PSFParam.NA/ndes, PSFParam.NA/nreal); % the aperture is defined by the minimum aperture
lambda=PSFParam.lambdaEm / PSFParam.n;
% propagate field
prop_pupil=propagator(lambda,scales,sz);
big_prop_pupil=prop_pupil; % default

if Options>=1 % to avoid wrap-around, define a region to damp. Use of Eq. 14 in the manuscript
    costhetamax=cos(asin(NA));
    z=squeeze(zz([1 1 ImageParam.Size(3)])).*ImageParam.Sampling(3);
    kxmax=2*NA/lambda;
    kzmax=(1-costhetamax)/lambda;
    a=kxmax/kzmax;
    W = floor(abs(4.*(z.*a + (1.3.*lambda/(2*NA)) )./ImageParam.Sampling(1))); % Eq. (14) in 1st manuscript
    if Options == 1
        Wmsk=W>sz(1); % requires window bigger than the given window to avoid wrap-around
        DampFac=newim(Wmsk);
        DampFac(Wmsk)=(W(Wmsk)-sz(1))/sz(1)/2; % percentage of window to dim
        DampFac(~Wmsk)=min(DampFac(Wmsk)); % set the DampEdge factor, where there should not be wrap-around dimmed, to be the smallest factor as per the previous line
        DampFac(DampFac>=0.15)=max(DampFac(DampFac<0.15)); % there might be wrap-around over the whole window if the given window is too small but only consider to dim 15% at the edge otherwise the signal itself will be removed
    end % --
end

DampFac(:) = 0.05; 

jpsf=jincPSF(ImageParam,PSFParam,ndes);
Aperture=real(ft(jpsf)); Aperture=Aperture./max(Aperture);

pupil_orig=squeeze(FPlane.*Aperture);
pupil=pupil_orig;
szslice=[sz(1:2) 3];
for j=zfoc:-1:0 % focal  position is included here
    slice=ift2d(pupil);
    amp3d(:,:,j,:)=extract(slice,szslice);
    if j==0
        break % end the loop as the following lines won't be needed anymore
    end
   
    if Options==1 % window size is kept the same but we apply DampEdge to dim the edges
        df=DampFac(j); % region to damp // window hanning
        pupil = ft2d(DampEdge(slice,df,2,0,1 )); % dim down to 1 0.05% of the border region in the rectangular 2D grid
    elseif Options==2 % zero-pad window size to account for the wrap-around
        sznew=szslice; sznew(1)=max([sz(1),W(j)]); sznew(2)=sznew(1); % it is assumed that window is a square
        if sznew(1)~=sz(1)
            slice=extract(slice,sznew);
            scalesnew=[scales(1:2).*(sznew(1:2)./sz(1:2)) scales(3)];
            big_prop_pupil=extract(propagator(lambda,scalesnew,sz),sznew(1:2));
        end
        pupil=ft2d(slice);
        
    end
    pupil=pupil.*big_prop_pupil;
%     pupil=pupil.* prop_pupil; 
end

% HALF PART OF THE PSF

if AddPhase==0 && isempty(AddParams)% there is no additional phase that can break the symmetry about the focal position
    for j=zfoc+1:1:sz(3)-1%-mod(sz(3)+1,2)
        jpos=sz(3)-j;
        amp3d(:,:,j,:)=amp3d(:,:,jpos,:);
    end
else % there is additional phase which can break the symmetry about the focal position
    pupil=pupil_orig;
    big_prop_pupil=conj(prop_pupil);
%     sznew=size(pupil);
    for j=zfoc+1:1:sz(3)-1%-mod(sz(3)+1,2)
%         big_prop_pupil=extract(prop_pupil,sznew(1:2));
        pupil=pupil.* big_prop_pupil;
        slice=ift2d(pupil);
        amp3d(:,:,j,:)=extract(slice,szslice);
        if j==sz(3)-1
            break % end the loop as the following two lines won't be needed anymore
        end
        if Options==1
            df=DampFac(j); % region to damp // window hanning
            pupil = ft2d(DampEdge(slice,df,2,0,1 )); % dim down to 1 0.05% of the border region in the rectangular 2D grid
        elseif Options==2
% %             szslice=size(slice);
%             sznew=szslice; sznew(1)=max([sz(1),W(j)]); sznew(2)=sznew(1); % it is assumed that window is a square
%             slice=extract(slice,sznew);
%             pupil=ft2d(slice);
%         end
        %**
            sznew=szslice; sznew(1)=max([sz(1),W(j)]); sznew(2)=sznew(1); % it is assumed that window is a square
            if sznew(1)~=sz(1)
                slice=extract(slice,sznew);
                scalesnew=[scales(1:2).*(sznew(1:2)./sz(1:2)) scales(3)];
                big_prop_pupil=extract(conj(propagator(lambda,scalesnew,sz)),sznew(1:2));
            end
            pupil=ft2d(slice);
        end
%**
    end
end
 
h=squeeze(abssqr(amp3d(:,:,:,0)) + abssqr(amp3d(:,:,:,1)) + abssqr(amp3d(:,:,:,2)) );
h=h/sum(h);
end

function prop_pupil=propagator(lambda,scales,sz)
    kxysqr=(abssqr(xx(sz(1:2),'freq')/scales(1))+abssqr(yy(sz(1:2),'freq')/scales(2)));
    kz=2*pi*sqrt(1/(lambda*lambda)-kxysqr);
    kz(1/(lambda*lambda)-kxysqr < 0)=0;
    prop_pupil=exp(1i.*kz.*scales(3));
end
