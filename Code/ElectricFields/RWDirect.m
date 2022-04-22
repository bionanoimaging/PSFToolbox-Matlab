function [h, amp3d, FPlane] = RWDirect(ImageParam,PSFParam,smoothAperture)

if nargin<3
    smoothAperture=1;
end

if ~isfield(PSFParam,'PosMidZ')
    PosMidZ=0;
else
    PosMidZ=PSFParam.PosMidZ;
end

sz=ImageParam.Size(1:2);
scales=ImageParam.Sampling;
lambda=PSFParam.lambdaEm./PSFParam.n;
NA=PSFParam.NA./PSFParam.n;

ktotal=2*pi/lambda;
ktotalsqr=ktotal*ktotal;
% kxysqr=(abssqr(xx(sz,'freq')/scales(1))+abssqr(yy(sz,'freq')/scales(2)));
% kzsqr=(ktotalsqr-kxysqr);
% kzsqr(kzsqr < 0)=0;
% kz=sqrt(kzsqr);

kxysqr=(abssqr(xx(sz,'freq')/scales(1))+abssqr(yy(sz,'freq')/scales(2)));
% kz=2*pi*sqrt(1/(lambda*lambda)-kxysqr);
kz=(1/(lambda*lambda)-kxysqr);
kz(kz < 0)=0;
kz=2*pi*sqrt(kz);
% Question: this was of calculating the kz seems ok for the psf but why calculating kz this way and not the above kz (in comment format)


cosalpha=kz./ktotal; % angle alpha is in the immersion medium in real condition
sinalpha=sqrt(kxysqr)./ktotal;

phi = atan2(ImageParam.Sampling(2)*yy(sz),ImageParam.Sampling(1)*xx(sz));
cosphi = cos(phi);
sinphi = sin(phi);

HardAperture= kxysqr/ktotalsqr < NA*NA;
if (smoothAperture)
        jpsf=jincPSF(ImageParam,PSFParam); % n defines in PSFParam corresponds to the ri of the medium through which light propagates in reality and ndes defines the ri of the medium through which the medium is supposed to propagate as described by the manufacturer
    if smoothAperture>1  % means: Do Damping
        jpsf=DampEdge(jpsf,0.05,2,1,3);
    end
    Aperture= real(ft(jpsf));  % force it to be real
else
    Aperture= HardAperture;
end

FPlane=newim([sz 1 3],'scomplex');

FPlane(:,:,:,0) = cosalpha + sinphi.^2.*(1-cosalpha);
FPlane(:,:,:,1) = (cosalpha-1).*cosphi.*sinphi;
FPlane(:,:,:,2) = sinalpha.*cosphi;

AMask=cosalpha>1e-6;
Aperture= dip_image(Aperture,'sfloat') .* AMask; 
Aperture(AMask) = Aperture(AMask) .* (cosalpha(AMask)).^(-1/2);  % this factor is a combination of AF=sqrt(cosalpha); % aplanatic factor and sz=cosalpha in Eq. 2.2 and 2.18a


Propagator=exp(-1i*kz*(ImageParam.Sampling(3)*zz([1 1 ImageParam.Size(3)])+PosMidZ)); % constructs the propagator Eq. 2.2

FPlane = FPlane.*Aperture;

constantfactor=-1i/lambda; % see Eq. 2.2

amp3d = constantfactor .* ift2d(FPlane.*Propagator);



h=abssqr(amp3d(:,:,:,0));
for fc=1:size(amp3d,4)-1
    h=h+abssqr(amp3d(:,:,:,fc)); % other vector field components (Y,Z) if existing
end


