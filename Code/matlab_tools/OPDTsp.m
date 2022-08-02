% function [res, Ts, Tp ] = OPDTsp(AddParams,PSFParam,ImageParam)
% Calculates the optical path difference (stored in res), transmission coefficient for s-polarized light (Ts) and for p-polarized light (Tp)
% This function takes into account a tilt under which the coverslip might be
% Parameters:
%         R0(x,y,z) : reference under which the rays in design condition are calculated
%         R1(u,v,n) : reference under which the rays in real condition are in due to the tilt in the coverslip
%         AddParams.tilt(1) : elevated angle of the coverslip from y-axis about x-axis (angle denoted by alpha)
%         AddParams.tilt(2) : azimuthal angle of the unit vector u of the axis of the coverslip from x-axis about the axial axis z (angle denoted by beta)
%         phi : azimuthal angle under which the OPD is calculated 
%         theta_j : elevation angle of the rays at medium j
% Example:
%         ImageParam=struct('Sampling',[65.5 65.5 160],'Size',[256 256 64]);
%         PSFParam=struct('NA',1.4,'n',1.518,'lambdaEm',542);
%         AddParams=struct('ns',1.33,'ng',1.518,'ni',1.515,'ng0',1.518,'ni0',1.518,'ts',1e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5,'tilt',[10 50]);
%         [res, Ts, Tp, angle] = OPDTsp(AddParams,PSFParam,ImageParam);

% Dina R - 02.10.21
% Update - 01.12.21 // tilt is included 
function [res, Tp, Ts, angle] = OPDTsp(AddParams,PSFParam,ImageParam)

if length(AddParams.ns)>1
    error('Type of system not allowed. This code can only support one layer of sample.')
end

% PARAMETERS
sz = ImageParam.Size(1:2);
scales = ImageParam.Sampling;
% zeroimg=newim([sz 3]);
myxx = xx(sz); myyy = yy(sz);
phi = atan2(myyy,myxx); sinphi = sin(phi); cosphi = cos(phi); % azimuthal angle 
tg0 = AddParams.tg0 ; wd = AddParams.wd; tg = AddParams.tg; ts = AddParams.ts ;
ni = AddParams.ni; ni0 = AddParams.ni0; ng = AddParams.ng; ng0 = AddParams.ng0; ns = AddParams.ns;

% TILT PROPERTIES
if ~isfield(AddParams,'tilt')
    AddParams.tilt=[0 0];
end
alpha = AddParams.tilt(1)*pi/180; cosa = cos(alpha); sina=sin(alpha); % elevation angle of the coverslip from y-axis about x-axis in counterclockwise direction and in degrees
beta = AddParams.tilt(2)*pi/180;   cosb = cos(beta);  sinb=sin(beta); % azimuth angle from the x-axis about z-axis in counterclockwise direction and in degrees
% M = [cosb, -cosa*sinb, sina*sinb; sinb, cosa*cosb, -cosb*sina; 0, sina, cosa]; % rotation matrix to go from reference R1 to reference R0

% REAL CONDITIONS
NA=min(PSFParam.NA/ni0,PSFParam.NA/ni);  
% NA=min([PSFParam.NA, ni, ni0, ng, ng0, ns]);
lambda=PSFParam.lambdaEm / ni;  % to obtain lambda in the medium
ktotal=2*pi/lambda; 
ktotalsqr=ktotal*ktotal;
kxysqr=(abssqr(2*pi*xx(sz,'freq')/scales(1))+abssqr(2*pi*yy(sz,'freq')/scales(2)));
kzsqr=(ktotalsqr-kxysqr);
kzsqr(kzsqr < 0)=0;
kz=sqrt(kzsqr);
costhetai=kz./ktotal; % angle in the immersion medium in real condition
sinthetai=sqrt(kxysqr)./ktotal;
sinbp=sinb.*cosphi - cosb.*sinphi ; % sin(beta-phi) % this term should be -sin(phi) if beta=0
cosgammai=cosa*costhetai + sina*sinthetai*sinbp; % = costhetai if alpha=0
singammai=real(sqrt(1-cosgammai.^2)); 

[singammas,cosgammas]=SnellLaw(ni,ns,singammai,cosgammai);
[singammag,cosgammag]=SnellLaw(ni,ng,singammai,cosgammai);

nis=ni/ns; 
nig=ni/ng;
OPDs=ns*ts*(cosgammas + nis*(-cosgammai + (1-nis)*costhetai));
OPDg=ng*tg*(cosgammag + nig*(-cosgammai + (1-nig)*costhetai));

% DESIGN CONDITIONS
[sinthetai0, costhetai0] = SnellLaw(ni, ni0, sinthetai, costhetai);
[sinthetag0, costhetag0] = SnellLaw(ni0, ng0, sinthetai0, costhetai0);
OPDg0 = ng0.*tg0.*(costhetag0 - (ni./ng0).^2.*costhetai); % path inside the coverslip
OPDi0 = ni0.*wd.*(costhetai0 - (ni./ni0).^2.*costhetai); % path inside the coverslip

% TOTAL OPD
Aperture=kxysqr/ktotalsqr < NA*NA;
res = OPDs + OPDg - OPDg0 - OPDi0;
res = res.*Aperture;

% CORRECT FOR RADIAL SHIFT 


% FRESNEL TRANSMISSION COEFFICIENTS 
if nargout>1
    listRI=[AddParams.ns AddParams.ng AddParams.ni];
    AnglesReal.cosalpha{1}=cosgammas.*Aperture; 
    AnglesReal.cosalpha{2}=cosgammag.*Aperture; 
    AnglesReal.cosalpha{3}=cosgammai.*Aperture;  
    Ts = transmCoeff(listRI,AnglesReal,'s');
    Tp = transmCoeff(listRI,AnglesReal,'p');
end

% ANGLES
angle.cos = costhetai;
angle.sin=sinthetai;

end

% functions we need
function res = R10(in,M) % rotation matrix from reference R1(u,v,n) to reference R0(x,y,z) i.e. expressing vector defined in R1 in R0
    res=newim(in);
    for k=0:2
        j=k+1;
        res(:,:,k)=M(j,1).*in(:,:,0) + M(j,2)*in(:,:,1) + M(j,3)*in(:,:,2);
    end
end

function res = R01(in,M) % rotation matrix from reference R0(x,y,z) to reference R1(u,v,n) i.e. expressing vector defined in R0 in R1
    res=newim(in);
    for k=0:2
        j=k+1;
        res(:,:,k)=M(1,j).*in(:,:,0) + M(2,j)*in(:,:,1) + M(3,j)*in(:,:,2);
    end
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
