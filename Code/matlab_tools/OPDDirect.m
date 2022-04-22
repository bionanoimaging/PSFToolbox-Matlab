% function OPD = OPDDirect(AddParams,ImageParam,PSFParam) 
% OPD direct calculation
function [res, tp, ts] = OPDDirect(AddParams,ImageParam,PSFParam) 

sz = ImageParam.Size(1:2);
scales = ImageParam.Sampling(1:2);
lambda=PSFParam.lambdaEm / AddParams.ni;  
ktotal=2*pi/lambda;
ktotalsqr=ktotal*ktotal;
kxysqr=(abssqr(2*pi*xx(sz,'freq')/scales(1))+abssqr(2*pi*yy(sz,'freq')/scales(2)));
kzsqr=(ktotalsqr-kxysqr);
kzsqr(kzsqr < 0)=0;
kz=sqrt(kzsqr);

cosalpha=kz./ktotal; % angle alpha is in the immersion medium in real condition
sinalpha=sqrt(kxysqr)./ktotal;

listRI=[AddParams.ns AddParams.ng AddParams.ni]; % length should be equal to : length_stratified_sample + 1 (coverslip) + 1 (immersion medium)
AnglesReal=angleGenerator(listRI,sinalpha,cosalpha,'backward'); % the third dimension of this should be the same as the length of the listRI
ktotDes=2*pi/(PSFParam.lambdaEm/AddParams.ni0); % design condition
kzsqrDes=(ktotDes*ktotDes-kxysqr);
kzsqrDes(kzsqrDes < 0)=0;
kzDes=sqrt(kzsqrDes);
sinImDes=sqrt(kxysqr)./ktotDes; % sin in immersion medium in design condition
cosImDes=kzDes./ktotDes; % cos in immersion medium in design condition

AnglesDesign=angleGenerator([AddParams.ng0 AddParams.ni0],sinImDes,cosImDes,'backward');

res = OPD(AddParams,AnglesReal,AnglesDesign); 

[IncoherentRadius,PupilRadius,EwaltRadii]=AbbeMaxRadiusFromPSF(PSFParam,ImageParam);
rho=rrscale(ImageParam.Size(1:2), 1.0 ./ PupilRadius);
mask=rho<=1.0 ;

if nargout>1
    listRI=[AddParams.ns AddParams.ng AddParams.ni];
    ts = transmCoeff(listRI,AnglesReal,'s');
    tp = transmCoeff(listRI,AnglesReal,'p');
end

res = res.*mask;