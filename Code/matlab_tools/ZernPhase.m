% function TheoPhase = ZernPhase(AddParams,NA,n_order)
% compute the phase due to refractive index mismatch in terms of Zernike
% polynomials with maximum order n_order. The Zernike coefficients is first calculated using
% ZernCoeff.m and is used to generate the theoretical phase 
% Source: 
%         [1] Lakshminarayanan V, Fleck A. Zernike polynomials: a guide. Journal of Modern Optics. 2011 Apr 10;58(7):545-61.
%         [2] P Török, P Varga, and G Nemeth. Analytical solution of the diffraction integrals and interpretation of wave-front distortion when
%              light is focused through a planar interface between materials of mismatched refractive indices. JOSA A, 12(12):2660–2671, 1995.
function [TheoPhase,res] = ZernPhase(AddParams,ImageParam,PSFParam,n_order)
if nargin<4
    n_order=6; % maximum order to compute
end

NA=PSFParam.NA; % NA argument can be PSFParam instead

ni0 = AddParams.ni0;
listri = [AddParams.ns AddParams.ng AddParams.ng0 AddParams.ni] ;
listthick = [AddParams.ts AddParams.tg AddParams.tg0 AddParams.wd] ;
res = dip_image(ZernCoeff(NA, listri, listthick, ni0, n_order)); % finding the Zernike coefficients

Zke = [];
for orderlist  = 2.*(0:1:n_order) % compute the Zernike polynomials of order ranging from 0 to 2*n_order in step of 2
    Zk0 = ZernikePoly({[orderlist 0]},ImageParam,PSFParam);
    Zke = cat(3,Zke,Zk0);
end

Zke = expanddim(Zke,4);
res = expanddim(res,4); 
res = permute(res,[4 3 1 2]); 

TheoPhase = sum(Zke.*res,[],3); % multiply Zernike with the corresponding coefficients for each set of ri and thickness
TheoPhase = squeeze(sum(TheoPhase(:,:,:,0:end-2),[],4) -  sum(TheoPhase(:,:,:,end-1:end) ,[],4)); % the last two terms corresponds to design condition for coverslip and immesion medium
TheoPhase = real(TheoPhase);