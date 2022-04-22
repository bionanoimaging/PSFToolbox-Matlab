% function res = ZernCoeff(NA, listri, listthick, ni, n_order)
% function generates the Zernike coefficients for a given refractive index and thickness in the corresponding medium of interest
%     listri : a vector of the refractive index at each layer
%     listthick : a vector of the thicknesses of each layer
%     ni : refractive index of the immersion medium in real condition
%     n_order : maximum order of the Zernike to compute
%
% Example:
%     listri = [1.33 1.52 1.518 1.518] ; % sample (real condition), coverslip in real and design condition respectively, immersion medium (design condition)
%     listthick = [2e3 1.7e5 1.7e5 1.5e5] ; % thicknesses of sample, coverslips real and design condition, working distance
%     ni = 1.515; % ri in real condition
%     n_order = 4; % Zernikes up to order 4*2 in step of 2
%     res = ZernCoeff(NA, listri, listthick, ni, n_order)
% 
% Source: P Török, P Varga, and G Nemeth. Analytical solution of the diffraction integrals and interpretation of wave-front distortion when
%              light is focused through a planar interface between materials of mismatched refractive indices. JOSA A, 12(12):2660–2671, 1995.

function res = ZernCoeff(NA, listri, listthick, ni, n_order)

% if length(listri) ~= length(listthick)
%     error("The length of RI vector and the thickness are not the same.")
% end
% res = zeros(length(listri)+1,n_order + 1); % +1 because order 0, 1, ..., n_order

k = reorient(0:1:n_order,4); % order of the Zernike to which we need coefficients is 2*k // along the column of the output matrix
sz = size(listri);
if length(sz)<3
    listri=reorient(listri,3); % else the data must be 3D and there is no need to reorient it
end

listri(listri<NA)=NA; % function is only valid and real if ri>=NA

sz = size(listthick);
if length(sz)<3
    listthick=reorient(listthick,3); % else the data must be 3D and there is no need to reorient it
end

nratio = ni./listri;
sinalpha2 = (NA/ni)^2;
c1 = (2.*listri.^2)./(NA^2) - 1; % argument of the first Kfun
c2 = 2/sinalpha2 - 1; % argument of the second Kfun

K1 = Kfun(k, c1);
K2 = nratio.*Kfun(k, c2);

cfac = listthick*NA*(2.*k+1)./4; % constant factor
ksqrt2=ones(size(k)); ksqrt2(k==0)=sqrt(2);
cfac = cfac.*ksqrt2; % order zero. See Eq. 35 for the definition of epsilon_nm 

res = cfac.*(K1 - K2);

res(abs(res)<1e-10) = 0;

res=squeeze(res);
% res = real(res);
end

function kres = Kfun(k,c) % Eq. 45
    c1 = sqrt(c.^2 - 1);
    c1 = c + c1;
    kk = 2.*k;
    fac1 = (kk - 1).*(kk + 1).*(c1.^(k - 1/2));
    fac2 = (kk - 1)./( (kk + 3) .* (c1.^2));
    kres = - (sqrt(2)./fac1).*(1 - fac2);
end