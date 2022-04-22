% res = OPD(AddParams,AnglesReal,AnglesDesign) calculate the otical difference between PSF with real parameters and PSF with design condition parameters
% Arguments
% """"""""""""
% AddParams: structure having as fields
%                       e.g. AddParams=struct('ns',1.33,'ng',1.518,'ni',1.516,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5); 
%                         'ns': ri of sample, can be in a vector format if the sample is in a stratified medium.  First element of ns (if a vector) is the ri of where the emitter is at, it is the ri of the medium which the furthest from the coverslip. The last element of the vector is the ri of the medium closest to the coverslip
%                         'ng': ri of coverslip in real condition
%                         'ni': ri of immersion medium in real condition
%                         'ni0': ri of immersion medium in design condition
%                         'ts': thickness of the sample or position where the emitter is at, this should be with the same length as the vector of ns and within the same order as ns.
%                         'tg': thickness of coverslip in real condition
%                         'tg0': thickness of coverslip in design condition
%                         'wd': working distance which is supposed to be the thickness of the immersion medium in design condition
% AnglesReal: structure containing cos(theta) at different layer of the stratified medium in real condition
% AnglesDesign: structure containing cos(theta) at different layer of the stratified medium in design condition
%
% Exemple:
% """"""""""
% AddParams=struct('RIReal',[1.33 1.518 1.516],'RIDesign',[1.518 PSFParam.n],'ThicknessReal',[2e3 1.7e5],'ThicknessDesign',[1.7e5 1.5e5]);
% 
% [1] Haeberlé O. Focusing of light through a stratified medium: a practical approach for computing microscope point spread functions. Part I: Conventional microscopy. Optics communications. 2003 Feb 1;216(1-3):55-63.
% [2] Eq (3) in Gibson SF, Lanni F. Experimental test of an analytical model of aberration in an oil-immersion objective lens used in three-dimensional light microscopy. JOSA A. 1991 Oct 1;8(10):1601-13.
% Copied from the last version in folder "Files RI Version7.0"
% Updated 16.08.21 - R. Dina
% Updated 16.08.21 : add of tilt
function res = OPD(AddParams,AnglesReal,AnglesDesign)

cosDes = cat(3,AnglesDesign.cosalpha{:}); % this should be of length 2 in the 3rd dimension
cosReal = cat(3,AnglesReal.cosalpha{:}); % this should be of length N+2 in the 3rd dimension with N being the number of layers in the sample region

clear AnglesDesign AnglesReal

nicosthival = (AddParams.ni).^2*cosReal(:,:,end); % all the OPD at each layer has this term in common 

% edit today 27.06.21    
% a = rireal(end).*zrange.*cosReal(:,:,end); % we add this in defocusing optical path factor in the PropagateFieldBlock for SP, FT(half sinc) in the SincR method and McCutchen from IFTA for VS
a=0; 
res = 0;
rireal=[AddParams.ns AddParams.ng];
thickreal=[AddParams.ts AddParams.tg];
for j=1:1:length(rireal) % values corresponding to the real conditions in a stratified medium sample and coverslip in real condition 
    opdj=PathDif(rireal(j),thickreal(j),cosReal(:,:,j-1),nicosthival,AddParams.ni);
    res = res + opdj;    
end 
opdCSDes=PathDif(AddParams.ng0,AddParams.tg0,cosDes(:,:,end-1),nicosthival,AddParams.ni); % coverslip in design condition, phase equal to this is corrected by the microscope 
opdImDes=PathDif(AddParams.ni0,AddParams.wd,cosDes(:,:,end),nicosthival,AddParams.ni); % immersion medium, phase equal to this is corrected by the microscope 
res = res - a - opdCSDes - opdImDes; % it is -a to get the same orientation as in PSFGenerator but it is +a in the GL article and in my calculation
res=squeeze(res);

end

function res = PathDif(ns,ts,costs,nicosthival,ri)
    if ns == ri
        res = 0;
    else
        res = ns*ts*(costs - nicosthival./(ns.^2)); 
    end
end