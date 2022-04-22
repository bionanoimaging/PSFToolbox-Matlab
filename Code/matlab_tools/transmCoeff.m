% function [res,tirfAngle,BrewsterAngle] = transmCoeffUpdated(listRI,AnglesReal,pol)
% Compute the transmission coefficients from 3 layers (sample - coverslip - immersion media) based on the theory of Török P and Varga P
% Source: Török P, Varga P. Electromagnetic diffraction of light focused through a stratified medium. Applied optics. 1997 Apr 10;36(11):2305-12.
% Edited version based on Example 4.1 page 92 of Peatross J, Ware M. Physics of Light and Optics (2015). Brigham Young University, Department of Physics.:101-19.
%-----------------------------------------------------------
% Input Arguments:
% listRI: list of RI (array vector)
% AnglesReal: list of cosines angles (Here it is a structure and the element AnglesReal.cosalpha of the structure is what we need)
% pol: polarization 's' or 'p'
%-----------------------------------------------------------
% Outupt arguments:
% res: transmission coefficients
% tirfAngle: angle at which there is tirff for each interface
% BrewsterAngle: Brewster Angle at each interface
%-----------------------------------------------------------
% Note: 
% The phase factor which is written in Török manuscript is not necessarily included in the 
% running function but can be accessed to by uncommenting the section at the bottom of this page.
% This phase factor is taken care of in the Propagator of the PSF simulator.
%-----------------------------------------------------------
% Exemple1: Glass-Air interface
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
% AddParamsDesign=struct('RIReal',[1.5 1],'RIDesign',[1.33 1.518 1.518],'ThicknessReal',[0 1.7e5 1.5e5],'ThicknessDesign',[2e3 1.7e5 1.5e5]);
% initial_angle=dip_image(0:0.1:pi/2); % it's better to have everything in dip_image otherwise double will have problems
% sinalphai=sin(initial_angle);
% cosalphai=cos(initial_angle);
% Angles = angleGenerator(AddParamsDesign.RIReal,sinalphai,cosalphai,'forward');
% [res_upd,tirfAngle_upd,BrewsterAngle_upd] = transmCoeff(AddParamsDesign.RIReal,PSFParam,Angles,'p');
% [res_usd,tirfAngle_usd,BrewsterAngle_usd] = transmCoeff(AddParamsDesign.RIReal,PSFParam,Angles,'s'); 
% figure; plot(double(initial_angle)*180/pi,double(res_upd))
% hold on
% plot(double(initial_angle)*180/pi,double(res_usd))
% hold off
% legend('p-polarized','s-polarized')
% title('Fresnel Amplitude Coeff Glass to Air')
% xlabel('Angle of incidence')
%-----------------------------------------------------------
% Exemple2: Water-Glass-Oil interface
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
% AddParamsDesign=struct('RIReal',[1.33 1.5 1.518],'RIDesign',[1.33 1.518 1.518],'ThicknessReal',[0 1.7e5 1.5e5],'ThicknessDesign',[2e3 1.7e5 1.5e5]);
% initial_angle=dip_image(0:0.1:pi/2); % it's better to have everything in dip_image otherwise double will have problems
% sinalphai=sin(initial_angle);
% cosalphai=cos(initial_angle);
% Angles = angleGenerator(AddParamsDesign.RIReal,sinalphai,cosalphai,'forward');
% [res_p,tirfAngle_upd,BrewsterAngle_upd] = transmCoeff(AddParamsDesign.RIReal,PSFParam,Angles,'p');
% [res_s,tirfAngle_usd,BrewsterAngle_usd] = transmCoeff(AddParamsDesign.RIReal,PSFParam,Angles,'s'); 
%-----------------------------------------------------------
% Update [28.06.21 - R. Dina]
%-----------------------------------------------------------
function [res,tirfAngle,BrewsterAngle] = transmCoeff(listRI,AnglesReal,pol)

if length(listRI)~=length(AnglesReal.cosalpha)
    nlen=min(length(listRI),length(AnglesReal.cosalpha));
    listRI=listRI(1:nlen);
    AnglesReal2.cosalpha=AnglesReal.cosalpha{1:nlen}; 
    clear AnglesReal
    AnglesReal=AnglesReal2;
    fprintf('Warning: the number of mediums in the RI list is not the same as the number of mediums in the angles.\nThe lenghts are adjusted to match so the transmission coefficient is calculated with less number of interfaces.\n');
end

tirfAngle=zeros(1,length(listRI)-1); %list of incident angles for different interfaces at which there is TIRF
BrewsterAngle=zeros(1,length(listRI)-1); %list of incident angles for different interfaces at which there is only transmitted light and no reflected 

res=newim(size2d(AnglesReal.cosalpha{1}));
refl=res+1;
transm=refl;
for k = 1:length(listRI)-1 % -1 because there are length(listRI)-1 interfaces
    if listRI(k)>listRI(k+1)
        theta=asin(listRI(k+1)/listRI(k)); % incident angle at which there is TIRF
    else
        theta=0;
    end
    tirfAngle(k)=theta;
    BrewsterAngle(k)=atan(listRI(k+1)/listRI(k));

    switch pol
        case 's'
            a0 = listRI(k)*AnglesReal.cosalpha{k};
            a1 = listRI(k+1)*AnglesReal.cosalpha{k+1}; 
            denomTr = a0+a1; 
            mymask = denomTr~=0;
            t = res; 
            t(mymask) = 2*a0(mymask)./denomTr(mymask);
%             r = res;
%             numR = a0-a1;
%             r(mymask) = numR(mymask)./denomTr(mymask);

        case 'p'
            a0 = listRI(k+1)*AnglesReal.cosalpha{k};
            a1 = listRI(k)*AnglesReal.cosalpha{k+1};
            denomTr = a0+a1;
            t = res;
            numT = 2*listRI(k)*AnglesReal.cosalpha{k};
            mymask = denomTr~=0;
            t(mymask) = numT(mymask)./denomTr(mymask);
%             r = res;
%             numR = a0-a1;
%             r(mymask) = numR(mymask)./denomTr(mymask);

    end
    
    transm=transm.*t;
%     refl=refl.*r;

end

res=transm; % this is valid when the the denominator having the reflection term is close to 1 and we include the phase term in the OPD function. Otherwise, uncomment the rest of the code! 


% if length(listRI)==2
%     res = transm;
% elseif length(listRI)==3
%     % beta=AnglesReal.cosalpha{2};
%     % denominator=1+(refl.*exp(2*1i.*beta));
%     % Start editing on 18.07.20
%     thick=AddParamsDesign.ThicknessReal(2); % thickness of the medium in the midddle which is the coverslip
%     k1=2*pi*listRI(2)/PSFParam.lambdaEm; % medium in the middle which is the coverslip 
%     beta=k1.*thick.*AnglesReal.cosalpha{2};
%     denominator=1-(refl.*exp(2*1i.*beta));
%     % End editing
%     validmask = abssqr(denominator)~=0 & ~isnan(abssqr(denominator));
%     numerator=transm.*exp(1i*beta);
%     clear res;
%     res=newim(numerator,'dcomplex');
%     res(validmask)=numerator(validmask)./denominator(validmask);
% end


    