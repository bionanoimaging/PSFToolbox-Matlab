function [sinalpha2,cosalpha2]=SnellLaw(n1,n2,sinalpha1,cosalpha1)

if n1==n2 
    sinalpha2=sinalpha1;
    cosalpha2=cosalpha1;
else
    sinalpha2=(n1/n2)*sinalpha1; 
    sinalpha2(sinalpha2>1)=1;%nan;
    cosalpha2=cos(asin(sinalpha2)); % no angle periodicity problem here because angles are between ]0,pi/2[    
end
cosalpha2(abs(cosalpha2)<1e-6)=0;
sinalpha2(abs(sinalpha2)<1e-6)=0;