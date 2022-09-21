function res = ZoomPhase(AddPhase,ImageParam,PSFParam, ZoomToApply) % zoom the input phase to match the pupil size in the CZT technique
% Size: output size
% ZoomToApply: self-explanatory parameter
% AddPhase: initial phase to zoom in
% res: output

bigSize=ImageParam.Size(1:2);
% lb=min(AddPhase);
% ub=max(AddPhase);
% InitPhase=AddPhase-lb; % to avoid negative values, this is needed as we need go back and forth in real and Fourier space
InitPhase=exp(1i.*AddPhase);
updatephase=extract(phase(ft(extract(ift(InitPhase), round(bigSize.*ZoomToApply)))),bigSize);
% ImageParam.Size(1:2)=round(sz.*ZoomToApply);
% jpsf=extract(jincPSF(ImageParam,PSFParam),bigSize);
% Aperture= real(ft(jpsf)); 
% Aperture=Aperture./max(Aperture); % normalize the max
res=updatephase;%SetRange(Aperture.*updatephase,lb,ub);

end



function res = SetRange(in,lb,ub,mask)
% this function set the range of the values of the input data in to be in
% between the lower bound lb and upper bound ub
% example:
% res=SetRange(1:10,-1.5,10)
if nargin<4
    mask=newim(in)+1==1;
end
classin=class(in);
indip = castType(in,'dip_image'); % we change the type becaus some of the function below work properly and specifically for the case of dip_image
a=min(indip(mask));
b=max(indip(mask));
s=(-1*sign(lb));
res=newim(indip);
res(mask) = (indip(mask)-a)./(b-a).*(s.*abs(lb)+abs(ub))-s.*abs(lb);
res = castType(res,classin); % convert it back to the same class as the input argument in 
end