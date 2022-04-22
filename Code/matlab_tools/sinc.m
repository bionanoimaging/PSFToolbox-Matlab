% out=sinc(in)
%
% sin(in)/in, sinc(0):=1

function out=sinc(in)
if isa(in,'double')
    if numel(in)~=1 || in ~=0
        out=sin(in)./in;
        out(in==0)=1.0;
    else
        out=1.0;
    end
else
    out=sin(in);
    out(in~=0)=out(in~=0)./in(in~=0);
    out(in==0)=1;
end