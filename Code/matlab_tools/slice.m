function xout = slice(xin,dim,nb)
sz=size(xin);
in=zeros(1,length(sz));
for k=1:length(sz)
    if k==dim
        in(k)=nb;
    else
        in(k)=':';
    end
end
xout = xin(in);