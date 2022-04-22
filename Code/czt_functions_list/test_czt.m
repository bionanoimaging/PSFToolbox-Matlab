% Comparison of the existing czt function 
s = 15:35;
res = [];
sz = [35 50 4];
for k = 1:length(s)
    xin=extract(readim,[s(k) 50]);
    tic; xout1 = czt_1d(xin, 1 ,1); toc % copy of the code from Julia // 
    tic; xout2 = czt_1d2(xin, 1 ,1); toc % Dina own code copied from Rabiner's manuscript directly
     tic; xoutGross = czt_dip(xin, 1,1,2,'forward'); toc  % Prof Gross code
    fref = ft1d(xin); % reference
    resk = cat(3,fref,xout1,xout2,xoutGross);
    res = cat(4, res, extract(resk,sz));
end