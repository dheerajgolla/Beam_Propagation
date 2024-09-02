function ic = IC(t,r)

la = 800.0e-9;
wr = 100.0e-06;
wt = 30.0e-15;
o0 =  2*pi*3.0e8/la;

ic = exp(-(r/wr)^2  - (t/wt)^2 - 1i*o0*t);

end
