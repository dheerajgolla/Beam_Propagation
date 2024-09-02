function H = dgDHT(Rmax,M);
J0zeros = dlmread('J0zeros.dat');
selectzeros = J0zeros(1:M);
H.rgrid = Rmax*(selectzeros/J0zeros(M+1));
H.kgrid = (1/Rmax)*(selectzeros);
num = 2* besselj(0,(selectzeros*selectzeros')./J0zeros(M+1));
den = selectzeros(M)*(besselj(1,ones(M,1)*selectzeros')).^2;
H.DHT = num./den;
H.iden = H.DHT*H.DHT;
H.norm = sum((diag(H.iden)).^2)/M;
%H.DHT = H.DHT./H.norm;

