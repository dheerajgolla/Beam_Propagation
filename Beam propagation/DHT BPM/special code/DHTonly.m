close all;
n=1;
Rmax = 1; %maximum spatial grid location in transverse plane in meters
c = 3*10^8; %speed of flight in m/s
lambda0 = 785 * 10^(-9); %carrier wavelength in meters
w = 2*pi*(c/lambda0); %in s^-1 units
M = 2000; %zeros of Bessel to consider - increase to increase resolution
z = 0.1; %propagation in metres
H = dgDHT(Rmax,M); % the H struct contains the DHT matrix and r grid and k grid
r=H.rgrid;
k=H.kgrid;
sigma = Rmax/100; %beam size(FWHM) in meters
E.space = initwave(r,sigma);
b = (n .* w)/c;
Kz = sqrt(b.^2 - k.^2);
prop = exp(1i.*(Kz - b).*z); % linear propagator
hank_wave = H.DHT*E.space; % DHT operation
BPM_wave = H.DHT * (prop.*hank_wave);
figure; plot(r,E.space)
figure;plot(r,BPM_wave)