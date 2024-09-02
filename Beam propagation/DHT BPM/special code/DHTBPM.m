%Beam propagation method - DHT for spatial profile propagation and FFT for
%temporal
close all;
%%%%%%%%temporal profile construction%%%%%%%%%%

tmax = 300; %total width of time domain in fs
c = 3*10^8 ; %speed of flight in m/s
N = 2048; %time resolution
dt=tmax/(N-1);
t = (-tmax/2:dt:tmax/2);
lambda0 = 785 * 10^(-9); %carrier wavelength in meters
w0 = 2*pi*(c/lambda0) * (10^-15);%  carrier frequency in fs^-1
phit = zeros(1,length(t)); %phase as a function of time
%phit = -(0.09)*(t.^1); %linear phase equals shift 
%phit = (0.01)*(t.^2); %experiment with the chirp (0.04 default)
sigsq = (3.5*10).^2; %sigma squared = width of the envelope squared in fs^2
I = exp((-t.^2)/(2*sigsq)); %intensity envelope
ewt = exp(1i*(w0*t - phit));
E=[];
E.time = (1).*sqrt(I).*((ewt)); %electric field in time domain
E.w = fftshift(fft(ifftshift(E.time))); %electric field in frequency domain through FFT
dw = 2*pi/tmax;
wmax = 2*pi/dt;
w = -wmax/2:dw:wmax/2; % frequency array in fs^-1 units
%n=1.0*ones(1,length(w)); % refractive index
x = ((2*pi*c)./(w*10^15)).*(10^6); %Wavelength in microns
n=sqrt(1+1.03961212./(1-0.00600069867./(x.^2))+0.231792344./(1-0.0200179144./(x.^2))+1.01046945./(1-103.560653./(x.^2)));
I_w = (abs(E.w)).^2; %intensity of electric field in frequency domain
A = (I_w/max(I_w))>0.0001; %phase blanking to avoid random phase values
I_w = A.*I_w;
phiw = A.*unwrap(angle(E.w)); %phase blanking to avoid random phase values
figure; plot(t,I);title('temporal profile before propagation')

%%%%%%%%spatial domain construction (transverse coordinates)%%%%%%%%%%

sigma = 100.0e-06; %spot size at z=0 in meters
z = pi.*(sigma.^2)/lambda0; %propagation length in meters (equal to Rayleigh range in this example)
steps = 10;
dz = z/steps;
Rmax = 5.0e-04; %maximum spatial grid location in transverse plane in meters
M = 200; %spatial resolution
T = (n.*z/c)*10^15; %group delay in fs
H = dgDHT(Rmax,M); % the H struct contains the DHT matrix and r grid and k grid
r=H.rgrid; %spatial cordinates
k=H.kgrid; %transverse cordinates
E.space = exp((-r.^2)./sigma.^2); %Electric field as a function of transverse coordinates
figure;
plot(r,abs(E.space).^2);
title('spatial profile before propagation')
z_w = zeros(M,N);
b = zeros(1,N);prop = zeros(M,N);Kz = zeros(M,N);hank_wave = zeros(M,N);step = zeros(M,N);
E_0 = (E.space)*E.time; % initial Electric field - 2d matrix - transverse coordinate and time
E.z_w = E.space*E.w; % initial Electric field - 2d matrix - transverse coordinate and frequency

for pp = 1:steps
b = zeros(1,N);prop = zeros(M,N);Kz = zeros(M,N);hank_wave = zeros(M,N);
for m=1:N
if (abs(w(m))<5 && abs(w(m))>0.5) %only propagate these frequencies - optional
b(m) = ((n(m) .* (w(m)))/c)*10^15;
Kz(:,m) = sqrt(b(m).^2 - k.^2); % k is not a function of m
prop(:,m) = exp(1i*(Kz(:,m) - b(m))*dz); % linear propagator
hank_wave(:,m) = H.DHT*E.z_w(:,m); % DHT operation
E.z_w(:,m) = H.DHT * (prop(:,m).*hank_wave(:,m)); %inverse DHT operarion after propagation
end
end
end
E_t_z = ifftshift(ifft(fftshift(E.z_w),N,2));
figure;
plot(r,abs(E_t_z(:,N/2)).^2);title('spatial profile after propagation')
figure;
plot(t,abs(E_t_z(1,:)).^2);title('temporal profile after propagation')
figure; 
surf(t,r,abs(E_t_z).^2);title('after propagation')
shading flat
view(2)
figure;
surfc(t,r,abs(E_0).^2);title('before propagation')
shading flat
view(2)