%%%%%%%%%pulse prop in time domain%%%%%%%%%
close all
tmax = 800; %in fs
N = 2000;
dt=tmax/(N-1);
t = (-tmax/2:dt:tmax/2);
c = 3*10^8; %speed of flight in m/s
lambda0 = 785 * 10^(-9); %carrier wavelength in meters
w0 = 2*pi*(c/lambda0) * 10^-15; %in fs^-1 units
phit = zeros(1,length(t)); %phase as a function of time
%phit = -(0.09)*(t.^1); %linear phase equals shift 
%phit = (0.04)*(t.^2); %experiment with the chirp (0.04 default)
sigsq = 4*10; %sigma squared = width of the envelope squared in fs^2
I = exp((-t.^2)/(2*sigsq)); %intensity envelope
ewt = exp(1i*(w0*t - phit));
E_real = (1).*sqrt(I).*(real(ewt));
E_comp = sqrt(I).*exp(-1i*phit);
 figure
 hold on 
% plot(t,E_comp)
% plot(t,I,'r')
 plot(t,I,'k')
% hold off
%%
dw = 2*pi/tmax;
wmax = 2*pi/dt;
w = -wmax/2:dw:wmax/2;
E_w = fftshift(fft(ifftshift(E_comp)));
S_w = (abs(E_w)).^2;
A = (S_w/max(S_w))>0.0001; %phase blanking to avoid random phase values
S_w = A.*S_w;
phiw = A.*unwrap(angle(E_w)); %phase blanking to avoid random phase values
%figure
%plotyy(w,S_w,w,phiw)
%%
E_w2 = sqrt(S_w) .* exp(-1i*phiw);
E_t2 = ifftshift(ifft(fftshift(E_w)));
%figure
plot(t,(abs(E_t2)).^2,'r')
%%
x = ((2*pi*c)./(w*10^15)).*(10^6); %Wavelength in microns
n=sqrt(1+1.03961212./(1-0.00600069867./(x.^2))+0.231792344./(1-0.0200179144./(x.^2))+1.01046945./(1-103.560653./(x.^2)));
% figure
% plot(x,n)

