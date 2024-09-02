tmax = 1000; %in fs
N = 10000;
totreps=10;
dt=tmax/(N-1);
t = (-tmax/2:dt:tmax/2);
c = 3*10^8; %speed of flight in m/s
lambda0 = 785 * 10^(-9); %carrier wavelength in meters
w0 = 2*pi*(c/lambda0) * 10^(-15); %in fs^-1 units
phit = zeros(1,length(t)); %phase as a function of time
%phit = -(0.09)*(t.^1); %linear phase equals shift 
%phit = (0.004)*(t.^2); %experiment with the chirp (0.04 default)
sigsq = 30; %sigma squared = width of the envelope squared in fs^2
I = exp((-t.^2)/(2*sigsq)); %intensity envelope
ewt = exp(1i*(w0*t - phit));
E_real = (1).*sqrt(I).*(ewt);E_train=[];
for rep=1:totreps
    E_train=[E_train E_real];
end
time =linspace((-totreps*tmax/2),(totreps*tmax/2),totreps*N);
figure;
plot(time,E_train)
%%%%%%%%%%%%%%%%%%%%
dw = 2*pi/(totreps*tmax);
wmax = 2*pi/dt;
%w = -wmax/2:dw:wmax/2;
w = linspace(-wmax/2,wmax/2,totreps*N);
E_w = fftshift(fft(ifftshift(E_train)))./(sqrt(totreps*N)); %FFT of E_real
S_w = (abs(E_w)).^2;
A = (S_w/max(S_w))>0.00001; %phase blanking to avoid random phase values
%S_w = A.*S_w;
phiw = A.*unwrap(angle(E_w)); %phase blanking to avoid random phase values
figure
plot(w,S_w)%,w,phiw)
xlim([2 3])
title('the frequency picture')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_w_recons = sqrt(S_w) .* exp(-1i*phiw);
E_t = (sqrt(totreps*N)).*ifftshift(ifft(fftshift(E_w))); %IFFT of E_w
I_t = (abs(E_t)).^2;
B = (I_t/max(I_t))>0.00001;
phi_t = B.*unwrap(angle(E_t));

figure
hold on
plot(time,abs(E_train).^2,'k')
plot(time,I_t,'r')
hold off