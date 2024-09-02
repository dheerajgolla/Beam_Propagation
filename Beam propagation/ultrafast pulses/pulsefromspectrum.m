% contruct pulse from spectrum and phase
close all
w0 = 2.4; %in fs^-1
w_max = 50;%in fs^-1

N = 1001;
dw = w_max/(N-1);

w = -w_max/2:dw:w_max/2;
sigfsq = 0.02; % sigma squared of frequency squared = width in freq (use 0.3 for chirp of 10)
Sw = (exp((-(w-w0).^2)./(2*sigfsq)));% + (exp((-(w+w0).^2)/(2*sigfsq)));
%phiw = 0.*w;
phiw = 1.*(0.*((w - w0).^2)); % use 10 as coeff of w^2 for chirp

E_w = sqrt(Sw) .* (exp(-1i.*phiw)); %+ exp(1i.*phiw));
E_t = sqrt(N).*fftshift(ifft(ifftshift(E_w)));
I_t = abs(E_t).^2;

E_w_w = fftshift(fft(ifftshift(E_t)))./sqrt(N);
S_w_w = abs(E_w_w).^2;

tmax = 2*pi/dw;
dt = 2*pi/w_max;
t = -tmax/2:dt:tmax/2;

energyspectrum = sum(Sw).*dw;
energytime = sum(I_t).*dt;
ratio = energytime/energyspectrum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plotyy(w,Sw,w,phiw)
% xlim([0 5])
% xlim([0 5])
title('Spectrum and phase')
xlabel('frequency')
ylabel('Intensity')

figure
plot(t,real(E_t))
xlim([-100 100])
title('Real field')
xlabel('time')
ylabel('Field strength')

figure
plot(t,abs((E_t).^2))
title('Intensity')
xlim([-80 80])
xlabel('time')
ylabel('Intensity')

figure
plot(w,S_w_w)
title('Intensity recalculated')
xlabel('frequency')
ylabel('Intensity')