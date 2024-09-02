% contruct pulse from spectrum and phase
clear all
close all
%%%%%%%%
%[I, K, lambda] = thinfilmreflectivity; % stack reflectivity;

w0 = 2.4; %in fs^-1
w_max = 12;
c = 3*10^8;
N = 1001;
dw = w_max/(N);
w = -w_max/2:dw:w_max/2; %frequency in fs^-1
lamb = ((2*pi*c)./(w*10^15)).*(10^6);
sigfsq = 0.02; % sigma squared of frequency squared = width in freq (use 0.3 for chirp)
Sw = (exp((-(w-w0).^2)./(2*sigfsq)));% + (exp((-(w+w0).^2)/(2*sigfsq)));
phiw = 0.*w; 
%phiw = 1.*(10.*((w - w0).^2)); % use 10 for chirp
A = (Sw/max(Sw))>0.0005; %phase blanking to avoid random phase values
E_w = A.*sqrt(Sw) .* exp(-1i*phiw);
E_t = fftshift(fft(ifftshift(E_w)))./sqrt(N);
I_t = abs(E_t).^2;

tmax = 2*pi/dw;
dt = 2*pi/w_max;
t = -tmax/2:dt:tmax/2;
%%%%%%%%%%%%%%%%%%%%%%%pulse shape%%%%%%%%%
figure
plotyy(w,Sw,w,phiw)

figure;plot(lamb(A),Sw(A));
% xlim([0 5])
% xlim([0 5])
title('Spectrum and phase')
xlabel('frequency')
ylabel('Intensity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Si medium information%%%%%%%%
S = dlmread('CRYSTALS_Si_Palik.csv.txt');
lam = S(:,1); %wavelength in nm
%w_r = ((2*pi*c)./(lam.*10^-6)) *10^-15;% frequency in fs^-1 units
%x =  ((2*pi*c)./(w*10^15)).*(10^6); %Wavelength in microns
%xx = (270.7:0.1:820.6)';

ww = w(A);
Sww = Sw(A);
%phiww = phiw(A);

lambda2 = ((2*pi*c)./(ww*10^15)).*(10^6); %Wavelength in microns
B = lambda2 < 0.8341;
lambda = lambda2(B);
%ww = ww(B);%((2*pi*c)./(lambda.*10^15)) .*(10^6);
yr = S(:,2);
yi= S(:,3);
yyr = interp1(lam,yr,lambda,'spline','extrap'); %real part ref index of Si
yyi = interp1(lam,yi,lambda,'spline','extrap'); %imag part of ref index of Si
x=lambda;
% xx = xx(1700:4500);  %limiting wavelength 480-660nm
% yyr = yyr(1700:4500);  
% yyi = yyi(1700:4500);
figure
plotyy(lambda,yyr,lambda,yyi)
figure
plotyy(lam,yr,lam,yi)

Sww = Sww(B);%phiww=phiw(B);

n0 = ones(size(x)); % refractive index of air
n1 = 2.7*ones(size(x)) - 1i*(5.446/2.7)*(lambda/1000);
n2 = (1.8).* ones(size(x)); %ref index of hb

n3 = 1*sqrt( 1 + 0.6961663*power(x,2)./(power(x,2)-power(0.0684043,2)) +(0.4079426*power(x,2))./(power(x,2)-power(0.1162414,2)) + (0.8974794*power(x,2))./(power(x,2)-power(9.896161,2)));  %refractive index of sio2
n4 = yyr + 1i.*yyi; %ref index of silicon

r(:,1) = (n0-n1)./(n0+n1); % reflection coefficient 
r(:,2) = (n1-n2)./(n1+n2);
r(:,3) = (n2-n3)./(n2+n3);
r(:,4) = (n3-n4)./(n3+n4);

d2 = 0.010; %size of hbn layer in micron
%l*0.0004; 
d1 = 0.00034; %size of graphene monolayer micron
d3 = .282; %size of sio2 layer in microns

p(:,1) = exp(-2i*(2*pi*n1*d1)./x);  %graphene phase
p(:,2) = exp(-2i*(2*pi*n2*d2)./x);  %hbn phase
p(:,3) = exp(-2i*(2*pi*n3*d3)./x);  %Sio2 phase

g(:,4) = r(:,4);

for j= 1:3
    g(:,4-j) = (r(:,4-j) + g(:,5-j).* p(:,4-j))./ (1 + r(:,4-j).*g(:,5-j).*p(:,4-j));
end

I = (abs(g(:,1))).^2; %full stack g+hbn+sio2

Sww = I'.*Sww;
figure;
plot(ww(B),Sww)


%%%%%%%%%%%%%%%%%%%%%%%%%%% SiO2 propagation %%%%%%%%%%%%%%%%%%%
%xx = 300:0.1:900;
% x =  ((2*pi*c)./(w*10^15)).*(10^6); %Wavelength in microns
% xx = A*.1000*x; %wavelength in nm
% A = 0.3 < x & x < 0.9;
% w = ((2*pi*c)./(xx.*10^-9)) *10^-15;% frequency in fs^-1 units
% % propagation through glass
% nSiO2 = 1*sqrt(1.0 + 0.81996*power(x,2)./(power(x,2)-power(0.10596,2)) - 0.01082*power(x,2));
% n = A.*nSiO2;
% n= 1.4*ones(1,length(w));
% figure
% plot(w,n,'r')
% k_w = n .* (w/c) .* 10^15; %wavevector in m^-1
% z = 0.01;%propagation length in metres
% 
% E_w_z = exp(1i.*k_w.*z).*(E_w);
% E_t_z = fftshift(fft(ifftshift(E_w_z)));
% I_t_z = abs(E_t_z).^2;
% 
% figure
% plot(t,I_t)
% title('Intensity')
% xlabel('time')
% xlim([-100 100])
% ylabel('Intensity')
% figure
% plot(t,I_t_z)
% title('Intensity after glass')
% xlabel('time')
% xlim([-100 100])
% ylabel('Intensity')

%%