%clear all;
function thinfilmreflectivity
S = dlmread('CRYSTALS_Si_Palik.csv.txt');
lam = 1000*S(:,1);
xx = (270.7:0.1:820.6)';
yr = S(:,2);
yi= S(:,3);
yyr = spline(lam,yr,xx); %real part ref index of Si
yyi = spline(lam,yi,xx); %imag part of ref index of Si

% xx = xx(1700:4500);  %limiting wavelength 480-660nm
% yyr = yyr(1700:4500);  
% yyi = yyi(1700:4500);

x=xx./1000; %wavelength in micron
lambda = xx; %wavelength in nanometers

l=26; %layers of hbn

n0 = ones(size(x)); % refractive index of air
n1 = 2.7*ones(size(x)) - 1i*(5.446/2.7)*(lambda/1000);
%n1 =  (2.3 - 1.6i) .* ones(size(x));%ref index of graphene/graphite
%n1 =  (ng + (kg)*1i) .* ones(size(x));
%n2 = 1.8 - (0.00069375 .* (lambda - 480));%ref index of hBN with 3% dispersion
n2 = (1.8).* ones(size(xx)); %ref index of hb

n3 = 1*sqrt( 1 + 0.6961663*power(x,2)./(power(x,2)-power(0.0684043,2)) +(0.4079426*power(x,2))./(power(x,2)-power(0.1162414,2)) + (0.8974794*power(x,2))./(power(x,2)-power(9.896161,2)));  %refractive index of sio2
%n2=n3;%n4=n3;
%n3=n0;
n4 = yyr + 1i.*yyi; %ref index of silicon
%n4 = n3; %on glass or cover slip or fused Si

r(:,1) = (n0-n1)./(n0+n1); % reflection coefficient 
r(:,2) = (n1-n2)./(n1+n2);
r(:,3) = (n2-n3)./(n2+n3);
r(:,4) = (n3-n4)./(n3+n4);

d2 = 0.093; %size of hbn layer in micron
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

I = (abs(g(:,1))).^2 %full stack g+hbn+sio2


rn(:,1) = (n0 - n2)./(n0 + n2); %no graphene
rn(:,2) = r(:,3);
rn(:,3) = r(:,4);
h(:,3) = rn(:,3);
q(:,2) = p(:,3);
q(:,1) = p(:,2);


for j=1:2   
    h(:,3-j) = (rn(:,3-j) + h(:,4-j).*q(:,3-j))./(1 + rn(:,3-j).*h(:,4-j).*q(:,3-j));
end
K = abs(h(:,1)).^2; %no graphene h+s


rnnn(:,1) = r(:,1); %air graphene
rnnn(:,2) = (n1-n3)./(n3+n1); %graphene sio2
rnnn(:,3) = r(:,4);
hh(:,3) = rnnn(:,3);
qq(:,2) = p(:,3);
qq(:,1) = p(:,1);


for j=1:2   
    hh(:,3-j) = (rnnn(:,3-j) + hh(:,4-j).*qq(:,3-j))./(1 + rnnn(:,3-j).*hh(:,4-j).*qq(:,3-j));
end
M = abs(hh(:,1)).^2 % g+s


rnn(:,1) = (n0-n3)./(n0+n3); %air sio2
rnn(:,2) = r(:,4); %sio2 si
t = p(:,3);
f(:,2) = rnn(:,2);

f(:,1) = (rnn(:,1) + f(:,2).*t)./(1 + rnn(:,1).*f(:,2).*t);
L = abs(f(:,1)).^2 %just s

%contrast calculation
C1 = (K - I)./K; %h+s - g+h+s
C2 = (L - K)./L; % s - h+s
C3 = (L - I)./L; %s - g+h+s
C4 = (M - I)./M; %g+s - g+h+s
C5 = (M - I - K + L)./M; %g+s - g+h+s - h+s + s
C6 = (L - M)./M; %s - g+s
%close all
%figure
 hold  on
 plot(lambda,C6)
% %xlim([410 740]) 
 %title('graphene on Si/SiO2');
% figure
% plot(lambda,C2)
% xlim([410 740]) 
% title('hBN on Si');
% grid on
% figure
% plot(lambda,C3,'b')
% plot(lambda,C4,'black')
% plot(lambda,C6,'c')
% xlim([480 760]) 
% title('graphene on Si');
plot(lambda,M)
%hold on 
%plot(lambda,L)