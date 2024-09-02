% 2-D radially symmetric pulsed beam propagator, FFT+DHT based method

% radial domain
LR     = 5.0e-04;
NR     = 100;

% temporal domain
LT     = 300.0e-15;
NT     = 1024;
omegamin = 0.5e+15;
omegamax = 5.0e+15;

% propagation along z
LZ     = pi*(100.0e-06)^2/800.0e-9; % set to Rayleigh range
stps   = 10;
dz     = LZ/stps;

% auxiliary parameters
cc = 3.0e8;
k0 = 2.0*pi/800.0e-9;

% prepare Fourier transform
FT = myFFT(LT,NT);

% prepare Hankel transform
HT = myDHT(LR,NR);

% amplitude holder, define an initial condition
am0 = zeros(NT,NR);
for t=1:NT
for x=1:NR
	am0(t,x) = IC(FT.ct(t), HT.cr(x));
end
end

am = am0;

% propagator
pr = zeros(NT,NR);
for o=1:NT
for k=1:NR

    % only propagate active frequencies
	if ((abs(FT.om(o)) > omegamin)&&(abs(FT.om(o)) < omegamax) )

	% paraxial propagator: no frequency dependence	
	pr(o,k) = exp(-1i*dz*HT.kt(k)^2/(2.0*k0));

	% paraxial propagator: with frequency dependence	
	pr(o,k) = exp(-1i*dz*HT.kt(k)^2/(2.0*FT.om(o)/cc));

	% non-paraxial propagator: with frequency dependence and moving frame	
	pr(o,k) = exp(+1i*dz*( sqrt(FT.om(o)^2/cc^2 - HT.kt(k)^2) - 1.00005*abs(FT.om(o)/cc) ));
    end


end
end


% evolve amplitude
for s=1:stps
	fprintf(1,'executing %d out of %d, distance = %d\n',s,stps,s*dz);

    % FFT
    for r=1:NR
        am(:,r) = fft(am(:,r));
    end
%    figure; plot(FT.om,am(:,r))
    % DHT
    for t=1:NT
        am(t,:) = HT.T*am(t,:)';
    end
    
    % propagator
    am = pr.*am;
    
    % inverse DHT
    for t=1:NT
        am(t,:) = HT.T*am(t,:)';
    end
    
    % inverse FFT
    for r=1:NR
        am(:,r) = ifft(am(:,r));
    end
    
end
		 
figure(1);
plot(10^15*FT.ct,real(am(:,1)));
hold on;
plot(10^15*FT.ct,abs(am( :,1)).^2,'b');
plot(10^15*FT.ct,abs(am0(:,1)).^2,'r');

figure(2);
imagesc(abs(am).^2); colorbar;

figure(3);
imagesc(abs(am0).^2);
figure(4);
plot(HT.cr,abs(am( 1,:)).^2,'b');title('spatial profile after')
