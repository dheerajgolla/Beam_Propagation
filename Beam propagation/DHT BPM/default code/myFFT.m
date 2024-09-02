
function FFT = myFFT(T, N);

FFT.tmax  = T;
FFT.Nt    = N;

% coordinates (temporal) in real space
deltat = T/N;
cr   = zeros(1,N);

size(cr)

for t=1:N
	cr(t) = deltat*(t-N/2);
end


% coordinates (omegas) in spectral space
deltaomega = 2*pi/T;
cs    = zeros(1,N);
for o=1:N/2
	cs(o) = deltaomega*(o - 1);
end
for o=N/2:N
	cs(o) = deltaomega*(o - N - 1);
end

FFT.ct = cr;
FFT.om = cs;
