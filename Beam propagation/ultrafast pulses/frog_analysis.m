function [  ] = frog_analysis( directory )
%frog_analysis Summary of this function goes here
%   Provide the directory of the FROGv3.exe outputs and this function will
%   save all the information to variables.

p1 = fopen(strcat(directory,'/arecon.dat')); % reconstructed spectrogram
p2 = fopen(strcat(directory,'/a.dat')); % experimental spectrogram

%%% retrieve data from arecon.dat

n1 = cell2mat(textscan(p1,'%f',2));
minmax1 = cell2mat(textscan(p1,'%f',2));
waves1 = cell2mat(textscan(p1,'%f',n1(1)));
fre1 = 2.*pi.*3e8./waves1;
delays1 = cell2mat(textscan(p1,'%f',n1(2)));
data1 = reshape(cell2mat(textscan(p1,'%f',n1(1)*n1(2))),[n1(1) n1(2)]);

%%%retrieve original frog trace

n2 = cell2mat(textscan(p2,'%f',2));
minmax2 = cell2mat(textscan(p2,'%f',2));
waves2 = cell2mat(textscan(p2,'%f',n2(1)));
fre2 = 2.*pi.*3e8./waves2;
delays2 = cell2mat(textscan(p2,'%f',n2(2)));
data2 = reshape(cell2mat(textscan(p2,'%f',n2(1)*n2(2))),[n2(1) n2(2)]);

%%% retrieve reconstructed electric field and spectrum

stufft = load(strcat(directory,'/Ek.dat'));
t = stufft(:,1);
intt = stufft(:,2);
pt = stufft(:,3);

stuffw = load(strcat(directory,'/speck.dat'));
lam = stuffw(:,1);
w = 2.*pi.*197.3./lam; % in eV
intw = stuffw(:,2);
pw = stuffw(:,3);

%%% plot all the things

figure;
surf(delays1,waves1,data1');
view(2);shading flat;
title('Reconstructed Spectrogram');
xlabel('Delay (fs)');
ylabel('Wavelength (nm)');
xlim([min(delays1) max(delays1)])
ylim([335 500])

figure;
surf(delays2,waves2,data2');
view(2);shading flat;
title('Experimental Spectrogram');
xlabel('Delay (fs)');
ylabel('Wavelength (nm)');
xlim([min(delays2) max(delays2)])
ylim([335 500])

figure;
plot(t,intt.*max(pt)./max(intt),t,pt);
title('Electric field intensity in Time')
xlabel('Time (fs)')
ylabel('Intensity (arb. units), phase (rad.s)')

figure
plot(lam,intw.*max(pw)./max(intw),lam,pw);
title('Electric field in Frequency');
xlabel('Wavelength (nm)')
ylabel('Intensity (arb. units), phase(rad.s)')
xlim([600 1000])

%%% recreate the full field

time = min(t):0.01:max(t);
field = zeros(length(time));

for i = 1:length(time)
field(i) = sum(intw .* exp(-1i .* ((2 * pi * 3e2 ./ lam) * time(i) + pw)));
end

figure
plot(time,real(field),time,imag(field),time,abs(field))
title('Field reconstructed from spectrum')
xlabel('Time, (fs)')
ylabel('Field (arb.)')

end