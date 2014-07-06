%WAVETEST Example Matlab script for WAVELET, using NINO3 SST dataset
%
% See "http://paos.colorado.edu/research/wavelets/"
% Written January 1998 by C. Torrence
%
% Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
%   changed all "log" to "log2", changed logarithmic axis on GWS to
%   a normal axis.
% Modified Dec 2013, by Dan Bernardo to incorporate 1/f noise significance
%   testing for application to neural ECoG data
% 2013.12.09 - replaced Red Noise w/Pink Noise Markov model for 1/f noise
% 
% 
% load 'sst_nino3.dat'   % input SST time series
% sst = sst_nino3;

EEG = pop_loadset('/Users/idaniboy/Documents/MATLAB/shazamEEG/fpData/fingerprintData1/JY_samples/3sec_677_680_UF_JY.set');

sst = EEG.data;

%------------------------------------------------------ Computation

% normalize by standard deviation (not necessary, but makes it easier
% to compare with plot on Interactive Wavelet page, at
% "http://paos.colorado.edu/research/wavelets/plot/"
variance = std(sst)^2;
sst = (sst - mean(sst))/sqrt(variance) ;

n = length(sst);
dt = 1 ;
time = [0:length(sst)-1]*dt ;  % construct time array
xlim = [0,2400];  % plotting range
pad = 1;      % pad the time series with zeroes (recommended)
% dj = 0.25;    % this will do 4 sub-octaves per octave 0.25 -->
% s0 = 3*dt;    % this says start at a scale of 6 months 2 -->
% j1 = 16/dj;    % this says do 7 powers-of-two with dj sub-octaves each

dj = 0.08;    % this will do 4 sub-octaves per octave 0.25 -->
s0 = 3*dt;    % this says start at a scale of 6 months 2 -->
j1 = 63;    % this says do 7 powers-of-two with dj sub-octaves each

lag1 = 0.0001;   % lag-1 autocorrelation for 1/f noise, varying 0.38-->rednoise
lag2 = 0.4;   % still varying...
lag3 = 0.8;   % still varying 
mother = 'Morlet';

% Wavelet transform:
tic
[wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
toc
power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Significance levels: (variance=1 for the normalized SST)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:
global_ws = variance*(sum(power')/n);   % time-average over all times
dof = n - scale;  % the -scale corrects for padding at edges
global_signif = wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);

% Scale-average between El Nino periods of 2--8 years
avg = find((scale >= 2) & (scale < 8));
Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scale')*(ones(1,n));  % expand scale --> (J+1)x(N) array
scale_avg = power ./ scale_avg;   % [Eqn(24)]
scale_avg = variance*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
scaleavg_signif = wave_signif(variance,dt,scale,2,lag1,-1,[2,7.9],mother);

whos

%------------------------------------------------------ Plotting

%--- Plot time series
subplot('position',[0.1 0.75 0.65 0.2])
plot(time,sst)
set(gca,'XLim',xlim(:))
xlabel('Time')
ylabel('Voltage')
title('a) Raw ECoG Signal')
hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.37 0.65 0.28])
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%contour(time,log2(period),log2(power),log2(levels));  %*** or use 'contourfill'
imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
xlabel('Time')
ylabel('Scale')
title('b) CWT w/95% significance contours')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(time,log2(period),sig95,[-99,1],'k');
hold on
% cone-of-influence, anything "below" is dubious
plot(time,log2(coi),'k')
hold off

%--- Plot global wavelet spectrum
subplot('position',[0.77 0.37 0.2 0.28])
plot(global_ws,log2(period))
hold on
plot(global_signif,log2(period),'--')
hold off
xlabel('Power (degC^2)')
title('c) Global Wavelet Spectrum')
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel','')
set(gca,'XLim',[0,1.25*max(global_ws)])

%--- Plot 2--8 yr scale-average time series
subplot('position',[0.1 0.07 0.65 0.2])
plot(time,scale_avg)
set(gca,'XLim',xlim(:))
xlabel('Time (year)')
ylabel('Avg variance (degC^2)')
title('d) 2-8 yr Scale-average Time Series')
hold on
plot(xlim,scaleavg_signif+[0,0],'--')
hold off


% Plot 1/f model
% figure
% syms u;ezplot(1/u);hold on
% plot(fft_theor,'--r');
% axis auto
%figure
%contour(time,log2(period),sig95,[-99,1],'k');
% end of code


%% OLD USING MATLAB'S CWT
figure
    SR = 800; % sampling rate
    dt =1/SR; % dt = recording period
    s0  = 3*dt; % s0 = smallest scale; default = 6
    ds = 0.08;  % ds = spacing between scales; default = 0.15
    NbSc = 64;  % nb = number of scales
    SCA = {s0,ds,NbSc}; 
tic
    cwtsig = cwtft({sst,dt},'scales',SCA,'wavelet','morl');
    toc
    cfs=abs(flipud(cwtsig.cfs));
    MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
    scales = cwtsig.scales.*MorletFourierFactor;
    freq = 1./scales;

    S = cfs;
    imagesc(flipud(S));hold on
    
    % Determine power spectrum
    power = (S).^2 ;        % compute wavelet power spectrum
    
% Determine significance levels    

tic 
[signif,fft_theor] = wave_signif(1.0,dt,scales,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant
toc
    contour(flipud(sig95),[-99,1],'k');
figure
    contour(time,log2(freq),flipud(sig95),[-99,1],'k');

figure
imagesc(sig95);
    
