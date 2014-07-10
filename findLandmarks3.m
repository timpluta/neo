function [L,S,T,maxes] = findLandmarks3(D,freq,ID,clipNum)
    % Includes:
    %   - weightedCentroid based landmarking
    %   - peakfinding based on fastPeakFinder
    %
    % 2013.1028
    % added two fields to L, for illustration purposes in illustrate match
    % ChA and ChB; A is originating Channel
    % 2013.1115
    %   DAN --> added multiChannel capability. 
    %           collapases all landmarks found across all channels into a single matrix. 
    % 2013.1208
    %   TIM --> implemented spatio-spectro-temporal hash
    % 2013.1210
    %   DAN --> cleaned up code
    %       --> added 95% confidence interval for spectra
    % 2014.0623 --> added back top band
    %           --> removed code for Circos near end (where it sorts L)
    %           --> Note that COI is for 2400 (assumes freq 800 and length 3 sec)
    %           --> Note that smallest scale can change with freq
    %           --> Note that Sig levels could be moved outside find_landmarks2
    %           -------to increase speed
    %           --> Note that there is still some Circos code that saves temporal
    %           matrix lTA! this is slowing stuff down a ton

    % SETTINGS 
    if nargin < 3
        ID=0;
        clipNum=1;
    end
    
    if nargin < 4
        error('this version of find_landmarks requires clip number')
    end

    % global totalFFTtime

    % %tim commented out info
    % tInfo = strcat('Start find_landmarks2','@',datestr(now));
    % disp(tInfo);

    % Set resolution of spectra:
    rezBoxes = 128; % How many boxes in spectra per sec. changed from 32-->128

    % Prevent landmark peaks within this many blocks (more for temporally dense single Ch
    % data)
    barrier=0;

    % Limit the number of pairs that we'll accept from each peak
    maxpairsperpeak=40;   % moved to front by DAn

    % Define target area (time-frequency)
    targetdf = 31;  % +/- 50 bins in freq (LIMITED TO -32..31 IN LANDMARK2HASH)

    % Time to look ahead
    targetdt = 3;  % (LIMITED TO <64 IN landmark2hash by current bits)  64 --> 2

    switch freq
        case 5000
            %downsample if freq is 5000
            dfreq = 5;
            Ds = zeros(size(D,1),freq/5);
            for j=1:size(D,1)
                Ds(j,:) = decimate(D(j,:),dfreq);
            end
            D=Ds;
            freq = freq/dfreq;

            %Adaptive threshold for centroid finder
            thresh = 1.2E5;

            % add zero buffer for total signal length of 2048
            zB = zeros(size(D,1),512);
            D = [zB,D,zB];

            % mask (1024 long)
            mask = ones(64,1024);
            a = flipud(triu((ones(50,50))));
            b = kron(a,[1 1 1]);
            br = fliplr(b);
            mask(1:3,:) = 0;
            mask(:,1:32) = 0;
            mask(:,992:end) = 0;
            mask(11:60,33:182) = b;
            mask(11:60,end-150-32+1:end-32) = br;  
            mask = logical(~mask);

        case 400
        % signal is == to 400 HZ or 500 Hz
            %Adaptive threshold for centroid finder
            thresh = 2E5; %best value ~ 2 for dog_1

            % add zero buffer for total signal length of 1024
            zB = zeros(size(D,1),312);
            D = [zB,D,zB];

            % mask
            mask = ones(64,512);
            a = flipud(triu((ones(40,40))));
            b = kron(a,[1 1 1]);
            br = fliplr(b);
            mask(1:3,:) = 0;
            mask(:,1:64) = 0;
            mask(:,448:end) = 0;
            mask(11:50,65:184) = b;
            mask(11:50,393-64:end-64) = br;
            mask = logical(~mask);

        case 500 % freq == 500 / this is patient 1
            thresh = 1E6;

            % add zero buffer for total signal length of 1024        
            zB = zeros(size(D,1),256);
            D = [zB,D,zB];  

            % mask
            mask = ones(64,512);
            a = flipud(triu((ones(40,40))));
            b = kron(a,[1 1 1]);
            br = fliplr(b);
            mask(1:3,:) = 0;
            mask(:,1:26) = 0;  
            mask(:,486:end) = 0; 
            mask(11:50,27:146) = b;
            mask(11:50,end-120-26+1:end-26) = br;
            mask = logical(~mask);
    end

    % Other settings:
    D=single(D);
    verbose = 0;
    nmaxes3 = 0;
    maxes3 = [];
    lenMask = size(mask,2); %this is goal length of signal 
    
    % Morlet parameters
%     dt =1/freq; % dt = recording period

    dt = (1/lenMask); % dt = recording period (neo_vers accounts for zero pads)
    s0  = 3*dt; % s0 = smallest scale; default = 6
    ds = 0.08;  % ds = spacing between scales; default = 0.15
    NbSc = 64;  % nb = number of scales
    SCA = {s0,ds,NbSc};

    % 95% sig params
    lag1 = 0.0001;
%    lengthD = size(D,2);

    % disp(strcat('Original size of file: ',num2str(length(D))));
        
    % disp('Channel:');
    %# Loop through channels and add maxes to nmaxes2
    for h=1:length(fieldnames(ID))
        % disp(h);
        
        %tim added 2014-06-26, use faster morf for 1sec segments
        
        [S,scales] = morFinger1sec(D(h,:),dt,SCA);
        %oS = S;%for debugging
        
%         %discriminates to only use for 1sec clips
%         if length(D)<=freq+20
%             [S,scales] = morFinger1sec(D(h,:),dt,SCA);
%         else
%             [S,scales] = morFinger(D(h,:),dt,ID,SCA);
%         end
%        figure;imagesc(S)
        
        % remove zero padding
        if size(mask,2)==1024
            S = S(:,513:1536); % for freq 1000 (2048 --> 1024 samps)
        else
            S = S(:,257:768); % 1024-->512samps 
        end
        
        % use mask to reflect COI
        S(mask) = 0;
               
        % lengthD = size(D,2); %length of sample
        % lengthSig = 1; Set to 1sec for neo, since always dealing w/1sec
        % clips
        
        % Determine power spectrum for significance levels
        power = (S).^2 ;        % compute wavelet power spectrum

        % Determine significance levels     
        [signif,fft_theor] = wave_signif(1.0,dt,scales,0,lag1,-1,-1,'morlet');
        sig95 = (signif')*(ones(1,lenMask));  % expand signif --> (J+1)x(N) array
        S = power ./ sig95;         % where ratio > 1, power is significant
        
        %Helps survival of smaller power oscillations through downsampling
        S = S.*10;
        
        %% DOWNSAMPLE MATRIX
        % Strategy For extending temporal landmark pairs efficiently  
        % newlen should be relative to size of the data sample!... so fingerprint and source will be different
        % make multiples of 100
        % for changing resolution 

        newmatrix=neointerp(S,1,rezBoxes); %to 1 b/c always set to 1 sec clips
        S=round(newmatrix);

        % 2d FFT Filter
        % For improved temporal localization
        %Determine good padding for Fourier transform
        PQ = 2*size(S);%replaced paddedsize, equivalent to 1 arg like itwas
        %Create a Gaussian Highpass filter 5% the width of the Fourier transform
        D0 = 0.10*PQ(1);
        H = hpfilter('gaussian', PQ(1), PQ(2), D0);
        % Calculate the discrete Fourier transform of the image
        F=fft2(double(S),size(H,1),size(H,2));
        % Apply the highpass filter to the Fourier spectrum of the image
        HPFS_S = H.*F;
        % convert the result to the spacial domain.
        HPF_S=real(ifft2(HPFS_S));
        % Crop the image to undo padding
        HPF_S=HPF_S(1:size(S,1), 1:size(S,2));
        % Store HPF_S
        S = HPF_S;
        % REMOVE values less than a threshold (<0)
        zeroS = S<=0;
        S(zeroS) = 0;
        % END 2d FFT FILTER CODE
        
        % BAND-PASS (already done by mask)
        %S(1:10,:) = 0; % REMOVE LOW FREQS < 10 Hz
        %S(54:64,:) = 0;  % REMOVE high FREQS > 100 Hz
 
        % Linear intensity filter (in future make 1/f)
        intensityMod = 1:64;
        S=diag(intensityMod')^1.7*S;

        %Centroid finder    
%         filt = (fspecial('gaussian', [3 3],0.5)); %note filter is turned off in centroidFinder
%         p=FastPeakFind(S',thresh,filt,0,2);  
        %dan added next 2014-06-27 with updated FPF for more sensitivity
        filt = (fspecial('gaussian', 5,1));
        [p,d]=FastPeakFind(S',thresh,filt,2,1); 
        p=floor(p);
        
        nmaxes2=size(p,1)/2;
        maxes2=zeros(2,nmaxes2);
        maxes2=reshape(p,2,nmaxes2)';  

        %Swap columns; took out when removed ' from S' in line 256
        maxes2 = fliplr(maxes2);
        
%         disp(maxes2);pause

% %        debugging landmark finding algo
%         figure;
%         subplot(2,2,1);imagesc(oS);hold on; scatter(maxes2(:,1)+rezBoxes+4,maxes2(:,2));
%         subplot(2,2,2);imagesc(S);hold on; scatter(maxes2(:,1),maxes2(:,2));
%         subplot(2,2,3);imagesc(d');hold on; scatter(maxes2(:,1),maxes2(:,2));drawnow
%         %subplot(2,2,4);plot(Ds(h,:));
%         pause
%         close
%         close
        
        %add column with channel info
        maxes2 = horzcat(maxes2,repmat(h,size(maxes2,1),1))';  
        
        %# concatenate maxes2 into maxes3
        maxes3 = horzcat(maxes3,maxes2);
        nmaxes3 = nmaxes3 + nmaxes2;

    end %end of multichannel for loop

% disp('Done creating matrix with peaks. Making peak pairs...');

%% Pack the maxes into nearby pairs = landmarks
  
% Limit the number of pairs that we'll accept from each peak
% maxpairsperpeak=3;   % moved to front by DAn

% Landmark is <starttime F1 endtime F2>
L = zeros(nmaxes3*maxpairsperpeak,7);

% Store maxes3 into maxes2;
maxes2 = maxes3;
clear maxes3

% Reset nlmarks marker
nlmarks = 0;

for i =1:nmaxes3
  startt = maxes2(1,i);     
  F1 = maxes2(2,i);  
  maxt = startt + targetdt;
  minf = F1 - targetdf;
  maxf = F1 + targetdf;
  %finds only maxes outside of barrier
  %matchmaxs =
  %find((maxes2(1,:)>startt)&(maxes2(1,:)<maxt)&(maxes2(2,:)>minf)&(maxes2(2,:)<maxf)&(abs(maxes2(2,:)-F1)>barrier|(maxes2(1,:)>startt+barrier)));
  % old code finds landmarks greater than dt of 1. 
  
  %tim modified 2014-06-27. conditions are (startt < peakTime < maxt)
  %%& (minf < peakFreq < maxf)) without including maxes2(:,i)
  tempmaxes2=maxes2;
  tempmaxes2(1,i)=-100;
  matchmaxs = find((tempmaxes2(1,:)>=startt)&(tempmaxes2(1,:)<maxt)&(tempmaxes2(2,:)>minf)&(tempmaxes2(2,:)<maxf));

  if length(matchmaxs) > maxpairsperpeak
    % limit the number of pairs we make; take first ones, as they
    % will be closest in time
    matchmaxs = matchmaxs(1:maxpairsperpeak);
  end
  
  
  
  for match = matchmaxs
      %OPTIONAL if maxes2(1,match)-startt < 10  %store only peaks within 10 ms 
      %this is also decided by targetDt
      
        nlmarks = nlmarks+1;
        L(nlmarks,1) = startt;
        L(nlmarks,2) = F1; 
        L(nlmarks,3) = maxes2(2,match);  % frequency row
        L(nlmarks,4) = maxes2(1,match)-startt;  % time column difference 
        L(nlmarks,5) = maxes2(3,i);  % originating Channel number
        L(nlmarks,6) = maxes2(3,match);  % ending Channel number  
        L(nlmarks,7) = clipNum; %added 2014-07-10 clip number for eliminating self matches

      %end
  end
  
  
end

% L = L(1:nlmarks,:);  % old code... added sort, in order to gen timed rec plot

% sort L, in order to make temporal matrix (FOR CIRCOS! and NEMObinary)
% L=sortrows(L(1:nlmarks,:),1);

% disp('Done creating matrix with peak pairs.');

if verbose
  disp(['find_landmarks: ',num2str(length(D)/targetfreq),' secs, ',...
      num2str(size(S,2)),' cols, ', ...
      num2str(nmaxes),' maxes, ', ...
      num2str(nmaxes3),' bwd-pruned maxes, ', ...
      num2str(nlmarks),' lmarks']);
end
  
%tInfo = strcat('finished find_landmarks2','@',datestr(now));
%disp(tInfo);

% for debug return, return the pruned set of maxes
maxes = maxes2;

%tim added 2014-06-26 to fix all 0 hashes but dont know why they happen
L( ~any(L,2), : ) = [];