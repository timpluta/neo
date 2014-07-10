function [R,L] = match_query2(clipLoc,clipName)
% [R,L] = match_query(D,SR)
% FOR FINDING THEN MATCHING FRESH LANDMARKS
%     Match landmarks from an audio query against the database.
%     Rows of R are potential maxes, in format
%      songID  modalDTcount modalDT
%     i.e. there were <modalDTcount> occurrences of hashes 
%     that occurred in the query and reference with a difference of 
%     <modalDT> frames.
%     L returns the actual landmarks that this implies.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

%tim 2014-06-27 clipLoc is the directory of the segment file with no
%filesep. clipName is the exact file name with no filesep.


% clear global Lmarks
% global Lmarks
%Rt = get_hash_hits(landmark2hash(find_landmarks(D,SR)));
%Lq = find_landmarks(D,SR,handles);

% designate that we are querying for 1 sec clips
isClip = 1; 

% strcat(clipLoc,filesep,clipName,'_test_segment_',num2str(nthSeizure),'.mat'))
% [data,freq,latency,channels]=sParLoad(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_concatenated_',num2str(nthSeizure),'.mat'));
[data,freq,latency,channels]=sParLoad(strcat(clipLoc,filesep,clipName));
D=data;
SR=freq;
%Lq = fuzzify_landmarks(Lq);

% Augment with landmarks calculated half-a-window advanced too
% landmarks_hopt = (handles.EEG.times(end)-handles.EEG.times(1))/totalFFTtime(1,2)/1000;
% landmarks_hopt=

%max number of matches to find per channel/and for recurrence plot
perChan=100;
rperChan=100;

%how many dt blocks to meld close matches
meldWindow=1;

%allow all landmarks in 2*lmBox+1 side length box
lmBox=0;

%input hash constraints (only used for landmark box, not needed otherwise)
ldT=0;
udT=2;
lF=11;
uF=62;

%check block around initial 2nd landmark point, except for fingerprint (already done)
% if fp==1 && size(fpLm,1)~=0
%     Lq=fpLm;
% else


% choose findLandmarks3 or findLandmarks2
% in future automate this off of autorandom sampling of # of high freq...
% landmarks found
if length(fieldnames(channels)) == 30
%     %non clip, segment search    
%     Lq= findLandmarks3(data,freq,channels);
    %clip search
    a=strfind(clipName,'_');
    b=strfind(clipName,'.');
    clipNum=clipName(a(4)+1:b(1)-1);
    Lq= findLandmarks3(data,freq,channels,str2double(clipNum));

else
%     %non clip, segment search
%     Lq= find_landmarks2(data,freq,channels);
    %clip search
    a=strfind(clipName,'_');
    b=strfind(clipName,'.');
    clipNum=clipName(a(4)+1:b(1)-1);
    Lq= find_landmarks2(data,freq,channels,str2double(clipNum));
end

if lmBox>0
    preLq=Lq;
    %preallocation
    Lq=zeros(((2*lmBox+1)^2)*length(Lq),6);
    %row counter
    i=1;
for h=1:size(preLq,1)
    for timeBox=-lmBox:lmBox
        for freqBox=-lmBox:lmBox
            Lq(i,1)=preLq(h,1);
            if preLq(h,2)+freqBox>=lF && preLq(h,2)+freqBox<=uF
                Lq(i,2)=preLq(h,2)+freqBox;
            else
                Lq(i,2)=preLq(h,2);
            end
            if preLq(h,3)+freqBox>=lF && preLq(h,3)+freqBox<=uF
                Lq(i,3)=preLq(h,3)+freqBox;
            else
                Lq(i,3)=preLq(h,3);
            end
            if preLq(h,4)+timeBox>=ldT && preLq(h,4)+timeBox<=udT
                Lq(i,4)=preLq(h,4)+timeBox;
            else
                Lq(i,4)=preLq(h,4);
            end
            Lq(i,5)=preLq(h,5);
            Lq(i,6)=preLq(h,6);
            i=i+1;
        end
    end
end
end
    %     preLq=[preLq;find_landmarks(D(i,round(landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(landmarks_hopt/2*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(3*landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(landmarks_hopt*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(5*landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(6*landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(7*landmarks_hopt/4*SR):end),SR,handles)];
% end
%Lq=unique(Lq,'rows');
%Lq = [Lq;find_landmarks(D(round(landmarks_hopt/4*SR):end),SR)];
%Lq = [Lq;find_landmarks(D(:,round(landmarks_hopt/2*SR):end),SR,handles)];
%Lq = [Lq;find_landmarks(D(round(3*landmarks_hopt/4*SR):end),SR)];
% add in quarter-hop offsets too for even better recall

if size(Lq,1)<=1
%     Rt=[];
    R=[0 0 0];
else
    Hq = landmark2hash(Lq);
%     Rt = get_hash_hits2(Hq,clipLoc);
    R=get_hash_hits2(Hq,clipLoc);
    if isempty(R)
        R=[0 0 0];
    end
end



% nr = size(Rt,1);
% 
% if nr > 0
% tkR=Rt;
% offsets=[];
% j=1;
% % %comment out next line to return no zeros
% % R=zeros(perChan,3);
% while j<=perChan
%     if j>1
%         rows_to_remove = any((tkR==dts(xx)), 2);
%         tkR(rows_to_remove,:) = [];
%     end
%     [dts,xx] = unique(sort(tkR(:,2)),'first');
%     dtcounts = 1+diff([xx',size(tkR,1)]);
%     [vv,xx] = max(dtcounts);
%     R(j,:) = [1,vv,dts(xx)];
%     doubles=0;
%     for k=1:j-1
%         if abs(R(j-doubles,3)-R(k,3))<meldWindow
%             offsets=[offsets; R(k,3) R(j-doubles,3)];
%             R(k,2)=R(k,2)+R(j-doubles,2);
%             R(j-doubles,:)=[];
%             doubles=doubles+1;
%         end
%     end
%     j=j+1-doubles;
% 
%     if length(tkR)<1
%         j=perChan+1;
%     end
% end
% 
%   % Sort by descending match count
%   [vv,xx] = sort(R(:,2),'descend');
%   R = R(xx,:);
%   %added by tim to fix unknown t=0 hash problem
% %   R(R(:,3)<41,:) = [];
% 
% %L = [];
%   % Extract the actual landmarks for plotting match...
%   % if needed, then uncomment L = []; 
% %   for m=1:perChan
% %       if size(R,1)>=m
% %           H = Rt(Rt(:,2)==R(m,3),:);
% %           try
% %             extras=offsets(offsets(:,1)==R(m,3),2);
% %           catch me
% %               extras=[];
% %           end
% %           for o=1:length(extras)
% %               H=[H;Rt(Rt(:,2)==extras(o),:)];
% %           end
% %           % Restore the original times
% %           for i = 1:size(H,1)
% %             hix = find(Hq(:,3)==H(i,3));
% %             hix = hix(1);  % if more than one...
% %             H(i,2) = H(i,2)+Hq(hix,2);
% %             L(i,:) = hash2landmark(H(i,:));
% %             if R(m,3)==0
% %                 R(m,3)=1;
% %             end
% %             Lmarks(i,:,m)=[L(i,:) R(m,3) H(i,1)];
% %            end
% 
% 
%       % Return no more than 10 hits, and only down to half the #hits in
%       % most popular
% %       if size(R,1) > 10
% %         R = R(1:10,:);
% %       end
% %           maxhits = R(m,2);
% %           nuffhits = R(:,2)>(maxhits/2);
% %       end
%       %R = R(nuffhits,:);
% %   end
%   
%   %added as inaccurate quick fix of lots of 1 landmarks. figure out t=0/1
%   %landmarks to speed up
% %   if R(1,3)==1
% %       R=R(2:end,:);
% %   end
  L=1;
% % else
%   R = [0 0 0];
%   L = [];
%   disp(['*** NO HITS FOUND in ' clipName]);
end
% end
