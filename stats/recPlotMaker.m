% For use with sparse Neo
% Based off of createSubmission to validate method, generate AUC from
% training data... 

%subject to run it on (do 1 at a time for now)
clipDirNames={'Dog_2'};%{a.name}; when we automate

%number of training segment randsamples to take (without replacement)
%to find threshold (can be higher than number available, just doesnt use them all)
% nSamples=178;

%threshold to use (found after running this program and looking at ictal/interictal mean/std
thresh=[8.5];

%*****************comment switch or change ictal/interictal leaving all *_ characters in
 typeI='_ictal_segment';
%typeI='_interictal_segment';
% typeI='*_test_segment*';

%your matlab directory
matdir='/Users/idaniboy/Documents/MATLAB/';
%matdir='C:\Users\paul\Documents\MATLAB\';

%clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
clipsdir='kaggleClips/ictalSegs/';

% %automation will require data to be in separate folders per patient with _ in folder name
% a=dir(strcat(matdir,clipsdir,'*_*'));
% clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});
clipLoc=strcat(matdir,clipsdir,clipDirNames{1});

% a=dir(strcat(clipLoc,filesep,'*',typeI,'*.mat'));
% srsOfSegments={a(randsample(length(a),min(nSamples,length(a)))).name};
% nTests=length(srsOfSegments);
% total=zeros(1,nTests);

% num ictal clips
nIClips=length(dir(strcat(matdir,clipsdir,clipDirNames{1},'/*_ictal_s*.mat')));

%     progress=floor(linspace(1,nTests,20));
%     npercent=zeros(length(progress));
rMatches=zeros(1E6,4);

tic
k=1;
for i=1:nIClips
    disp(i) %progress
    clipName = strcat(clipDirNames{1},'_ictal_segment_',num2str(i));
	R=match_query2(clipLoc,clipName);
    Rs = size(R,1);
    if i == 20
        disp(i);
    end
    
    if Rs > 0
        matchRow=round((R(:,3)))+1; % had to add '/32' to fix X positions
        % matchRow = R(:,1);
        rMatches(k:k-1+Rs,:)=[matchRow,repmat(i,[Rs,1]),R(:,2),R(:,3)]; %R2 = intensity?
        k = k+Rs;
    end

    if rem(i,10) == 0    % progress update, plot, or save after every Nth chunk
        disp(i)
%       save('shazExp.mat','rMatches','rMarks','sampLengthSecs','numWindows');
    end
end


disp('Total time ran for recurrence plot:');
toc




%% Plot recurrence plot
rM = rMatches;
rM(rM(:,1)<0) = 0;
rM = rM(all(rM,2),:); % remove zeros from rM


%for neo in current version as on 6/30/2014, where it makes funky recurrence plots
% numWindows = 4973;
% spM = sparse(rM(:,1),rM(:,2),rM(:,3),numWindows,numWindows);

%for adjusted recurrence plot
%numWindows = max(ceil((rM(:,1))/29));
% numWindows = 180;
% spM = sparse(ceil((rM(:,1))/28),rM(:,2),rM(:,3),numWindows,numWindows);

% for accumulated recurrence plot
a = ceil((rM(:,1))/49);
numWindows = 180;
spM = sparse(a,rM(:,2),rM(:,2),numWindows,numWindows);


figure;spy(spM,'.',1); axis tight
set(gcf,'color','w');
title('Sparsity Recurrence Plot');


% %% Transpose
% disp(strcat('Size of rM before pruning:',num2str(numel(rM))));
% 
% % Transpose
% spM = spM + spM';
% 
% temprM = find(spM);
% size(temprM,1)
% disp(strcat('Size of rM after pruning:',num2str(numel(temprM))));
% clear temprM
% 
% figure;spy(spM,'.',1); axis tight
% set(gcf,'color','w');
% title('Sparsity Recurrence Plot');

%% Thresholded heatmap

[a,b,c] = find(spM);
trM = [a,b,c];
clearvars a b c

% Some Stats
U=mean(trM(:,3));
B=std(trM(:,3));
disp('Overall mean intensity and STD');
disp(U);
disp(B);

trM(trM(:,3)<U,:)=0; %find diags where intensity is > than U+1B (2 std above mean)
trM(all(trM==0,2),:)=[];
%rM = rM(logical(spM(:,3)),:);


% Reconsitute sparse matrix after thresholding
fM = sparse(trM(:,1),trM(:,2),trM(:,3),numWindows,numWindows);
fM = full(fM);

figure;imagesc(fM);axis tight;axis square
set(gcf,'color','w');
title('Recurrence Plot Heatmap');