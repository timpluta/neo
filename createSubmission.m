function createSubmission(  )

%checks both ictal and interictal segments using match_query on the
%ictal-only hash table. total number of hits within top perchan (set in match_query)
%matches are averaged. the mean and std of this total are used to find a 
%threshold. test segments above this threshold are considered ictal for the
%submission.


tic

%number of training segment randsamples to take (without replacement)
%to find threshold (can be higher than number available, just doesnt use them all)
nSamples=200;

% %manual thresholds to use (if threshold is run below, these will be replaced
% %before running submission)
% searchScore=[0];

%your matlab directory
    if ~ispc
        matdir='/Users/idaniboy/Documents/MATLAB/';
        clipsdir='kaggleClips/ictalSegs/';
    else
        matdir='C:\Users\paul\Documents\MATLAB\';
        clipsdir='kaggleShazam\clips\';
    end

% %automation requires data to be in separate folders per patient with _ in folder name
% %automated
a=dir(strcat(matdir,clipsdir,'*_*'));
clipDirNames={a.name};
% clipDirNames={'Patient_3','Patient_6'};
% clipDirNames={'Dog_2'};%,'Dog_2','Dog_3'};
nClips=length(clipDirNames);
total=zeros(nClips,nSamples,2);



%find searchScore

for ii=1:nClips
    
    %without this line, hash tables from other clips might get used, comment out
    %for faster 1-subject testing only
    clear get_hash_hits
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});    
    a=dir(strcat(clipLoc,filesep,'*_ictal_segment*.mat'));
    srsOfSegments={a(randsample(length(a),min(nSamples,length(a)))).name};
    nTests=length(srsOfSegments);
    
    %ghetto progress 
    disp(['finding searchScore for ' clipLoc])
    
    %roundabout way to allow for mean/std in parfor
    parfor i=1:nTests
        
        R=match_query2(clipLoc,srsOfSegments{i});
        if ~isempty(R)
            total(ii,i,1)=length(R);
        end
%             total(i)=mean(R(:,2));
        
    end
    
%     sSI=total;
%     sSII=total;
%     ictalMean=mean(total)
%     ictalStd=std(total)
    
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});    
    a=dir(strcat(clipLoc,filesep,'*_interictal_segment*.mat'));
    srsOfSegments={a(randsample(length(a),min(nSamples,length(a)))).name};
    nTests=length(srsOfSegments);
   
    parfor i=1:nTests

        R=match_query2(clipLoc,srsOfSegments{i});
        if ~isempty(R)
            total(ii,i,2)=length(R);
        end
%         total(i)=mean(R(:,2));
        
    end
%     sSII=total;
%     sSI=total;
%     interMean=mean(total)
%     interStd=std(total)
    
    
    %ghetto progress 
    disp(['finding searchScore for ' clipLoc])

    sSI=total(ii,total(ii,:,1)>0,1);
    sSII=total(ii,total(ii,:,2)>0,2);
    
    % Create true labels vector
    labels = [ones(1,size(sSI,2)),zeros(1,size(sSII,2))];

    % Tabulate scores
    
    scores = [sSI,sSII];
    
    % Define positive class
    posclass = 1;
    
    % Generate AUC
    [~,~,~,~,OPTROCPT] = perfcurve(labels,scores,posclass);
    optTPR = OPTROCPT(1,2); %optimal TPR
    
    [TPR,~,T,~] = perfcurve(labels,scores,posclass,'Xcrit','TPR');
    ind = find(TPR==optTPR);
    searchScore(ii) = mean(T(ind))
    
end

toc

save('mrscore.mat','searchScore','total')


%make submission
tic
resultsCell=cell(1,nClips);
% searchScore=load('mrscore.mat', 'searchScore');
for ii=1:nClips

    %without this line, hash tables from other clips might get used, comment out
    %for faster 1-file testing only
    clear get_hash_hits
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    seizure=zeros(nTests,1);
    clip = cell(nTests,1);
    
     %ghetto progress
    disp(['finding submission values' clipLoc])
    parfor i=1:nTests
        
%         searchScore=sParLoad('mrscore.mat');
        R=match_query2(clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        if i==7 || i==10
            seizure(i,1)=mean(R(:,2))<searchScore(ii);
        else
            seizure(i,1)=mean(R(:,2))>searchScore(ii);
        end
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
        
        
    end
    
    %assume early probability is same as seizure
    early=seizure;
    resultsCell{ii} = table(clip,seizure,early);
    save('sResults.mat','resultsCell');
    
end

submissionTable = vertcat(resultsCell{:});
writetable(submissionTable);

toc
end


