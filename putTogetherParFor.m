
%your matlab directory
matdir='/Users/idaniboy/Documents/MATLAB/';
%matdir='C:\Users\paul\Documents\MATLAB\';
%clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
clipsdir='kaggleClips/ictalSegs/';

% %upper limit for concatenated matrix length, only use if memory problems
% cMaxSize=1e7;

% patients={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};

patients={'Patient_8'};


%not useful, just initializing
offset=zeros(1,4);

for ii=1:length(patients) 

    %concatenate ictal
    %find info
    nIctalClips=length(dir(strcat(matdir,clipsdir,patients{ii},'/*_ictal_s*.mat')));
    [data,freq,latency,channels]=sParLoad(strcat(patients{ii},'_ictal_segment_1'));
    nSamplesPerSegment=size(data,2);
    nChannels=size(data,1);
    nthSeizure=0;
    D=zeros(nChannels,nIctalClips*nSamplesPerSegment);

    for i=1:nIctalClips
        [data,freq,latency]=sParLoad(strcat(patients{ii},'_ictal_segment_',num2str(i)));
        
        %group individual seizures together. assumes ordered by time.
        if latency==0 || i==nIctalClips
            
            %skip first clip to allow automation
            if i>1 
                %delete trailing zeros from preallocated D
                D=D(:,1:(i-offset(ii)-1)*nSamplesPerSegment);
                saveTo=fullfile(matdir,clipsdir,patients{ii},strcat(patients{ii},'_ictal_concatenated_',num2str(nthSeizure),'.mat'));
                sParSave(saveTo,channels,D,freq);
            end
            offset(ii)=i-1;
            nthSeizure=nthSeizure+1;
            D=zeros(nChannels,nIctalClips*nSamplesPerSegment);
        end
        
        
        
        %write current seizure data into D starting at D(:,1)
        Dindex=i-offset(ii);
        D(:,(Dindex-1)*nSamplesPerSegment+1:Dindex*nSamplesPerSegment)=single(data);
        
    end
    
%     %break long segments up, only use if memory errors
%     newLength=nSamplesPerSegment*nIctalClips;
%     shorter=min(newLength,cMaxSize/nChannels);
%     for j=1:ceil(newLength/shorter)
%         D=zeros(nChannels,min(shorter,nIctalClips-shorter*j));
%         %concatenate
%         for i=(j-1)*shorter+1:min(shorter*j,nIctalClips)
%             load(strcat(subject{1},'_ictal_segment_',num2str(i)))
%             %part of dropped frames fix, assumes all segments are 1sec
%             if length(data)<freq
%                 data(:,end+1:freq)=data(:,end);
%             end
%             D(:,(i-1)*nSamplesPerSegment+1:i*nSamplesPerSegment)=single(data);
%         end
%         save(strcat(subject{1},'_ictal_concatenated_',num2str(j),'.mat'),'channels','D','freq');
%     end
    
    
%     %concatenate interictal
%     %find info
%     nInter=length(dir(strcat(matdir,clipsdir,subject{1},'\*rictal_s*.mat')));
%     load(strcat(subject{1},'_interictal_segment_1'))
%     nSamplesPerSegment=size(data,2);
%     nChannels=size(data,1);
%     sr=freq;
%     %break long segments up to avoid memory problems
%     newLength=nSamplesPerSegment*nInter;
%     shorter=min(newLength,cMaxSize/nChannels);
%     for j=1:ceil(newLength/shorter)
%         D=zeros(nChannels,min(shorter,nInter-shorter*j));
%         %concatenate
%         for i=(j-1)*shorter+1:min(shorter*j,nInter)
%             load(strcat(subject{1},'_interictal_segment_',num2str(i)))
%             %part of dropped frames fix, assumes all segments are 1sec
%             if length(data)<freq
%                 data(:,end+1:freq)=data(:,end);
%             end
%             D(:,(i-1)*nSamplesPerSegment+1:i*nSamplesPerSegment)=single(data);
%         end
%         save(strcat(subject{1},'_interictal_concatenated_',num2str(j),'.mat'),'channels','D','freq');
%     end
    
end