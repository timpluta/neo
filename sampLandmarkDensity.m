% show landmark density across channels

%your matlab directory
matdir='/Users/idaniboy/Documents/MATLAB/';
%matdir='C:\Users\paul\Documents\MATLAB\';

%clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
clipsdir='kaggleClips/ictalSegs/';

%subject to run it on (do 1 at a time for now)
%clipDirNames={'Dog_2'};
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};


for ii = 1:length(clipDirNames)
    % %automation will require data to be in separate folders per patient with _ in folder name
    % a=dir(strcat(matdir,clipsdir,'*_*'));
    % clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    
    %type of segment to view
    typeI='_ictal_segment';
    % typeI='*_test_segment*';

    %number of training segment randsamples to view 
    nSamples=30;
    a=dir(strcat(clipLoc,filesep,'*',typeI,'*.mat'));
    srsOfSegments={a(randsample(length(a),min(nSamples,length(a)))).name};
    nTests=length(srsOfSegments);

    for i=1:nTests

        [data,freq,latency,channels]=sParLoad(strcat(clipLoc,filesep,srsOfSegments{i}));
        
        if ii == 1
            preLq= findLandmarks3(data,freq,channels);
        else
            preLq= find_landmarks2(data,freq,channels);
        end
        
        Ld(ii,i) = size(preLq,1)/length(fieldnames(channels)); % landmark density = #landmarks in clip / # channels

    end
    disp(strcat('current subject:',num2str(ii)));
end

%% Plot landmark densities
% hist(Ld,1:10:100)
linespec = {'b-.o', 'r-', 'g--o','m--p','c-v','b-','g-','r-^', 'g--o','m-o','c-','y-','rv','r--o', 'b--o'};
figure
for i=1:size(Ld,1)
    plot(Ld(i,:),linespec{i})
    hold on
end
title('landmark density: low-mid neo for p6');
hleg = legend(clipDirNames,...
              'Location','NorthEastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])