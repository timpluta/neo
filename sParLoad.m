function [data,freq,latency,channels] = sParLoad( fname )
    load(fname);
    
    %D is name for data in files concatenated by putTogether
    if exist('D','var')
        data=D;
    end
    
    %for general load
    if ~exist('data','var')
        data=searchScore;
        freq=0;
        channels=0;
    end
    
    %they said sometimes some frames got dropped, screws up shazam, this fixes, it assumes 1sec prints
    freq=ceil(freq);
    if length(data)<freq
        data(:,end+1:freq)=data(:,end);
    end
    
    %for test data with no latency,easier for automation to just =0 it
    if ~exist('latency','var')
        latency=0;
    end
    
    
end

