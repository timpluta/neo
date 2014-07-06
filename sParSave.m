function [a]=sParSave( fname,channels,D,freq )
    %easier if concatenations have arbitrary latency value
    latency=0;
    save(fname,'channels','D','freq','latency');
end

