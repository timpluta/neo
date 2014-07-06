function [cfs scales] = morFinger1sec(D,dt,SCA)
%% 

if size(D,1)>1
    error('morFinger1sec must be passed to 1 channel at a time')
end
    
MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));  

cwtsig = cwtft({D,dt},'scales',SCA,'wavelet','morl');
MM=abs(cwtsig.cfs);

cfs=flipud(MM);
scales = cwtsig.scales.*MorletFourierFactor;
end




