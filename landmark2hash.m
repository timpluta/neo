function H = landmark2hash(L,S)
% H = landmark2hash(L,S)
%  Convert a set of 4-entry landmarks <t1 f1 f2 dt> 
%  into a set of <songid time hash> triples ready to store.
%  S is a scalar songid, or one per landmark (defaults to 0)
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

% Hash value is 20 bits: 8 bits of F1, 6 bits of delta-F, 6 bits of delta-T
%number of time boxes used
targetdt=3;

if nargin < 2
  S = 0;
end
if length(S) == 1 && S~=0
  S = repmat(S, size(L,1), 1);
end

% if S==0
%     chan=1;
%     S=[];
%     index = true(1, size(L, 1));
%     for i=1:size(L,1)
%         if sum(L(i,:))==0
%             chan=chan+1;
%             index(i)=false;
%         else
%             S(i)=chan;
%         end
%     end
%     L = L(index, :);
%     S=S(index)';
%     CH=rem(S-1,2^4);
% else
%     CH = rem(S-1,2^4);
% end
DF = round(L(:,3)-L(:,2)); 
%replace DF with Ch1 - Ch2?

H = uint32(L(:,1));
% Make sure F1 is 0..255, not 1..256
F1 = rem(round(L(:,2)-1),2^6);
if DF < 0
  DF = DF + 2^6;
end

DF = rem(DF,2^6);
DT = rem(abs(round(L(:,4))), targetdt);
CH1= L(:,5);
CH2= L(:,6);
S=CH1;

try
   H = [S,H,uint32(CH1*(2^19*targetdt)+CH2*(2^12*targetdt)+F1*(2^6*targetdt)+DF*(targetdt)+DT)];
catch err
   disp(L);
   disp(S);
end
%tim added 2014-06-26 to fix all 0 hashes but dont know why they happen
% H( ~any(H,2), : ) = [];


