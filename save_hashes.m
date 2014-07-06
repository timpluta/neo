function save_hashes(H,clipLoc,clipName)
% save_hashes(H)
%   Record the set of hashes that are rows of H in persistent
%   database.
%   Format of H rows are 3 columns:
%   <song id> <start time index> <hash>
% song ID is 32 bit
% time index is 32 bit
%   (32 ms basic resolution = 30/sec, so 600 sec song has time indices
%    up to 18,000 = 14 bits?)
% Hash is 20 bit = 1M slots
%
% 2008-12-24 Dan Ellis dpwe@ee.columbia.edu

% This version uses an in-memory global with one row per hash
% value, and a series of song ID / time ID entries per hash


% if exist('HashTable','var') == 0 || length(HashTable) == 0
%    clear_hashtable;
% end


% tic
% %maxnentries = size(HashTable,1)/2;
% persistent HashTable
% 
% %workaround initialization of map structure
% if isempty(HashTable)
%     HashTable=containers.Map(uint32(1), [1 1]);
%     remove(HashTable,1);
% end
% 
% nhash = size(H,1);
% for i=1:nhash
%   time = H(i,2);
%   %if time< line added by tim 2014-06-28 to fix unknown t=0 hash problem
% %   if time<.1
% %   else
%       song = H(i,1);
% 
%       hash = H(i,3)+1;  % avoid hash == 0
%     %  htcol = HashTable(:,hash);
%     %  nentries = max([0;find(htcol ~= 0)])/2;
%     %  if nentries < maxnentries
%     % HashTable(i+bounds-1).channel=song;
%     % HashTable(i+bounds-1).time=time;
%     % HashTable(i+bounds-1).hash=hash;
%         try 
%             HashTable(hash) = [HashTable(hash); [song time]];
%         catch me
%             HashTable(hash)= [song time];
%         end
%     %  end
% %   end
% end
%  save(strcat(clipLoc,filesep,clipName,'_HashTable.mat'),'HashTable','-v7.3')
% toc


tic

i=1;
j=2;
H=sortrows(H,3);
a(1,:)=uint32(find(diff(H(:,3))==0));
rs=ones(1,size(H,1));
na=a(i,:);

while sum(a(i,:))~=0
% H(a(1):a(1)+1,:)
    if i>1
        a(i-1,:)=a(i-1,:)+1;
        na=a(i-1,:);
        na=na(a(i,a(i,:)~=0));
        a(i,1:length(na))=na;
        a=a(2:end,1:size(na, 2));
    else
        i=i+1;
    end
    
    
    b=uint32(find(diff(a(i-1,:))==1));
    a(i,1:length(b))=b;
    
    rs(na+1)=j;
    j=j+1;
end

z=sparse(rs',double(H(:,3)),double(H(:,2)),j,double(max(H(:,3))));
save(strcat(clipLoc,filesep,clipName,'_HushTable.mat'),'z','-v7.3');
toc



