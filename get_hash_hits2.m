function R = get_hash_hits2(H,clipLoc)
% R = get_hash_hits(H)
%    Return values from song hash table for particular hashes
%    Each element of H is a <(20 bit) hash value>
%    Each row of R is a hit in format:
%    <song id> <start time index> <hash>
%    If H is a 2 column matrix, the first element is taken as a
%    time base which is subtracted from the start time index for
%    the retrieved hashes.
%    If H is a 3 column matrix, the first element is taken as a
%    songID and discarded.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

persistent z

%load first found hashtable in clipLoc
a=dir(strcat(clipLoc,filesep,'*Hush*.mat'));

% %load first clip table
% a=dir(strcat(clipLoc,filesep,'*By*.mat'));

if size(a,1)==0
    error(strcat('no hash table in ',clipLoc))
end
tableName=a(1).name;

if isempty(z)
    load(strcat(clipLoc,filesep,tableName));
end


if size(H,2) == 3
    clips=H(:,1);
  H = H(:,[2 3]);
end

if min(size(H))==1
  H = [zeros(length(H),1),H(:)];
end


Rsize = 100;
% R = zeros(Rsize,3);
R=zeros(1,3);
Rmax = 0;

for i = 1:size(H,1)
  hash = H(i,2);
%   htime = double(H(i,1));

%not by clip
try
  R(1,2)=R(1,2)+nnz(z(:,hash));
catch me
    hash
end

% %by clip
% try
%   R(1,2)=R(1,2)+length(find(nonzeros(z(:,hash))~=str2double(clips(i))));
% catch me
%     hash
% end

%   try
% %       htcol = nonzeros(z(:,hash));
%       htcol = z(z(:,hash)~=0,hash);
%   catch me
%       htcol=[];
%   end
%   nentries = length(htcol);
%   
%   for j = 1:nentries
%     time = htcol(j);
%       %if time... line tim added to not count t=0 hashes
% %     if time==htime
% %     else
%         Rmax = Rmax + 1;
%         if Rmax > Rsize
%           R = [R;zeros(Rsize,3)];
%           Rsize = size(R,1);
%         end
%         dtime = time-htime;
%         R(Rmax,:) = [1, dtime, double(hash)];
% %     end
%   end
end

% R = unique(R(1:Rmax,:),'rows');





% 
% 
% disp('old')
% 
% 
% 
% 
% 
% 
% 
% 
% tic
% persistent HashTable
% 
% if size(H,2) == 3
%   H = H(:,[2 3]);
% end
% 
% if min(size(H))==1
%   H = [zeros(length(H),1),H(:)];
% end
% 
% 
% %load first found hashtable in clipLoc
% a=dir(strcat(clipLoc,filesep,'*Hash*.mat'));
% if size(a,1)==0
%     error(strcat('no hash table in ',clipLoc))
% end
% tableName=a(1).name;
% 
% if isempty(HashTable)
%     load(strcat(clipLoc,filesep,tableName));
% end
% 
% 
% 
% Rsize = 100;
% R = zeros(Rsize,3);
% Rmax = 0;
% 
% for i = 1:length(H)
%   hash = H(i,2);
%   htime = double(H(i,1));
%   try
%       htcol = HashTable(hash+1);
%   catch me
%       htcol=[];
%   end
%   nentries = size(htcol,1);
% 
%   if size(htcol,2)==0
%       nentries=0;
%   end
%   for j = 1:nentries
%     song = htcol(j,1);
%     time = htcol(j,2);
%       %if time... line tim added to not count t=0 hashes
% %     if time==htime
% %     else
%         Rmax = Rmax + 1;
%         if Rmax > Rsize
%           R = [R;zeros(Rsize,3)];
%           Rsize = size(R,1);
%         end
%         dtime = time-htime;
%         R(Rmax,:) = [double(song), dtime, double(hash)];
% %     end
%   end
% end
% 
% R = unique(R(1:Rmax,:),'rows');
% 
% 
% 
% 
% toc




