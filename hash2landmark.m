function L = hash2landmark(H)
% L = hash2landmark(H)
%  Convert a set of <time hash> pairs ready from store 
%  into a set of 4-entry landmarks <t1 f1 f2 dt>.
%  If H is 3 cols, first col (song ID) is discarded.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

% Hash value is 20 bits: 8 bits of F1, 6 bits of F2-F1, 6 bits of delta-T

dthash=3;
dfhash=2^6;
f1hash=2^6;
ch2hash=2^6;
ch1hash=2^6;
%H = uint32(L(:,1));
%F1 = rem(round(L(:,2)),2^7);
%F2 = rem(round(L(:,4)),2^7);
%DT = rem(abs(round(L(:,3) - L(:,1))), 2^6);
%H = [H,uint32(F1*(2^13)+F2*(2^6)+DT)];

if size(H,2) == 3
  H = H(:,[2 3]);
end

H1 = H(:,1);
H2 = double(H(:,2));
CH1 = floor(H2/(ch2hash*f1hash*dfhash*dthash));
H2 = H2 - ch2hash*f1hash*dfhash*dthash * CH1;
CH2 = floor(H2/(f1hash*dfhash*dthash));
H2 = H2 - f1hash*dfhash*dthash * CH2;
F1 = floor(H2/(dfhash*dthash));
H2 = H2 - dfhash*dthash * F1;
F1 = F1 + 1;
DF = floor(H2/(dthash));
H2 = H2 - dthash * DF;
if DF > 2^5
  DF = DF-dfhash;
end
F2 = F1+DF;

DT = H2;

L = [H1,F1,F2,DT,CH1,CH2];
