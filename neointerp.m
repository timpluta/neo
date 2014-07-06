function newmatrix=neointerp(matrix,lengthSig,rezBoxes)



%averages y values within a horizontal window of size len/newlen centered
%around points in xx. each original x value contained in that window
%contributes its y value weighted by the percentage of the total window
%occupied by a window of length 1 centered at the original x value.

clear neointerp
persistent intMatrix

len=size(matrix,2);
if len~=size(intMatrix,1)
    newlen = lengthSig*rezBoxes; 

    xx=1:(len-1)/(newlen-1):len;
    rr=zeros(1,ceil(len^2/newlen));
    cc=zeros(1,ceil(len^2/newlen));
    proportions=zeros(1,ceil(len^2/newlen));
    entryCount=1;

    for nn=1:newlen
        try 
            lbound = 3*xx(nn-1)/4 + xx(nn+1)/4;
        catch me
            switch nn
                case 1
                    lbound = xx(1);
                    ubound = (xx(1)+xx(2))/2;
                case newlen
                    lbound = (xx(end-1)+xx(end))/2;
                    ubound = xx(end);
            end
        end
        leftXold=round(lbound);
        try
            ubound = xx(nn-1)/4 + 3*xx(nn+1)/4;
        catch me
        end
        rightXold=round(ubound);
        lfraction = .5 + (leftXold - lbound);
        ufraction = .5 - (rightXold - ubound);
        baseContribution=1/(rightXold-leftXold+lfraction+ufraction-1);

        for nnn=leftXold:rightXold
            rr(entryCount)=nnn;
            cc(entryCount)=nn;
            if nnn==leftXold
                proportions(entryCount)=lfraction*baseContribution;
            elseif nnn==rightXold
                proportions(entryCount)=ufraction*baseContribution;
            else
                proportions(entryCount)=baseContribution;
            end
            entryCount=entryCount+1;
        end
    end
    rr=rr(rr~=0);
    cc=cc(cc~=0);
    proportions=proportions(1:length(rr));
    intMatrix=sparse(rr,cc,proportions);
end

newmatrix=double(matrix)*intMatrix;
% toc
end
