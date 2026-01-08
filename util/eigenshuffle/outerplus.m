function [minval,rIdx,cIdx]=outerplus(M,x,y)
[nx,ny]=size(M);
minval=inf;
for r=1:nx
    x1=x(r);
    for c=1:ny
        M(r,c)=M(r,c)-(x1+y(c));
        if minval>M(r,c)
            minval=M(r,c);
        end
    end
end
[rIdx,cIdx]=find(M==minval);
end