function [NNLt] = createNNLtriangle(V, Fb, sn, nsn, hs, bw, mw)

NNLt = cell(1,nsn);
mx = max(1,int8(bw/mw));
head(1,1:mx*mx*mx) = -1;
list = zeros(1,size(Fb,1));
for i = 1:size(Fb,1)
    cog = (V(Fb(i,1),:)+V(Fb(i,2),:)+V(Fb(i,3),:))/3.0;
    xa = int8((cog(1)+0.5*bw)/bw*mx);
    ya = int8((cog(2)+0.5*bw)/bw*mx);
    za = int8((cog(3)+0.5*bw)/bw*mx);
    tmp = mx*mx*za+mx*ya+xa;
    list(1,i) = head(1,tmp);
    head(1,tmp) = i;
end;

parfor i = 1:nsn
    pt = sn(i);
    xa = int8((V(pt,1)+0.5*bw)/bw*mx);
    ya = int8((V(pt,2)+0.5*bw)/bw*mx);
    za = int8((V(pt,3)+0.5*bw)/bw*mx);
    for xi = max(0,xa-1):min(mx-1,xa+1)
    for yi = max(0,ya-1):min(mx-1,ya+1)
    for zi = max(0,za-1):min(mx-1,za+1)
                tri = head(1,mx*mx*zi + mx*yi + xi);
                while tri~=-1
                    if pt~=Fb(tri,1)&&pt~=Fb(tri,2)&&pt~=Fb(tri,3)
                        [clp,u,v,w] = closestPointTriangle(V(pt,:), V(Fb(tri,1),:), V(Fb(tri,2),:), V(Fb(tri,3),:));
                        if length(clp - V(pt,:)) < hs
                            NNLt{i}=[NNLt{i},tri];
                        end;
                    end;
                    tri = list(1,tri);
                end;
            end;
        end;
    end;
end;
          
          
