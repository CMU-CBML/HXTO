function [Cx,Cy] = CurvePoint(n,p,U,CP,u)
uspan = FindSpan(n,p,u,U);
Nu = Der1BasisFuns(uspan,u,p,U);
uind = uspan - p;
Cx = zeros(size(u,2),1);
Cy = zeros(size(u,2),1);
for i = 0:p
    tempCP = CP(uind+i+1,:);
    Cx = Cx + Nu(i+1,:)' .* tempCP(:,1);
    Cy = Cy + Nu(i+1,:)' .* tempCP(:,2);
end