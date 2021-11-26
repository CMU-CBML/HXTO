function [VBCurveX,VBCurveY] = VoidBoundaryCurve(NV,n,p,U,CP,u)
sVBU = size(u,2);
VBCurveX = zeros(NV,sVBU); 
VBCurveY = zeros(NV,sVBU);
for vi = 1:NV
    viCP = squeeze(CP(vi,:,:));
    [tempCx,tempCy] = CurvePoint(n,p,U,viCP,u);
    VBCurveX(vi,:) = tempCx;
    VBCurveY(vi,:) = tempCy;
end