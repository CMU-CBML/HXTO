%   Moving morphable voids and isogeometric analysis:
%% INPUT
% DW = Width of the design domain
% DH = Height of the design domain
% nelx = Number of elements of the initial coarse NURBS patch in x-direction
% nely = Number of elements of the initial coarse NURBS patch in y-direction
% numRefine = Number of refinement of the coarse NURBS patch
% Nvx = Number of voids in x-direction
% Nvy = Number of voids in y-direction
% Nvicp = Number of independent comtrol points in every void
%
function [test] = top_iga_mmv_new(update,cpt_matrix,DW,DH,nelx,nely,Nvx,Nvy,Nvicp)
%% Information of Material and Voids
% E = 1000; 
% nu = 0.3; 
% th = 0.1; 
% D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2.0];
% update = update';
pvb = 2; 
%% Initialization of NURBS basis functions of the Isogeometric mesh
px = 2; py = 2;
numknotX = nelx + 1; numknotY = nely + 1;

% CP = NurG.controlPoints;
CP = cpt_matrix;
% Weights = NurG.weights;
% numCPX = length(U) - px - 1;
numCPX = numknotX;
% numCPY = length(V) - py - 1;
numCPY = numknotY;
CPgrid.x = reshape(CP(:,1),numCPY,numCPX);
CPgrid.y = reshape(CP(:,2),numCPY,numCPX);

%% Preparation for Isogeometric analysis
%% Initialization of design variables 
NV = Nvx*Nvy;
dvx = DW/Nvx; 
dvy = DH/Nvy; 

% NV = 5;
VoidCenterx = dvx/2:dvx:DW;    VoidCenterx = repmat(VoidCenterx,1,Nvy)';
VoidCentery = dvy/2:dvy:DH;    VoidCentery = kron(VoidCentery,ones(1,Nvx))';
VoidCenter = [VoidCenterx,VoidCentery];

cent_update = zeros(NV,2);
% for ic = 1:NV
for ic = 1:NV/2
    cent_update(ic,1) = update(2*ic-1,1)*1.0;
    cent_update(ic,2) = update(2*ic,1)*1.0;
end
% cent_update(3,1) = update(5,1)*1.0;
VoidCenter = VoidCenter + cent_update;

d0 = 0.0; %% min(dvx,dvy)/8;
%% FOR circle
dd = repmat(d0,NV,Nvicp);
dd_update = zeros(NV,Nvicp);
n_start = 0;
n_end = 0;
for idd = 1:NV/2
    n_start = NV +(idd-1)*Nvicp + 1;
    n_end = NV +idd*Nvicp;
    dd_update(idd,:) = update(n_start:n_end,1)*1.0;
end
dd = dd + dd_update;

angle = 2*pi/Nvicp;
theta = linspace(angle/2,2*pi-angle/2,Nvicp);
VoidCPax = sin(repmat(theta,NV,1)).*dd;
VoidCPay = cos(repmat(theta,NV,1)).*dd;
VoidCP(:,:,1) = (repmat(VoidCenter(:,1),1,Nvicp) + VoidCPax);
VoidCP(:,:,2) = (repmat(VoidCenter(:,2),1,Nvicp) + VoidCPay);
% DESIGNVARIABLES = [VoidCenter,dd]';
% xval = DESIGNVARIABLES(:); 
%% Forming enclosed void boundary curves
VBKnotVector = [zeros(1,pvb),linspace(0,1,Nvicp+1),ones(1,pvb)];
plotVBKnots = 0:0.005:VBKnotVector(end);
VoidAllCP = zeros(NV,Nvicp+2,2);
VoidAllCP(:,2:end-1,:) = VoidCP;
VoidAllCP(:,1,:) = (VoidAllCP(:,2,:) + VoidAllCP(:,end-1,:))/2;
VoidAllCP(:,end,:) = VoidAllCP(:,1,:);
[VBCurveX,VBCurveY] = VoidBoundaryCurve(NV,Nvicp+1,pvb,VBKnotVector,VoidAllCP,plotVBKnots);

%% For symmetry constraint 
VBCurveX(NV/2+1:NV,:) = VBCurveX(1:NV/2,:);
VBCurveY(NV/2+1:NV,:) = DH - VBCurveY(1:NV/2,:);

%% Assemble the voids
[CPgridm,CPgridn] = size(CPgrid.x);
Phitemp = zeros(CPgridm,CPgridn,NV);
% [Subsizem,Subsizen] = size(NURBSSubFigSurfaceMesh.x);
% FigPhitemp = zeros(Subsizem,Subsizen,NV);
for vi=1:NV
    in = inpolygon(CPgrid.x,CPgrid.y,VBCurveX(vi,:),VBCurveY(vi,:));
    Phitemp(:,:,vi) = -2*in + 1;
end
Phi = min(Phitemp,[],3);
% Phi = (fliplr(Phi)+Phi)/2;
Phi_label = Phi(:);
test = Phi_label;
