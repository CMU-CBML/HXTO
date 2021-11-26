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
% CPx0 = linspace(0,DW,numknotX+px-1);
% CPy0 = linspace(0,DH,numknotY+py-1);
% CPx0 = repmat(CPx0,1,numknotY+py-1)';
% CPy0 = kron(CPy0,ones(1,numknotX+px-1))';
% CP_First = [CPx0,CPy0]; 
% Weights_First = ones((numknotX+px-1)*(numknotY+py-1),1); 
% U_First = [zeros(1,px),linspace(0,1,numknotX),ones(1,px)]; 
% V_First = [zeros(1,py),linspace(0,1,numknotY),ones(1,py)];
% %% Refine the initial Isogeometric mesh, and distribute Control Points nonuniformly
% NurG = isogat(px,py,U_First,V_First,CP_First,Weights_First);
% if numRefine < 0
%     error('the Parameter numRefine must be Non-Negative!!!') 
% elseif numRefine > 0
%     for i = 1:numRefine
%         NurG = refine(NurG);
%     end
% end
% U = NurG.knotVectorU;
% V = NurG.knotVectorV;
% CP = NurG.controlPoints;
CP = cpt_matrix;
% Weights = NurG.weights;
% numCPX = length(U) - px - 1;
numCPX = numknotX;
% numCPY = length(V) - py - 1;
numCPY = numknotY;
% KnotVectorU = unique(U);
% KnotVectorV = unique(V);
% IntervalU = KnotVectorU(2:end) - KnotVectorU(1:end-1); 
% IntervalV = KnotVectorV(2:end) - KnotVectorV(1:end-1); 
CPgrid.x = reshape(CP(:,1),numCPY,numCPX);
CPgrid.y = reshape(CP(:,2),numCPY,numCPX);
% eledof = 2*(px+1)*(py+1); 
% sysdof = 2*numCPX*numCPY; 
%% Preparation for Isogeometric analysis
% Gauss Points and the corresponding weights of Gauss Points
% [GaussIPU,GaussWeightU] = gquad(KnotVectorU,px+1); 
% [GaussIPV,GaussWeightV] = gquad(KnotVectorV,py+1); 
% [KE,IndexCP,RBasisGP] = IGABasicKe_GS(D,eledof,px,py,U,V,KnotVectorU,KnotVectorV,CP,numCPX,numCPY,Weights,th,GaussIPU,GaussWeightU,GaussIPV,GaussWeightV);
% IndexCPdof(1:2:eledof,:) = 2 * IndexCP - 1; IndexCPdof(2:2:eledof,:) = 2 * IndexCP; 
% iK = kron(IndexCPdof,ones(1,eledof));
% jK = kron(IndexCPdof,ones(eledof,1));

%% Initialization of design variables 
NV = Nvx*Nvy;
dvx = DW/Nvx; 
dvy = DH/Nvy; 

% if mod(Nvy, 2) == 1
%     NV = NV + Nvx;
% end 

% NV = 5;
VoidCenterx = dvx/2:dvx:DW;    VoidCenterx = repmat(VoidCenterx,1,Nvy)';
VoidCentery = dvy/2:dvy:DH;    VoidCentery = kron(VoidCentery,ones(1,Nvx))';
% VoidCenterx = DW*[1/4,3/4,1/2,1/4,3/4];
% VoidCenterx = VoidCenterx';
% VoidCentery = DH*[1/4,1/4,1/2,3/4,3/4];
% VoidCentery = VoidCentery';
VoidCenter = [VoidCenterx,VoidCentery];
% if mod(Nvy, 2) == 1
%     VoidCenter = [VoidCenter;zeros(Nvx,2)];
% end

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
% for idd = 1:NV
% for idd = 1:NV/2
% %    n_start = 2*NV+(idd-1)*Nvicp + 1;
% %    n_end = 2*NV+idd*Nvicp;
%     n_start = NV+(idd-1)*Nvicp/4 + 1;
%     n_end = NV+idd*Nvicp/4;
%     dd_update(idd,1:Nvicp/4) = update(n_start:n_end,1)*1.0;
%     dd_update(idd,Nvicp/2:-1:Nvicp/4+1) = update(n_start:n_end,1)*1.0;
%     n_start = NV+idd*Nvicp/4 + 1;
%     n_end = NV+idd*Nvicp/2;
%     dd_update(idd,Nvicp/2+1:Nvicp*3/4) = update(n_start:n_end,1)*1.0;
%     dd_update(idd,Nvicp:-1:Nvicp*3/4+1) = update(n_start:n_end,1)*1.0;
% end

n_start = 0;
n_end = 0;
for idd = 1:NV/2
%    n_start = 2*NV+(idd-1)*Nvicp + 1;
%    n_end = 2*NV+idd*Nvicp;
    n_start = NV +(idd-1)*Nvicp + 1;
    n_end = NV +idd*Nvicp;
    dd_update(idd,:) = update(n_start:n_end,1)*1.0;
end
% dd_update(3,1:Nvicp/4) = update(n_end+1:n_end+Nvicp/4,1)*1.0;
% dd_update(3,Nvicp/2:-1:Nvicp/4+1) = update(n_end+1:n_end+Nvicp/4,1)*1.0;
% dd_update(3,Nvicp/2+1:3*Nvicp/4) = update((n_end+Nvicp/4+1):n_end+Nvicp/2,1)*1.0;
% dd_update(3,Nvicp:-1:(3*Nvicp/4+1)) = update((n_end+Nvicp/4+1):n_end+Nvicp/2,1)*1.0;
dd = dd + dd_update;
% 
dd_tmp = dd;
% % for i_nv = 1:NV/2
% for i_nv = 1:2
%     for iter = 1:Nvicp
%         if iter == 1
%             dd(i_nv, iter) = (dd_tmp(i_nv, Nvicp) + 4*dd_tmp(i_nv, iter) + dd_tmp(i_nv, iter+1))/6;
%         elseif iter == Nvicp
%             dd(i_nv, iter) = (dd_tmp(i_nv, iter-1) + 4*dd_tmp(i_nv, iter) + dd_tmp(i_nv, 1))/6;
%         else
%             dd(i_nv, iter) = (dd_tmp(i_nv, iter-1) + 4*dd_tmp(i_nv, iter) + dd_tmp(i_nv, iter+1))/6;
%         end
%     end
% end

% dd_tmp = dd;
% for i_nv = 1:NV/2
%     for iter = 1:Nvicp
%         if iter == 1
%             dd(i_nv, iter) = (dd_tmp(i_nv, iter) + dd_tmp(i_nv, iter+1))/2;
%         elseif iter == Nvicp
%             dd(i_nv, iter) = (dd_tmp(i_nv, iter) + dd_tmp(i_nv, 1))/2;
%         else
%             dd(i_nv, iter) = (dd_tmp(i_nv, iter) + dd_tmp(i_nv, iter+1))/2;
%         end
%     end
% end

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
%% For boundary gap constraint
% gap = 0.05;
% for id = 1:NV/2
%     if min(VBCurveX(id,:)) < gap
%          VBCurveX(id,:) = VBCurveX(id,:) + abs(min(VBCurveX(id,:))-gap);
%     end
%     if min(VBCurveY(id,:)) < gap
%          VBCurveY(id,:) = VBCurveY(id,:) + abs(min(VBCurveY(id,:))-gap);
%     end
%     if max(VBCurveX(id,:)) > DW-gap
%          VBCurveX(id,:) = VBCurveX(id,:) - abs(max(VBCurveX(id,:))-DW + gap);
%     end
%     if max(VBCurveY(id,:)) > DH-gap
%          VBCurveY(id,:) = VBCurveY(id,:) - abs(max(VBCurveY(id,:))-DH + gap);
%     end
% end
%% For symmetry constraint 
VBCurveX(NV/2+1:NV,:) = VBCurveX(1:NV/2,:);
VBCurveY(NV/2+1:NV,:) = DH - VBCurveY(1:NV/2,:);
% VBCurveX(3+1:NV,:) = VBCurveX(1:2,:);
% VBCurveY(3+1:NV,:) = DH - VBCurveY(1:2,:);
% %% Parameters of MMA
% xold1 = xval;
% xold2 = xval;
% xmin = [0,0,linspace(0,0,Nvicp)];
% xmin = repmat(xmin',NV,1);
% low = xmin;
% xmax  = [DW,DH,linspace(min(DW,DH),min(DW,DH),Nvicp)];
% xmax = repmat(xmax',NV,1);
% upp = xmax;
% nn = NV*(Nvicp+2); 
% m = 1; 
% c = 1000 * ones(m,1);
% d = zeros(m,1);
% a0 = 1;
% a = zeros(m,1);
%% Defination of loads and supports (For example, Short Beam A)
% alldofs = 1:sysdof;
% fixedCPs = find(CP(:,1) == 0);
% fixeddofs = [2*fixedCPs-1;2*fixedCPs];
% loadCPs =numCPX * [numCPY/2,numCPY/2+1]; 
% loaddofs = 2 * loadCPs;
% freedofs = setdiff(alldofs,fixeddofs);
% F = sparse(loaddofs,1,-0.5,sysdof,1);
%% Assemble the voids
[CPgridm,CPgridn] = size(CPgrid.x);
Phitemp = zeros(CPgridm,CPgridn,NV);
% [Subsizem,Subsizen] = size(NURBSSubFigSurfaceMesh.x);
% FigPhitemp = zeros(Subsizem,Subsizen,NV);
for vi=1:NV
    in = inpolygon(CPgrid.x,CPgrid.y,VBCurveX(vi,:),VBCurveY(vi,:));
    Phitemp(:,:,vi) = -2*in + 1;
%     Figin = inpolygon(NURBSSubFigSurfaceMesh.x,NURBSSubFigSurfaceMesh.y,VBCurveX(vi,:),VBCurveY(vi,:));
%     FigPhitemp(:,:,vi) = -2*Figin + 1;
end
Phi = min(Phitemp,[],3);
% Phi = (fliplr(Phi)+Phi)/2;
Phi_label = Phi(:);
test = Phi_label;
% FigPhi = min(FigPhitemp,[],3);
% SignFigPhi = SignDist(NURBSSubFigSurfaceMesh.x,NURBSSubFigSurfaceMesh.y,VBCurveX,VBCurveY,sign(FigPhi));
% 
% CountPhi = Phi;
% CountPhi(CountPhi == -1) =1e-6;
% %% Draw initial figures
% parent_dir_name1 = 'Topology Results\';
% mkdir(parent_dir_name1);
% h1 = figure(1);
% hold on
% pl1=ones(DW*10,1);
% pl2=zeros(10,1);
% pl3=linspace(0,DH,10);
% pl4=linspace(0,DW,DW*10);
% plot(pl2      , pl3    , 'k','linewidth',2.0);
% plot(pl2 + DW , pl3    , 'k','linewidth',2.0);
% plot(pl4      , pl1-1  , 'k','linewidth',2.0);
% plot(pl4      , pl1*DH , 'k','linewidth',2.0);
% for vi = 1:NV
%     plot(VBCurveX(vi,:),VBCurveY(vi,:),'k','linewidth',2.0);
% end
% axis manual
% axis image;
% % axis off % axis equal;
% axis([0 DW 0 DH])
% set(gca,'color','none')
% set(gcf,'color','white')
% xticks(0:0.5:DW);
% yticks(0:0.5:DH);
% hold off
% 
