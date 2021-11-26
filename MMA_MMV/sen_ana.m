function [d_var] = sen_ana(delta,update,cpt_matrix,DW,DH,nelx,nely,nvx,nvy,nvicp)
% Summary of this function goes here
% Detailed explanation goes here
[m,n] = size(update);
d_var = zeros(m,2);

[HXP_0,PD_0] = HXP_sen(update,cpt_matrix,DW,DH,nelx,nely,nvx,nvy,nvicp);

HXP = zeros(m,n);
PD = zeros(m,n);

for i = 1:m
    update_iter = update;
    update_iter(i) = update(i) + delta;
    [HXP_i,PD_i] = HXP_sen(update_iter,cpt_matrix,DW,DH,nelx,nely,nvx,nvy,nvicp);
    HXP(i) = HXP_i;
    PD(i) = PD_i;
end

d_var(:,1) = (HXP - HXP_0)./delta;
d_var(:,2) = (PD - PD_0)./delta;

end

