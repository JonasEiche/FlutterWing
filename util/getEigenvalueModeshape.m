function [EV, MS] = getEigenvalueModeshape(CL,num_modes)

num_vel = length(CL);
orderCL=order(CL);
assert(all(orderCL == orderCL(1)),'All Cl must have the same order')
num_X = orderCL(1);
PHIV = zeros(num_modes,num_X,num_vel); 
EV   = zeros(num_X,num_vel);
for i_v = 1:num_vel
    CL_i = CL(:,:,i_v);
    [phi,e]        =  eig(full(CL_i.A));
    PHIV(:,:,i_v)  =  phi(num_modes+1:num_modes*2,:);
    EV(:,i_v)      =  diag(e);
    % track modes
    if i_v>=2
	    idx = eigenshuffle(EV(:,i_v-1),PHIV(:,:,i_v-1),EV(:,i_v),PHIV(:,:,i_v));
    else
        idx = 1:size(EV,1);
    end
    PHIV(:,:,i_v) = PHIV(:,idx,i_v);
    EV(:,i_v)     = EV(idx,i_v);
end

[~, idx] = sort(EV(:,1),'descend','ComparisonMethod','real');
% [~, idx] = sort(EV(:,num_vel),'descend','ComparisonMethod','real'); % use the last velocity to determine pole order
EV = EV(idx,:); 
MS = PHIV(:,idx,:);

end