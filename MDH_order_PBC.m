function Jmaxpos_vec = MDH_order_PBC(L,J,Jz)
%Jmaxpos_vec = MDH_order_PBC(L,J,Jz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - MDH_order_PBC
% find singlet order using MDH
% 
% Andrew Goldsborough - 30/09/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jmaxpos_vec = zeros(1,(L/2)-1);

for iteration = 1:(L/2)-1
    %find the maximum J
    [~,Jmaxpos] = max(abs(J));
    Jmaxpos_vec(iteration) = Jmaxpos;
    
    %renormalise coupling
    J(PBC_pos(Jmaxpos-1,size(J,2))) = J(PBC_pos(Jmaxpos-1,size(J,2)))*J(PBC_pos(Jmaxpos+1,size(J,2)))/((1+Jz)*J(Jmaxpos));
    
    %remove two sites
    J(Jmaxpos) = [];
    J(PBC_pos(Jmaxpos,size(J,2))) = [];
end
end