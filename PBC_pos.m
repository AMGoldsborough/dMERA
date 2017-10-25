function x = PBC_pos(x,L)
%x = PBC_pos(x,L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - PBC_pos
% gives position (x) with PBCs with length L
% 
% Andrew Goldsborough - 07/10/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = mod(x-1,L)+1;
end