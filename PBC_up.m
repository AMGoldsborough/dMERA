function x = PBC_up(x,Jmaxpos,L)
%x = PBC_up(x,Jmaxpos,L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - PBC_up
% gives position of the site one level up
% 
% Andrew Goldsborough - 15/11/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = PBC_pos(x,L);
if Jmaxpos == L
    if any(x == [Jmaxpos-1,Jmaxpos,1]) == 1
        x = Jmaxpos - 2;
    else
        x = x - 1;
    end
else
    if any(x == [PBC_pos(Jmaxpos-1,L),Jmaxpos,PBC_pos(Jmaxpos+1,L)]) == 1
        x = PBC_pos(Jmaxpos - 1,L-2);
    elseif x > Jmaxpos
        x = x - 2;
    end
end

%slower but neater method
% if min(mod(x-Jmaxpos,L),mod(Jmaxpos-x,L)) <= 1
%     x = PBC_pos(Jmaxpos - 1,L);
% end
% 
% lattice = 1:L;
% lattice(Jmaxpos) = [];
% lattice(PBC_pos(Jmaxpos,L-1)) = [];
% x = find(lattice == x);
end