function chi_inc = random_start_PBC(L,chi_w,chi_sing,Jmaxpos_vec,chi_inc,at_on)
%chi_inc = random_start_PBC(L,chi_w,chi_sing,Jmaxpos_vec,chi_inc,at_on)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - random_start_PBC
% create random start tensors
% crossed legs on u
% 
% Andrew Goldsborough - 24/11/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global w u h;

for iteration = 1:(L/2)-2
    Jmaxpos = Jmaxpos_vec(iteration);
    
    Lcur = size(h,1);
    
    if size(h{PBC_pos(Jmaxpos-1,Lcur)},1) * size(h{PBC_pos(Jmaxpos-1,Lcur)},3) <= chi_w && size(h{PBC_pos(Jmaxpos+1,Lcur)},1) * size(h{PBC_pos(Jmaxpos+1,Lcur)},3) <= chi_w
        %increase chi, store
        chi_inc(iteration) = 1;
        
        %generate random tensors
        [uleg1,uleg2,uleg3,uleg4] = deal(size(h{Jmaxpos},3),size(h{Jmaxpos},1),size(h{Jmaxpos},1),size(h{Jmaxpos},3));        
        u{iteration,1} = random('unif',-1,1,[uleg1,uleg2,uleg3,uleg4]);
        
        [w1leg1,w1leg2,w1leg3] = deal(min(size(h{PBC_pos(Jmaxpos-1,Lcur)},1)*size(u{iteration,1},1),chi_w),size(h{PBC_pos(Jmaxpos-1,Lcur)},1),size(u{iteration,1},1));
        w{iteration,1} = random('unif',-1,1,[w1leg1,w1leg2,w1leg3]);
        
        [w2leg1,w2leg2,w2leg3] = deal(min(size(u{iteration,1},3)*size(h{PBC_pos(Jmaxpos+1,Lcur)},3),chi_w),size(u{iteration,1},3),size(h{PBC_pos(Jmaxpos+1,Lcur)},3));
        w{iteration,2} = random('unif',-1,1,[w2leg1,w2leg2,w2leg3]);
        
        %update Hamiltonian
        raiseHam_PBC(Jmaxpos,iteration,chi_inc(iteration));
    else
        %normal structure
        %generate random tensors
        [uleg1,uleg2,uleg3,uleg4] = deal(size(h{Jmaxpos},3),size(h{Jmaxpos},1),size(h{Jmaxpos},1),size(h{Jmaxpos},3));
        u{iteration,1} = random('unif',-1,1,[uleg1,uleg2,uleg3,uleg4]);
        
        [w1leg1,w1leg2,w1leg3,w1leg4] = deal(min(chi_w,ceil(size(h{PBC_pos(Jmaxpos-1,Lcur)},1)*uleg1/chi_sing)),size(h{PBC_pos(Jmaxpos-1,Lcur)},1),chi_sing,uleg1);
        w{iteration,1} = random('unif',-1,1,[w1leg1,w1leg2,w1leg3,w1leg4]);
        
        [w2leg1,w2leg2,w2leg3,w2leg4] = deal(chi_sing,uleg3,min(chi_w,ceil(uleg3*size(h{PBC_pos(Jmaxpos+1,Lcur)},3)/chi_sing)),size(h{PBC_pos(Jmaxpos+1,Lcur)},3));
        w{iteration,2} = random('unif',-1,1,[w2leg1,w2leg2,w2leg3,w2leg4]);
        
        %update Hamiltonian
        raiseHam_PBC(Jmaxpos,iteration,chi_inc(iteration));
    end
end

%top set
iteration = (L/2)-1;
Jmaxpos = Jmaxpos_vec(iteration);
Lcur = size(h,1);

[uleg1,uleg2,uleg3,uleg4] = deal(max(size(h{Jmaxpos},1),size(h{Jmaxpos},3)),size(h{Jmaxpos},1),max(size(h{Jmaxpos},1),size(h{Jmaxpos},3)),size(h{Jmaxpos},3));
u{iteration,1} = random('unif',-1,1,[uleg1,uleg2,uleg3,uleg4]);

[w1leg1,w1leg2,w1leg3,w1leg4] = deal(chi_sing,size(h{PBC_pos(Jmaxpos-1,Lcur)},1),chi_sing,uleg1);
w{iteration,1} = random('unif',-1,1,[w1leg1,w1leg2,w1leg3,w1leg4]);

[w2leg1,w2leg2,w2leg3,w2leg4] = deal(chi_sing,uleg3,chi_sing,size(h{PBC_pos(Jmaxpos+1,Lcur)},3));
w{iteration,2} = random('unif',-1,1,[w2leg1,w2leg2,w2leg3,w2leg4]);

%alternative top on/off
if at_on == 1
    %random without leg swap
    [uleg1,uleg2,uleg3,uleg4] = deal(size(w{(L/2)-1,2},4),size(w{(L/2)-1,2},4),size(w{(L/2)-1,1},2),size(w{(L/2)-1,1},2));
    u{L/2,1} = random('unif',-1,1,[uleg1,uleg2,uleg3,uleg4]);
else %at_on == 0
    %identity
    [uleg1,uleg2,uleg3,uleg4] = deal(size(w{(L/2)-1,2},4),size(w{(L/2)-1,2},4),size(w{(L/2)-1,1},2),size(w{(L/2)-1,1},2));
    u{L/2,1} = zeros([uleg1,uleg2,uleg3,uleg4]);
    u{L/2,1}(1:min(uleg1,uleg2),1:min(uleg1,uleg2),1:min(uleg3,uleg4),1:min(uleg3,uleg4)) = ncon({eye(min(uleg1,uleg2)),eye(min(uleg3,uleg4))},{[-1,-2],[-3,-4]});
end

end