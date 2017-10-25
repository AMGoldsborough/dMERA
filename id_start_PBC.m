function chi_inc = id_start_PBC(L,chi_w,chi_sing,Jmaxpos_vec,chi_inc)
%chi_inc = id_start_PBC(L,chi_w,chi_sing,Jmaxpos_vec,chi_inc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - id_start_PBC
% create identity start tensors
% crossed legs on u
% 
% Andrew Goldsborough - 24/11/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global w u h VMDH;

for iteration = 1:(L/2)-2
    Jmaxpos = Jmaxpos_vec(iteration);
    
    Lcur = size(h,1);
    
    if size(h{PBC_pos(Jmaxpos-1,Lcur)},1) * size(h{PBC_pos(Jmaxpos-1,Lcur)},3) <= chi_w && size(h{PBC_pos(Jmaxpos+1,Lcur)},1) * size(h{PBC_pos(Jmaxpos+1,Lcur)},3) <= chi_w
        %increase chi. VMDH in u, id in w
        chi_inc(iteration) = 1;
        
        [uleg1,uleg2,uleg3,uleg4] = deal(size(h{Jmaxpos},3),size(h{Jmaxpos},1),size(h{Jmaxpos},1),size(h{Jmaxpos},3));
        u{iteration,1} = zeros([uleg1,uleg2,uleg3,uleg4]);
        u{iteration,1}(1,1:2,1,1:2) = VMDH;
        
        [w1leg1,w1leg2,w1leg3] = deal(min(size(h{PBC_pos(Jmaxpos-1,Lcur)},1)*size(u{iteration,1},1),chi_w),size(h{PBC_pos(Jmaxpos-1,Lcur)},1),size(u{iteration,1},1));
        w{iteration,1} = zeros([w1leg1,w1leg2,w1leg3]);
        w{iteration,1}(1:w1leg2,1:w1leg2,1) = eye(w1leg2);
        
        [w2leg1,w2leg2,w2leg3] = deal(min(size(u{iteration,1},3)*size(h{PBC_pos(Jmaxpos+1,Lcur)},3),chi_w),size(u{iteration,1},3),size(h{PBC_pos(Jmaxpos+1,Lcur)},3));
        w{iteration,2} = zeros([w2leg1,w2leg2,w2leg3]);
        w{iteration,2}(1:w2leg3,1,1:w2leg3) = eye(w2leg3);

        %update Hamiltonian
        raiseHam_PBC(Jmaxpos,iteration,chi_inc(iteration));
    else
        %normal structure
        [uleg1,uleg2,uleg3,uleg4] = deal(size(h{Jmaxpos},3),size(h{Jmaxpos},1),size(h{Jmaxpos},1),size(h{Jmaxpos},3));
        u{iteration,1} = zeros([uleg1,uleg2,uleg3,uleg4]);
        u{iteration,1}(1:min(uleg1,uleg2),1:min(uleg1,uleg2),1:min(uleg3,uleg4),1:min(uleg3,uleg4)) = ncon({eye(min(uleg1,uleg2)),eye(min(uleg3,uleg4))},{[-1,-2],[-3,-4]});
        
        [w1leg1,w1leg2,w1leg3,w1leg4] = deal(min(chi_w,ceil(size(h{PBC_pos(Jmaxpos-1,Lcur)},1)*uleg1/chi_sing)),size(h{PBC_pos(Jmaxpos-1,Lcur)},1),chi_sing,uleg1);
        w{iteration,1} = zeros([w1leg1,w1leg2,w1leg3,w1leg4]);
        w{iteration,1}(1:min(w1leg1,w1leg2),1:min(w1leg1,w1leg2),1:min(w1leg3,w1leg4),1:min(w1leg3,w1leg4)) = ncon({eye(min(w1leg1,w1leg2)),eye(min(w1leg3,w1leg4))},{[-1,-2],[-3,-4]});

        [w2leg1,w2leg2,w2leg3,w2leg4] = deal(chi_sing,uleg3,min(chi_w,ceil(uleg3*size(h{PBC_pos(Jmaxpos+1,Lcur)},3)/chi_sing)),size(h{PBC_pos(Jmaxpos+1,Lcur)},3));
        w{iteration,2} = zeros([w2leg1,w2leg2,w2leg3,w2leg4]);
        w{iteration,2}(1:min(w2leg1,w2leg2),1:min(w2leg1,w2leg2),1:min(w2leg3,w2leg4),1:min(w2leg3,w2leg4)) = ncon({eye(min(w2leg1,w2leg2)),eye(min(w2leg3,w2leg4))},{[-1,-2],[-3,-4]});

        %update Hamiltonian
        raiseHam_PBC(Jmaxpos,iteration,chi_inc(iteration));
    end
end

%top set
iteration = (L/2)-1;
Jmaxpos = Jmaxpos_vec(iteration);
Lcur = size(h,1);

%normal structure
[uleg1,uleg2,uleg3,uleg4] = deal(max(size(h{Jmaxpos},1),size(h{Jmaxpos},3)),size(h{Jmaxpos},1),max(size(h{Jmaxpos},1),size(h{Jmaxpos},3)),size(h{Jmaxpos},3));
u{iteration,1} = zeros([uleg1,uleg2,uleg3,uleg4]);
u{iteration,1}(1:min(uleg1,uleg2),1:min(uleg1,uleg2),1:min(uleg3,uleg4),1:min(uleg3,uleg4)) = ncon({eye(min(uleg1,uleg2)),eye(min(uleg3,uleg4))},{[-1,-2],[-3,-4]});

[w1leg1,w1leg2,w1leg3,w1leg4] = deal(chi_sing,size(h{PBC_pos(Jmaxpos-1,Lcur)},1),chi_sing,uleg1);
w{iteration,1} = zeros([w1leg1,w1leg2,w1leg3,w1leg4]);
w{iteration,1}(1:min(w1leg1,w1leg2),1:min(w1leg1,w1leg2),1:min(w1leg3,w1leg4),1:min(w1leg3,w1leg4)) = ncon({eye(min(w1leg1,w1leg2)),eye(min(w1leg3,w1leg4))},{[-1,-2],[-3,-4]});

[w2leg1,w2leg2,w2leg3,w2leg4] = deal(chi_sing,uleg3,chi_sing,size(h{PBC_pos(Jmaxpos+1,Lcur)},3));
w{iteration,2} = zeros([w2leg1,w2leg2,w2leg3,w2leg4]);
w{iteration,2}(1:min(w2leg1,w2leg2),1:min(w2leg1,w2leg2),1:min(w2leg3,w2leg4),1:min(w2leg3,w2leg4)) = ncon({eye(min(w2leg1,w2leg2)),eye(min(w2leg3,w2leg4))},{[-1,-2],[-3,-4]});

%alternative top is standard identity
[uleg1,uleg2,uleg3,uleg4] = deal(size(w{(L/2)-1,2},4),size(w{(L/2)-1,2},4),size(w{(L/2)-1,1},2),size(w{(L/2)-1,1},2));
u{L/2,1} = zeros([uleg1,uleg2,uleg3,uleg4]);
u{L/2,1}(1:min(uleg1,uleg2),1:min(uleg1,uleg2),1:min(uleg3,uleg4),1:min(uleg3,uleg4)) = ncon({eye(min(uleg1,uleg2)),eye(min(uleg3,uleg4))},{[-1,-2],[-3,-4]});

end