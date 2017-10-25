function update_PBC(Jmaxpos,iteration,optmax,ULmax,slow,epsilon,sweep,sweepmax,chi_incs)
%update_PBC(Jmaxpos,iteration,optmax,ULmax,slow,epsilon,sweep,sweepmax,chi_incs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - update_PBC
% performs update of coarse graining block
% 
% Andrew Goldsborough - 29/09/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global w u rho h VMDH;
global envu envw1 envw2;

%store current length
Lcur = size(h,1);

if chi_incs == 0
    %normal update
    
    %contraction orders from netcon
    envu_2_order = {[6,12,5,4],[12,14,13,-2],[13,1],[1,-4,4,3],[14,10,-1,11],[8,11,9,-3],[6,10,7,8],[7,2],[2,9,5,3]};
    envu_3_order = {[6,8,5,4],[8,10,9,-2],[9,1],[1,-4,4,3],[-1,13,-3,14],[11,13,12,14],[6,10,7,11],[7,2],[2,12,5,3]};
    envu_4_order = {[8,12,9,10],[12,14,13,-2],[13,1],[1,-4,10,11],[-3,6,11,7],[4,-1,5,6],[8,14,3,4],[3,2],[2,5,9,7]};
    envw1_1_order = {[3,2,5,-2],[2,3,-1,4],[-4,1],[5,4,6,-3],[6,1]};
    envw1_2_order = {[13,-2,14,12],[-4,1],[1,10,12,11],[-3,8,10,9],[-1,6,8,7],[4,7,5,9],[13,6,3,4],[3,2],[2,5,14,11]};
    envw1_3_order = {[13,-2,14,12],[-4,1],[1,8,12,9],[-3,5,8,6],[5,3,6,4],[10,3,7,4],[13,-1,11,10],[11,2],[2,7,14,9]};
    envw1_4_order = {[13,-2,14,12],[-4,1],[1,8,12,9],[-3,6,8,7],[7,4,9,5],[10,6,3,4],[13,-1,11,10],[11,2],[2,3,14,5]};
    envw2_2_order = {[13,12,14,-4],[12,8,1,9],[1,-2],[9,6,-1,7],[8,4,6,5],[3,5,10,7],[13,4,2,3],[2,11],[11,10,14,-3]};
    envw2_3_order = {[13,12,14,-4],[12,8,1,9],[1,-2],[9,5,-1,6],[5,3,6,4],[7,3,10,4],[13,8,2,7],[2,11],[11,10,14,-3]};
    envw2_4_order = {[13,12,14,-4],[12,10,1,11],[1,-2],[11,8,-1,9],[9,6,-3,7],[4,8,5,6],[13,10,3,4],[3,2],[2,5,14,7]};
    envw2_5_order = {[1,-2],[5,-4,3,2],[-3,4,2,3],[1,6],[6,-1,5,4]};
    
    %useful values
    [uleg1,uleg2,uleg3,uleg4] = size(u{iteration,1});
    [w1leg1,w1leg2,w1leg3,w1leg4] = size(w{iteration,1});
    [w2leg1,w2leg2,w2leg3,w2leg4] = size(w{iteration,2});
    
    if slow == 3 || slow == 4
        %shifts
        if sweep < sweepmax/10
            shift_frac = 0;
        else
            shift_frac = 1;
        end
        h1 = h{PBC_pos(Jmaxpos-2,Lcur)} - ...
            shift_frac*max(eig(0.5*(tfuse(tfuse(h{PBC_pos(Jmaxpos-2,Lcur)},[1,-2,1,-3]),[-1,2,2])...
            +tfuse(tfuse(h{PBC_pos(Jmaxpos-2,Lcur)},[1,-2,1,-3]),[-1,2,2])')))...
            *ncon({eye(size(h{PBC_pos(Jmaxpos-2,Lcur)},1)),eye(size(h{PBC_pos(Jmaxpos-2,Lcur)},3))},{[-1,-2],[-3,-4]});
        h2 = h{PBC_pos(Jmaxpos-1,Lcur)} - ...
            shift_frac*max(eig(0.5*(tfuse(tfuse(h{PBC_pos(Jmaxpos-1,Lcur)},[1,-2,1,-3]),[-1,2,2])...
            +tfuse(tfuse(h{PBC_pos(Jmaxpos-1,Lcur)},[1,-2,1,-3]),[-1,2,2])')))...
            *ncon({eye(size(h{PBC_pos(Jmaxpos-1,Lcur)},1)),eye(size(h{PBC_pos(Jmaxpos-1,Lcur)},3))},{[-1,-2],[-3,-4]});
        h3 = h{PBC_pos(Jmaxpos,Lcur)} - ...
            shift_frac*max(eig(0.5*(tfuse(tfuse(h{PBC_pos(Jmaxpos,Lcur)},[1,-2,1,-3]),[-1,2,2])...
            +tfuse(tfuse(h{PBC_pos(Jmaxpos,Lcur)},[1,-2,1,-3]),[-1,2,2])')))...
            *ncon({eye(size(h{PBC_pos(Jmaxpos,Lcur)},1)),eye(size(h{PBC_pos(Jmaxpos,Lcur)},3))},{[-1,-2],[-3,-4]});
        h4 = h{PBC_pos(Jmaxpos+1,Lcur)} - ...
            shift_frac*max(eig(0.5*(tfuse(tfuse(h{PBC_pos(Jmaxpos+1,Lcur)},[1,-2,1,-3]),[-1,2,2])...
            +tfuse(tfuse(h{PBC_pos(Jmaxpos+1,Lcur)},[1,-2,1,-3]),[-1,2,2])')))...
            *ncon({eye(size(h{PBC_pos(Jmaxpos+1,Lcur)},1)),eye(size(h{PBC_pos(Jmaxpos+1,Lcur)},3))},{[-1,-2],[-3,-4]});
        h5 = h{PBC_pos(Jmaxpos+2,Lcur)} - ...
            shift_frac*max(eig(0.5*(tfuse(tfuse(h{PBC_pos(Jmaxpos+2,Lcur)},[1,-2,1,-3]),[-1,2,2])...
            +tfuse(tfuse(h{PBC_pos(Jmaxpos+2,Lcur)},[1,-2,1,-3]),[-1,2,2])')))...
            *ncon({eye(size(h{PBC_pos(Jmaxpos+2,Lcur)},1)),eye(size(h{PBC_pos(Jmaxpos+2,Lcur)},3))},{[-1,-2],[-3,-4]});
    else
        %no shift
        h1 = h{PBC_pos(Jmaxpos-2,Lcur)};
        h2 = h{PBC_pos(Jmaxpos-1,Lcur)};
        h3 = h{PBC_pos(Jmaxpos,Lcur)};
        h4 = h{PBC_pos(Jmaxpos+1,Lcur)};
        h5 = h{PBC_pos(Jmaxpos+2,Lcur)};
    end
    
    for optcount = 1:optmax
        %u
        for ULcount = 1:ULmax
            
            %diagram 1 - left most
            envu_1 = 0;
            
            %diagram 2 - left of centre
            envu_2 = ncon({rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},w{iteration,1},VMDH,w{iteration,2},h2,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envu_2_order);
            
            %diagram 3 - centre
            envu_3 = ncon({rho{iteration,PBC_up(Jmaxpos,Jmaxpos,Lcur)},w{iteration,1},VMDH,w{iteration,2},h3,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envu_3_order);
            
            %diagram 4 - right of centre
            envu_4 = ncon({rho{iteration,PBC_up(Jmaxpos+1,Jmaxpos,Lcur)},w{iteration,1},VMDH,w{iteration,2},h4,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envu_4_order);
            
            %diagram 5 - right most
            envu_5 = 0;
            
            %optimise u by UL
            if slow == 1  && isempty(envu{iteration}) ~= 1
                if sweep <= sweepmax/2
                    %slower
                    envu{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * envu{iteration}...
                        + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envu_1 + envu_2 + envu_3 + envu_4 + envu_5);
                else
                    envu{iteration} = envu_1 + envu_2 + envu_3 + envu_4 + envu_5;
                end
            elseif (slow == 2 || slow == 4) && isempty(envu{iteration}) ~= 1
                if sweep <= sweepmax/2
                    %random perturbation
                    envu{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * random('unif',-1,1,size(envu_2))...
                        + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envu_1 + envu_2 + envu_3 + envu_4 + envu_5);
                else
                    envu{iteration} = envu_1 + envu_2 + envu_3 + envu_4 + envu_5;
                end
            else
                envu{iteration} = envu_1 + envu_2 + envu_3 + envu_4 + envu_5;
            end
            
            if uleg1 == 1
                envum = tfuse(envu{iteration},[2,-1,2]);
            else
                envum = tfuse(envu{iteration},[1,-2,1,-3]);
                envum = tfuse(envum,[-1,2,2]);
            end
            
            [U,~,V] = svd(envum,'econ');
            u{iteration,1} = -V*U';
            
            %split into a 4-index tensor
            if uleg1 == 1
                u{iteration,1} = tsplit(u{iteration,1},1,[uleg2,uleg4]);
                u{iteration,1} = permute(u{iteration,1},[4,1,3,2]);
            else
                u{iteration,1} = tsplit(u{iteration,1},1,[uleg1,uleg3]);
                u{iteration,1} = tsplit(u{iteration,1},3,[uleg2,uleg4]);
                u{iteration,1} = permute(u{iteration,1},[1,3,2,4]);
            end
        end
        
        %w1
        for ULcount = 1:ULmax
            
            %diagram 1 - left most
            envw1_1 = ncon({rho{iteration,PBC_up(Jmaxpos-2,Jmaxpos,Lcur)},h1,...
                VMDH,conj(w{iteration,1}),conj(VMDH)},envw1_1_order);
            
            %diagram 2 - left of centre
            envw1_2 = ncon({rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},VMDH,w{iteration,2},u{iteration,1},h2,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envw1_2_order);
            
            %diagram 3 - centre
            envw1_3 = ncon({rho{iteration,PBC_up(Jmaxpos,Jmaxpos,Lcur)},VMDH,w{iteration,2},u{iteration,1},h3,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envw1_3_order);
            
            %diagram 4 - right of centre
            envw1_4 = ncon({rho{iteration,PBC_up(Jmaxpos+1,Jmaxpos,Lcur)},VMDH,w{iteration,2},u{iteration,1},h4,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envw1_4_order);
            
            %diagram 5 - right most
            envw1_5 = 0;
                
            %optimise w1 by SVD
            if slow == 1 && isempty(envw1{iteration}) ~= 1
                if sweep <= sweepmax/2
                    %slower
                    envw1{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * envw1{iteration}...
                        + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5);
                else
                    envw1{iteration} = envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5;
                end
            elseif (slow == 2 || slow == 4) && isempty(envw1{iteration}) ~= 1
                if sweep <= sweepmax/2
                    %random perturbation
                    envw1{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * random('unif',-1,1,size(envw1_1))...
                        + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5);
                else
                    envw1{iteration} = envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5;
                end
            else
                envw1{iteration} = envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5;
            end
            
            if w1leg1 == 1
                envw1m = tfuse(envw1{iteration},[2,-1,2]);
            else
                envw1m = tfuse(envw1{iteration},[1,-2,1,-3]);
                envw1m = tfuse(envw1m,[-1,2,2]);
            end
            
            [U,~,V] = svd(envw1m,'econ');
            w{iteration,1} = -V*U';
            
            %split into a 4-index tensor
            if w1leg1 == 1
                w{iteration,1} = tsplit(w{iteration,1},2,[w1leg2,w1leg4]);
                w{iteration,1} = permute(w{iteration,1},[4,2,1,3]);
            else
                w{iteration,1} = tsplit(w{iteration,1},1,[w1leg1,w1leg3]);
                w{iteration,1} = tsplit(w{iteration,1},3,[w1leg2,w1leg4]);
                w{iteration,1} = permute(w{iteration,1},[1,3,2,4]);
            end
        end
        
        %w2
        for ULcount = 1:ULmax
            
            %diagram 1 - left most
            envw2_1 = 0;
            
            %diagram 2 - left of centre
            envw2_2 = ncon({rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},w{iteration,1},VMDH,u{iteration,1},h2,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envw2_2_order);
            
            %diagram 3 - centre
            envw2_3 = ncon({rho{iteration,PBC_up(Jmaxpos,Jmaxpos,Lcur)},w{iteration,1},VMDH,u{iteration,1},h3,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envw2_3_order);
            
            %diagram 4 - right of centre
            envw2_4 = ncon({rho{iteration,PBC_up(Jmaxpos+1,Jmaxpos,Lcur)},w{iteration,1},VMDH,u{iteration,1},h4,...
                conj(u{iteration,1}),conj(w{iteration,1}),conj(VMDH),conj(w{iteration,2})},envw2_4_order);
            
            %diagram 5 - right most
            envw2_5 = ncon({VMDH,rho{iteration,PBC_up(Jmaxpos+2,Jmaxpos,Lcur)},h5,...
                conj(VMDH),conj(w{iteration,2})},envw2_5_order);
            
            %optimise w2 by SVD
            if slow == 1 && isempty(envw2{iteration}) ~= 1
                if sweep <= sweepmax/2
                    %slower
                    envw2{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * envw2{iteration}...
                        + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5);
                else
                    envw2{iteration} = envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5;
                end
            elseif (slow == 2 || slow == 4) && isempty(envw2{iteration}) ~= 1
                if sweep <= sweepmax/2
                    %random perturbation
                    envw2{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * random('unif',-1,1,size(envw2_2))...
                        + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5);
                else
                    envw2{iteration} = envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5;
                end
            else
                envw2{iteration} = envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5;
            end
            
            if w2leg1 == 3
                envw2m = tfuse(envw2{iteration},[-1,2,2]);
            else
                envw2m = tfuse(envw2{iteration},[1,-2,1,-3]);
                envw2m = tfuse(envw2m,[-1,2,2]);
            end
            
            [U,~,V] = svd(envw2m,'econ');
            w{iteration,2} = -V*U';
            
            %split into a 4-index tensor
            if w2leg1 == 3
                w{iteration,2} = tsplit(w{iteration,2},2,[w2leg2,w2leg4]);
                w{iteration,2} = permute(w{iteration,2},[1,2,4,3]);
            else
                w{iteration,2} = tsplit(w{iteration,2},1,[w2leg1,w2leg3]);
                w{iteration,2} = tsplit(w{iteration,2},3,[w2leg2,w2leg4]);
                w{iteration,2} = permute(w{iteration,2},[1,3,2,4]);
            end
        end
    end
else %chi_incs == 1
    %chi increases
    if sweep == 2
        %set max ul/block sweeps
        optmax_inc = 2;
        ULmax_inc = 2;
        criterion_inc = 1e-20;
        
        %contraction orders from netcon
        envu_2_inc_order = {[9,7,-1,8],[5,8,6,-3],[4,7,5],[2,6,1],[4,10,2,3],[10,9,-2],[3,-4,1]};
        envu_3_inc_order = {[-1,9,-3,10],[7,9,8,10],[5,1,7],[3,8,2],[5,6,3,4],[6,1,-2],[4,-4,2]};
        envu_4_inc_order = {[-3,7,9,8],[5,-1,6,7],[2,1,5],[4,6,8],[2,3,4,10],[3,1,-2],[10,-4,9]};
        envw1_1_inc_order = {[1,2,4,-3],[2,1,-1,3],[4,3,-2]};
        envw1_2_inc_order = {[-2,9,10,11],[-1,7,9,8],[5,8,6,11],[4,7,5],[2,6,1],[4,-3,2,3],[3,10,1]};
        envw1_3_inc_order = {[-2,4,6,5],[4,2,5,3],[10,2,7,3],[11,-1,10],[8,7,1],[11,-3,8,9],[9,6,1]};
        envw1_4_inc_order = {[-2,4,2,3],[3,5,1,6],[10,4,7,5],[11,-1,10],[8,7,6],[11,-3,8,9],[9,2,1]};
        envw2_2_inc_order = {[2,3,-1,4],[1,6,3,5],[7,5,10,4],[8,6,7],[11,10,-2],[8,9,11,-3],[9,1,2]};
        envw2_3_inc_order = {[6,4,-1,5],[4,2,5,3],[7,2,10,3],[8,1,7],[11,10,-2],[8,9,11,-3],[9,1,6]};
        envw2_4_inc_order = {[9,10,-1,11],[11,7,-2,8],[5,10,6,7],[2,1,5],[4,6,8],[2,3,4,-3],[3,1,9]};
        envw2_5_inc_order = {[3,-3,2,1],[-2,4,1,2],[3,-1,4]};
        
        %useful values
        [~,w1leg2,w1leg3] = size(w{iteration,1});
        [~,w2leg2,w2leg3] = size(w{iteration,2});
        [uleg1,uleg2,uleg3,uleg4] = size(u{iteration,1});
        
        %no shift
        h1 = h{PBC_pos(Jmaxpos-2,Lcur)};
        h2 = h{PBC_pos(Jmaxpos-1,Lcur)};
        h3 = h{PBC_pos(Jmaxpos,Lcur)};
        h4 = h{PBC_pos(Jmaxpos+1,Lcur)};
        h5 = h{PBC_pos(Jmaxpos+2,Lcur)};
        
        %id rho to locally diagonalise
        rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)} = ncon({eye(size(rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},1)),...
            eye(size(rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},3))},{[-1,-2],[-3,-4]});
        rho{iteration,PBC_up(Jmaxpos-2,Jmaxpos,Lcur)} = ncon({eye(size(rho{iteration,PBC_up(Jmaxpos-2,Jmaxpos,Lcur)},1)),...
            eye(size(rho{iteration,PBC_up(Jmaxpos-2,Jmaxpos,Lcur)},3))},{[-1,-2],[-3,-4]});
        rho{iteration,PBC_up(Jmaxpos+2,Jmaxpos,Lcur)} = ncon({eye(size(rho{iteration,PBC_up(Jmaxpos+2,Jmaxpos,Lcur)},1)),...
            eye(size(rho{iteration,PBC_up(Jmaxpos+2,Jmaxpos,Lcur)},3))},{[-1,-2],[-3,-4]});
        
        for optcount = 1:optmax_inc
            
            %u
            for ULcount = 1:ULmax_inc
                
                %store u to check convergence
                old = u{iteration,1};
                
                %diagram 1 - left most
                envu_1 = 0;
                
                %diagram 2 - left of centre
                envu_2 = ncon({h2,conj(u{iteration,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
                    rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},w{iteration,1},w{iteration,2}},envu_2_inc_order);
                
                %diagram 3 - centre
                envu_3 = ncon({h3,conj(u{iteration,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
                    rho{iteration,PBC_up(Jmaxpos,Jmaxpos,Lcur)},w{iteration,1},w{iteration,2}},envu_3_inc_order);
                
                %diagram 4 - right of centre
                envu_4 = ncon({h4,conj(u{iteration,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
                    rho{iteration,PBC_up(Jmaxpos+1,Jmaxpos,Lcur)},w{iteration,1},w{iteration,2}},envu_4_inc_order);
                
                %diagram 5 - right most
                envu_5 = 0;
                
                %optimise u by UL
                envu{iteration} = envu_1 + envu_2 + envu_3 + envu_4 + envu_5;
                
                if uleg1 == 1
                    envum = tfuse(envu{iteration},[2,-1,2]);
                else
                    envum = tfuse(envu{iteration},[1,-2,1,-3]);
                    envum = tfuse(envum,[-1,2,2]);
                end
                
                [U,~,V] = svd(envum,'econ');
                u{iteration,1} = -V*U';
                
                %split into a 4-index tensor
                if uleg1 == 1
                    u{iteration,1} = tsplit(u{iteration,1},1,[uleg2,uleg4]);
                    u{iteration,1} = permute(u{iteration,1},[4,1,3,2]);
                else
                    u{iteration,1} = tsplit(u{iteration,1},1,[uleg1,uleg3]);
                    u{iteration,1} = tsplit(u{iteration,1},3,[uleg2,uleg4]);
                    u{iteration,1} = permute(u{iteration,1},[1,3,2,4]);
                end
                
                %check convergence
                change = ncon({u{iteration,1},conj(u{iteration,1})},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({u{iteration,1},conj(old)},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({old,conj(u{iteration,1})},{[1,2,3,4],[1,2,3,4]})...
                    + ncon({old,conj(old)},{[1,2,3,4],[1,2,3,4]});
                
%                 fprintf('iteration: %d, change in u = %.5e.\n',iteration,change);

                if change < criterion_inc
                    break
                end

            end
            
            %w1
            for ULcount = 1:ULmax_inc
                
                %store w1 to check convergence
                old = w{iteration,1};
                
                %diagram 1 - left most
                envw1_1 = ncon({rho{iteration,PBC_up(Jmaxpos-2,Jmaxpos,Lcur)},h1,...
                    conj(w{iteration,1})},envw1_1_inc_order);
                
                %diagram 2 - left of centre
                envw1_2 = ncon({u{iteration,1},h2,conj(u{iteration,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},w{iteration,2}},envw1_2_inc_order);
                
                %diagram 3 - centre
                envw1_3 = ncon({u{iteration,1},h3,conj(u{iteration,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),rho{iteration,PBC_up(Jmaxpos,Jmaxpos,Lcur)},w{iteration,2}},envw1_3_inc_order);
                
                %diagram 4 - right of centre
                envw1_4 = ncon({u{iteration,1},h4,conj(u{iteration,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),rho{iteration,PBC_up(Jmaxpos+1,Jmaxpos,Lcur)},w{iteration,2}},envw1_4_inc_order);
                
                %diagram 5 - right most
                envw1_5 = 0;
                
                %optimise w1 by SVD
                envw1{iteration} = envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5;
                envw1m = tfuse(envw1{iteration},[1,1,-2]);
                
                [U,~,V] = svd(envw1m,'econ');
                w{iteration,1} = -V*U';
                w{iteration,1} = tsplit(w{iteration,1},2,[w1leg2,w1leg3]);
                
                %check convergence
                change = ncon({w{iteration,1},conj(w{iteration,1})},{[1,2,3],[1,2,3]})...
                    - ncon({w{iteration,1},conj(old)},{[1,2,3],[1,2,3]})...
                    - ncon({old,conj(w{iteration,1})},{[1,2,3],[1,2,3]})...
                    + ncon({old,conj(old)},{[1,2,3],[1,2,3]});
                
%                 fprintf('iteration: %d, change in w1 = %.5e.\n',iteration,change);

                if change < criterion_inc
                    break
                end
            end
            
            %w2
            for ULcount = 1:ULmax_inc
                
                %store w2 to check convergence
                old = w{iteration,2};
                
                %diagram 1 - left most
                envw2_1 = 0;
                
                %diagram 2 - left of centre
                envw2_2 = ncon({u{iteration,1},h2,conj(u{iteration,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),rho{iteration,PBC_up(Jmaxpos-1,Jmaxpos,Lcur)},w{iteration,1}},envw2_2_inc_order);
                
                %diagram 3 - centre
                envw2_3 = ncon({u{iteration,1},h3,conj(u{iteration,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),rho{iteration,PBC_up(Jmaxpos,Jmaxpos,Lcur)},w{iteration,1}},envw2_3_inc_order);
                
                %diagram 4 - right of centre
                envw2_4 = ncon({u{iteration,1},h4,conj(u{iteration,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),rho{iteration,PBC_up(Jmaxpos+1,Jmaxpos,Lcur)},w{iteration,1}},envw2_4_inc_order);
                
                %diagram 5 - right most
                envw2_5 = ncon({rho{iteration,PBC_up(Jmaxpos+2,Jmaxpos,Lcur)},h5,...
                    conj(w{iteration,2})},envw2_5_inc_order);
                
                %optimise w2 by SVD
                envw2{iteration} = envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5;
                envw2m = tfuse(envw2{iteration},[1,1,-2]);
                
                [U,~,V] = svd(envw2m,'econ');
                w{iteration,2} = -V*U';
                w{iteration,2} = tsplit(w{iteration,2},2,[w2leg2,w2leg3]);
                
                %check convergence
                change = ncon({w{iteration,2},conj(w{iteration,2})},{[1,2,3],[1,2,3]})...
                    - ncon({w{iteration,2},conj(old)},{[1,2,3],[1,2,3]})...
                    - ncon({old,conj(w{iteration,2})},{[1,2,3],[1,2,3]})...
                    + ncon({old,conj(old)},{[1,2,3],[1,2,3]});
                
%                 fprintf('iteration: %d, UL = %d, opt = %d, change in w2 = %.5e.\n',iteration,ULcount,optcount,change);
                
                 if change < criterion_inc
                    break
                end
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find optimal contraction orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% envu_2_netcon(netcon({[11,1,14,4],[1,5,2,-2],[2,3],[3,-4,4,6],[5,7,-1,8],[9,8,10,-3],[11,7,12,9],[12,13],[13,10,14,6]},0,2,1,1)) = 1:14;
% envu_2_order = {[envu_2_netcon(11),envu_2_netcon(1),envu_2_netcon(14),envu_2_netcon(4)],...
%     [envu_2_netcon(1),envu_2_netcon(5),envu_2_netcon(2),-2],...
%     [envu_2_netcon(2),envu_2_netcon(3)],...
%     [envu_2_netcon(3),-4,envu_2_netcon(4),envu_2_netcon(6)],...
%     [envu_2_netcon(5),envu_2_netcon(7),-1,envu_2_netcon(8)],...
%     [envu_2_netcon(9),envu_2_netcon(8),envu_2_netcon(10),-3],...
%     [envu_2_netcon(11),envu_2_netcon(7),envu_2_netcon(12),envu_2_netcon(9)],...
%     [envu_2_netcon(12),envu_2_netcon(13)],...
%     [envu_2_netcon(13),envu_2_netcon(10),envu_2_netcon(14),envu_2_netcon(6)]};
% 
% envu_3_netcon(netcon({[11,1,14,4],[1,5,2,-2],[2,3],[3,-4,4,6],[-1,7,-3,8],[9,7,10,8],[11,5,12,9],[12,13],[13,10,14,6]},0,2,1,1)) = 1:14;
% envu_3_order = {[envu_3_netcon(11),envu_3_netcon(1),envu_3_netcon(14),envu_3_netcon(4)],...
%     [envu_3_netcon(1),envu_3_netcon(5),envu_3_netcon(2),-2],...
%     [envu_3_netcon(2),envu_3_netcon(3)],...
%     [envu_3_netcon(3),-4,envu_3_netcon(4),envu_3_netcon(6)],...
%     [-1,envu_3_netcon(7),-3,envu_3_netcon(8)],...
%     [envu_3_netcon(9),envu_3_netcon(7),envu_3_netcon(10),envu_3_netcon(8)],...
%     [envu_3_netcon(11),envu_3_netcon(5),envu_3_netcon(12),envu_3_netcon(9)],...
%     [envu_3_netcon(12),envu_3_netcon(13)],...
%     [envu_3_netcon(13),envu_3_netcon(10),envu_3_netcon(14),envu_3_netcon(6)]};
% 
% envu_4_netcon(netcon({[11,1,14,4],[1,5,2,-2],[2,3],[3,-4,4,6],[-3,7,6,8],[9,-1,10,7],[11,5,12,9],[12,13],[13,10,14,8]},0,2,1,1)) = 1:14;
% envu_4_order = {[envu_4_netcon(11),envu_4_netcon(1),envu_4_netcon(14),envu_4_netcon(4)],...
%     [envu_4_netcon(1),envu_4_netcon(5),envu_4_netcon(2),-2],...
%     [envu_4_netcon(2),envu_4_netcon(3)],...
%     [envu_4_netcon(3),-4,envu_4_netcon(4),envu_4_netcon(6)],...
%     [-3,envu_4_netcon(7),envu_4_netcon(6),envu_4_netcon(8)],...
%     [envu_4_netcon(9),-1,envu_4_netcon(10),envu_4_netcon(7)],...
%     [envu_4_netcon(11),envu_4_netcon(5),envu_4_netcon(12),envu_4_netcon(9)],...
%     [envu_4_netcon(12),envu_4_netcon(13)],...
%     [envu_4_netcon(13),envu_4_netcon(10),envu_4_netcon(14),envu_4_netcon(8)]};
% 
% envw1_1_netcon(netcon({[3,1,5,-2],[1,3,-1,4],[-4,2],[5,4,6,-3],[6,2]},0,2,1,1)) = 1:6;
% envw1_1_order = {[envw1_1_netcon(3),envw1_1_netcon(1),envw1_1_netcon(5),-2],...
%     [envw1_1_netcon(1),envw1_1_netcon(3),-1,envw1_1_netcon(4)],...
%     [-4,envw1_1_netcon(2)],...
%     [envw1_1_netcon(5),envw1_1_netcon(4),envw1_1_netcon(6),-3],...
%     [envw1_1_netcon(6),envw1_1_netcon(2)]};
% 
% envw1_2_netcon(netcon({[11,-2,14,2],[-4,1],[1,3,2,4],[-3,5,3,6],[-1,7,5,8],[9,8,10,6],[11,7,12,9],[12,13],[13,10,14,4]},0,2,1,1)) = 1:14;
% envw1_2_order = {[envw1_2_netcon(11),-2,envw1_2_netcon(14),envw1_2_netcon(2)],...
%     [-4,envw1_2_netcon(1)],...
%     [envw1_2_netcon(1),envw1_2_netcon(3),envw1_2_netcon(2),envw1_2_netcon(4)],...
%     [-3,envw1_2_netcon(5),envw1_2_netcon(3),envw1_2_netcon(6)],...
%     [-1,envw1_2_netcon(7),envw1_2_netcon(5),envw1_2_netcon(8)],...
%     [envw1_2_netcon(9),envw1_2_netcon(8),envw1_2_netcon(10),envw1_2_netcon(6)],...
%     [envw1_2_netcon(11),envw1_2_netcon(7),envw1_2_netcon(12),envw1_2_netcon(9)],...
%     [envw1_2_netcon(12),envw1_2_netcon(13)],...
%     [envw1_2_netcon(13),envw1_2_netcon(10),envw1_2_netcon(14),envw1_2_netcon(4)]};
% 
% envw1_3_netcon(netcon({[11,-2,14,2],[-4,1],[1,3,2,4],[-3,5,3,6],[5,7,6,8],[9,7,10,8],[11,-1,12,9],[12,13],[13,10,14,4]},0,2,1,1)) = 1:14;
% envw1_3_order = {[envw1_3_netcon(11),-2,envw1_3_netcon(14),envw1_3_netcon(2)],...
%     [-4,envw1_3_netcon(1)],...
%     [envw1_3_netcon(1),envw1_3_netcon(3),envw1_3_netcon(2),envw1_3_netcon(4)],...
%     [-3,envw1_3_netcon(5),envw1_3_netcon(3),envw1_3_netcon(6)],...
%     [envw1_3_netcon(5),envw1_3_netcon(7),envw1_3_netcon(6),envw1_3_netcon(8)],...
%     [envw1_3_netcon(9),envw1_3_netcon(7),envw1_3_netcon(10),envw1_3_netcon(8)],...
%     [envw1_3_netcon(11),-1,envw1_3_netcon(12),envw1_3_netcon(9)],...
%     [envw1_3_netcon(12),envw1_3_netcon(13)],...
%     [envw1_3_netcon(13),envw1_3_netcon(10),envw1_3_netcon(14),envw1_3_netcon(4)]};
% 
% envw1_4_netcon(netcon({[11,-2,14,2],[-4,1],[1,3,2,4],[-3,5,3,6],[6,7,4,8],[9,5,10,7],[11,-1,12,9],[12,13],[13,10,14,8]},0,2,1,1)) = 1:14;
% envw1_4_order = {[envw1_4_netcon(11),-2,envw1_4_netcon(14),envw1_4_netcon(2)],...
%     [-4,envw1_4_netcon(1)],...
%     [envw1_4_netcon(1),envw1_4_netcon(3),envw1_4_netcon(2),envw1_4_netcon(4)],...
%     [-3,envw1_4_netcon(5),envw1_4_netcon(3),envw1_4_netcon(6)],...
%     [envw1_4_netcon(6),envw1_4_netcon(7),envw1_4_netcon(4),envw1_4_netcon(8)],...
%     [envw1_4_netcon(9),envw1_4_netcon(5),envw1_4_netcon(10),envw1_4_netcon(7)],...
%     [envw1_4_netcon(11),-1,envw1_4_netcon(12),envw1_4_netcon(9)],...
%     [envw1_4_netcon(12),envw1_4_netcon(13)],...
%     [envw1_4_netcon(13),envw1_4_netcon(10),envw1_4_netcon(14),envw1_4_netcon(8)]};
% 
% envw2_2_netcon(netcon({[11,1,14,-4],[1,3,2,4],[2,-2],[4,5,-1,6],[3,7,5,8],[9,8,10,6],[11,7,12,9],[12,13],[13,10,14,-3]},0,2,1,1)) = 1:14;
% envw2_2_order = {[envw2_2_netcon(11),envw2_2_netcon(1),envw2_2_netcon(14),-4],...
%     [envw2_2_netcon(1),envw2_2_netcon(3),envw2_2_netcon(2),envw2_2_netcon(4)],...
%     [envw2_2_netcon(2),-2],...
%     [envw2_2_netcon(4),envw2_2_netcon(5),-1,envw2_2_netcon(6)],...
%     [envw2_2_netcon(3),envw2_2_netcon(7),envw2_2_netcon(5),envw2_2_netcon(8)],...
%     [envw2_2_netcon(9),envw2_2_netcon(8),envw2_2_netcon(10),envw2_2_netcon(6)],...
%     [envw2_2_netcon(11),envw2_2_netcon(7),envw2_2_netcon(12),envw2_2_netcon(9)],...
%     [envw2_2_netcon(12),envw2_2_netcon(13)],...
%     [envw2_2_netcon(13),envw2_2_netcon(10),envw2_2_netcon(14),-3]};
% 
% envw2_3_netcon(netcon({[11,1,14,-4],[1,3,2,4],[2,-2],[4,5,-1,6],[5,7,6,8],[9,7,10,8],[11,3,12,9],[12,13],[13,10,14,-3]},0,2,1,1)) = 1:14;
% envw2_3_order = {[envw2_3_netcon(11),envw2_3_netcon(1),envw2_3_netcon(14),-4],...
%     [envw2_3_netcon(1),envw2_3_netcon(3),envw2_3_netcon(2),envw2_3_netcon(4)],...
%     [envw2_3_netcon(2),-2],...
%     [envw2_3_netcon(4),envw2_3_netcon(5),-1,envw2_3_netcon(6)],...
%     [envw2_3_netcon(5),envw2_3_netcon(7),envw2_3_netcon(6),envw2_3_netcon(8)],...
%     [envw2_3_netcon(9),envw2_3_netcon(7),envw2_3_netcon(10),envw2_3_netcon(8)],...
%     [envw2_3_netcon(11),envw2_3_netcon(3),envw2_3_netcon(12),envw2_3_netcon(9)],...
%     [envw2_3_netcon(12),envw2_3_netcon(13)],...
%     [envw2_3_netcon(13),envw2_3_netcon(10),envw2_3_netcon(14),-3]};
% 
% envw2_4_netcon(netcon({[11,1,14,-4],[1,3,2,4],[2,-2],[4,5,-1,6],[6,7,-3,8],[9,5,10,7],[11,3,12,9],[12,13],[13,10,14,8]},0,2,1,1)) = 1:14;
% envw2_4_order = {[envw2_4_netcon(11),envw2_4_netcon(1),envw2_4_netcon(14),-4],...
%     [envw2_4_netcon(1),envw2_4_netcon(3),envw2_4_netcon(2),envw2_4_netcon(4)],...
%     [envw2_4_netcon(2),-2],...
%     [envw2_4_netcon(4),envw2_4_netcon(5),-1,envw2_4_netcon(6)],...
%     [envw2_4_netcon(6),envw2_4_netcon(7),-3,envw2_4_netcon(8)],...
%     [envw2_4_netcon(9),envw2_4_netcon(5),envw2_4_netcon(10),envw2_4_netcon(7)],...
%     [envw2_4_netcon(11),envw2_4_netcon(3),envw2_4_netcon(12),envw2_4_netcon(9)],...
%     [envw2_4_netcon(12),envw2_4_netcon(13)],...
%     [envw2_4_netcon(13),envw2_4_netcon(10),envw2_4_netcon(14),envw2_4_netcon(8)]};
% 
% envw2_5_netcon(netcon({[1,-2],[6,-4,4,2],[-3,3,2,4],[1,5],[5,-1,6,3]},0,2,1,1)) = 1:6;
% envw2_5_order = {[envw2_5_netcon(1),-2],...
%     [envw2_5_netcon(6),-4,envw2_5_netcon(4),envw2_5_netcon(2)],...
%     [-3,envw2_5_netcon(3),envw2_5_netcon(2),envw2_5_netcon(4)],...
%     [envw2_5_netcon(1),envw2_5_netcon(5)],...
%     [envw2_5_netcon(5),-1,envw2_5_netcon(6),envw2_5_netcon(3)]};
%
% envu_2_inc_netcon(netcon({[1,2,-1,3],[4,3,5,-3],[6,2,4],[7,5,8],[6,9,7,10],[9,1,-2],[10,-4,8]},0,2,1,1)) = 1:10;
% envu_2_inc_order = {[envu_2_inc_netcon(1),envu_2_inc_netcon(2),-1,envu_2_inc_netcon(3)],...
%     [envu_2_inc_netcon(4),envu_2_inc_netcon(3),envu_2_inc_netcon(5),-3],...
%     [envu_2_inc_netcon(6),envu_2_inc_netcon(2),envu_2_inc_netcon(4)],...
%     [envu_2_inc_netcon(7),envu_2_inc_netcon(5),envu_2_inc_netcon(8)],...
%     [envu_2_inc_netcon(6),envu_2_inc_netcon(9),envu_2_inc_netcon(7),envu_2_inc_netcon(10)],...
%     [envu_2_inc_netcon(9),envu_2_inc_netcon(1),-2],...
%     [envu_2_inc_netcon(10),-4,envu_2_inc_netcon(8)]};
% 
% envu_3_inc_netcon(netcon({[-1,1,-3,2],[3,1,4,2],[5,6,3],[7,4,8],[5,9,7,10],[9,6,-2],[10,-4,8]},0,2,1,1)) = 1:10;
% envu_3_inc_order = {[-1,envu_3_inc_netcon(1),-3,envu_3_inc_netcon(2)],...
%     [envu_3_inc_netcon(3),envu_3_inc_netcon(1),envu_3_inc_netcon(4),envu_3_inc_netcon(2)],...
%     [envu_3_inc_netcon(5),envu_3_inc_netcon(6),envu_3_inc_netcon(3)],...
%     [envu_3_inc_netcon(7),envu_3_inc_netcon(4),envu_3_inc_netcon(8)],...
%     [envu_3_inc_netcon(5),envu_3_inc_netcon(9),envu_3_inc_netcon(7),envu_3_inc_netcon(10)],...
%     [envu_3_inc_netcon(9),envu_3_inc_netcon(6),-2],...
%     [envu_3_inc_netcon(10),-4,envu_3_inc_netcon(8)]};
% 
% envu_4_inc_netcon(netcon({[-3,1,2,3],[4,-1,5,1],[7,8,4],[6,5,3],[7,9,6,10],[9,8,-2],[10,-4,2]},0,2,1,1)) = 1:10;
% envu_4_inc_order = {[-3,envu_4_inc_netcon(1),envu_4_inc_netcon(2),envu_4_inc_netcon(3)],...
%     [envu_4_inc_netcon(4),-1,envu_4_inc_netcon(5),envu_4_inc_netcon(1)],...
%     [envu_4_inc_netcon(7),envu_4_inc_netcon(8),envu_4_inc_netcon(4)],...
%     [envu_4_inc_netcon(6),envu_4_inc_netcon(5),envu_4_inc_netcon(3)],...
%     [envu_4_inc_netcon(7),envu_4_inc_netcon(9),envu_4_inc_netcon(6),envu_4_inc_netcon(10)],...
%     [envu_4_inc_netcon(9),envu_4_inc_netcon(8),-2],...
%     [envu_4_inc_netcon(10),-4,envu_4_inc_netcon(2)]};
% 
% envw1_1_inc_netcon(netcon({[1,2,4,-3],[2,1,-1,3],[4,3,-2]},0,2,1,1)) = 1:4;
% envw1_1_inc_order = {[envw1_1_inc_netcon(1),envw1_1_inc_netcon(2),envw1_1_inc_netcon(4),-3],...
%     [envw1_1_inc_netcon(2),envw1_1_inc_netcon(1),-1,envw1_1_inc_netcon(3)],...
%     [envw1_1_inc_netcon(4),envw1_1_inc_netcon(3),-2]};
% 
% envw1_2_inc_netcon(netcon({[-2,1,2,3],[-1,4,1,5],[6,5,7,3],[8,4,6],[9,7,10],[8,-3,9,11],[11,2,10]},0,2,1,1)) = 1:11;
% envw1_2_inc_order = {[-2,envw1_2_inc_netcon(1),envw1_2_inc_netcon(2),envw1_2_inc_netcon(3)],...
%     [-1,envw1_2_inc_netcon(4),envw1_2_inc_netcon(1),envw1_2_inc_netcon(5)],...
%     [envw1_2_inc_netcon(6),envw1_2_inc_netcon(5),envw1_2_inc_netcon(7),envw1_2_inc_netcon(3)],...
%     [envw1_2_inc_netcon(8),envw1_2_inc_netcon(4),envw1_2_inc_netcon(6)],...
%     [envw1_2_inc_netcon(9),envw1_2_inc_netcon(7),envw1_2_inc_netcon(10)],...
%     [envw1_2_inc_netcon(8),-3,envw1_2_inc_netcon(9),envw1_2_inc_netcon(11)],...
%     [envw1_2_inc_netcon(11),envw1_2_inc_netcon(2),envw1_2_inc_netcon(10)]};
% 
% envw1_3_inc_netcon(netcon({[-2,1,2,3],[1,4,3,5],[6,4,7,5],[8,-1,6],[9,7,10],[8,-3,9,11],[11,2,10]},0,2,1,1)) = 1:11;
% envw1_3_inc_order = {[-2,envw1_3_inc_netcon(1),envw1_3_inc_netcon(2),envw1_3_inc_netcon(3)],...
%     [envw1_3_inc_netcon(1),envw1_3_inc_netcon(4),envw1_3_inc_netcon(3),envw1_3_inc_netcon(5)],...
%     [envw1_3_inc_netcon(6),envw1_3_inc_netcon(4),envw1_3_inc_netcon(7),envw1_3_inc_netcon(5)],...
%     [envw1_3_inc_netcon(8),-1,envw1_3_inc_netcon(6)],[envw1_3_inc_netcon(9),envw1_3_inc_netcon(7),envw1_3_inc_netcon(10)],...
%     [envw1_3_inc_netcon(8),-3,envw1_3_inc_netcon(9),envw1_3_inc_netcon(11)],...
%     [envw1_3_inc_netcon(11),envw1_3_inc_netcon(2),envw1_3_inc_netcon(10)]};
% 
% envw1_4_inc_netcon(netcon({[-2,1,2,3],[3,4,5,6],[7,1,8,4],[9,-1,7],[10,8,6],[9,-3,10,11],[11,2,5]},0,2,1,1)) = 1:11;
% envw1_4_inc_order = {[-2,envw1_4_inc_netcon(1),envw1_4_inc_netcon(2),envw1_4_inc_netcon(3)],...
%     [envw1_4_inc_netcon(3),envw1_4_inc_netcon(4),envw1_4_inc_netcon(5),envw1_4_inc_netcon(6)],...
%     [envw1_4_inc_netcon(7),envw1_4_inc_netcon(1),envw1_4_inc_netcon(8),envw1_4_inc_netcon(4)],...
%     [envw1_4_inc_netcon(9),-1,envw1_4_inc_netcon(7)],...
%     [envw1_4_inc_netcon(10),envw1_4_inc_netcon(8),envw1_4_inc_netcon(6)],...
%     [envw1_4_inc_netcon(9),-3,envw1_4_inc_netcon(10),envw1_4_inc_netcon(11)],...
%     [envw1_4_inc_netcon(11),envw1_4_inc_netcon(2),envw1_4_inc_netcon(5)]};
%  
% envw2_2_inc_netcon(netcon({[1,2,-1,3],[4,5,2,6],[7,6,8,3],[9,5,7],[10,8,-2],[9,11,10,-3],[11,4,1]},0,2,1,1)) = 1:11;
% envw2_2_inc_order = {[envw2_2_inc_netcon(1),envw2_2_inc_netcon(2),-1,envw2_2_inc_netcon(3)],...
%     [envw2_2_inc_netcon(4),envw2_2_inc_netcon(5),envw2_2_inc_netcon(2),envw2_2_inc_netcon(6)],...
%     [envw2_2_inc_netcon(7),envw2_2_inc_netcon(6),envw2_2_inc_netcon(8),envw2_2_inc_netcon(3)],...
%     [envw2_2_inc_netcon(9),envw2_2_inc_netcon(5),envw2_2_inc_netcon(7)],...
%     [envw2_2_inc_netcon(10),envw2_2_inc_netcon(8),-2],...
%     [envw2_2_inc_netcon(9),envw2_2_inc_netcon(11),envw2_2_inc_netcon(10),-3],...
%     [envw2_2_inc_netcon(11),envw2_2_inc_netcon(4),envw2_2_inc_netcon(1)]};
% 
% envw2_3_inc_netcon(netcon({[1,2,-1,3],[2,4,3,5],[6,4,7,5],[8,9,6],[10,7,-2],[8,11,10,-3],[11,9,1]},0,2,1,1)) = 1:11;
% envw2_3_inc_order = {[envw2_3_inc_netcon(1),envw2_3_inc_netcon(2),-1,envw2_3_inc_netcon(3)],...
%     [envw2_3_inc_netcon(2),envw2_3_inc_netcon(4),envw2_3_inc_netcon(3),envw2_3_inc_netcon(5)],...
%     [envw2_3_inc_netcon(6),envw2_3_inc_netcon(4),envw2_3_inc_netcon(7),envw2_3_inc_netcon(5)],...
%     [envw2_3_inc_netcon(8),envw2_3_inc_netcon(9),envw2_3_inc_netcon(6)],...
%     [envw2_3_inc_netcon(10),envw2_3_inc_netcon(7),-2],...
%     [envw2_3_inc_netcon(8),envw2_3_inc_netcon(11),envw2_3_inc_netcon(10),-3],...
%     [envw2_3_inc_netcon(11),envw2_3_inc_netcon(9),envw2_3_inc_netcon(1)]};
% 
% envw2_4_inc_netcon(netcon({[1,2,-1,3],[3,4,-2,5],[6,2,7,4],[8,9,6],[10,7,5],[8,11,10,-3],[11,9,1]},0,2,1,1)) = 1:11;
% envw2_4_inc_order = {[envw2_4_inc_netcon(1),envw2_4_inc_netcon(2),-1,envw2_4_inc_netcon(3)],...
%     [envw2_4_inc_netcon(3),envw2_4_inc_netcon(4),-2,envw2_4_inc_netcon(5)],...
%     [envw2_4_inc_netcon(6),envw2_4_inc_netcon(2),envw2_4_inc_netcon(7),envw2_4_inc_netcon(4)],...
%     [envw2_4_inc_netcon(8),envw2_4_inc_netcon(9),envw2_4_inc_netcon(6)],...
%     [envw2_4_inc_netcon(10),envw2_4_inc_netcon(7),envw2_4_inc_netcon(5)],...
%     [envw2_4_inc_netcon(8),envw2_4_inc_netcon(11),envw2_4_inc_netcon(10),-3],...
%     [envw2_4_inc_netcon(11),envw2_4_inc_netcon(9),envw2_4_inc_netcon(1)]};
% 
% envw2_5_inc_netcon(netcon({[1,-3,4,3],[-2,2,3,4],[1,-1,2]},0,2,1,1)) = 1:4;
% envw2_5_inc_order = {[envw2_5_inc_netcon(1),-3,envw2_5_inc_netcon(4),envw2_5_inc_netcon(3)],...
%     [-2,envw2_5_inc_netcon(2),envw2_5_inc_netcon(3),envw2_5_inc_netcon(4)],...
%     [envw2_5_inc_netcon(1),-1,envw2_5_inc_netcon(2)]};