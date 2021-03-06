function energy = update_top_PBC(Jmaxpos,iteration,slow,epsilon,sweep,sweepmax,chi_incs,shift,at_on)
%energy = update_top_PBC(Jmaxpos,iteration,slow,epsilon,sweep,sweepmax,chi_incs,shift,at_on)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - update_top_PBC
% performs an upadate of the top coarse-graining block
% 
% Andrew Goldsborough - 16/12/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global w u h VMDH envu envw1 envw2;

%contraction orders from netcon
envu_at2_order = {[3,2],[1,4],[4,17,3,-2],[2,-4,1,18],[18,16,17,15],[15,13,-1,14],[11,14,12,-3],[10,16,9,13],[8,9,7,11],[6,12,5,10],[7,6],[5,8]};
envu_at3_order = {[3,2],[1,4],[4,9,3,-2],[2,-4,1,10],[-1,13,-3,14],[11,13,12,14],[8,9,7,11],[6,12,5,10],[7,6],[5,8]};
envu_at4_order = {[3,2],[1,4],[4,17,3,-2],[2,-4,1,18],[18,16,17,15],[-3,13,16,14],[11,-1,12,13],[10,14,9,15],[8,9,7,11],[6,12,5,10],[7,6],[5,8]};
envw1_at2_order ={[-4,18],[17,-2],[18,15,17,16],[-3,13,15,14],[16,12,-1,11],[11,9,13,10],[7,10,8,14],[6,12,5,9],[4,5,3,7],[2,8,1,6],[3,2],[1,4]};
envw1_at3_order = {[-4,14],[13,-2],[14,11,13,12],[-3,9,11,10],[9,7,10,8],[5,7,6,8],[4,-1,3,5],[2,6,1,12],[3,2],[1,4]};
envw1_at4_order = {[-4,18],[17,-2],[18,15,17,16],[-3,13,15,14],[16,12,-1,11],[14,9,12,10],[7,13,8,9],[6,10,5,11],[4,5,3,7],[2,8,1,6],[3,2],[1,4]};
envw1_at5_order = {[-4,14],[13,-2],[14,11,13,12],[12,10,-1,9],[10,8,9,7],[6,8,5,7],[4,5,3,-3],[2,11,1,6],[3,2],[1,4]};
energy_at2_order = {[21,22],[17,18],[18,19,21,20],[22,15,17,16],[20,13,15,14],[16,12,19,11],[11,9,13,10],[7,10,8,14],[6,12,5,9],[4,5,3,7],[2,8,1,6],[3,2],[1,4]};
energy_at3_order = {[17,18],[13,14],[14,15,17,16],[18,11,13,12],[16,9,11,10],[9,7,10,8],[5,7,6,8],[4,15,3,5],[2,6,1,12],[3,2],[1,4]};
energy_at4_order = {[21,22],[17,18],[18,19,21,20],[22,15,17,16],[20,13,15,14],[16,12,19,11],[14,9,12,10],[7,13,8,9],[6,10,5,11],[4,5,3,7],[2,8,1,6],[3,2],[1,4]};
energy_at5_order = {[17,18],[13,14],[14,15,17,16],[18,11,13,12],[12,10,15,9],[10,8,9,7],[6,8,5,7],[4,5,3,16],[2,11,1,6],[3,2],[1,4]};

%update top tensor
Lcur = 4;
optmax = 2;
ULmax = 2;

if sweepmax ~= 1
    if chi_incs == 0
        %normal update
        
        [uleg1,uleg2,uleg3,uleg4] = size(u{iteration,1});
        [u2leg1,u2leg2,u2leg3,u2leg4] = size(u{iteration+1,1});
        [w1leg1,w1leg2,w1leg3,w1leg4] = size(w{iteration,1});
        [w2leg1,w2leg2,w2leg3,w2leg4] = size(w{iteration,2});
        
        if slow == 3 || slow == 4
            %shifts
            if sweep < sweepmax/10
                shift_frac = 0;
            else
                shift_frac = 1;
            end
            
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
            h2 = h{PBC_pos(Jmaxpos-1,Lcur)};
            h3 = h{PBC_pos(Jmaxpos,Lcur)};
            h4 = h{PBC_pos(Jmaxpos+1,Lcur)};
            h5 = h{PBC_pos(Jmaxpos+2,Lcur)};
        end
        
        for optcount = 1:optmax
            %u1
            for ULcount = 1:ULmax
                old = u{iteration,1};
                
                %diagram 1 - left most
                envu_1 = 0;
                
                %diagram 2 - left of centre
                envu_2 = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration+1,1},...
                    h2,conj(u{iteration,1}),conj(u{iteration+1,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),conj(VMDH),conj(VMDH)},envu_at2_order);
                
                %diagram 3 - centre
                envu_3 = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},h3,...
                    conj(u{iteration,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
                    conj(VMDH),conj(VMDH)},envu_at3_order);
                
                %diagram 4 - right of centre
                envu_4 = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration+1,1},...
                    h4,conj(u{iteration,1}),conj(u{iteration+1,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),conj(VMDH),conj(VMDH)},envu_at4_order);
                
                %diagram 5 - right most
                envu_5 = 0;
                
                %optimise u by UL
                if slow == 1 && isempty(envu{iteration}) ~= 1
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
                
                %check convergence
                change = ncon({u{iteration,1},conj(u{iteration,1})},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({u{iteration,1},conj(old)},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({old,conj(u{iteration,1})},{[1,2,3,4],[1,2,3,4]})...
                    + ncon({old,conj(old)},{[1,2,3,4],[1,2,3,4]});
                
                if change < 1e-12
                    break
                end
            end
            
            %alternative top
            if at_on == 1
                %u2 - permute by 2
                for ULcount = 1:ULmax
                    old = u{iteration+1,1};
                    
                    %diagram 1 - left most
                    envu_1 = 0;
                    
                    %diagram 2 - left of centre
                    envu_2 = ncon({VMDH,VMDH,w{iteration,2},w{iteration,1},u{iteration,1},...
                        h4,conj(u{iteration+1,1}),conj(u{iteration,1}),conj(w{iteration,2}),...
                        conj(w{iteration,1}),conj(VMDH),conj(VMDH)},envu_at2_order);
                    
                    %diagram 3 - centre
                    envu_3 = ncon({VMDH,VMDH,w{iteration,2},w{iteration,1},h5,...
                        conj(u{iteration+1,1}),conj(w{iteration,2}),conj(w{iteration,1}),...
                        conj(VMDH),conj(VMDH)},envu_at3_order);
                    
                    %diagram 4 - right of centre
                    envu_4 = ncon({VMDH,VMDH,w{iteration,2},w{iteration,1},u{iteration,1},...
                        h2,conj(u{iteration+1,1}),conj(u{iteration,1}),conj(w{iteration,2}),...
                        conj(w{iteration,1}),conj(VMDH),conj(VMDH)},envu_at4_order);
                    
                    %diagram 5 - right most
                    envu_5 = 0;
                    
                    %optimise u by UL
                    if slow == 1 && isempty(envu{iteration+1}) ~= 1
                        if sweep <= sweepmax/2
                            %slower
                            envu{iteration+1} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * envu{iteration+1}...
                                + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envu_1 + envu_2 + envu_3 + envu_4 + envu_5);
                        else
                            envu{iteration+1} = envu_1 + envu_2 + envu_3 + envu_4 + envu_5;
                        end
                    elseif (slow == 2 || slow == 4) && isempty(envu{iteration+1}) ~= 1
                        if sweep <= sweepmax/2
                            %random perturbation
                            envu{iteration+1} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * random('unif',-1,1,size(envu_2))...
                                + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envu_1 + envu_2 + envu_3 + envu_4 + envu_5);
                        else
                            envu{iteration+1} = envu_1 + envu_2 + envu_3 + envu_4 + envu_5;
                        end
                    else
                        envu{iteration+1} = envu_1 + envu_2 + envu_3 + envu_4 + envu_5;
                    end
                    
                    if u2leg1 == 1
                        envum = tfuse(envu{iteration+1},[2,-1,2]);
                    else
                        envum = tfuse(envu{iteration+1},[1,-2,1,-3]);
                        envum = tfuse(envum,[-1,2,2]);
                    end
                    
                    [U,~,V] = svd(envum,'econ');
                    u{iteration+1,1} = -V*U';
                    
                    %split into a 4-index tensor
                    if uleg1 == 1
                        u{iteration+1,1} = tsplit(u{iteration+1,1},1,[u2leg2,u2leg4]);
                        u{iteration+1,1} = permute(u{iteration+1,1},[4,1,3,2]);
                    else
                        u{iteration+1,1} = tsplit(u{iteration+1,1},1,[u2leg1,u2leg3]);
                        u{iteration+1,1} = tsplit(u{iteration+1,1},3,[u2leg2,u2leg4]);
                        u{iteration+1,1} = permute(u{iteration+1,1},[1,3,2,4]);
                    end
                    
                    %check convergence
                    change = ncon({u{iteration+1,1},conj(u{iteration+1,1})},{[1,2,3,4],[1,2,3,4]})...
                        - ncon({u{iteration+1,1},conj(old)},{[1,2,3,4],[1,2,3,4]})...
                        - ncon({old,conj(u{iteration+1,1})},{[1,2,3,4],[1,2,3,4]})...
                        + ncon({old,conj(old)},{[1,2,3,4],[1,2,3,4]});
                    
                    if change < 1e-12
                        break
                    end
                end
            end
            
            %w1
            for ULcount = 1:ULmax
                old = w{iteration,1};
                
                %diagram 1 - left most
                envw1_1 = 0;
                
                %diagram 2 - left of centre
                envw1_2 = ncon({VMDH,VMDH,w{iteration,2},u{iteration,1},u{iteration+1,1},...
                    h2,conj(u{iteration,1}),conj(u{iteration+1,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),conj(VMDH),conj(VMDH)},envw1_at2_order);
                
                %diagram 3 - centre
                envw1_3 = ncon({VMDH,VMDH,w{iteration,2},u{iteration,1},h3,...
                    conj(u{iteration,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
                    conj(VMDH),conj(VMDH)},envw1_at3_order);
                
                %diagram 4 - right of centre
                envw1_4 = ncon({VMDH,VMDH,w{iteration,2},u{iteration,1},u{iteration+1,1},...
                    h4,conj(u{iteration,1}),conj(u{iteration+1,1}),conj(w{iteration,1}),...
                    conj(w{iteration,2}),conj(VMDH),conj(VMDH)},envw1_at4_order);
                
                %diagram 5 - right most
                envw1_5 = ncon({VMDH,VMDH,w{iteration,2},u{iteration+1,1},h5,...
                    conj(u{iteration+1,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
                    conj(VMDH),conj(VMDH)},envw1_at5_order);
                
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
                
                %check convergence
                change = ncon({w{iteration,1},conj(w{iteration,1})},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({w{iteration,1},conj(old)},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({old,conj(w{iteration,1})},{[1,2,3,4],[1,2,3,4]})...
                    + ncon({old,conj(old)},{[1,2,3,4],[1,2,3,4]});
                
                if change < 1e-12
                    break
                end
            end
            
            %w2 - permute by 2
            for ULcount = 1:ULmax
                old = w{iteration,2};
                
                %diagram 1 - left most
                envw1_1 = 0;
                
                %diagram 2 - left of centre
                envw1_2 = ncon({VMDH,VMDH,w{iteration,1},u{iteration+1,1},u{iteration,1},...
                    h4,conj(u{iteration+1,1}),conj(u{iteration,1}),conj(w{iteration,2}),...
                    conj(w{iteration,1}),conj(VMDH),conj(VMDH)},envw1_at2_order);
                
                %diagram 3 - centre
                envw1_3 = ncon({VMDH,VMDH,w{iteration,1},u{iteration+1,1},h5,...
                    conj(u{iteration+1,1}),conj(w{iteration,2}),conj(w{iteration,1}),...
                    conj(VMDH),conj(VMDH)},envw1_at3_order);
                
                %diagram 4 - right of centre
                envw1_4 = ncon({VMDH,VMDH,w{iteration,1},u{iteration+1,1},u{iteration,1},...
                    h2,conj(u{iteration+1,1}),conj(u{iteration,1}),conj(w{iteration,2}),...
                    conj(w{iteration,1}),conj(VMDH),conj(VMDH)},envw1_at4_order);
                
                %diagram 5 - right most
                envw1_5 = ncon({VMDH,VMDH,w{iteration,1},u{iteration,1},h3,...
                    conj(u{iteration,1}),conj(w{iteration,2}),conj(w{iteration,1}),...
                    conj(VMDH),conj(VMDH)},envw1_at5_order);
                
                %optimise w1 by SVD
                if slow == 1 && isempty(envw2{iteration}) ~= 1
                    if sweep <= sweepmax/2
                        %slower
                        envw2{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * envw2{iteration}...
                            + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5);
                    else
                        envw2{iteration} = envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5;
                    end
                elseif (slow == 2 || slow == 4) && isempty(envw2{iteration}) ~= 1
                    if sweep <= sweepmax/2
                        %random perturbation
                        envw2{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * random('unif',-1,1,size(envw1_1))...
                            + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5);
                    else
                        envw2{iteration} = envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5;
                    end
                else
                    envw2{iteration} = envw1_1 + envw1_2 + envw1_3 + envw1_4 + envw1_5;
                end
                
                if w2leg1 == 1
                    envw2m = tfuse(envw2{iteration},[2,-1,2]);
                else
                    envw2m = tfuse(envw2{iteration},[1,-2,1,-3]);
                    envw2m = tfuse(envw2m,[-1,2,2]);
                end
                
                [U,~,V] = svd(envw2m,'econ');
                w{iteration,2} = -V*U';
                
                %split into a 4-index tensor
                if w2leg1 == 1
                    w{iteration,2} = tsplit(w{iteration,2},2,[w2leg2,w2leg4]);
                    w{iteration,2} = permute(w{iteration,2},[4,2,1,3]);
                else
                    w{iteration,2} = tsplit(w{iteration,2},1,[w2leg1,w2leg3]);
                    w{iteration,2} = tsplit(w{iteration,2},3,[w2leg2,w2leg4]);
                    w{iteration,2} = permute(w{iteration,2},[1,3,2,4]);
                end
                
                %check convergence
                change = ncon({w{iteration,2},conj(w{iteration,2})},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({w{iteration,2},conj(old)},{[1,2,3,4],[1,2,3,4]})...
                    - ncon({old,conj(w{iteration,2})},{[1,2,3,4],[1,2,3,4]})...
                    + ncon({old,conj(old)},{[1,2,3,4],[1,2,3,4]});
                
                if change < 1e-12
                    break
                end
            end
        end
    else %chi_incs == 1
        error('chi increase top not implemented, use ED');
    end
end
    
%print energy
energy2 = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration,1},u{iteration+1,1},...
    h{PBC_pos(Jmaxpos-1,Lcur)},conj(u{iteration,1}),conj(u{iteration+1,1}),conj(w{iteration,1}),...
    conj(w{iteration,2}),conj(VMDH),conj(VMDH)},energy_at2_order);

energy3 = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration,1},...
    h{PBC_pos(Jmaxpos,Lcur)},conj(u{iteration,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
    conj(VMDH),conj(VMDH)},energy_at3_order);

energy4 = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration,1},u{iteration+1,1},...
    h{PBC_pos(Jmaxpos+1,Lcur)},conj(u{iteration,1}),conj(u{iteration+1,1}),conj(w{iteration,1}),...
    conj(w{iteration,2}),conj(VMDH),conj(VMDH)},energy_at4_order);

energy5 = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration+1,1},...
    h{PBC_pos(Jmaxpos+2,Lcur)},conj(u{iteration+1,1}),conj(w{iteration,1}),conj(w{iteration,2}),...
    conj(VMDH),conj(VMDH)},energy_at5_order);

energy = energy2 + energy3 + energy4 + energy5 + shift;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find optimal contraction orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% envu_at2_netcon(netcon({[1,2],[3,4],[4,5,1,-2],[2,-4,3,6],[6,8,5,7],[7,9,-1,10],[12,10,13,-3],[14,8,11,9],[18,11,15,12],[16,13,17,14],[15,16],[17,18]},0,2,1,1)) = 1:18;
% envu_at2_order = {[envu_at2_netcon(1),envu_at2_netcon(2)],...
%     [envu_at2_netcon(3),envu_at2_netcon(4)],...
%     [envu_at2_netcon(4),envu_at2_netcon(5),envu_at2_netcon(1),-2],...
%     [envu_at2_netcon(2),-4,envu_at2_netcon(3),envu_at2_netcon(6)],...
%     [envu_at2_netcon(6),envu_at2_netcon(8),envu_at2_netcon(5),envu_at2_netcon(7)],...
%     [envu_at2_netcon(7),envu_at2_netcon(9),-1,envu_at2_netcon(10)],...
%     [envu_at2_netcon(12),envu_at2_netcon(10),envu_at2_netcon(13),-3],...
%     [envu_at2_netcon(14),envu_at2_netcon(8),envu_at2_netcon(11),envu_at2_netcon(9)],...
%     [envu_at2_netcon(18),envu_at2_netcon(11),envu_at2_netcon(15),envu_at2_netcon(12)],...
%     [envu_at2_netcon(16),envu_at2_netcon(13),envu_at2_netcon(17),envu_at2_netcon(14)],...
%     [envu_at2_netcon(15),envu_at2_netcon(16)],...
%     [envu_at2_netcon(17),envu_at2_netcon(18)]};
% 
% envu_at3_netcon(netcon({[1,2],[3,4],[4,5,1,-2],[2,-4,3,6],[-1,7,-3,8],[9,7,10,8],[14,5,11,9],[12,10,13,6],[11,12],[13,14]},0,2,1,1)) = 1:14;
% envu_at3_order = {[envu_at3_netcon(1),envu_at3_netcon(2)],...
%     [envu_at3_netcon(3),envu_at3_netcon(4)],...
%     [envu_at3_netcon(4),envu_at3_netcon(5),envu_at3_netcon(1),-2],...
%     [envu_at3_netcon(2),-4,envu_at3_netcon(3),envu_at3_netcon(6)],...
%     [-1,envu_at3_netcon(7),-3,envu_at3_netcon(8)],...
%     [envu_at3_netcon(9),envu_at3_netcon(7),envu_at3_netcon(10),envu_at3_netcon(8)],...
%     [envu_at3_netcon(14),envu_at3_netcon(5),envu_at3_netcon(11),envu_at3_netcon(9)],...
%     [envu_at3_netcon(12),envu_at3_netcon(10),envu_at3_netcon(13),envu_at3_netcon(6)],...
%     [envu_at3_netcon(11),envu_at3_netcon(12)],...
%     [envu_at3_netcon(13),envu_at3_netcon(14)]};
% 
% envu_at4_netcon(netcon({[1,2],[3,4],[4,5,1,-2],[2,-4,3,6],[6,8,5,7],[-3,9,8,10],[12,-1,13,9],[14,10,11,7],[18,11,15,12],[16,13,17,14],[15,16],[17,18]},0,2,1,1)) = 1:18;
% envu_at4_order = {[envu_at4_netcon(1),envu_at4_netcon(2)],...
%     [envu_at4_netcon(3),envu_at4_netcon(4)],...
%     [envu_at4_netcon(4),envu_at4_netcon(5),envu_at4_netcon(1),-2],...
%     [envu_at4_netcon(2),-4,envu_at4_netcon(3),envu_at4_netcon(6)],...
%     [envu_at4_netcon(6),envu_at4_netcon(8),envu_at4_netcon(5),envu_at4_netcon(7)],...
%     [-3,envu_at4_netcon(9),envu_at4_netcon(8),envu_at4_netcon(10)],...
%     [envu_at4_netcon(12),-1,envu_at4_netcon(13),envu_at4_netcon(9)],...
%     [envu_at4_netcon(14),envu_at4_netcon(10),envu_at4_netcon(11),envu_at4_netcon(7)],...
%     [envu_at4_netcon(18),envu_at4_netcon(11),envu_at4_netcon(15),envu_at4_netcon(12)],...
%     [envu_at4_netcon(16),envu_at4_netcon(13),envu_at4_netcon(17),envu_at4_netcon(14)],...
%     [envu_at4_netcon(15),envu_at4_netcon(16)],...
%     [envu_at4_netcon(17),envu_at4_netcon(18)]};
% 
% envw1_at2_netcon(netcon({[-4,1],[2,-2],[1,3,2,4],[-3,6,3,7],[4,8,-1,5],[5,9,6,10],[12,10,13,7],[14,8,11,9],[18,11,15,12],[16,13,17,14],[15,16],[17,18]},0,2,1,1)) = 1:18;
% envw1_at2_order = {[-4,envw1_at2_netcon(1)],[envw1_at2_netcon(2),-2],...
%     [envw1_at2_netcon(1),envw1_at2_netcon(3),envw1_at2_netcon(2),envw1_at2_netcon(4)],...
%     [-3,envw1_at2_netcon(6),envw1_at2_netcon(3),envw1_at2_netcon(7)],...
%     [envw1_at2_netcon(4),envw1_at2_netcon(8),-1,envw1_at2_netcon(5)],...
%     [envw1_at2_netcon(5),envw1_at2_netcon(9),envw1_at2_netcon(6),envw1_at2_netcon(10)],...
%     [envw1_at2_netcon(12),envw1_at2_netcon(10),envw1_at2_netcon(13),envw1_at2_netcon(7)],...
%     [envw1_at2_netcon(14),envw1_at2_netcon(8),envw1_at2_netcon(11),envw1_at2_netcon(9)],...
%     [envw1_at2_netcon(18),envw1_at2_netcon(11),envw1_at2_netcon(15),envw1_at2_netcon(12)],...
%     [envw1_at2_netcon(16),envw1_at2_netcon(13),envw1_at2_netcon(17),envw1_at2_netcon(14)],...
%     [envw1_at2_netcon(15),envw1_at2_netcon(16)],...
%     [envw1_at2_netcon(17),envw1_at2_netcon(18)]};
% 
% envw1_at3_netcon(netcon({[-4,1],[2,-2],[1,3,2,4],[-3,5,3,6],[5,7,6,8],[9,7,10,8],[14,-1,11,9],[12,10,13,4],[11,12],[13,14]},0,2,1,1)) = 1:14;
% envw1_at3_order = {[-4,envw1_at3_netcon(1)],[envw1_at3_netcon(2),-2],...
%     [envw1_at3_netcon(1),envw1_at3_netcon(3),envw1_at3_netcon(2),envw1_at3_netcon(4)],...
%     [-3,envw1_at3_netcon(5),envw1_at3_netcon(3),envw1_at3_netcon(6)],...
%     [envw1_at3_netcon(5),envw1_at3_netcon(7),envw1_at3_netcon(6),envw1_at3_netcon(8)],...
%     [envw1_at3_netcon(9),envw1_at3_netcon(7),envw1_at3_netcon(10),envw1_at3_netcon(8)],...
%     [envw1_at3_netcon(14),-1,envw1_at3_netcon(11),envw1_at3_netcon(9)],...
%     [envw1_at3_netcon(12),envw1_at3_netcon(10),envw1_at3_netcon(13),envw1_at3_netcon(4)],...
%     [envw1_at3_netcon(11),envw1_at3_netcon(12)],...
%     [envw1_at3_netcon(13),envw1_at3_netcon(14)]};
% 
% envw1_at4_netcon(netcon({[-4,1],[2,-2],[1,3,2,4],[-3,6,3,7],[4,8,-1,5],[7,9,8,10],[12,6,13,9],[14,10,11,5],[18,11,15,12],[16,13,17,14],[15,16],[17,18]},0,2,1,1)) = 1:18;
% envw1_at4_order = {[-4,envw1_at4_netcon(1)],[envw1_at4_netcon(2),-2],...
%     [envw1_at4_netcon(1),envw1_at4_netcon(3),envw1_at4_netcon(2),envw1_at4_netcon(4)],...
%     [-3,envw1_at4_netcon(6),envw1_at4_netcon(3),envw1_at4_netcon(7)],...
%     [envw1_at4_netcon(4),envw1_at4_netcon(8),-1,envw1_at4_netcon(5)],...
%     [envw1_at4_netcon(7),envw1_at4_netcon(9),envw1_at4_netcon(8),envw1_at4_netcon(10)],...
%     [envw1_at4_netcon(12),envw1_at4_netcon(6),envw1_at4_netcon(13),envw1_at4_netcon(9)],...
%     [envw1_at4_netcon(14),envw1_at4_netcon(10),envw1_at4_netcon(11),envw1_at4_netcon(5)],...
%     [envw1_at4_netcon(18),envw1_at4_netcon(11),envw1_at4_netcon(15),envw1_at4_netcon(12)],...
%     [envw1_at4_netcon(16),envw1_at4_netcon(13),envw1_at4_netcon(17),envw1_at4_netcon(14)],...
%     [envw1_at4_netcon(15),envw1_at4_netcon(16)],...
%     [envw1_at4_netcon(17),envw1_at4_netcon(18)]};
%
% envw1_at5_netcon(netcon({[-4,1],[2,-2],[1,3,2,4],[4,6,-1,5],[6,8,5,7],[10,8,9,7],[14,9,11,-3],[12,3,13,10],[11,12],[13,14]},0,2,1,1)) = 1:14;
% envw1_at5_order = {[-4,envw1_at5_netcon(1)],[envw1_at5_netcon(2),-2],...
%     [envw1_at5_netcon(1),envw1_at5_netcon(3),envw1_at5_netcon(2),envw1_at5_netcon(4)],...
%     [envw1_at5_netcon(4),envw1_at5_netcon(6),-1,envw1_at5_netcon(5)],...
%     [envw1_at5_netcon(6),envw1_at5_netcon(8),envw1_at5_netcon(5),envw1_at5_netcon(7)],...
%     [envw1_at5_netcon(10),envw1_at5_netcon(8),envw1_at5_netcon(9),envw1_at5_netcon(7)],...
%     [envw1_at5_netcon(14),envw1_at5_netcon(9),envw1_at5_netcon(11),-3],...
%     [envw1_at5_netcon(12),envw1_at5_netcon(3),envw1_at5_netcon(13),envw1_at5_netcon(10)],...
%     [envw1_at5_netcon(11),envw1_at5_netcon(12)],...
%     [envw1_at5_netcon(13),envw1_at5_netcon(14)]};
% 
% energy_at2_netcon(netcon({[1,2],[3,4],[4,5,1,6],[2,7,3,8],[6,10,7,11],[8,12,5,9],[9,13,10,14],[16,14,17,11],[18,12,15,13],[22,15,19,16],[20,17,21,18],[19,20],[21,22]},0,2,1,1)) = 1:22;
% energy_at2_order = {[energy_at2_netcon(1),energy_at2_netcon(2)],...
%     [energy_at2_netcon(3),energy_at2_netcon(4)],...
%     [energy_at2_netcon(4),energy_at2_netcon(5),energy_at2_netcon(1),energy_at2_netcon(6)],...
%     [energy_at2_netcon(2),energy_at2_netcon(7),energy_at2_netcon(3),energy_at2_netcon(8)],...
%     [energy_at2_netcon(6),energy_at2_netcon(10),energy_at2_netcon(7),energy_at2_netcon(11)],...
%     [energy_at2_netcon(8),energy_at2_netcon(12),energy_at2_netcon(5),energy_at2_netcon(9)],...
%     [energy_at2_netcon(9),energy_at2_netcon(13),energy_at2_netcon(10),energy_at2_netcon(14)],...
%     [energy_at2_netcon(16),energy_at2_netcon(14),energy_at2_netcon(17),energy_at2_netcon(11)],...
%     [energy_at2_netcon(18),energy_at2_netcon(12),energy_at2_netcon(15),energy_at2_netcon(13)],...
%     [energy_at2_netcon(22),energy_at2_netcon(15),energy_at2_netcon(19),energy_at2_netcon(16)],...
%     [energy_at2_netcon(20),energy_at2_netcon(17),energy_at2_netcon(21),energy_at2_netcon(18)],...
%     [energy_at2_netcon(19),energy_at2_netcon(20)],...
%     [energy_at2_netcon(21),energy_at2_netcon(22)]};
% 
% energy_at3_netcon(netcon({[1,2],[3,4],[4,5,1,6],[2,7,3,8],[6,9,7,10],[9,11,10,12],[13,11,14,12],[18,5,15,13],[16,14,17,8],[15,16],[17,18]},0,2,1,1)) = 1:18;
% energy_at3_order = {[energy_at3_netcon(1),energy_at3_netcon(2)],...
%     [energy_at3_netcon(3),energy_at3_netcon(4)],...
%     [energy_at3_netcon(4),energy_at3_netcon(5),energy_at3_netcon(1),energy_at3_netcon(6)],...
%     [energy_at3_netcon(2),energy_at3_netcon(7),energy_at3_netcon(3),energy_at3_netcon(8)],...
%     [energy_at3_netcon(6),energy_at3_netcon(9),energy_at3_netcon(7),energy_at3_netcon(10)],...
%     [energy_at3_netcon(9),energy_at3_netcon(11),energy_at3_netcon(10),energy_at3_netcon(12)],...
%     [energy_at3_netcon(13),energy_at3_netcon(11),energy_at3_netcon(14),energy_at3_netcon(12)],...
%     [energy_at3_netcon(18),energy_at3_netcon(5),energy_at3_netcon(15),energy_at3_netcon(13)],...
%     [energy_at3_netcon(16),energy_at3_netcon(14),energy_at3_netcon(17),energy_at3_netcon(8)],...
%     [energy_at3_netcon(15),energy_at3_netcon(16)],...
%     [energy_at3_netcon(17),energy_at3_netcon(18)]};
% 
% energy_at4_netcon(netcon({[1,2],[3,4],[4,5,1,6],[2,7,3,8],[6,10,7,11],[8,12,5,9],[11,13,12,14],[16,10,17,13],[18,14,15,9],[22,15,19,16],[20,17,21,18],[19,20],[21,22]},0,2,1,1)) = 1:22;
% energy_at4_order = {[energy_at4_netcon(1),energy_at4_netcon(2)],...
%     [energy_at4_netcon(3),energy_at4_netcon(4)],...
%     [energy_at4_netcon(4),energy_at4_netcon(5),energy_at4_netcon(1),energy_at4_netcon(6)],...
%     [energy_at4_netcon(2),energy_at4_netcon(7),energy_at4_netcon(3),energy_at4_netcon(8)],...
%     [energy_at4_netcon(6),energy_at4_netcon(10),energy_at4_netcon(7),energy_at4_netcon(11)],...
%     [energy_at4_netcon(8),energy_at4_netcon(12),energy_at4_netcon(5),energy_at4_netcon(9)],...
%     [energy_at4_netcon(11),energy_at4_netcon(13),energy_at4_netcon(12),energy_at4_netcon(14)],...
%     [energy_at4_netcon(16),energy_at4_netcon(10),energy_at4_netcon(17),energy_at4_netcon(13)],...
%     [energy_at4_netcon(18),energy_at4_netcon(14),energy_at4_netcon(15),energy_at4_netcon(9)],...
%     [energy_at4_netcon(22),energy_at4_netcon(15),energy_at4_netcon(19),energy_at4_netcon(16)],...
%     [energy_at4_netcon(20),energy_at4_netcon(17),energy_at4_netcon(21),energy_at4_netcon(18)],...
%     [energy_at4_netcon(19),energy_at4_netcon(20)],...
%     [energy_at4_netcon(21),energy_at4_netcon(22)]};
% 
% energy_at5_netcon(netcon({[1,2],[3,4],[4,5,1,6],[2,7,3,8],[8,10,5,9],[10,12,9,11],[14,12,13,11],[18,13,15,6],[16,7,17,14],[15,16],[17,18]},0,2,1,1)) = 1:18;
% energy_at5_order = {[energy_at5_netcon(1),energy_at5_netcon(2)],...
%     [energy_at5_netcon(3),energy_at5_netcon(4)],...
%     [energy_at5_netcon(4),energy_at5_netcon(5),energy_at5_netcon(1),energy_at5_netcon(6)],...
%     [energy_at5_netcon(2),energy_at5_netcon(7),energy_at5_netcon(3),energy_at5_netcon(8)],...
%     [energy_at5_netcon(8),energy_at5_netcon(10),energy_at5_netcon(5),energy_at5_netcon(9)],...
%     [energy_at5_netcon(10),energy_at5_netcon(12),energy_at5_netcon(9),energy_at5_netcon(11)],...
%     [energy_at5_netcon(14),energy_at5_netcon(12),energy_at5_netcon(13),energy_at5_netcon(11)],...
%     [energy_at5_netcon(18),energy_at5_netcon(13),energy_at5_netcon(15),energy_at5_netcon(6)],...
%     [energy_at5_netcon(16),energy_at5_netcon(7),energy_at5_netcon(17),energy_at5_netcon(14)],...
%     [energy_at5_netcon(15),energy_at5_netcon(16)],...
%     [energy_at5_netcon(17),energy_at5_netcon(18)]};