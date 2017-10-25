function update_PBC_XX(Jmaxpos,iteration,optmax,ULmax,slow,epsilon,sweep,sweepmax,chi_inc)
%update_PBC_XX(Jmaxpos,iteration,optmax,ULmax,slow,epsilon,sweep,sweepmax,chi_inc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - update_PBC_XX
% performs update using XX parameterisation
% 
% Andrew Goldsborough - 29/09/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%update using the theta parameterisation
global w u rho h VMDH;
global envu envw1 envw2;

%only works for chi=2
if all(chi_inc == 0) ~= 1
    error('chi must be 2');
end

%store current length
Lcur = size(h,1);

%contraction orders from netcon
envw1_1_order = {[3,2,5,-2],[2,3,-1,4],[-4,1],[5,4,6,-3],[6,1]};
envw1_2_order = {[13,-2,14,12],[-4,1],[1,10,12,11],[-3,8,10,9],[-1,6,8,7],[4,7,5,9],[13,6,3,4],[3,2],[2,5,14,11]};
envw1_3_order = {[13,-2,14,12],[-4,1],[1,8,12,9],[-3,5,8,6],[5,3,6,4],[10,3,7,4],[13,-1,11,10],[11,2],[2,7,14,9]};
envw1_4_order = {[13,-2,14,12],[-4,1],[1,8,12,9],[-3,6,8,7],[7,4,9,5],[10,6,3,4],[13,-1,11,10],[11,2],[2,3,14,5]};
envw2_2_order = {[13,12,14,-4],[12,8,1,9],[1,-2],[9,6,-1,7],[8,4,6,5],[3,5,10,7],[13,4,2,3],[2,11],[11,10,14,-3]};
envw2_3_order = {[13,12,14,-4],[12,8,1,9],[1,-2],[9,5,-1,6],[5,3,6,4],[7,3,10,4],[13,8,2,7],[2,11],[11,10,14,-3]};
envw2_4_order = {[13,12,14,-4],[12,10,1,11],[1,-2],[11,8,-1,9],[9,6,-3,7],[4,8,5,6],[13,10,3,4],[3,2],[2,5,14,7]};
envw2_5_order = {[1,-2],[5,-4,3,2],[-3,4,2,3],[1,6],[6,-1,5,4]};

%useful values
[w1leg1,~,~,~] = size(w{iteration,1});
[w2leg1,~,~,~] = size(w{iteration,2});

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
    %u - no update
    u{iteration,1} = zeros([2,2,2,2]);
    u{iteration,1}(1:2,1:2,1:2,1:2) = -ncon({eye(2),eye(2)},{[-1,-4],[-3,-2]});
    u{iteration,1}(1,1,1,1) = 1;
    
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
        
        %optimise w1
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
        
        %find extrema
        theta = atan((envw1m(2,3) - envw1m(3,2))/(envw1m(2,2) + envw1m(3,3)));
        
        if (-(envw1m(2,3) - envw1m(3,2)) * sin(theta) - (envw1m(2,2) + envw1m(3,3)) * cos(theta)) > 0
            %theta is a min
        elseif (-(envw1m(2,3) - envw1m(3,2)) * sin(theta+pi) - (envw1m(2,2) + envw1m(3,3)) * cos(theta+pi)) > 0
            %theta + pi is a min
            theta = theta+pi;
        else
            %no minimum
            error('MERAupdate:no_min','no minimum theta');
        end
        
        %build w
        w{iteration,1} = zeros(2,2,2,2);
        w{iteration,1}(1,1,1,1) = 1;
        w{iteration,1}(2,2,2,2) = 1;
        w{iteration,1}(2,2,1,1) = cos(theta);
        w{iteration,1}(1,1,2,2) = cos(theta);
        w{iteration,1}(2,1,1,2) = sin(theta);
        w{iteration,1}(1,2,2,1) = -sin(theta);
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
        
        %optimise w2
        if slow == 1 && isempty(envw2{iteration}) ~= 1
            if sweep <= sweepmax/2
                %slower
                envw2{iteration} = (((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * envw2{iteration}...
                    + (1 - ((0.5*sweepmax - sweep + 1)/(0.5*sweepmax))*epsilon) * (envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5);
            else
                envw2{iteration} = envw2_1 + envw2_2 + envw2_3 + envw2_4 + envw2_5;
            end
        elseif slow == 2 && isempty(envu{iteration}) ~= 1
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
        
        %find extrema
        theta = atan((envw2m(2,3) - envw2m(3,2))/(envw2m(2,2) + envw2m(3,3)));
        
        if (-(envw2m(2,3) - envw2m(3,2)) * sin(theta) - (envw2m(2,2) + envw2m(3,3)) * cos(theta)) > 0
            %theta is a min
        elseif (-(envw2m(2,3) - envw2m(3,2)) * sin(theta+pi) - (envw2m(2,2) + envw2m(3,3)) * cos(theta+pi)) > 0
            %theta + pi is a min
            theta = theta+pi;
        else
            %no minimum
            error('MERAupdate:no_min','no minimum theta');
        end
        
        %build w
        w{iteration,2} = zeros(2,2,2,2);
        w{iteration,2}(1,1,1,1) = 1;
        w{iteration,2}(2,2,2,2) = 1;
        w{iteration,2}(2,2,1,1) = cos(theta);
        w{iteration,2}(1,1,2,2) = cos(theta);
        w{iteration,2}(2,1,1,2) = sin(theta);
        w{iteration,2}(1,2,2,1) = -sin(theta);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find optimal contraction orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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