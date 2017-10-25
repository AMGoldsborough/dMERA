function [ee,n_a,chi] = eeDM_PBC(L,Jmaxpos_vec,chi_inc,blocks)
%[ee,n_a,chi] = eeDM_PBC(L,Jmaxpos_vec,chi_inc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - eeDM_PBC
% calculates entanglement entropy by density matrix
% 
% Andrew Goldsborough - 03/02/2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global w u VMDH;

%contraction orders from netcon
at_order = {[3,2],[1,4],[4,5,3,7],[2,8,1,6],[7,-2,8,-3],[6,-4,5,-1]};
block_order = {[3,1],[-1,-3,3,4],[1,2,-2,-6],[4,-4,2,-5]};
    
%find which way is cheapest (tracing over A or B)
chi = zeros(1,2);
n_a = zeros(1,2);

for ver=1:2
    %Use original blocks for ver 1, inverse for ver 2
    if ver == 1
        blocks1 = blocks;
    end
    
    if ver == 2
        blocks1 = -blocks + 2;
    end

    %find blocks below each leg => bbelow
    %all A = 2, all B = 0, mix = 1
    %note which tensors are connected => tbelow
    bbelow = cell(size(Jmaxpos_vec,2),1);
    tbelow = cell(size(Jmaxpos_vec,2),1);
    tlegs = zeros(1,L);
    
    for i = 1:size(Jmaxpos_vec,2)
        
        %propagate blocks up
        Lcur = size(blocks1,2);
        bbelow{i} = blocks1(PBC_pos(Jmaxpos_vec(i)-1:Jmaxpos_vec(i)+2,Lcur));
        
        if bbelow{i} == 2*ones(1,4)
            %all in A
            blocks1(PBC_pos(Jmaxpos_vec(i)-1,Lcur)) = 2;
            blocks1(PBC_pos(Jmaxpos_vec(i)+2,Lcur)) = 2;
            
        elseif bbelow{i} == zeros(1,4)
            %all in B
            blocks1(PBC_pos(Jmaxpos_vec(i)-1,Lcur)) = 0;
            blocks1(PBC_pos(Jmaxpos_vec(i)+2,Lcur)) = 0;
            
        else
            %in both A and B
            blocks1(PBC_pos(Jmaxpos_vec(i)-1,Lcur)) = 1;
            blocks1(PBC_pos(Jmaxpos_vec(i)+2,Lcur)) = 1;
        end
        
        %propagate tensor number up
        tbelow{i} = tlegs(PBC_pos(Jmaxpos_vec(i)-1:Jmaxpos_vec(i)+2,Lcur));
        tlegs(PBC_pos(Jmaxpos_vec(i)-1,Lcur)) = i;
        tlegs(PBC_pos(Jmaxpos_vec(i)+2,Lcur)) = i;
        
        %remove two sites
        blocks1(Jmaxpos_vec(i)) = [];
        blocks1(PBC_pos(Jmaxpos_vec(i),Lcur-1)) = [];
        tlegs(Jmaxpos_vec(i)) = [];
        tlegs(PBC_pos(Jmaxpos_vec(i),Lcur-1)) = [];
    end
    
    %save bbelow and tbelow for ver 1
    if ver == 1
        bbelow_1 = bbelow;
        tbelow_1 = tbelow;
    end
    
    %keep track of the size of the indices
    numele = zeros(1,4);
    [~,numele(4),~,numele(1)] = size(u{L/2,1});
    [~,numele(2),~,numele(3)] = size(u{(L/2)-1,1});
    
    %note the full blocks below of the density matrix
    bbelow_dm = bbelow{(L/2)-1};
    tbelow_dm = tbelow{(L/2)-1};
    
    %start with top tensor, trace over any zeros
    tbelow_dm(bbelow_dm==0) = [];
    numele(bbelow_dm==0) = [];
    bbelow_dm(bbelow_dm==0) = [];
    
    for iteration = (L/2)-2:-1:1
        %"contract" and trace over zeros
        idx = find(tbelow_dm==iteration);

        if isempty(idx) == 1
            %tensor not to be contracted
            continue
        end
        
        if idx(1)==1 && idx(2)==size(tbelow_dm,2)
            %PBC term
            tbelow_dm = [tbelow{iteration}(3:4),tbelow_dm(idx(1)+1:idx(2)-1),tbelow{iteration}(1:2)];
            numele = [size(u{iteration,1},4),size(w{iteration,2},4),numele(idx(1)+1:idx(2)-1),size(w{iteration,1},2),size(u{iteration,1},2)];
            bbelow_dm = [bbelow{iteration}(3:4),bbelow_dm(idx(1)+1:idx(2)-1),bbelow{iteration}(1:2)];
        else
            %normal
            tbelow_dm = [tbelow_dm(1:idx(1)-1),tbelow{iteration},tbelow_dm(idx(2)+1:end)];
            numele = [numele(1:idx(1)-1),size(w{iteration,1},2),size(u{iteration,1},2),size(u{iteration,1},4),size(w{iteration,2},4),numele(idx(2)+1:end)];
            bbelow_dm = [bbelow_dm(1:idx(1)-1),bbelow{iteration},bbelow_dm(idx(2)+1:end)];
        end
        
        %trace over any zeros
        tbelow_dm(bbelow_dm==0) = [];
        numele(bbelow_dm==0) = [];
        bbelow_dm(bbelow_dm==0) = [];
    end
    
    chi(ver) = prod(numele);
    n_a(ver) = size(numele,2);
end

%choose cheapest option
if chi(1) > chi(2)
    n_a = n_a(2);
    chi = chi(2);
else
    n_a = n_a(1);
    chi = chi(1);
    bbelow = bbelow_1;
    tbelow = tbelow_1;
end

%start with top
iteration = (L/2)-1;

%note the full blocks below of the density matrix
bbelow_dm = bbelow{iteration};
tbelow_dm = tbelow{iteration};

rho_v = zeros(1,8);

%find B blocks
idx = find(bbelow_dm==0);
v_idx = 1;

%set to be traced over
for i = 1:size(idx,2)
    rho_v(2*idx(i)-1) = v_idx;
    rho_v(2*idx(i)) = v_idx;
    v_idx = v_idx + 1;
end

%set rest to -ve
rho_v(rho_v==0) = -(1:sum(rho_v==0));

%start by contracting the top
rdm = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration,1},u{iteration+1,1}},at_order);
rdm = ncon({rdm,conj(rdm)},{[rho_v(2),rho_v(4),rho_v(6),rho_v(8)],[rho_v(1),rho_v(3),rho_v(5),rho_v(7)]});

%remove traced over indices
tbelow_dm(bbelow_dm==0) = [];
bbelow_dm(bbelow_dm==0) = [];

%contract down
for iteration = (L/2)-2:-1:1
    
    %check to see if complete
    if sum(bbelow_dm==1)==0
        %no more contractions
        break
    end
    
    %find where the next block connects
    t_idx = find(tbelow_dm==iteration);
    
    if isempty(t_idx) == 1
        %tensor not to be contracted
        continue
    end
    
    %contract block
    if chi_inc(iteration) == 0
        %normal
        block = ncon({VMDH,w{iteration,1},w{iteration,2},u{iteration,1}},block_order);
    else
        %increase chi
        block = ncon({w{iteration,1},w{iteration,2},u{iteration,1}},...
            {[-1,-3,1],[-2,2,-6],[1,-4,2,-5]});
    end
        
    %make coarse graining blocks
    block_v = zeros(1,6);
    block_v([1,2]) = [1,2];
    block_d_v = zeros(1,6);
    block_d_v([1,2]) = [3,4];
    v_idx = 5;
    
    %find B blocks
    idx = find(bbelow{iteration}==0);
    
    %set to be traced over
    for i = 1:size(idx,2)
        block_v(idx(i)+2) = v_idx;
        block_d_v(idx(i)+2) = v_idx;
        v_idx = v_idx + 1;
    end
    
    %remove traced over indices
    tbelow{iteration}(bbelow{iteration}==0) = [];
    bbelow{iteration}(bbelow{iteration}==0) = [];
    
    %add new block in
    if t_idx(1)==1 && t_idx(2)==size(tbelow_dm,2)
        %PBC term
        block_v(block_v==0) = -(2:2:2*sum(block_v==0)) - 2*(t_idx(2)-2);
        block_d_v(block_d_v==0) = -(1:2:2*sum(block_d_v==0)) - 2*(t_idx(2)-2);
        rho_v = [[4,2],-(1:2*(t_idx(2)-2)),[3,1]];
        
        tbelow_dm = [tbelow_dm(t_idx(1)+1:t_idx(2)-1),tbelow{iteration}];
        bbelow_dm = [bbelow_dm(t_idx(1)+1:t_idx(2)-1),bbelow{iteration}];
    else
        %normal
        block_v(block_v==0) = -(2:2:2*sum(block_v==0)) - 2*(t_idx(1)-1);
        block_d_v(block_d_v==0) = -(1:2:2*sum(block_d_v==0)) -2*(t_idx(1)-1);
        rho_v = [-(1:2*(t_idx(1)-1)),[3,1,4,2],-(1:2*(size(tbelow_dm,2)-t_idx(2))) + min(block_v)];
        
        tbelow_dm = [tbelow_dm(1:t_idx(1)-1),tbelow{iteration},tbelow_dm(t_idx(2)+1:end)];
        bbelow_dm = [bbelow_dm(1:t_idx(1)-1),bbelow{iteration},bbelow_dm(t_idx(2)+1:end)];
    end
    
    %contract
    rdm = ncon({rdm,block,conj(block)},{rho_v,block_v,block_d_v});
end

%fuse to get a matrix
fuse_v = ones(1,ndims(rdm));
fuse_v(2:2:end) = -(2:(ndims(rdm)/2)+1);
ee = eentropydm(tfuse(tfuse(rdm,fuse_v),[-1,2*ones(1,ndims(rdm)/2)]));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netcon commands for orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at_netcon(netcon({[1,2],[3,4],[4,5,1,6],[2,7,3,8],[6,-2,7,-3],[8,-4,5,-1]},0,2,1,1)) = 1:8;
% at_order = {[at_netcon(1),at_netcon(2)],...
%     [at_netcon(3),at_netcon(4)],...
%     [at_netcon(4),at_netcon(5),at_netcon(1),at_netcon(6)],...
%     [at_netcon(2),at_netcon(7),at_netcon(3),at_netcon(8)],...
%     [at_netcon(6),-2,at_netcon(7),-3],...
%     [at_netcon(8),-4,at_netcon(5),-1]};
%
% block_netcon(netcon({[1,2],[-1,-3,1,3],[2,4,-2,-6],[3,-4,4,-5]},0,2,1,1)) = 1:4;
% block_order = {[block_netcon(1),block_netcon(2)],...
%     [-1,-3,block_netcon(1),block_netcon(3)],...
%     [block_netcon(2),block_netcon(4),-2,-6],...
%     [block_netcon(3),-4,block_netcon(4),-5]};