function ee = eeSVD_PBC(L,Jmaxpos_vec,chi_inc,blocks)
%ee = eeSVD_PBC(L,Jmaxpos_vec,chi_inc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - eeSVD_PBC
% calculates entanglement entropy by SVD
% 
% Andrew Goldsborough - 08/02/2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global w u VMDH;

blocks1 = blocks;
    
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
    
%contract top
iteration = (L/2)-1;
psi = ncon({VMDH,VMDH,w{iteration,1},w{iteration,2},u{iteration,1},...
    u{iteration+1,1}},{[3,2],[1,4],[4,5,3,7],[2,8,1,6],...
    [7,-2,8,-3],[6,-4,5,-1]});

%note the full blocks below
bbelow_psi = bbelow{(L/2)-1};
tbelow_psi = tbelow{(L/2)-1};

for iteration = (L/2)-2:-1:1
    
    %check to see if complete
    if sum(bbelow_psi==1)==0
        %no more contractions
        break
    end
    
    %contract down
    idx = find(tbelow_psi==iteration);
    
    if isempty(idx) == 1
        %tensor not to be contracted
        continue
    end
    
    %check if contraction is needed
    if bbelow_psi(idx(1)) == 1 && bbelow_psi(idx(2)) == 1
        
        %normal
        if chi_inc(iteration) == 0
            block = ncon({w{iteration,1},VMDH,w{iteration,2},u{iteration,1}},...
                {[-1,-3,1,3],[1,2],[2,4,-2,-6],[3,-4,4,-5]});
        else %chi_inc(iteration) == 1
            block = ncon({w{iteration,1},w{iteration,2},u{iteration,1}},...
                {[-1,-3,1],[-2,2,-6],[1,-4,2,-5]});
        end
        
        if idx(1)==1 && idx(2)==size(tbelow_psi,2)
            %PBC term
            psi = ncon({psi,block},{[2,-(1:(idx(2)-2)),1],...
                [1,2,-((idx(2)-1):idx(2)+2)]});
            
            tbelow_psi = [tbelow_psi(idx(1)+1:idx(2)-1),tbelow{iteration}];
            bbelow_psi = [bbelow_psi(idx(1)+1:idx(2)-1),bbelow{iteration}];
        else 
            psi = ncon({psi,block},{[-(1:(idx(1)-1)),[1,2],-(idx(2)+3:(size(tbelow_psi,2)+2))],...
                -[-1,-2,idx(1):idx(1)+3]});
            
            tbelow_psi = [tbelow_psi(1:idx(1)-1),tbelow{iteration},tbelow_psi(idx(2)+1:end)];
            bbelow_psi = [bbelow_psi(1:idx(1)-1),bbelow{iteration},bbelow_psi(idx(2)+1:end)];
        end
    end
end

%svd for ee

%create indices for which parts are A(2) and which are B(0)
svd_v = bbelow_psi;
svd_v(bbelow_psi==0) = 1;
svd_v(bbelow_psi==2) = -(2:(sum(bbelow_psi==2)+1));
eeM = tfuse(tfuse(psi,svd_v),[-1,2*ones(1,sum(bbelow_psi==2))]);
[~,S,~] = svd(eeM,'econ');
ee = eentropy(S);
end