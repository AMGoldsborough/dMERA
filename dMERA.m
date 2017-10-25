function dMERA(L,Jstr,Jdis,Jz,chi_w,Pdist,Jseed,ULmax,sweepmax,slow)
%dMERA(L,Jstr,Jdis,Jz,chi_w,Pdist,Jseed,ULmax,sweepmax,slow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA PBC
% code for paper by Goldsborough and Evenbly
%
% Andrew Goldsborough - 27/02/2017
%
% index convention
%   -1     -1  -3  
%    |      |___|  
%   /_\     |___|  
%   | |     |   |  
% -2  -3   -2  -4  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%inputs
%when compiled the command line inputs are strings, convert to numbers
if ischar(L)==1
  L = str2double(L);
end
if ischar(Jstr)==1
  Jstr = str2double(Jstr);
end
if ischar(Jdis)==1
  Jdis = str2double(Jdis);
end
if ischar(Jz)==1
  Jz = str2double(Jz);
end
if ischar(chi_w)==1
  chi_w = str2double(chi_w);
end
if ischar(Pdist)==1
  Pdist = str2double(Pdist);
end
if ischar(Jseed)==1
  Jseed = str2double(Jseed);
end
if ischar(ULmax)==1
  ULmax = str2double(ULmax);
end
if ischar(sweepmax)==1
  sweepmax = str2double(sweepmax);
end
if ischar(slow)==1
  slow = str2double(slow);
end

% L = 10;         %chain length
% Jstr = 1;       %overall J strength
% Jdis = 2;       %disorder strength
% Jz = 1;         %anisotropy
% chi_w = 8;      %max chi w
% Pdist = 6;      %coupling distribution
% Jseed = 2;      %seed for rng, 0 => shuffle
% ULmax = 1;      %number of updates on single tensor
optmax = ULmax;   %number of updates on tensor block
% sweepmax = 500; %number of full network sweeps
% slow = 0;       %update options

%coupling distribution
%0 => manual
%1 => 2 theta(K-1/2)
%2 => 1
%3 => uniform around Jstr normalised by Jstr
%4 => uniform around Jstr un-normalised Jdis
%5 => box distribution of Hikihara AF (10.1103/PhysRevB.60.12116)
%6 => Laflorencie's infinite disorder distribution (10.1103/PhysRevB.72.140408)

%update options (slow) (0 = off, 1 = slow, 2 = random, 3 = energy shift, 4 = 3 and 2)
epsilon = 1e-10;

%convergence
energy_old = 0;
criterion = 1e-12;

%start type (0 = random start, 1 = singlet start)
start_type = 1;

%alternative top on/off
at_on = 0;

%calculate things
if L > 30
    ee_on = 0;
else
    ee_on = 1;
end
corr_on = 1;
exp_on = 1;

%turn off warnings
warning('off','MATLAB:eigs:SigmaChangedToSA');
warning('off','ncon:suboptimalsequence');

%storage for tensors
global w u rho h VMDH;
global envu envw1 envw2;
w = cell((L/2)-1,2);
u = cell((L/2),1);
rho = cell((L/2)-2,L-2);
envu = cell((L/2),1);
envw1 = cell((L/2)-1,1);
envw2 = cell((L/2)-1,1);
h = cell(L,1);
h_full = cell((L/2),L);

%define MDH singlet
VMDH = 0.5*sqrt(2)*[0 1; -1 0];
chi_sing = 2; %dimension of singlet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate couplings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set seed, set by clock if zero
if Jseed > 0
    rng(Jseed);
else
    rng('shuffle');
end

%set the probability distribution
if Pdist==0
    %custom set
    J = [0.1,Jdis,0.1,1,Jdis,1,0.1,Jdis,0.1,0.01];
elseif Pdist==1
    %P(K) = 2 theta(K-1/2)
    J = zeros(1,L) + Jstr*random('unif',0.5,1,[1,L]);
elseif Pdist==2
    %P(K) = 1
    J = zeros(1,L) + Jstr*(rand(L,1));
elseif Pdist==3
    %uniform around Jstr normalised by Jstr
    J = zeros(1,L) + Jstr + Jstr*Jdis*(rand(1,L) - 0.5);
elseif Pdist==4
    %uniform around Jstr un-normalised Jdis
    J = zeros(1,L) + Jstr + Jdis*(rand(1,L) - 0.5);
elseif Pdist==5
    %box distribution of Hikihara AF
    J = zeros(1,L) + Jdis*random('unif',0,1,[1,L]);
elseif Pdist==6
    %Laflorencie's infinite disorder distribution
    J = rand(1,L).^Jdis;
end

rng('shuffle');

%print J to file
fprintf('printing interaction strengths\n');
fname = strcat('./J/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_J.txt');
fidJ = fopen(fname, 'w');
for i = 1:L
    fprintf(fidJ,'%.15e\n',J(i));
end
fclose(fidJ);

%create hamiltonian
shift = 0;
for i=1:L
    [h{i},~,shift2] = Heisham_twosite(J(i),Jz);
    shift = shift + shift2;
end

h0 = h;
for i = 1:size(h,1)
    h_full{1,i} = h{i};
end

%store which iterations increase chi
chi_inc = zeros(1,(L/2)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate start tensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find order by MDH
Jmaxpos_vec = MDH_order_PBC(L,J,Jz);

%initialise tensors    
if start_type == 0
    %random tensors
    chi_inc = random_start_PBC(L,chi_w,chi_sing,Jmaxpos_vec,chi_inc,at_on);
elseif start_type == 1
    %MDH start => identities 4 - crossed u
    chi_inc = id_start_PBC(L,chi_w,chi_sing,Jmaxpos_vec,chi_inc);
else
    error('start_type = 0 (random) or 1 (identity)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sweep over the blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%print output to see convergence
fname = strcat('./output/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi_w),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_output_dMERA.txt');
fidout = fopen(fname, 'w');

if sweepmax == 1
    %fill tensors and calculate energy
    
    %update rho
    makerho_PBC(L,Jmaxpos_vec,chi_inc);
    
    %reset h
    h = h0;

    for iteration = 1:(L/2)-2
        %update Hamiltonian
        raiseHam_PBC(Jmaxpos_vec(iteration),iteration,chi_inc(iteration));
    end
    
    %top tensor
    energy = update_top_PBC(Jmaxpos_vec((L/2)-1),(L/2)-1,slow,epsilon,sweepmax,sweepmax,chi_inc((L/2)-1),shift,at_on);
    
else
    for sweep = 2:sweepmax
        
        if sweep ~= 2
            %reupdate top
            [~] = update_top_PBC(Jmaxpos_vec((L/2)-1),(L/2)-1,slow,epsilon,sweep,sweepmax,chi_inc((L/2)-1),shift,at_on);
        end
        
        %update rho
        makerho_update_PBC(L,Jmaxpos_vec,chi_inc,sweep,optmax,ULmax,slow,epsilon,sweepmax,h_full);
        
        %reset h
        h = h0;
        
        %update tensors
        for iteration = 1:(L/2)-2
            
            %update tensors
            update_PBC(Jmaxpos_vec(iteration),iteration,optmax,ULmax,slow,epsilon,sweep,sweepmax,chi_inc(iteration));
            
            %update Hamiltonian
            raiseHam_PBC(Jmaxpos_vec(iteration),iteration,chi_inc(iteration));
            
            for i = 1:size(h,1)
                h_full{iteration+1,i} = h{i};
            end
        end
        
        %top tensor
        energy = update_top_PBC(Jmaxpos_vec((L/2)-1),(L/2)-1,slow,epsilon,sweep,sweepmax,chi_inc((L/2)-1),shift,at_on);
        
        fprintf('%.15e\n',energy);
        
        %print energy to file
        fprintf(fidout,'%.15e\n',energy);
        
        %check convergence
        if abs(energy - energy_old) < criterion
            fprintf('converged!\n');
            fprintf(fidout,'converged!\n');
            break
        else
            energy_old = energy;
        end
    end
end

fclose(fidout);

%print energy to file
fprintf('printing energy\n');
fname = strcat('./energy/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi_w),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_energy_dMERA.txt');
fidenergy = fopen(fname, 'w');
fprintf(fidenergy,'%.15e\n',energy);
fclose(fidenergy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sz correlation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corr_on == 1    
    %make operators
    op_L = [0.5 0;0 -0.5];
    op_R = [0.5 0;0 -0.5];
    
    %open file to print to
    fprintf('printing Sz correlation functions\n');
    fname = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi_w),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_Szcorr_dMERA.txt');
    fidSzcorr = fopen(fname, 'w');
    
    %loop over all pairs of sites
    for SL = 1:(L-1)
        for SR = (SL+1):L            
            %calculate separation wrt PBC
            PBC_r = min(PBC_pos(SL-SR,L),PBC_pos(SR-SL,L));
            
            %nn
            if PBC_r == 1
                %put both spins in one operator
                twosite_L = ncon({op_L,op_R},{[-1,-2],[-3,-4]});
                twosite_R = [];
                
                if SL == 1 && SR == L
                    Sj = SL;
                    Si = SR;
                else
                    Sj = SR;
                    Si = SL;
                end
            elseif PBC_r == 2
                %nnn
                %check costs and choose cheapest
                cost = zeros(1,2);
                cost(1) = corrcheck_PBC(L,Jmaxpos_vec,SL,SR,chi_inc);
                cost(2) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),PBC_pos(SR-1,L),chi_inc);
                [~,corr_setup] = min(cost);
                
                if corr_setup == 1
                    %LL
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = SR;
                else %corr_setup == 2
                    %RR
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = PBC_pos(SR - 1,L);
                end
            else %PBC_r >= 3
                
                %check costs and choose cheapest
                cost = zeros(1,4);
                cost(1) = corrcheck_PBC(L,Jmaxpos_vec,SL,SR,chi_inc);
                cost(2) = corrcheck_PBC(L,Jmaxpos_vec,SL,PBC_pos(SR-1,L),chi_inc);
                cost(3) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),SR,chi_inc);
                cost(4) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),PBC_pos(SR-1,L),chi_inc);
                [~,corr_setup] = min(cost);
                
                if corr_setup == 1
                    %LL
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = SR;
                elseif corr_setup == 2
                    %LR
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = SR - 1;
                elseif corr_setup == 3
                    %RL
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = SR;
                else %corr_setup == 4
                    %RR
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = PBC_pos(SR - 1,L);
                end
            end
            
            %calculate correlation
            try
                corr = corr_PBC(L,Jmaxpos_vec,Si,Sj,u,w,VMDH,twosite_L,twosite_R,chi_inc);
                fprintf(fidSzcorr,'%d %d %.15e\n',SL,SR,corr);
            catch
                fprintf('Szcorr: Skipping %d %d\n',SL,SR);
                continue
            end
        end
    end
    
    fclose(fidSzcorr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SpSm correlation functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %make operators
    op_L = [0 1;0 0];
    op_R = [0 0;1 0];
    
    %open file to print to
    fprintf('printing SpSm correlation functions\n');
    fname = strcat('./SpSmcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi_w),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_SpSmcorr_dMERA.txt');
    fidSzcorr = fopen(fname, 'w');
    
    %loop over all pairs of sites
    for SL = 1:(L-1)
        for SR = (SL+1):L
            %calculate separation wrt PBC
            PBC_r = min(PBC_pos(SL-SR,L),PBC_pos(SR-SL,L));
            
            %nn
            if PBC_r == 1
                %put both spins in one operator
                twosite_L = ncon({op_L,op_R},{[-1,-2],[-3,-4]});
                twosite_R = [];
                
                if SL == 1 && SR == L
                    Sj = SL;
                    Si = SR;
                else
                    Sj = SR;
                    Si = SL;
                end
            elseif PBC_r == 2
                %nnn
                %check costs and choose cheapest
                cost = zeros(1,2);
                cost(1) = corrcheck_PBC(L,Jmaxpos_vec,SL,SR,chi_inc);
                cost(2) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),PBC_pos(SR-1,L),chi_inc);
                [~,corr_setup] = min(cost);
                
                if corr_setup == 1
                    %LL
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = SR;
                else %corr_setup == 2
                    %RR
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = PBC_pos(SR - 1,L);
                end
            else %PBC_r >= 3
                
                %check costs and choose cheapest
                cost = zeros(1,4);
                cost(1) = corrcheck_PBC(L,Jmaxpos_vec,SL,SR,chi_inc);
                cost(2) = corrcheck_PBC(L,Jmaxpos_vec,SL,PBC_pos(SR-1,L),chi_inc);
                cost(3) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),SR,chi_inc);
                cost(4) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),PBC_pos(SR-1,L),chi_inc);
                [~,corr_setup] = min(cost);
                
                if corr_setup == 1
                    %LL
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = SR;
                elseif corr_setup == 2
                    %LR
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = PBC_pos(SR - 1,L);
                elseif corr_setup == 3
                    %RL
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = SR;
                else %corr_setup == 4
                    %RR
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = PBC_pos(SR - 1,L);
                end
            end
            
            %calculate correlation
            try
                corr = corr_PBC(L,Jmaxpos_vec,Si,Sj,u,w,VMDH,twosite_L,twosite_R,chi_inc);
                fprintf(fidSzcorr,'%d %d %.15e\n',SL,SR,corr);
            catch
                fprintf('SpSmcorr: Skipping %d %d\n',SL,SR);
                continue
            end
        end
    end
    
    fclose(fidSzcorr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SmSp correlation functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %make operators
    op_L = [0 0;1 0];
    op_R = [0 1;0 0];
    
    %open file to print to
    fprintf('printing SmSp correlation functions\n');
    fname = strcat('./SmSpcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi_w),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_SmSpcorr_dMERA.txt');
    fidSzcorr = fopen(fname, 'w');
    
    %loop over all pairs of sites
    for SL = 1:(L-1)
        for SR = (SL+1):L
            %calculate separation wrt PBC
            PBC_r = min(PBC_pos(SL-SR,L),PBC_pos(SR-SL,L));
            
            %nn
            if PBC_r == 1
                %put both spins in one operator
                twosite_L = ncon({op_L,op_R},{[-1,-2],[-3,-4]});
                twosite_R = [];
                
                if SL == 1 && SR == L
                    Sj = SL;
                    Si = SR;
                else
                    Sj = SR;
                    Si = SL;
                end
            elseif PBC_r == 2
                %nnn
                %check costs and choose cheapest
                cost = zeros(1,2);
                cost(1) = corrcheck_PBC(L,Jmaxpos_vec,SL,SR,chi_inc);
                cost(2) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),PBC_pos(SR-1,L),chi_inc);
                [~,corr_setup] = min(cost);
                
                if corr_setup == 1
                    %LL
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = SR;
                else %corr_setup == 2
                    %RR
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = PBC_pos(SR - 1,L);
                end
            else %PBC_r >= 3
                
                %check costs and choose cheapest
                cost = zeros(1,4);
                cost(1) = corrcheck_PBC(L,Jmaxpos_vec,SL,SR,chi_inc);
                cost(2) = corrcheck_PBC(L,Jmaxpos_vec,SL,PBC_pos(SR-1,L),chi_inc);
                cost(3) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),SR,chi_inc);
                cost(4) = corrcheck_PBC(L,Jmaxpos_vec,PBC_pos(SL-1,L),PBC_pos(SR-1,L),chi_inc);
                [~,corr_setup] = min(cost);
                
                if corr_setup == 1
                    %LL
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = SR;
                elseif corr_setup == 2
                    %LR
                    twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = SL;
                    Sj = PBC_pos(SR - 1,L);
                elseif corr_setup == 3
                    %RL
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({op_R,eye(2)},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = SR;
                else %corr_setup == 4
                    %RR
                    twosite_L = ncon({eye(2),op_L},{[-1,-2],[-3,-4]});
                    twosite_R = ncon({eye(2),op_R},{[-1,-2],[-3,-4]});
                    
                    Si = PBC_pos(SL - 1,L);
                    Sj = PBC_pos(SR - 1,L);
                end
            end
            
            %calculate correlation
            try
                corr = corr_PBC(L,Jmaxpos_vec,Si,Sj,u,w,VMDH,twosite_L,twosite_R,chi_inc);
                fprintf(fidSzcorr,'%d %d %.15e\n',SL,SR,corr);
            catch
                fprintf('SmSpcorr: Skipping %d %d\n',SL,SR);
                continue
            end
        end
    end
    
    fclose(fidSzcorr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sz expectation values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exp_on == 1
    
    %make operators
    op_L = [0.5 0;0 -0.5];
    
    %open file to print to
    fprintf('printing Sz expectation values\n');
    fname = strcat('./Szexp/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi_w),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_Szexp_dMERA.txt');
    fidSzexp = fopen(fname, 'w');
    
    %loop over all pairs of sites
    for SL = 1:L
        SR = PBC_pos(SL+1,L);
        
        %put both spins in one operator
        twosite_L = ncon({op_L,eye(2)},{[-1,-2],[-3,-4]});
        twosite_R = [];
        
        %calculate correlation
        try
            corr = corr_PBC(L,Jmaxpos_vec,SL,SR,u,w,VMDH,twosite_L,twosite_R,chi_inc);
            fprintf(fidSzexp,'%d %.15e\n',SL,corr);
        catch
            fprintf('Szexp: Skipping %d\n',SL);
            continue
        end
        
    end
    
    fclose(fidSzexp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%entanglement entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ee_on == 1
    
    %print to file
    fprintf('printing ee\n');
    fname = strcat('./ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi_w),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_ee_dMERA.txt');
    fidee = fopen(fname, 'w');
    
    %loop over the bipartition positions (splitting chain into two)
    % AAAA/BB/AAAA
    for i=1:L-1
        for j=i+1:L
            %create indices for which parts are A(2) and which are B(0)
            blocks = [zeros(1,i),2*ones(1,j-i),zeros(1,L-j)];
            
            try
                ee = eeSVD_PBC(L,Jmaxpos_vec,chi_inc,blocks);
%                 ee = eeDM_PBC(L,Jmaxpos_vec,chi_inc,blocks);
                fprintf(fidee,'%d %d %.15e\n',i,j,ee);
            catch 
                fprintf('ee failed: %d %d\n',i,j);
                continue
            end
        end
    end
    
    fclose(fidee);
end
toc
end
