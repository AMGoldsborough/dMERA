function [h,d,shift] = HeisHam_twosite(J,Jz)
%[h,d,shift] = HeisHam_twosite(J,Jz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HeisHam_twosite
% creates the two site Heisenberg XXZ Hamiltonian
% note energy shift for use with MERA
% 
% Andrew Goldsborough - 11/11/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define pauli matrices
pauli_p = [0 1; 0 0];
pauli_m = [0 0; 1 0];
pauli_z = [1 0; 0 -1];

%h = Jz*ZZ + 0.5*(PM + MP)
h = J*(Jz*ncon({0.5*pauli_z,0.5*pauli_z},{[-1,-2],[-3,-4]})...
    + 0.5*(ncon({pauli_p,pauli_m},{[-1,-2],[-3,-4]})...
    + ncon({pauli_m,pauli_p},{[-1,-2],[-3,-4]})));

%shift energy spectrum by largest eigenvalue
shift = max(eig(tfuse(tfuse(h,[1,-2,1,-3]),[-1,2,2])));
h = h - shift * ncon({eye(2),eye(2)},{[-1,-2],[-3,-4]});

%dimension of spin = 2
d = 2;
end