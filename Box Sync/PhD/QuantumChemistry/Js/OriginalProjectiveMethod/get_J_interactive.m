disp ("Assumes the presence of the following files: MOs_1.txt, MOs_2.txt, MOs_pair.txt, S_1.txt, S_2.txt, S_pair.txt, Evls_pair.txt")
disp ("!!! Octave counts from 1 !!!")

nbasis_A = input("Number of MOs (of molA)")		%... of monomer
nhomo_mon_A = input("HOMO (of molA)")	%... of monomer (!! Octave counts from 1 !!)
nbasis_B = input("Number of MOs of molB")
nhomo_mon_B = input("HOMO (of molB)")
%nbasis=446
%nhomo_mon=73
%nbasis=698
%nhomo_mon=165
%nbasis=46
%nhomo_mon=8

%sprintf('read data')
load('MOs_1.txt');
load('MOs_2.txt');
load('MOs_pair.txt');
load('S_1.txt');
load('S_2.txt');
load('S_pair.txt');
load('Evls_pair.txt');

zero_mat_A=zeros(nbasis_A,nbasis_B);
zero_mat_B=zeros(nbasis_B,nbasis_A);
MOs_mons = [MOs_1, zero_mat_A; zero_mat_B, MOs_2];
S_mons   = [S_1,   zero_mat_A; zero_mat_B, S_2 ];
D_mons   = chol(S_mons);
D_pair   = chol(S_pair);

D_mons_sqrt = sqrtm(S_mons);
D_pair_sqrt = sqrtm(S_pair);

%D=D_mons'*D_mons;

%C=D(1:10,1:10)
%D_orig=S_mons(1:10,1:10)

A=D_mons_sqrt(1:10,1:10)
B=D_mons(1:10,1:10)

%b3lyp - have to orthogonalise
MOs_mons_orth = D_mons * MOs_mons;
MOs_pair_orth = D_pair * MOs_pair;

B = MOs_mons_orth' * MOs_pair_orth;
Evls=diag(Evls_pair');
H_eff = B * Evls * B';

%test
%B_test = MOs_mons' * MOs_pair;
%H_eff_test = B_test * Evls * B_test';
%H00_test=H_eff_test(nhomo_mon+nbasis, nhomo_mon)
%           -0.5*(H_eff_test(nhomo_mon+nbasis,nhomo_mon+nbasis)+H_eff_test(nhomo_mon,nhomo_mon))*S_pair(nhomo_mon+nbasis,nhomo_mon)
%	  )/(1-S_pair(nhomo_mon+nbasis,nhomo_mon)*S_pair(nhomo_mon+nbasis,nhomo_mon))

%<HOMO | F | HOMO > 
%<LUMO | F | LUMO > 
%Output.  N.B. Energies already converted to eV by rewrite_S_phi_E.cpp

%non-degenerate
%H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%E1_E2=(H_eff(nhomo_mon_A,nhomo_mon_A)-H_eff(nhomo_mon_B+nbasis_A,nhomo_mon_B+nbasis_A))
%half_HOMO_diff=(Evls(2*nhomo_mon,2*nhomo_mon)-Evls(2*nhomo_mon-1,2*nhomo_mon-1))/2
Non_degenerate_HOMO_coupling=H_eff(nhomo_mon_B+nbasis_A,   nhomo_mon_A  )
Non_degenerate_LUMO_coupling=H_eff(nhomo_mon_B+1+nbasis_A, nhomo_mon_A+1)
%Non_degenerate_HOMO_LUMO_coupling=H_eff(nhomo_mon+nbasis, nhomo_mon+1)
%Non_degenerate_LUMO_HOMO_coupling=H_eff(nhomo_mon+1+nbasis, nhomo_mon)
%L00=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%Non_degenerate_LUMO_coupling=abs(L00)

%doubly degenerate
%H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%H01=H_eff(nhomo_mon+nbasis,   nhomo_mon-1);
%H10=H_eff(nhomo_mon+nbasis-1, nhomo_mon  );
%H11=H_eff(nhomo_mon+nbasis-1, nhomo_mon-1);
%L00=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%L01=H_eff(nhomo_mon+nbasis+1, nhomo_mon);
%L10=H_eff(nhomo_mon+nbasis,   nhomo_mon+1);
%L11=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%Doubly_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H10)+abs(H11))/4
%Doubly_degenerate_LUMO_coupling=(abs(L00)+abs(L01)+abs(L10)+abs(L11))/4

%triply degenerate
%H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%H01=H_eff(nhomo_mon+nbasis,   nhomo_mon-1);
%H02=H_eff(nhomo_mon+nbasis,   nhomo_mon-2);
%H10=H_eff(nhomo_mon+nbasis-1, nhomo_mon  );
%H11=H_eff(nhomo_mon+nbasis-1, nhomo_mon-1);
%H12=H_eff(nhomo_mon+nbasis-1, nhomo_mon-2);
%H20=H_eff(nhomo_mon+nbasis-2, nhomo_mon  );
%H21=H_eff(nhomo_mon+nbasis-2, nhomo_mon-1);
%H22=H_eff(nhomo_mon+nbasis-2, nhomo_mon-2);
%L00=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%L01=H_eff(nhomo_mon+nbasis+1, nhomo_mon);
%L02=H_eff(nhomo_mon+nbasis+1, nhomo_mon-1);
%L10=H_eff(nhomo_mon+nbasis,   nhomo_mon+1);
%L11=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%L12=H_eff(nhomo_mon+nbasis,   nhomo_mon-1);
%L20=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%L21=H_eff(nhomo_mon+nbasis+1, nhomo_mon  );
%L22=H_eff(nhomo_mon+nbasis+1, nhomo_mon-1);
%Triply_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H02)+abs(H10)+abs(H11)+abs(H12)+abs(H20)+abs(H21)+abs(H22))/9
%Triply_degenerate_LUMO_coupling=(abs(L00)+abs(L01)+abs(L02)+abs(L10)+abs(L11)+abs(L12)+abs(L20)+abs(L21)+abs(L22))/9

%H03=H_eff(nhomo_mon+nbasis,   nhomo_mon-3);
%H04=H_eff(nhomo_mon+nbasis,   nhomo_mon-4);
%H13=H_eff(nhomo_mon+nbasis-1, nhomo_mon-3);
%H14=H_eff(nhomo_mon+nbasis-1, nhomo_mon-4);
%H22=H_eff(nhomo_mon+nbasis-2, nhomo_mon-3);
%H22=H_eff(nhomo_mon+nbasis-2, nhomo_mon-4);


%Four_fold_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H02)+abs(H03)
%				  + abs(H10)+abs(H11)+abs(H12)+abs(H13)
%				  + abs(H20)+abs(H21)+abs(H22)+abs(H23)
%				  + abs(H30)+abs(H31)+abs(H32)+abs(H33))/16
%Five_fold_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H02)+abs(H03)+abs(H04)
%				  + abs(H10)+abs(H11)+abs(H12)+abs(H13)+abs(H14)
%				  + abs(H20)+abs(H21)+abs(H22)+abs(H23)+abs(H24)
%				  + abs(H30)+abs(H31)+abs(H32)+abs(H33)+abs(H34)
%				  + abs(H40)+abs(H41)+abs(H42)+abs(H43)+abs(H44))/25
