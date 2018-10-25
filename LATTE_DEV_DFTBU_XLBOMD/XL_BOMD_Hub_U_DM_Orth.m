%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Extended Lagrangian Born-Oppenheimer Molecular Dynamics       %%%
%%         for H2 with restricted Hartree-Fock theory                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        A.M.N. Niklasson J. Chem. Phys. 147, 054102 (2017)         %%%
%% A.M.N. Niklasson & M. Cawkwell J. Chem. Phys. 141, 164123 (2014)  %%%
%% P. Souvatzis & A.M.N. Niklasson J. Chem. Phys. 140, 044117 (2014) %%%
%%        P. Steneteg et al. Phys. Rev. B 82, 075110 (2010)          %%%
%%       A.M.N Niklasson Phys. Rev. Lett. 100, 123001 (2008)         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HF Potential: Fockian = T + A + J_K(g,D) + H_Hub(D)              %%%
%%  U = 2*trace((T+A)*D) + trace(J_K*D) + EEnt + NNr                 %%%
%%            + 2*EHub + EEnt + NNr                                  %%%
%%  EHub = 0.5*trace(U*(DS - PS*DS - DS*PS + PS*PS))                 %%%
%%  H_Hub = 0.25*(SU-SUPS-SPSU+US-SPUS-USPS)                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
r = 3.7; v = 0.0;              % Initial position and velocity
%r = 1.45; v = 0.0;              % Initial position and velocity
c = 0.6; dr = 0.001;          % Mixing and finite difference step
c = 0.9; dr = 0.001;          % Mixing and finite difference step
Nocc = 1;
Temp = 10000;
%Temp = 100;
U = 1.0*diag([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]); % HUBBARD U VECTOR
%U = 0.0*diag([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]); % HUBBARD U VECTOR

%% Charge independent integrals and potentials
S = Overlap(r); T = Kinetic(r); A = External(r); g = TwoElInt(r); NNr = 1/r;
Z = S^(-0.5);                 % Inverse overlap factor
I = eye(8); PO = I/8;         % Initial guess of density matrix

%% SCF optimization P = SCF[P,r,n,c]
for i = 1:100
  P = Z*PO*Z';
  F = T+A+J_K(g,P) + H_Hub(U,P,S);  % Fockian, H_Hub = 0.25*(X + X'), X = S*U-S*U*P*S-S*P*S*U
  [DO,S_Ent,QQ,e,Fe_occ,mu0] = Fermi_Q(Z'*F*Z,Temp,Nocc);
  dDO = DO-PO; d_DO = dDO/norm(dDO,'fro');
  d_D = Z*d_DO*Z';
  F1 = J_K(g,d_D) + H_Hub(U,d_D,S); 
  [PO1] = Fermi_PRT(Z'*F*Z,Z'*F1*Z,Nocc,Temp,QQ,e,Fe_occ,mu0);  % PO1 = dD[F(PO)+c*F(d_D)]/dc at c=0 in orth. repr.
  VO = PO1 + ((1-c)/c)*d_DO;
  PO = PO + c*dDO + 1*(c^2*trace(d_DO'*dDO))/(1-c*trace(d_DO'*VO))*VO; % Rank-1 updated quantum response mixing
  SCF = norm(Z*dDO*Z');
end
P = Z*PO*Z';
SCF
pause
EEnt = -2.D0*Temp*S_Ent;

%% Calculate initial potential energy U[P,r]
JK = J_K(g,P); F = T+A+JK+ H_Hub(U,P,S);    % Fockian F[P] including Hubbard term
[DO,S_Ent,QQ,e,Fe_occ,mu0] = Fermi_Q(Z'*F*Z,Temp,Nocc);
D = Z*DO*Z';
[EHub] = E_Hub(U,D,S,P);  % EHub = 0.5*trace(U*(DS - PS*DS - DS*PS + PS*PS))
UPot = 2*trace((T+A)*D) + trace(JK*(2*D-P)) + 2*EHub + EEnt + NNr;  % Shadow potential

%% and force Force[P,r]
dS = dOverlap(r,dr); dT = dKinetic(r,dr); dA = dExternal(r,dr);
dg = dTwoElInt(r,dr); dJK = J_K(dg,P); dNNr = -1/r^2;
Pulay = -2*trace(Z*Z'*F*D*dS); % Pulay term for fractional occupation
[dEHub] = dE_Hub(U,D,S,P,dS);  % dEHub = 0.5*trace(U*(DdS-PdS*DS-PS*DdS-DdS*PS-DS*PdS+PdS*PS+PS*PdS))
dUdR = 2*trace((dT+dA)*D) + trace(dJK*(2*D-P)) + dNNr + Pulay + 2*dEHub;

%% Initial BC for P_n and coefficients for the modified Verlet scheme
% Kmax = 6 (or 5) in Niklasson, JCP 147, 054103 (2017) Integration coefficients
PO = DO;   P_0 = DO; P_1 = DO; P_2 = DO; P_3 = DO; P_4 = DO; P_5 = DO; P_6 = DO;
C0 = -14; C1 = 36; C2 = -27; C3 = -2; C4 = 12;  C5 = -6; C6 = 1;  
kappa = 1.84; alpha = 0.00055; 
%C0 = -6; C1 = 14; C2 = -8; C3 = -3; C4 = 4;  C5 = -1; C6 = 0.0;
%kappa = 1.82; alpha = 0.018;

dt = 1.0;    % Integration time step, unknown unit, minimum 15-20 per period in R(t)
%dt = 0.25;  % Integration time step, unknown unit, minimum 15-20 per period in R(t)
M_init = 6;  % Extra SFC convergence the first ~10 steps for smooth initialization

for i = 1:100  %%%%%%%% MAIN MD LOOP  %%%%%%%%%% 

  %% Check total energy, i.e. shadow Hamiltonian constant of motion
  Energy(i) = 0.5*v^2 + UPot; time(i) = (i-1)*dt; RR(i) =  r; vv(i) = v; UU(i) = UPot;
  EGap(i) = e(Nocc+1)-e(Nocc);
  [i,0.5*v^2 + UPot]

  v = v - dt*dUdR/2;          % Update velocity
  r = r + dt*v;               % Update position

  dDO = DO-PO; d_DO = dDO/norm(dDO,'fro');
  d_D = Z*d_DO*Z';
  F1 = J_K(g,d_D) + H_Hub(U,d_D,S);
  [PO1] = Fermi_PRT(Z'*F*Z,Z'*F1*Z,Nocc,Temp,QQ,e,Fe_occ,mu0); % DM response from change d_D
  VO = PO1 + ((1-c)/c)*d_DO;
  if i > 1
    dPO2dt2 = c*dDO + (c^2*norm(dDO,'fro')/(1-c*trace(d_DO'*VO)))*VO; % K as rank-1 updated scaled delta function
  else
    dPO2dt2 = c*dDO;  % K as scaled delta function
  end

  %% Integrate the XL-BOMD equations of motion, P'' = w^2*K(D-P), with modified Verlet, i.e.  
  %% P(t+dt) = 2P(t)-P(t-dt)+dt^2*P''(t) + Weak Dissipation
  if i <= M_init P_0 = DO; end % Use modified initial integration with high dissipation
  PO = 2*P_0 - P_1 + kappa*dPO2dt2 + alpha*(C0*P_0 + C1*P_1 + C2*P_2 + C3*P_3 + C4*P_4 + C5*P_5 + C6*P_6);
  P_6 = P_5; P_5 = P_4; P_4 = P_3; P_3 = P_2; P_2 = P_1; P_1 = P_0; P_0 = PO; % Shift values

  %% Update charge-independent integrals and potentials for new r
  S = Overlap(r); T = Kinetic(r); A = External(r); g = TwoElInt(r); Z = S^(-0.5); NNr = 1/r;

  Nr_SCF_It = 0;                    % Set Nr_SCF_It > 0 for "exact" BOMD
  if i < M_init Nr_SCF_It = 18; end % Use SCF initially 

  PO_tmp = PO;
  %% SCF optimization P = SCF[P,n,r,c], but only initially for i < M_init
  for k = 1:Nr_SCF_It 
    F = T+A+J_K(g,Z*PO*Z') + H_Hub(U,Z*PO*Z',S);     % Fockian
    [DO,S_Ent,QQ,e,Fe_occ,mu0] = Fermi_Q(Z'*F*Z,Temp,Nocc);
    dDO = DO-PO; d_DO = dDO/norm(dDO,'fro');
    F1 = J_K(g,Z*d_DO*Z') + H_Hub(U,Z*d_DO*Z',S);
    [PO1] = Fermi_PRT(Z'*F*Z,Z'*F1*Z,Nocc,Temp,QQ,e,Fe_occ,mu0);  % PO1 = dD[F(PO)+c*F(d_D)]/dc at c=0 in orth. repr.
    VO = PO1 + ((1-c)/c)*d_DO;
    PO = PO + c*dDO + 1*(c^2*norm(dDO,'fro')/(1-c*trace(d_DO'*VO)))*VO; % Rank-1 updated quantum response mixing
    SCF = norm(Z*dDO*Z');
  end
  P = Z*PO*Z';

  %% Shadow Hamiltonian U of XL-BOMD in 0-SCF limit {U[P,r]}
  JK = J_K(g,P); F = T+A+JK + H_Hub(U,P,S);  % Fockian F[P]
  [DO,S_Ent,QQ,e,Fe_occ,mu0] = Fermi_Q(Z'*F*Z,Temp,Nocc);
  D = Z*DO*Z';
  EEnt = -2.D0*Temp*S_Ent; % Entropy
  [EHub] = E_Hub(U,D,S,P); 
  UPot = 2*trace((T+A)*D) + trace(JK*(2*D-P)) + 2*EHub + EEnt + NNr;   % Shadow potential

  %% with exact force dU/dr {Force[P,r]}
  dS = dOverlap(r,dr); dT = dKinetic(r,dr); dA = dExternal(r,dr); 
  dg = dTwoElInt(r,dr); dJK = J_K(dg,P); dNNr = -1/r^2;
  Pulay = -2*trace(Z*Z'*F*D*dS); % Pulay term for fractional occupation
  [dEHub] = dE_Hub(U,D,S,P,dS); %  dEHub = 0.5*trace(U*(DdS-PdS*DS-PS*DdS-DdS*PS-DS*PdS+PdS*PS+PS*PdS))
  dUdR = 2*trace((dT+dA)*D) + trace(dJK*(2*D-P)) + dNNr + Pulay + 2*dEHub; % Total force

  v = v - dt*dUdR/2;          % Update velocity

end

plot(time,Energy-Energy(1),'-g*')
hold on
