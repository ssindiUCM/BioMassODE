function [time,y] = CheckAgain/TestMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py CheckAgain/Test.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
CheckAgain/TestParams

% Set the Initial Conditions
CheckAgain/TestIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs,flow), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
CheckAgain/TestRename
%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, A, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('A');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p,nbs,flow)

dy = zeros(23,1);


% Rename Variables 

A   = y(1); 
B   = y(2); 
K   = y(3); 
H   = y(4); 
N   = y(5); 
L   = y(6); 
C   = y(7); 
D   = y(8); 
E   = y(9); 
F   = y(10); 
L_TF   = y(11); 
V   = y(12); 
V_s   = y(13); 
X   = y(14); 
X_s   = y(15); 
P10Avail   = y(16); 
X_b   = y(17); 
II   = y(18); 
X_sbV_s   = y(19); 
P_SUB   = y(20); 
PL   = y(21); 
PL_S   = y(22); 
PL_V   = y(23); 


% Rename Kinetic Parameters 
k_0 = p(1);  
d = p(2);  
r = p(3);  
p1 = p(4);  
k_1 = p(5);  
k_6 = p(6);  
k_2 = p(7);  
k_3 = p(8);  
kon1 = p(9);  
koff1 = p(10);  
kon2 = p(11);  
koff2 = p(12);  
kon3 = p(13);  
koff3 = p(14);  
kcat = p(15);  
kon4 = p(16);  
koff4 = p(17);  
k1 = p(18);  
kflow = p(19);  
k_pla_plus = p(20);  
k_pla_minus = p(21);  
k_pla_act = p(22);  
A(II,e2P) = p(23);  
kact_e2Big = p(24);  


% Rename Binding Site 
nbs_v = nbs(1);  
nbs_x = nbs(2);  
np_x = nbs(3);  


% ODEs from reaction equations 

% A
 dy(1)  =  +  k_0  -  k_1 * A * B^2  +   0   -  k_1 * A * B^3  +  k_2 * D  -  k_1 * A  -  k1 * A * B  +  kflow * Aup  -  kflow * A;

% B
 dy(2)  =  -  d * B  -  2 * k_1 * A * B^2  -  k_6 * A * B  -  3 * k_1 * A * B^3  +  3 * k_2 * D  -  k1 * A * B;

% K
 dy(3)  =  -  r * K;

% H
 dy(4)  =  +  2 * r * K;

% N
 dy(5)  =  -  p1 * N^2;

% L
 dy(6)  =  +  p1 * N^2;

% C
 dy(7)  =  +  k_1 * A * B^2  +  k_6 * A * B  +  2 * k1 * A * B  -  kflow * C  +  kflow * Cup;

% D
 dy(8)  =  +  k_1 * A * B^3  -  k_2 * D  -  k_3 * D * E;

% E
 dy(9)  =  -  k_3 * D * E;

% F
 dy(10)  =  +  k_3 * D * E  +  k_1 * A;

% L_TF
 dy(11)  =  -  kon1 * L_TF * V  +  koff1*nbs_v * V_s  -  kon2 * L_TF * X  +  koff2*nbs_x * X_s;

% V
 dy(12)  =  -  kon1/nbs_v * L_TF * V  +  koff1 * V_s;

% V_s
 dy(13)  =  +  kon1/nbs_v * L_TF * V  -  koff1 * V_s  -  kon4 * V_s * X_s  +  koff4 * X_sbV_s;

% X
 dy(14)  =  -  kon2/nbs_x * L_TF * X  +  koff2 * X_s  -  kon3 * X * P10Avail  +  koff3 * X_b;

% X_s
 dy(15)  =  +  kon2/nbs_x * L_TF * X  -  koff2 * X_s  -  kon4 * V_s * X_s  +  koff4 * X_sbV_s;

% P10Avail
 dy(16)  =  -  kon3 * X * P10Avail  +  koff3 * X_b;

% X_b
 dy(17)  =  +  kon3 * X * P10Avail  -  koff3 * X_b  +   0 ;

% II
 dy(18)  =  +  kcat * X_b;

% X_sbV_s
 dy(19)  =  +  kon4 * V_s * X_s  -  koff4 * X_sbV_s;

% P_SUB
 dy(20)  =  -  k_pla_plus * P_SUB * PL  -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S  -  kact_e2Big * P_SUB * PL;

% PL
 dy(21)  =  -  k_pla_plus * P_SUB * PL  -  k_pla_act * PL * PL_S  -  k_pla_act * PL * PL_V  -  A(II,e2P) * PL  -  kact_e2Big * P_SUB * PL  -  kflow * PL  +  kflow * PLup;

% PL_S
 dy(22)  =  +  k_pla_plus * P_SUB * PL  +  k_pla_plus * P_SUB * PL_V  -  k_pla_minus * PL_S  +   0   +  kact_e2Big * P_SUB * PL;

% PL_V
 dy(23)  =  -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S  +  k_pla_act * PL * PL_S  +  k_pla_act * PL * PL_V  +  A(II,e2P) * PL;





end