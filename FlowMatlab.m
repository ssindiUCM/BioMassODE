function [time,y] = FlowMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py Flow.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
FlowParams

% Set the Initial Conditions
FlowIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs,flowUp,otherArgs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
FlowRename
%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, L_TF, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('L_TF');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p,nbs,flowUp,otherArgs)

dy = zeros(19,1);


% Rename Variables 

L_TF   = y(1); 
X   = y(2); 
X_st   = y(3); 
IIa   = y(4); 
V   = y(5); 
PL   = y(6); 
P_SUB   = y(7); 
PL_S   = y(8); 
PL_V   = y(9); 
p2   = y(10); 
IIa_sp   = y(11); 
p5   = y(12); 
V_sp   = y(13); 
TF   = y(14); 
VII   = y(15); 
TFbVII_se   = y(16); 
XI   = y(17); 
TFPI   = y(18); 
Vh   = y(19); 


% Rename Kinetic Parameters 
kon_x = p(1);  
koff_x = p(2);  
kflow = p(3);  
k_pla_plus = p(4);  
k_pla_minus = p(5);  
k_pla_act = p(6);  
kact_e2 = p(7);  
kon_IIa_sp = p(8);  
koff_IIa_sp = p(9);  
kon_v_sp = p(10);  
koff_v_sp = p(11);  
kon_se = p(12);  


% Rename Binding Site 
nbs_x = nbs(1);  
np2 = nbs(2);  
np5 = nbs(3);  
nv = nbs(4);  
nhv = nbs(5);  


% Rename Function Arguments Site 
e2P = otherArgs(1);  
VolP = otherArgs(2);  


% Rename Flow Up 
IIa_up = flowUp(1);  
V_up = flowUp(2);  
PL_up = flowUp(3);  
X_up = flowUp(4);  
XI_up = flowUp(5);  
TFPI_up = flowUp(6);  


% ODEs from reaction equations 

% L_TF
 dL_TF =   -  kon_x/nbs_x * L_TF * X  +  koff_x * X_st - L_TF * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% X
 dX =   -  kon_x/nbs_x * L_TF * X  +  koff_x * X_st  -  kflow * X  +  kflow * X_up - X * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% X_st
 dX_st =   +  kon_x/nbs_x * L_TF * X  -  koff_x * X_st - X_st * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% IIa
 dIIa =   -  kflow * IIa  +  kflow * IIa_up  -  kon_IIa_sp * IIa * p2  +  koff_IIa_sp * IIa_sp - IIa * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% PL
 dPL =   +  kflow * PL_up  -  kflow * PL  -  k_pla_plus * PL * P_SUB  -  k_pla_act * PL * PL_S  -  k_pla_act * PL * PL_V  -  kact_e2 * B(IIa,e2P) * PL - PL * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% P_SUB
 dP_SUB =   -  k_pla_plus * PL * P_SUB  -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S - P_SUB * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% PL_S
 dPL_S =   +  k_pla_plus * PL * P_SUB  +  k_pla_plus * P_SUB * PL_V  -  k_pla_minus * PL_S  +   0  - PL_S * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% PL_V
 dPL_V =   -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S  +  k_pla_act * PL * PL_S  +  k_pla_act * PL * PL_V  +  kact_e2 * B(IIa,e2P) * PL - PL_V * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% IIa_sp
 dIIa_sp =   +  kon_IIa_sp * IIa * p2  -  koff_IIa_sp * IIa_sp - IIa_sp * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% V_sp
 dV_sp =   +  kon_v_sp * V * p5  -  koff_v_sp * V_sp - V_sp * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% TF
 dTF =   -  kon_se * TF * VII - TF * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% VII
 dVII =   -  kon_se * TF * VII - VII * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% TFbVII_se
 dTFbVII_se =   +  kon_se * TF * VII - TFbVII_se * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% XI
 dXI =   -  kflow * XI  +  kflow * XI_up - XI * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% TFPI
 dTFPI =   -  kflow * TFPI  +  kflow * TFPI_up - TFPI * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% V
 dV =   -  kflow * V  +  kflow * V_up  -  kon_v_sp * V * p5  +  koff_v_sp * V_sp + nv * dPL_S  + nv * dPL_V  - V * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% p2
 dp2 =   -  kon_IIa_sp * IIa * p2  +  koff_IIa_sp * IIa_sp + np2 * dPL_S  + np2 * dPL_V  - p2 * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% p5
 dp5 =   -  kon_v_sp * V * p5  +  koff_v_sp * V_sp + np5 * dPL_S  + np5 * dPL_V  - p5 * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% Vh
 dVh =  + nhv * dPL_S  + nhv * dPL_V  - Vh * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

 dy = [ dL_TF, dX, dX_st, dIIa, dV, dPL, dP_SUB, dPL_S, dPL_V, dp2, dIIa_sp, dp5, dV_sp, dTF, dVII, dTFbVII_se, dXI, dTFPI, dVh ]';


end

%Beginning of Helper Functions
function output = A(IIa, e2P)
    % Function: A
    % Arguments: IIa, e2P
    % Body: IIa/(e2P + IIa)
    output = IIa/(e2P + IIa);
end


function output = B(x, e2P)
    % Function: B
    % Arguments: x, e2P
    % Body: x/(e2P+x)
    output = x/(e2P+x);
end


function output = Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P)
    % Function: Dilution
    % Arguments: VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P
    % Body: (VolP)/((1-VolP)*(PL_S + PL_V))*dPdt(P_SUB,PL, PL_S, PL_V,IIa,k_pla_act,k_pla_plus,kact_e2,e2P)
    output = (VolP)/((1-VolP)*(PL_S + PL_V))*dPdt(P_SUB,PL, PL_S, PL_V,IIa,k_pla_act,k_pla_plus,kact_e2,e2P);
end


function output = dPdt(P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P)
    % Function: dPdt
    % Arguments: P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P
    % Body: +  k_pla_act * PL * PL_S  + k_pla_act * PL * PL_V  +  kact_e2 * A(IIa,e2P) * PL +  k_pla_plus * PL * P_SUB
    output = +  k_pla_act * PL * PL_S  + k_pla_act * PL * PL_V  +  kact_e2 * A(IIa,e2P) * PL +  k_pla_plus * PL * P_SUB;
end


%End of Helper Functions
