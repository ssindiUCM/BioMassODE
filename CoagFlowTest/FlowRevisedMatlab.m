function [time,y] = FlowRevisedMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py CoagFlowTest/FlowRevised.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
FlowRevisedParams

% Set the Initial Conditions
FlowRevisedIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs,flowUp,otherArgs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
FlowRevisedRename
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

dy = zeros(13,1);


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
IIa_p   = y(11); 
p5   = y(12); 
V_p   = y(13); 


% Rename Kinetic Parameters 
kon_x = p(1);  
koff_x = p(2);  
kflow = p(3);  
k_pla_plus = p(4);  
k_pla_minus = p(5);  
k_pla_act = p(6);  
kact_e2 = p(7);  
kon_IIa_p = p(8);  
koff_IIa_p = p(9);  
kon_v_p = p(10);  
koff_v_p = p(11);  


% Rename Binding Site 
nbs_x = nbs(1);  


% Rename Function Arguments Site 
e2P = otherArgs(1);  
VolP = otherArgs(2);  


% Rename Flow Up 
IIa_up = flowUp(1);  
V_up = flowUp(2);  
PL_up = flowUp(3);  


% ODEs from reaction equations 

% L_TF
 dL_TF =   -  kon_x/nbs_x * L_TF * X  +  koff_x * X_st - L_TF * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% X
 dX =   -  kon_x/nbs_x * L_TF * X  +  koff_x * X_st - X * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% X_st
 dX_st =   +  kon_x/nbs_x * L_TF * X  -  koff_x * X_st - X_st * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% IIa
 dIIa =   -  kflow * IIa  +  kflow * IIa_up  -  kon_IIa_p * IIa * p2  +  koff_IIa_p * IIa_p - IIa * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% PL
 dPL =   +  kflow * PL_up  -  kflow * PL  -  k_pla_plus * PL * P_SUB  -  k_pla_act * PL * PL_S  -  k_pla_act * PL * PL_V  -  kact_e2 * A(IIa,e2P) * PL - PL * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% P_SUB
 dP_SUB =   -  k_pla_plus * PL * P_SUB  -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S - P_SUB * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% PL_S
 dPL_S =   +  k_pla_plus * PL * P_SUB  +  k_pla_plus * P_SUB * PL_V  -  k_pla_minus * PL_S  +   0  - PL_S * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% PL_V
 dPL_V =   -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S  +  k_pla_act * PL * PL_S  +  k_pla_act * PL * PL_V  +  kact_e2 * A(IIa,e2P) * PL - PL_V * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% IIa_p
 dIIa_p =   +  kon_IIa_p * IIa * p2  -  koff_IIa_p * IIa_p - IIa_p * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% V_p
 dV_p =   +  kon_v_p * V * p5  -  koff_v_p * V_p - V_p * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% V
 dV =  + dPL_S  + dPL_V  - V * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% p2
 dp2 =  + dPL_S  + dPL_V  - p2 * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

% p5
 dp5 =  + dPL_S  + dPL_V  - p5 * Dilution(VolP, P_SUB, PL, PL_S, PL_V, IIa, k_pla_act, k_pla_plus, kact_e2, e2P);

 dy = [ dL_TF, dX, dX_st, dIIa, dV, dPL, dP_SUB, dPL_S, dPL_V, dp2, dIIa_p, dp5, dV_p ]';


end

%Beginning of Helper Functions
function output = A(IIa, e2P)
    % Function: A
    % Arguments: IIa, e2P
    % Body: IIa/(e2P + IIa)
    output = IIa/(e2P + IIa);
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
