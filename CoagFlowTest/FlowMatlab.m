function [time,y] = CoagFlowTest/FlowMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py CoagFlowTest/Flow.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
CoagFlowTest/FlowParams

% Set the Initial Conditions
CoagFlowTest/FlowIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs,flowUp,otherArgs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
CoagFlowTest/FlowRename
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
p5avail   = y(9); 
p2avail   = y(10); 
PL_V   = y(11); 
IIa_p   = y(12); 
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
np_ii = nbs(2);  
np_v = nbs(3);  


% Rename Function Arguments Site 
e2P = otherArgs(1);  


% Rename Flow Up 
IIa_up = flowUp(1);  
V_up = flowUp(2);  
PL_up = flowUp(3);  


% ODEs from reaction equations 

% L_TF
 dy(1)  =  -  kon_x/nbs_x * L_TF * X  +  koff_x * X_st;

% X
 dy(2)  =  -  kon_x/nbs_x * L_TF * X  +  koff_x * X_st;

% X_st
 dy(3)  =  +  kon_x/nbs_x * L_TF * X  -  koff_x * X_st;

% IIa
 dy(4)  =  -  kflow * IIa  +  kflow * IIa_up  -  kon_IIa_p * IIa * p2avail  +  koff_IIa_p * IIa_p;

% V
 dy(5)  =  -  kflow * V  +  kflow * V_up  -  kon_v_p * V * p5avail  +  koff_v_p * V_p;

% PL
 dy(6)  =  +  kflow * PL_up  -  kflow * PL  -  k_pla_plus * PL * P_SUB  -  k_pla_act * PL * PL_S  -  k_pla_act * PL * PL_V  -  kact_e2 * A(IIa,e2P) * PL;

% P_SUB
 dy(7)  =  -  k_pla_plus * PL * P_SUB  -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S;

% PL_S
 dy(8)  =  +  k_pla_plus * PL * P_SUB  +  k_pla_plus * P_SUB * PL_V  -  k_pla_minus * PL_S  +   0 ;

% p5avail
 dy(9)  =  +  k_pla_plus * np_v * PL * P_SUB  +  k_pla_act * np_v * PL * PL_S  +  k_pla_act * np_v * PL * PL_V  +  kact_e2 * A(IIa,e2P) * np_v * PL  -  kon_v_p * V * p5avail  +  koff_v_p * V_p;

% p2avail
 dy(10)  =  +  k_pla_plus * np_ii * PL * P_SUB  +  k_pla_act * np_v * PL * PL_S  +  k_pla_act * np_v * PL * PL_V  +  kact_e2 * A(IIa,e2P) * np_ii * PL  -  kon_IIa_p * IIa * p2avail  +  koff_IIa_p * IIa_p;

% PL_V
 dy(11)  =  -  k_pla_plus * P_SUB * PL_V  +  k_pla_minus * PL_S  +  k_pla_act * PL * PL_S  +  k_pla_act * PL * PL_V  +  kact_e2 * A(IIa,e2P) * PL;

% IIa_p
 dy(12)  =  +  kon_IIa_p * IIa * p2avail  -  koff_IIa_p * IIa_p;

% V_p
 dy(13)  =  +  kon_v_p * V * p5avail  -  koff_v_p * V_p;





end

%Beginning of Helper Functions
function output = A(x, e2P)
    % Function: A
    % Arguments: x, e2P
    % Body: x/(e2P + x)
    output = x/(e2P + x);
end


%End of Helper Functions
