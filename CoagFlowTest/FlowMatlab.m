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
[time,y] = ode15s(@(t,y)RHS(t,y,p,flowUp), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
CoagFlowTest/FlowRename
%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, IIa, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('IIa');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p,flowUp)

dy = zeros(10,1);


% Rename Variables 

IIa   = y(1); 
V   = y(2); 
PL   = y(3); 
P_SUB   = y(4); 
PL_S   = y(5); 
p5avail   = y(6); 
p2avail   = y(7); 
PL_V   = y(8); 
IIa_p   = y(9); 
V_p   = y(10); 


% Rename Kinetic Parameters 
kflow = p(1);  
k_pla_plus = p(2);  
k_pla_minus = p(3);  
k_pla_act = p(4);  
A(IIa,e2P) = p(5);  
kon_IIa_p = p(6);  
koff_IIa_p = p(7);  
kon_v_p = p(8);  
koff_v_p = p(9);  


% ODEs from reaction equations 

