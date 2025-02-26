function [time,y] = Example/Test2Matlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py Example/Test2.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
Example/Test2Params

% Set the Initial Conditions
Example/Test2IC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
Example/Test2Rename
%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, L_TF, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('L_TF');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p,nbs)

dy = zeros(6,1);


% Rename Variables 

L_TF   = y(1); 
V   = y(2); 
V_s   = y(3); 
A   = y(4); 
B   = y(5); 
C   = y(6); 


% Rename Kinetic Parameters 
kon = p(1);  
koff = p(2);  
k1 = p(3);  


% Rename Binding Site 
nbs_v = nbs(1);  


% ODEs from reaction equations 

% L_TF
 dy(1)  =  -  kon * L_TF * V  +  koff*nbs_v * V_s;

% V
 dy(2)  =  -  kon/nbs_v * L_TF * V  +  koff * V_s;

% V_s
 dy(3)  =  +  kon/nbs_v * L_TF * V  -  koff * V_s;

% A
 dy(4)  =  -  k1 * A * B;

% B
 dy(5)  =  -  k1 * A * B;

% C
 dy(6)  =  +  k1 * A * B;





end
