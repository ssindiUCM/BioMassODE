function [time,y] = Example/TestCoagMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py Example/TestCoag.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
Example/TestCoagParams

% Set the Initial Conditions
Example/TestCoagIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
Example/TestCoagRename
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

dy = zeros(11,1);


% Rename Variables 

L_TF   = y(1); 
V   = y(2); 
V_s   = y(3); 
X   = y(4); 
X_s   = y(5); 
P10Avail   = y(6); 
X_b   = y(7); 
X_sbV_s   = y(8); 
A   = y(9); 
B   = y(10); 
C   = y(11); 


% Rename Kinetic Parameters 
kon1 = p(1);  
kon2 = p(2);  
kon3 = p(3);  
kon4 = p(4);  
k1 = p(5);  
kflow = p(6);  


% Rename Binding Site 
nbs_v = nbs(1);  
nbs_x = nbs(2);  
np_x = nbs(3);  


% ODEs from reaction equations 

% L_TF
 dy(1)  =  -  kon1 * L_TF * V  -  kon2 * L_TF * X;

% V
 dy(2)  =  -  kon1/nbs_v * L_TF * V;

% V_s
 dy(3)  =  +  kon1/nbs_v * L_TF * V  -  kon4 * V_s * X_s;

% X
 dy(4)  =  -  kon2/nbs_x * L_TF * X  -  kon3 * X^2 * P10Avail;

% X_s
 dy(5)  =  +  kon2/nbs_x * L_TF * X  -  kon4 * V_s * X_s;

% P10Avail
 dy(6)  =  -  kon3 * X^2 * P10Avail;

% X_b
 dy(7)  =  +  kon3 * 2 * X^2 * 2 * P10Avail;

% X_sbV_s
 dy(8)  =  +  kon4 * V_s * X_s;

% A
 dy(9)  =  -  k1 * A * B  +  kflow;

% B
 dy(10)  =  -  k1 * A * B;

% C
 dy(11)  =  +  k1 * 2 * A * 2 * B  -  kflow * C;





end