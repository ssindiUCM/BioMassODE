function [time,y] = Test2Matlab(t_start,t_final)
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
Test2Params

% Set the Initial Conditions
Test2IC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
Test2Rename
%  
% Place plots or other calculations here
%   
% Example: 
plot(time, L_TF, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('L_TF');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p,nbs)

dy = zeros(9,1);


% Rename Variables 

L_TF   = y(1); 
V   = y(2); 
V_s   = y(3); 
X   = y(4); 
X_s   = y(5); 
X_sbV_s   = y(6); 
A   = y(7); 
B   = y(8); 
C   = y(9); 


% Rename Kinetic Parameters 
kon1 = p(1);  
koff1 = p(2);  
kon2 = p(3);  
koff2 = p(4);  
kon3 = p(5);  
koff3 = p(6);  
k1 = p(7);  


% Rename Binding Site 
nbs_v = nbs(1);  
nbs_x = nbs(2);  


% ODEs from reaction M

% L_TF
 dy(1)  =  -  kon1 * L_TF * V  +  koff1*nbs_v * V_s  -  kon2 * L_TF * X  +  koff2*nbs_x * X_s;

% V
 dy(2)  =  -  kon1/nbs_v * L_TF * V  +  koff1 * V_s;

% V_s
 dy(3)  =  +  kon1/nbs_v * L_TF * V  -  koff1 * V_s  -  kon3 * V_s * X_s  +  koff3 * X_sbV_s;

% X
 dy(4)  =  -  kon2/nbs_x * L_TF * X  +  koff2 * X_s;

% X_s
 dy(5)  =  +  kon2/nbs_x * L_TF * X  -  koff2 * X_s  -  kon3 * V_s * X_s  +  koff3 * X_sbV_s;

% X_sbV_s
 dy(6)  =  +  kon3 * V_s * X_s  -  koff3 * X_sbV_s;

% A
 dy(7)  =  -  k1 * A * B;

% B
 dy(8)  =  -  k1 * A * B;

% C
 dy(9)  =  +  k1 * A * B;





end
