function [time,y] = Example/TestCoagShortMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py Example/TestCoagShort.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
Example/TestCoagShortParams

% Set the Initial Conditions
Example/TestCoagShortIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
Example/TestCoagShortRename
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

dy = zeros(3,1);


% Rename Variables 

L_TF   = y(1); 
A   = y(2); 
A_st   = y(3); 


% Rename Kinetic Parameters 
kon = p(1);  
koff = p(2);  
kflow = p(3);  


% Rename Binding Site 
nbs_a = nbs(1);  


% ODEs from reaction equations 

% L_TF
 dy(1)  =  -  kon * L_TF * A  +  koff*nbs_a * A_st;

% A
 dy(2)  =  -  kon/nbs_a * L_TF * A  +  koff * A_st  +  kflow;

% A_st
 dy(3)  =  +  kon/nbs_a * L_TF * A  -  koff * A_st;





end