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
[time,y] = ode15s(@(t,y)RHS(t,y,p), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
Example/TestCoagShortRename
%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, A, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('A');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p)

dy = zeros(2,1);


% Rename Variables 

A   = y(1); 
B   = y(2); 


% Rename Kinetic Parameters 
kflow = p(1);  
k1 = p(2);  


% ODEs from reaction equations 

% A
 dy(1)  =  +  3 * kflow  +   0 ;

% B
 dy(2)  =  -  2 * k1 * A * B^2;





end