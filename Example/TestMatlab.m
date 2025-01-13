function [time,y] = Example/TestMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py Example/Test.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
Example/TestParams

% Set the Initial Conditions
Example/TestIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
Example/TestRename
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

dy = zeros(6,1);


% Rename Variables 

A   = y(1); 
B   = y(2); 
C   = y(3); 
D   = y(4); 
E   = y(5); 
F   = y(6); 


% Rename Kinetic Parameters 
k_0 = p(1);  
d = p(2);  
k_1 = p(3);  
k_6 = p(4);  
k_2 = p(5);  
k_3 = p(6);  


% ODEs from reaction equations 

% A
 dy(1)  =  +  k_0  -  k_1 * A * B^2  +   0   -  k_1 * A * B^2  +  k_2 * D  -  k_1 * A;

% B
 dy(2)  =  -  d * B  -  k_1 * A * B^2  -  k_6 * A * B  -  k_1 * A * B^2  +  k_2 * 2 * D;

% C
 dy(3)  =  +  k_1 * A * B^2  +  k_6 * A * B;

% D
 dy(4)  =  +  k_1 * A * B^2  -  k_2 * D  -  k_3 * D * E;

% E
 dy(5)  =  -  k_3 * D * E;

% F
 dy(6)  =  +  k_3 * D * E  +  k_1 * A;





end