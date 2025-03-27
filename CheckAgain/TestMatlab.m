function [time,y] = CheckAgain/TestMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py CheckAgain/Test.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
CheckAgain/TestParams

% Set the Initial Conditions
CheckAgain/TestIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
CheckAgain/TestRename
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

dy = zeros(10,1);


% Rename Variables 

A   = y(1); 
B   = y(2); 
K   = y(3); 
H   = y(4); 
N   = y(5); 
L   = y(6); 
C   = y(7); 
D   = y(8); 
E   = y(9); 
F   = y(10); 


% Rename Kinetic Parameters 
k_0 = p(1);  
d = p(2);  
r = p(3);  
p1 = p(4);  
k_1 = p(5);  
k_6 = p(6);  
k_2 = p(7);  
k_3 = p(8);  


% ODEs from reaction equations 

% A
 dy(1)  =  +  k_0  -  k_1 * A * B^2  +   0   -  k_1 * A * B^3  +  k_2 * D  -  k_1 * A;

% B
 dy(2)  =  -  d * B  -  2 * k_1 * A * B^2  -  k_6 * A * B  -  3 * k_1 * A * B^3  +  3 * k_2 * D;

% K
 dy(3)  =  -  r * K;

% H
 dy(4)  =  +  2 * r * K;

% N
 dy(5)  =  -  p1 * N^2;

% L
 dy(6)  =  +  p1 * N^2;

% C
 dy(7)  =  +  k_1 * A * B^2  +  k_6 * A * B;

% D
 dy(8)  =  +  k_1 * A * B^3  -  k_2 * D  -  k_3 * D * E;

% E
 dy(9)  =  -  k_3 * D * E;

% F
 dy(10)  =  +  k_3 * D * E  +  k_1 * A;





end