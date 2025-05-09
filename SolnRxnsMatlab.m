function [time,y] = SolnRxnsMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py SolnRxns.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
SolnRxnsParams

% Set the Initial Conditions
SolnRxnsIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
SolnRxnsRename
%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, IXa, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('IXa');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p)

dy = zeros(9,1);


% Rename Variables 

IXa   = y(1); 
M   = y(2); 
IXabM   = y(3); 
X   = y(4); 
XbM   = y(5); 
Xa   = y(6); 
XabM   = y(7); 
IXabMbX   = y(8); 
IXabMbXa   = y(9); 


% Rename Kinetic Parameters 
kpixam = p(1);  
kmixam = p(2);  
kpxm = p(3);  
kmxm = p(4);  
kpxam = p(5);  
kmxam = p(6);  
kpixamx = p(7);  
kmixamx = p(8);  
kcat = p(9);  


% ODEs from reaction equations 

% IXa
 dIXa =   -  kpixam * IXa * M  +  kmixam * IXabM  -  kpixam * IXa * XbM  +  kmixam * IXabMbX  +  kpixam * IXabMbXa  -  kmixam * IXa * XabM;

% M
 dM =   -  kpixam * IXa * M  +  kmixam * IXabM  -  kpxm * M * X  +  kmxm * XbM  -  kpxam * M * Xa  +  kmxam * XabM;

% IXabM
 dIXabM =   +  kpixam * IXa * M  -  kmixam * IXabM  -  kpixamx * IXabM * X  +  kmixamx * IXabMbX  +  kpxam * IXabMbXa  -  kmxam * IXabM * Xa;

% X
 dX =   -  kpxm * M * X  +  kmxm * XbM  -  kpixamx * IXabM * X  +  kmixamx * IXabMbX;

% XbM
 dXbM =   +  kpxm * M * X  -  kmxm * XbM  -  kpixam * IXa * XbM  +  kmixam * IXabMbX;

% Xa
 dXa =   -  kpxam * M * Xa  +  kmxam * XabM  +  kpxam * IXabMbXa  -  kmxam * IXabM * Xa;

% XabM
 dXabM =   +  kpxam * M * Xa  -  kmxam * XabM  +  kpixam * IXabMbXa  -  kmixam * IXa * XabM;

% IXabMbX
 dIXabMbX =   +  kpixamx * IXabM * X  -  kmixamx * IXabMbX  +  kpixam * IXa * XbM  -  kmixam * IXabMbX  -  kcat * IXabMbX;

% IXabMbXa
 dIXabMbXa =   +  kcat * IXabMbX  -  kpxam * IXabMbXa  +  kmxam * IXabM * Xa  -  kpixam * IXabMbXa  +  kmixam * IXa * XabM;

 dy = [ dIXa, dM, dIXabM, dX, dXbM, dXa, dXabM, dIXabMbX, dIXabMbXa ]';


end

%Beginning of Helper Functions
%End of Helper Functions
