function [time,y] = Example/LakshmananMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py Example/Lakshmanan.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
Example/LakshmananParams

% Set the Initial Conditions
Example/LakshmananIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
Example/LakshmananRename
%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, TF, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('TF');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p)

dy = zeros(43,1);


% Rename Variables 

TF   = y(1); 
VII   = y(2); 
TFbVII   = y(3); 
VIIa   = y(4); 
TFbVIIa   = y(5); 
Xa   = y(6); 
IIa   = y(7); 
X   = y(8); 
TFbVIIabX   = y(9); 
TFbVIIabXa   = y(10); 
IX   = y(11); 
TFbVIIabIX   = y(12); 
TFbVIIabiX   = y(13); 
IXa   = y(14); 
II   = y(15); 
VIII   = y(16); 
VIIIa   = y(17); 
IXabVIIIa   = y(18); 
IXabVIIIabX   = y(19); 
VIIIa1   = y(20); 
VIIIa2   = y(21); 
V   = y(22); 
Va   = y(23); 
XabVa   = y(24); 
XabVabII   = y(25); 
mIIa   = y(26); 
XI   = y(27); 
IIabXI   = y(28); 
XIa   = y(29); 
XIabIX   = y(30); 
IXabX   = y(31); 
TFPI   = y(32); 
XabTFPI   = y(33); 
TFbVIIabXabTFPI   = y(34); 
AT   = y(35); 
XabAT   = y(36); 
mIIabAT   = y(37); 
IXabAT   = y(38); 
IIabAT   = y(39); 
TFbVIIabAT   = y(40); 
XIabAT   = y(41); 
C1INH   = y(42); 
XIbC1INH   = y(43); 


% Rename Kinetic Parameters 
ka_tf_vii = p(1);  
kd_tfvii = p(2);  
ka_tf_viia = p(3);  
kd_tfviia = p(4);  
ka_tfviia_vii = p(5);  
ka_tfviia_tfvii = p(6);  
ka_tfvii_xa = p(7);  
ka_tfvii_iia = p(8);  
ka_tfviia_x = p(9);  
kd_tfviiax = p(10);  
ka_tfviiax = p(11);  
ka_tfviia_xa = p(12);  
kd_tfviiaxa = p(13);  
ka_tfviia_ix = p(14);  
kd_tfviiaix = p(15);  
ka_tfviiaix = p(16);  
ka_xa_vii = p(17);  
ka_iia_vii = p(18);  
ka_xa_ii = p(19);  
ka_iia_viii = p(20);  
ka_viiia_ixa = p(21);  
kd_ixaviiia = p(22);  
ka_ixaviiia_x = p(23);  
kd_ixaviiiax = p(24);  
ka_ixaviiiax = p(25);  
ka_viiia = p(26);  
kd_viiia1_viiia2 = p(27);  
ka_ixaviiiax_2 = p(28);  
ka_ixaviiia = p(29);  
ka_iia_v = p(30);  
ka_xa_va = p(31);  
kd_xava = p(32);  
ka_xava_ii = p(33);  
kd_xavaii = p(34);  
ka_xavaii = p(35);  
ka_miia_xava = p(36);  
ka_iia_xi = p(37);  
kd_iiaxi = p(38);  
ka_iiaxi = p(39);  
ka_xia_ix = p(40);  
kd_xiaix = p(41);  
ka_xiaix = p(42);  
ka_ixa_x = p(43);  
kd_ixax = p(44);  
ka_ixax = p(45);  
ka_xa_tfpi = p(46);  
kd_xatfpi = p(47);  
ka_tfviiaxa_tfpi = p(48);  
kd_tfviiaxatfpi = p(49);  
ka_tfviia_xatfpi = p(50);  
ka_xa_at = p(51);  
ka_miia_at = p(52);  
ka_ixa_at = p(53);  
ka_iia_at = p(54);  
ka_tfviia_at = p(55);  
ka_xia_at = p(56);  
ka_xia_c1inh = p(57);  


% ODEs from reaction equations 

% TF
 dy(1)  =  -  ka_tf_vii * TF * VII  +  kd_tfvii * TFbVII  -  ka_tf_viia * TF * VIIa  +  kd_tfviia * TFbVIIa;

% VII
 dy(2)  =  -  ka_tf_vii * TF * VII  +  kd_tfvii * TFbVII  -  ka_tfviia_vii * VII * TFbVIIa  -  ka_xa_vii * VII * Xa  -  ka_iia_vii * VII * IIa;

% TFbVII
 dy(3)  =  +  ka_tf_vii * TF * VII  -  kd_tfvii * TFbVII  -  ka_tfviia_tfvii * TFbVII * TFbVIIa  -  ka_tfvii_xa * TFbVII * Xa  -  ka_tfvii_iia * TFbVII * IIa;

% VIIa
 dy(4)  =  -  ka_tf_viia * TF * VIIa  +  kd_tfviia * TFbVIIa  +  ka_tfviia_vii * VII * TFbVIIa  +  ka_xa_vii * VII * Xa  +  ka_iia_vii * VII * IIa;

% TFbVIIa
 dy(5)  =  +  ka_tf_viia * TF * VIIa  -  kd_tfviia * TFbVIIa  +   0   +  ka_tfviia_tfvii * TFbVII  +  ka_tfvii_xa * TFbVII * Xa  +  ka_tfvii_iia * TFbVII * IIa  -  ka_tfviia_x * TFbVIIa * X  +  kd_tfviiax * TFbVIIabX  -  ka_tfviia_xa * TFbVIIa * Xa  +  kd_tfviiaxa * TFbVIIabXa  -  ka_tfviia_ix * TFbVIIa * IX  +  kd_tfviiaix * TFbVIIabIX  +  ka_tfviiaix * TFbVIIabiX  -  ka_tfviia_xatfpi * TFbVIIa * XabTFPI  -  ka_tfviia_at * TFbVIIa * AT;

% Xa
 dy(6)  =  +   0   -  ka_tfviia_xa * TFbVIIa * Xa  +  kd_tfviiaxa * TFbVIIabXa  +   0   +   0   +  ka_ixaviiiax * IXabVIIIabX  -  ka_xa_va * Xa * Va  +  kd_xava * XabVa  +  ka_ixax * IXabX  -  ka_xa_tfpi * Xa * TFPI  +  kd_xatfpi * XabTFPI  -  ka_xa_at * Xa * AT;

% IIa
 dy(7)  =  +   0   +   0   +  ka_xa_ii * Xa * II  +   0   +   0   +  ka_miia_xava * XabVa * mIIa  -  ka_iia_xi * IIa * XI  +  kd_iiaxi * IIabXI  +  ka_iiaxi * IIabXI  -  ka_iia_at * IIa * AT;

% X
 dy(8)  =  -  ka_tfviia_x * TFbVIIa * X  +  kd_tfviiax * TFbVIIabX  -  ka_ixaviiia_x * X * IXabVIIIa  +  kd_ixaviiiax * IXabVIIIabX  +  ka_ixaviiiax_2 * IXabVIIIabX  -  ka_ixa_x * X * IXa  +  kd_ixax * IXabX;

% TFbVIIabX
 dy(9)  =  +  ka_tfviia_x * TFbVIIa * X  -  kd_tfviiax * TFbVIIabX  -  ka_tfviiax * TFbVIIabX;

% TFbVIIabXa
 dy(10)  =  +  ka_tfviiax * TFbVIIabX  +  ka_tfviia_xa * TFbVIIa * Xa  -  kd_tfviiaxa * TFbVIIabXa  -  ka_tfviiaxa_tfpi * TFbVIIabXa * TFPI  +  kd_tfviiaxatfpi * TFbVIIabXabTFPI;

% IX
 dy(11)  =  -  ka_tfviia_ix * TFbVIIa * IX  +  kd_tfviiaix * TFbVIIabIX  +  ka_ixaviiiax_2 * IXabVIIIabX  -  ka_xia_ix * IX * XIa  +  kd_xiaix * XIabIX;

% TFbVIIabIX
 dy(12)  =  +  ka_tfviia_ix * TFbVIIa * IX  -  kd_tfviiaix * TFbVIIabIX;

% TFbVIIabiX
 dy(13)  =  -  ka_tfviiaix * TFbVIIabiX;

% IXa
 dy(14)  =  +  ka_tfviiaix * TFbVIIabiX  -  ka_viiia_ixa * IXa * VIIIa  +  kd_ixaviiia * IXabVIIIa  +  ka_ixaviiia * IXabVIIIa  +  ka_xiaix * XIabIX  -  ka_ixa_x * X * IXa  +  kd_ixax * IXabX  +  ka_ixax * IXabX  -  ka_ixa_at * IXa * AT;

% II
 dy(15)  =  -  ka_xa_ii * Xa * II  -  ka_xava_ii * II * XabVa  +  kd_xavaii * XabVabII;

% VIII
 dy(16)  =  -  ka_iia_viii * IIa * VIII;

% VIIIa
 dy(17)  =  +  ka_iia_viii * IIa * VIII  -  ka_viiia_ixa * IXa * VIIIa  +  kd_ixaviiia * IXabVIIIa  -  ka_viiia * VIIIa  +  kd_viiia1_viiia2 * VIIIa1 * VIIIa2;

% IXabVIIIa
 dy(18)  =  +  ka_viiia_ixa * IXa * VIIIa  -  kd_ixaviiia * IXabVIIIa  -  ka_ixaviiia_x * X * IXabVIIIa  +  kd_ixaviiiax * IXabVIIIabX  +  ka_ixaviiiax * IXabVIIIabX  -  ka_ixaviiia * IXabVIIIa;

% IXabVIIIabX
 dy(19)  =  +  ka_ixaviiia_x * X * IXabVIIIa  -  kd_ixaviiiax * IXabVIIIabX  -  ka_ixaviiiax * IXabVIIIabX  -  ka_ixaviiiax_2 * IXabVIIIabX;

% VIIIa1
 dy(20)  =  +  ka_viiia * VIIIa  -  kd_viiia1_viiia2 * VIIIa1 * VIIIa2  +  ka_ixaviiiax_2 * IXabVIIIabX  +  ka_ixaviiia * IXabVIIIa;

% VIIIa2
 dy(21)  =  +  ka_viiia * VIIIa  -  kd_viiia1_viiia2 * VIIIa1 * VIIIa2  +  ka_ixaviiiax_2 * IXabVIIIabX  +  ka_ixaviiia * IXabVIIIa;

% V
 dy(22)  =  -  ka_iia_v * IIa * V;

% Va
 dy(23)  =  +  ka_iia_v * IIa * V  -  ka_xa_va * Xa * Va  +  kd_xava * XabVa;

% XabVa
 dy(24)  =  +  ka_xa_va * Xa * Va  -  kd_xava * XabVa  -  ka_xava_ii * II * XabVa  +  kd_xavaii * XabVabII  +  ka_xavaii * XabVabII  +   0 ;

% XabVabII
 dy(25)  =  +  ka_xava_ii * II * XabVa  -  kd_xavaii * XabVabII  -  ka_xavaii * XabVabII;

% mIIa
 dy(26)  =  +  ka_xavaii * XabVabII  -  ka_miia_xava * XabVa * mIIa  -  ka_miia_at * mIIa * AT;

% XI
 dy(27)  =  -  ka_iia_xi * IIa * XI  +  kd_iiaxi * IIabXI;

% IIabXI
 dy(28)  =  +  ka_iia_xi * IIa * XI  -  kd_iiaxi * IIabXI  -  ka_iiaxi * IIabXI;

% XIa
 dy(29)  =  +  ka_iiaxi * IIabXI  -  ka_xia_ix * IX * XIa  +  kd_xiaix * XIabIX  +  ka_xiaix * XIabIX  -  ka_xia_at * XIa * AT  -  ka_xia_c1inh * XIa * C1INH;

% XIabIX
 dy(30)  =  +  ka_xia_ix * IX * XIa  -  kd_xiaix * XIabIX  -  ka_xiaix * XIabIX;

% IXabX
 dy(31)  =  +  ka_ixa_x * X * IXa  -  kd_ixax * IXabX  -  ka_ixax * IXabX;

% TFPI
 dy(32)  =  -  ka_xa_tfpi * Xa * TFPI  +  kd_xatfpi * XabTFPI  -  ka_tfviiaxa_tfpi * TFbVIIabXa * TFPI  +  kd_tfviiaxatfpi * TFbVIIabXabTFPI;

% XabTFPI
 dy(33)  =  +  ka_xa_tfpi * Xa * TFPI  -  kd_xatfpi * XabTFPI  -  ka_tfviia_xatfpi * TFbVIIa * XabTFPI;

% TFbVIIabXabTFPI
 dy(34)  =  +  ka_tfviiaxa_tfpi * TFbVIIabXa * TFPI  -  kd_tfviiaxatfpi * TFbVIIabXabTFPI  +  ka_tfviia_xatfpi * TFbVIIa * XabTFPI;

% AT
 dy(35)  =  -  ka_xa_at * Xa * AT  -  ka_miia_at * mIIa * AT  -  ka_ixa_at * IXa * AT  -  ka_iia_at * IIa * AT  -  ka_tfviia_at * TFbVIIa * AT  -  ka_xia_at * XIa * AT;

% XabAT
 dy(36)  =  +  ka_xa_at * Xa * AT;

% mIIabAT
 dy(37)  =  +  ka_miia_at * mIIa * AT;

% IXabAT
 dy(38)  =  +  ka_ixa_at * IXa * AT;

% IIabAT
 dy(39)  =  +  ka_iia_at * IIa * AT;

% TFbVIIabAT
 dy(40)  =  +  ka_tfviia_at * TFbVIIa * AT;

% XIabAT
 dy(41)  =  +  ka_xia_at * XIa * AT;

% C1INH
 dy(42)  =  -  ka_xia_c1inh * XIa * C1INH;

% XIbC1INH
 dy(43)  =  +  ka_xia_c1inh * XIa * C1INH;





end