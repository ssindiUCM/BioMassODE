function [time,y] = StaticCoag/StaticCoagWithBSMatlab(t_start,t_final)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 30*60 
%
%
% This file was created by issuing command: 
%     python createMatlabFile.py StaticCoag/StaticCoagWithBS.txt
%

if nargin == 0
     t_start = 0;     % Default start time is 0
     t_final = 30*60; % Default final time is 30*60
elseif nargin~=2
   disp('Need to Specify t_start, t_end')
   return
end


% Set the Kinetic Parameters
StaticCoag/StaticCoagWithBSParams

% Set the Initial Conditions
StaticCoag/StaticCoagWithBSIC

options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------- Main Solve ----------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs), t_start:1:t_final, init_cond, options);
%-----------------------------------------------------------%


% Rename solution components
StaticCoag/StaticCoagWithBSRename
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

dy = zeros(153,1);


% Rename Variables 

L_TF   = y(1); 
II   = y(2); 
II_st   = y(3); 
L_noTF   = y(4); 
II_s   = y(5); 
IIf   = y(6); 
IIf_st   = y(7); 
IIf_s   = y(8); 
V   = y(9); 
V_st   = y(10); 
V_s   = y(11); 
Vh   = y(12); 
Vh_st   = y(13); 
Vh_s   = y(14); 
Va   = y(15); 
Va_st   = y(16); 
Va_s   = y(17); 
VII   = y(18); 
VII_st   = y(19); 
VII_s   = y(20); 
VIIa   = y(21); 
VIIa_st   = y(22); 
VIIa_s   = y(23); 
VIII   = y(24); 
VIII_st   = y(25); 
VIII_s   = y(26); 
VIIIa   = y(27); 
VIIIa_st   = y(28); 
VIIIa_s   = y(29); 
IX   = y(30); 
IX_st   = y(31); 
IX_s   = y(32); 
IXa   = y(33); 
IXa_st   = y(34); 
IXa_s   = y(35); 
X   = y(36); 
X_st   = y(37); 
X_s   = y(38); 
Xa   = y(39); 
Xa_st   = y(40); 
Xa_s   = y(41); 
XabTFPI   = y(42); 
Xa_stbTFPI   = y(43); 
Xa_sbTFPI   = y(44); 
TFPIbVh   = y(45); 
TFPIbVh_st   = y(46); 
TFPIbVh_s   = y(47); 
IXabAT   = y(48); 
IXa_sbAT   = y(49); 
IXa_stbAT   = y(50); 
XabAT   = y(51); 
Xa_sbAT   = y(52); 
Xa_stbAT   = y(53); 
PS   = y(54); 
PS_st   = y(55); 
PS_s   = y(56); 
TF   = y(57); 
TFbVII_st   = y(58); 
TFbVIIa_st   = y(59); 
IX_stbTFbVIIa_st   = y(60); 
TFbVIIa_stbX_st   = y(61); 
TFbVIIa_stbXa_st   = y(62); 
TFbVII_stbXa_st   = y(63); 
TFPI   = y(64); 
TFbVIIa_stbXa_stbTFPI   = y(65); 
Xa_stbTFPIbTFbVIIa_st   = y(66); 
Xa_stbV_st   = y(67); 
Xa_stbVII_st   = y(68); 
Xa_stbVIII_st   = y(69); 
VIIIa_stbIXa_st   = y(70); 
X_stbVIIIa_stbIXa_st   = y(71); 
Xa_stbVh_st   = y(72); 
II_stbXa_stbVh_st   = y(73); 
IIa   = y(74); 
IIabV_st   = y(75); 
IIabVh_st   = y(76); 
IIabVIII_st   = y(77); 
XIa   = y(78); 
XIabIX_st   = y(79); 
Xa_stbVa_st   = y(80); 
II_stbXa_stbVa_st   = y(81); 
IXa_stbX_st   = y(82); 
Xa_sbV_s   = y(83); 
Xa_sbVII_s   = y(84); 
Xa_sbVIII_s   = y(85); 
VIIIa_sbIXa_s   = y(86); 
X_sbVIIIa_sbIXa_s   = y(87); 
Xa_sbVh_s   = y(88); 
II_sbXa_sbVh_s   = y(89); 
IIabV_s   = y(90); 
IIabVh_s   = y(91); 
IIabVIII_s   = y(92); 
XIabIX_s   = y(93); 
Xa_sbVa_s   = y(94); 
II_sbXa_sbVa_s   = y(95); 
IXa_sbX_s   = y(96); 
XabV   = y(97); 
XabVII   = y(98); 
XabVIII   = y(99); 
IIabV   = y(100); 
IIabVh   = y(101); 
IIabVIII   = y(102); 
XIabIX   = y(103); 
XI   = y(104); 
IIabXI   = y(105); 
IXabX   = y(106); 
TFPIbXa_sbVh_s   = y(107); 
TFPIbXa_stbVh_st   = y(108); 
Xa_sbVh_sbTFPI   = y(109); 
Xa_stbVh_stbTFPI   = y(110); 
VhbTFPIbXa   = y(111); 
Vh_sbTFPIbXa   = y(112); 
Vh_stbTFPIbXa   = y(113); 
Vh_sbTFPIbXa_s   = y(114); 
Vh_stbTFPIbXa_st   = y(115); 
VhbTFPIbXa_s   = y(116); 
VhbTFPIbXa_st   = y(117); 
AT   = y(118); 
IIabAT   = y(119); 
XIabAT   = y(120); 
XII   = y(121); 
SN   = y(122); 
XII_sn   = y(123); 
XIIa   = y(124); 
XIIa_sn   = y(125); 
XIIa_snbXII_sn   = y(126); 
PK   = y(127); 
XIIabPK   = y(128); 
PKa   = y(129); 
PKabXII   = y(130); 
XIIabXI   = y(131); 
PKabIX   = y(132); 
XIIa_snbAT   = y(133); 
XIIabAT   = y(134); 
C1INH   = y(135); 
XIabC1INH   = y(136); 
XIIabC1INH   = y(137); 
XIIa_snbC1INH   = y(138); 
PKabC1INH   = y(139); 
PKabIX_s   = y(140); 
HK   = y(141); 
XIbHK   = y(142); 
XIbHK_sn   = y(143); 
XIbHK_snbXIIa_sn   = y(144); 
XIabHK_sn   = y(145); 
PKbHK   = y(146); 
PKbHK_sn   = y(147); 
PKbHK_snbXIIa_sn   = y(148); 
PKabHK_sn   = y(149); 
HKa   = y(150); 
BK   = y(151); 
PKabHK   = y(152); 
PKabHK_snbXII_sn   = y(153); 


% Rename Kinetic Parameters 
kon_ii = p(1);  
koff_ii = p(2);  
kon_iif = p(3);  
koff_iif = p(4);  
kon_v = p(5);  
koff_v = p(6);  
kon_vh = p(7);  
koff_vh = p(8);  
kon_va = p(9);  
koff_va = p(10);  
kon_vii = p(11);  
koff_vii = p(12);  
kon_viia = p(13);  
koff_viia = p(14);  
kon_viii = p(15);  
koff_viii = p(16);  
kon_viiia = p(17);  
koff_viiia = p(18);  
kon_ix = p(19);  
koff_ix = p(20);  
kon_ixa = p(21);  
koff_ixa = p(22);  
kon_x = p(23);  
koff_x = p(24);  
kon_xa = p(25);  
koff_xa = p(26);  
kon_tfpixa = p(27);  
koff_tfpixa = p(28);  
kon_tfpivh = p(29);  
koff_tfpivh = p(30);  
kon_ixat = p(31);  
koff_ixat = p(32);  
kon_xaat = p(33);  
koff_xaat = p(34);  
kon_ps = p(35);  
koff_ps = p(36);  
ka_s_tf_vii = p(37);  
kd_s_tfvii = p(38);  
ka_s_tf_viia = p(39);  
kd_s_tfviia = p(40);  
ka_s_ix_tfviia = p(41);  
kd_s_ixtfviia = p(42);  
kc_s_ixtfviia = p(43);  
ka_s_x_tfviia = p(44);  
kd_s_xtfviia = p(45);  
kc_s_xtfviia = p(46);  
ka_s_xa_tfviia = p(47);  
kd_s_xatfviia = p(48);  
ka_s_xa_tfvii = p(49);  
kd_s_xatfvii = p(50);  
kc_s_xatfvii = p(51);  
ka_on_s_tfpi_xa = p(52);  
kd_off_s_tfpixa = p(53);  
ka_tfpi_xa = p(54);  
kd_tfpixa = p(55);  
ka_s_tfpixa_tfviia = p(56);  
kd_s_tfviiaxatfpi = p(57);  
k_p1_tfviiaxatfpi = p(58);  
k_m1_xatfpitfviia = p(59);  
ka_on_s_tfpi_tfviiaxa = p(60);  
kd_off_s_tfviiaxatfpi = p(61);  
ka_s_xa_v = p(62);  
kd_s_xav = p(63);  
kc_s_xav = p(64);  
ka_s_xa_vii = p(65);  
kd_s_xavii = p(66);  
kc_s_xavii = p(67);  
ka_s_xa_viii = p(68);  
kd_s_xaviii = p(69);  
kc_s_xaviii = p(70);  
ka_s_viiia_ixa = p(71);  
kd_s_viiiaixa = p(72);  
ka_s_x_viiiaixa = p(73);  
kd_s_xviiiaixa = p(74);  
kc_s_xviiiaixa = p(75);  
ka_s_xa_vh = p(76);  
kd_s_xavh = p(77);  
ka_s_ii_xavh = p(78);  
kd_s_iixavh = p(79);  
kc_s_iixavh = p(80);  
ka_on_s_iia_v = p(81);  
kd_off_s_iiav = p(82);  
kc_s_iiav = p(83);  
ka_on_s_iia_vh = p(84);  
kd_off_s_iiavh = p(85);  
kc_s_iiavh = p(86);  
ka_on_s_iia_viii = p(87);  
kd_off_s_iiaviii = p(88);  
kc_s_iiaviii = p(89);  
ka_on_s_xia_ix = p(90);  
kd_off_s_xiaix = p(91);  
kc_s_xiaix = p(92);  
ka_s_xa_va = p(93);  
kd_s_xava = p(94);  
ka_s_ii_xava = p(95);  
kd_s_iixava = p(96);  
kc_s_iixava = p(97);  
ka_s_ixa_x = p(98);  
kd_s_ixax = p(99);  
kc_ixax = p(100);  
kc_xav = p(101);  
kc_xaviii = p(102);  
kc_s_ixax = p(103);  
ka_xa_v = p(104);  
kd_xav = p(105);  
ka_xa_vii = p(106);  
kd_xavii = p(107);  
kc_xavii = p(108);  
ka_xa_viii = p(109);  
kd_xaviii = p(110);  
ka_iia_v = p(111);  
kd_iiav = p(112);  
kc_iiav = p(113);  
ka_iia_vh = p(114);  
kd_iiavh = p(115);  
kc_iiavh = p(116);  
ka_iia_viii = p(117);  
kd_iiaviii = p(118);  
kc_iiaviii = p(119);  
ka_xia_ix = p(120);  
kd_xiaix = p(121);  
kc_xiaix = p(122);  
ka_iia_xi = p(123);  
kd_iiaxi = p(124);  
kc_iiaxi = p(125);  
ka_ixa_x = p(126);  
kd_ixax = p(127);  
ka_tfpi_vh = p(128);  
kd_tfpivh = p(129);  
ka_on_s_tfpi_vh = p(130);  
kd_off_s_tfpivh = p(131);  
ka_on_s_tfpixa_vh = p(132);  
kd_off_s_tfpixavh = p(133);  
ka_s_tfpixa_vh = p(134);  
kd_s_tfpixavh = p(135);  
ka_s_tfpivh_xa = p(136);  
kd_s_tfpivhxa = p(137);  
ka_iia_at = p(138);  
kd_iiaat = p(139);  
ka_ixa_at = p(140);  
kd_ixaat = p(141);  
ka_on_s_ixa_at = p(142);  
kd_off_s_ixaat = p(143);  
ka_xa_at = p(144);  
kd_xaat = p(145);  
ka_on_s_xa_at = p(146);  
kd_off_s_xaat = p(147);  
ka_xia_at = p(148);  
kd_xiaat = p(149);  
kon_sn_xii = p(150);  
koff_sn_xii = p(151);  
kon_sn_xiia = p(152);  
koff_sn_xiia = p(153);  
ka_sn_xiia_xii = p(154);  
kd_sn_xiiaxii = p(155);  
kc_sn_xiiaxii = p(156);  
kc_sn_xii = p(157);  
ka_xiia_pk = p(158);  
kd_xiiapk = p(159);  
kc_xiiapk = p(160);  
ka_pka_xii = p(161);  
kd_pkaxii = p(162);  
kc_pkaxii = p(163);  
ka_xiia_xi = p(164);  
kd_xiiaxi = p(165);  
kc_xiiaxi = p(166);  
ka_pka_ix = p(167);  
kd_pkaix = p(168);  
kc_pkaix = p(169);  
ka_sn_xiia_at = p(170);  
ka_xiia_at = p(171);  
ka_xia_c1inh = p(172);  
ka_xiia_c1inh = p(173);  
ka_sn_xiia_c1inh = p(174);  
ka_pka_c1inh = p(175);  
ka_s_pka_ix = p(176);  
kd_s_pkaix = p(177);  
kc_s_pkaix = p(178);  
ka_xi_hk = p(179);  
kd_xihk = p(180);  
kon_sn_xihk = p(181);  
koff_sn_xihk = p(182);  
ka_sn_xihk_xiia = p(183);  
kd_sn_xihkxiia = p(184);  
kc_sn_xihkxiia = p(185);  
koff_sn_xiahk = p(186);  
ka_pk_hk = p(187);  
kd_pkhk = p(188);  
kon_sn_pkhk = p(189);  
koff_sn_pkhk = p(190);  
ka_sn_pkhk_xiia = p(191);  
kd_sn_pkhkxiia = p(192);  
kc_sn_pkhkxiia = p(193);  
koff_sn_pkahk_split = p(194);  
ka_pka_hk = p(195);  
kd_pkahk = p(196);  
kon_sn_pkahk = p(197);  
koff_sn_pkahk = p(198);  
ka_sn_pkahk_xii = p(199);  
kd_sn_pkahkxii = p(200);  
kc_sn_pkahkxii = p(201);  


% Rename Binding Site 
nbs_ii = nbs(1);  
nbs_v = nbs(2);  
nbs_vii = nbs(3);  
nbs_viii = nbs(4);  
nbs_ix = nbs(5);  
nbs_ixa = nbs(6);  
nbs_x = nbs(7);  
nbs_tfpixa = nbs(8);  
nbs_tfpivh = nbs(9);  
nbs_ixaAT = nbs(10);  
nbs_xaAT = nbs(11);  
nbs_PS = nbs(12);  


% ODEs from reaction equations 

% L_TF
 dy(1)  =  -  kon_ii/nbs_ii * L_TF * II  +  koff_ii*nbs_ii * II_st  -  kon_iif/nbs_ii * L_TF * IIf  +  koff_iif*nbs_ii * IIf_st  -  kon_v/nbs_v * L_TF * V  +  koff_v*nbs_v * V_st  -  kon_vh/nbs_v * L_TF * Vh  +  koff_vh*nbs_v * Vh_st  -  kon_va/nbs_v * L_TF * Va  +  koff_va*nbs_v * Va_st  -  kon_vii/nbs_vii * L_TF * VII  +  koff_vii*nbs_vii * VII_st  -  kon_viia/nbs_vii * L_TF * VIIa  +  koff_viia*nbs_vii * VIIa_st  -  kon_viii/nbs_viii * L_TF * VIII  +  koff_viii*nbs_viii * VIII_st  -  kon_viiia/nbs_viii * L_TF * VIIIa  +  koff_viiia*nbs_viii * VIIIa_st  -  kon_ix/nbs_ix * L_TF * IX  +  koff_ix*nbs_ix * IX_st  -  kon_ixa/nbs_ixa * L_TF * IXa  +  koff_ixa*nbs_ixa * IXa_st  -  kon_x/nbs_x * L_TF * X  +  koff_x*nbs_x * X_st  -  kon_xa/nbs_x * L_TF * Xa  +  koff_xa*nbs_x * Xa_st  -  kon_tfpixa/nbs_tfpixa * L_TF * XabTFPI  +  koff_tfpixa*nbs_tfpixa * Xa_stbTFPI  -  kon_tfpivh/nbs_tfpivh * L_TF * TFPIbVh  +  koff_tfpivh*nbs_tfpivh * TFPIbVh_st  -  kon_ixat/nbs_ixaAT * L_TF * IXabAT  +  koff_ixat*nbs_ixaAT * IXa_stbAT  -  kon_xaat/nbs_xaAT * L_TF * XabAT  +  koff_xaat*nbs_xaAT * Xa_stbAT  -  kon_ps/nbs_PS * L_TF * PS  +  koff_ps*nbs_PS * PS_st;

% II
 dy(2)  =  -  kon_ii/nbs_ii * L_TF * II  +  koff_ii * II_st  -  kon_ii/nbs_ii * II * L_noTF  +  koff_ii * II_s;

% II_st
 dy(3)  =  +  kon_ii * L_TF * II  -  koff_ii/nbs_ii * II_st  -  ka_s_ii_xavh * II_st * Xa_stbVh_st  +  kd_s_iixavh * II_stbXa_stbVh_st  -  ka_s_ii_xava * II_st * Xa_stbVa_st  +  kd_s_iixava * II_stbXa_stbVa_st;

% L_noTF
 dy(4)  =  -  kon_ii/nbs_ii * II * L_noTF  +  koff_ii*nbs_ii * II_s  -  kon_iif/nbs_ii * L_noTF * IIf  +  koff_iif*nbs_ii * IIf_s  -  kon_v/nbs_v * L_noTF * V  +  koff_v*nbs_v * V_s  -  kon_vh/nbs_v * L_noTF * Vh  +  koff_vh*nbs_v * Vh_s  -  kon_va/nbs_v * L_noTF * Va  +  koff_va*nbs_v * Va_s  -  kon_vii/nbs_vii * L_noTF * VII  +  koff_vii*nbs_vii * VII_s  -  kon_viia/nbs_vii * L_noTF * VIIa  +  koff_viia*nbs_vii * VIIa_s  -  kon_viii/nbs_viii * L_noTF * VIII  +  koff_viii*nbs_viii * VIII_s  -  kon_viiia/nbs_viii * L_noTF * VIIIa  +  koff_viiia*nbs_viii * VIIIa_s  -  kon_ix/nbs_ix * L_noTF * IX  +  koff_ix*nbs_ix * IX_s  -  kon_ixa/nbs_ixa * L_noTF * IXa  +  koff_ixa*nbs_ixa * IXa_s  -  kon_x/nbs_x * L_noTF * X  +  koff_x*nbs_x * X_s  -  kon_xa/nbs_x * L_noTF * Xa  +  koff_xa*nbs_x * Xa_s  -  kon_tfpixa/nbs_tfpixa * L_noTF * XabTFPI  +  koff_tfpixa*nbs_tfpixa * Xa_sbTFPI  -  kon_tfpivh/nbs_tfpivh * L_noTF * TFPIbVh  +  koff_tfpivh*nbs_tfpivh * TFPIbVh_s  -  kon_ixat/nbs_ixaAT * L_noTF * IXabAT  +  koff_ixat*nbs_ixaAT * IXa_sbAT  -  kon_xaat/nbs_xaAT * L_noTF * XabAT  +  koff_xaat*nbs_xaAT * Xa_sbAT  -  kon_ps/nbs_PS * L_noTF * PS  +  koff_ps*nbs_PS * PS_s;

% II_s
 dy(5)  =  +  kon_ii * II * L_noTF  -  koff_ii/nbs_ii * II_s  -  ka_s_ii_xavh * II_s * Xa_sbVh_s  +  kd_s_iixavh * II_sbXa_sbVh_s  -  ka_s_ii_xava * II_s * Xa_sbVa_s  +  kd_s_iixava * II_sbXa_sbVa_s;

% IIf
 dy(6)  =  -  kon_iif/nbs_ii * L_TF * IIf  +  koff_iif * IIf_st  -  kon_iif/nbs_ii * L_noTF * IIf  +  koff_iif * IIf_s;

% IIf_st
 dy(7)  =  +  kon_iif * L_TF * IIf  -  koff_iif/nbs_ii * IIf_st  +  kc_s_iixavh * II_stbXa_stbVh_st  +  kc_s_iixava * II_stbXa_stbVa_st;

% IIf_s
 dy(8)  =  +  kon_iif * L_noTF * IIf  -  koff_iif/nbs_ii * IIf_s  +  kc_s_iixavh * II_sbXa_sbVh_s  +  kc_s_iixava * II_sbXa_sbVa_s;

% V
 dy(9)  =  -  kon_v/nbs_v * L_TF * V  +  koff_v * V_st  -  kon_v/nbs_v * L_noTF * V  +  koff_v * V_s  -  ka_xa_v * V * Xa  +  kd_xav * XabV  -  ka_iia_v * V * IIa  +  kd_iiav * IIabV;

% V_st
 dy(10)  =  +  kon_v * L_TF * V  -  koff_v/nbs_v * V_st  -  ka_s_xa_v * V_st * Xa_st  +  kd_s_xav * Xa_stbV_st  -  ka_on_s_iia_v * V_st * IIa  +  kd_off_s_iiav * IIabV_st;

% V_s
 dy(11)  =  +  kon_v * L_noTF * V  -  koff_v/nbs_v * V_s  -  ka_s_xa_v * V_s * Xa_s  +  kd_s_xav * Xa_sbV_s  -  ka_on_s_iia_v * V_s * IIa  +  kd_off_s_iiav * IIabV_s;

% Vh
 dy(12)  =  -  kon_vh/nbs_v * L_TF * Vh  +  koff_vh * Vh_st  -  kon_vh/nbs_v * L_noTF * Vh  +  koff_vh * Vh_s  +  kc_xav * XabV  -  ka_iia_vh * Vh * IIa  +  kd_iiavh * IIabVh  -  ka_tfpi_vh * Vh * TFPI  +  kd_tfpivh * TFPIbVh  -  ka_tfpi_vh * Vh * XabTFPI  +  kd_tfpivh * VhbTFPIbXa;

% Vh_st
 dy(13)  =  +  kon_vh * L_TF * Vh  -  koff_vh/nbs_v * Vh_st  +  kc_s_xav * Xa_stbV_st  -  ka_s_xa_vh * Vh_st * Xa_st  +  kd_s_xavh * Xa_stbVh_st  -  ka_on_s_iia_vh * Vh_st * IIa  +  kd_off_s_iiavh * IIabVh_st  -  ka_on_s_tfpi_vh * Vh_st * TFPI  +  kd_off_s_tfpivh * TFPIbVh_st  -  ka_on_s_tfpixa_vh * Vh_st * XabTFPI  +  kd_off_s_tfpixavh * Vh_stbTFPIbXa  -  ka_s_tfpixa_vh * Vh_st * Xa_stbTFPI  +  kd_s_tfpixavh * Vh_stbTFPIbXa_st;

% Vh_s
 dy(14)  =  +  kon_vh * L_noTF * Vh  -  koff_vh/nbs_v * Vh_s  +  kc_xav * Xa_sbV_s  -  ka_s_xa_vh * Vh_s * Xa_s  +  kd_s_xavh * Xa_sbVh_s  -  ka_on_s_iia_vh * Vh_s * IIa  +  kd_off_s_iiavh * IIabVh_s  -  ka_on_s_tfpi_vh * Vh_s * TFPI  +  kd_off_s_tfpivh * TFPIbVh_s  -  ka_on_s_tfpixa_vh * Vh_s * XabTFPI  +  kd_off_s_tfpixavh * Vh_sbTFPIbXa  -  ka_s_tfpixa_vh * Vh_s * Xa_sbTFPI  +  kd_s_tfpixavh * Vh_sbTFPIbXa_s;

% Va
 dy(15)  =  -  kon_va/nbs_v * L_TF * Va  +  koff_va * Va_st  -  kon_va/nbs_v * L_noTF * Va  +  koff_va * Va_s  +  kc_iiav * IIabV  +  kc_iiavh * IIabVh;

% Va_st
 dy(16)  =  +  kon_va * L_TF * Va  -  koff_va/nbs_v * Va_st  +  kc_s_iiav * IIabV_st  +  kc_s_iiavh * IIabVh_st  -  ka_s_xa_va * Va_st * Xa_st  +  kd_s_xava * Xa_stbVa_st;

% Va_s
 dy(17)  =  +  kon_va * L_noTF * Va  -  koff_va/nbs_v * Va_s  +  kc_s_iiav * IIabV_s  +  kc_s_iiavh * IIabVh_s  -  ka_s_xa_va * Va_s * Xa_s  +  kd_s_xava * Xa_sbVa_s;

% VII
 dy(18)  =  -  kon_vii/nbs_vii * L_TF * VII  +  koff_vii * VII_st  -  kon_vii/nbs_vii * L_noTF * VII  +  koff_vii * VII_s  -  ka_xa_vii * VII * Xa  +  kd_xavii * XabVII;

% VII_st
 dy(19)  =  +  kon_vii * L_TF * VII  -  koff_vii/nbs_vii * VII_st  -  ka_s_tf_vii * VII_st * TF  +  kd_s_tfvii * TFbVII_st  -  ka_s_xa_vii * VII_st * Xa_st  +  kd_s_xavii * Xa_stbVII_st;

% VII_s
 dy(20)  =  +  kon_vii * L_noTF * VII  -  koff_vii/nbs_vii * VII_s  -  ka_s_xa_vii * VII_s * Xa_s  +  kd_s_xavii * Xa_sbVII_s;

% VIIa
 dy(21)  =  -  kon_viia/nbs_vii * L_TF * VIIa  +  koff_viia * VIIa_st  -  kon_viia/nbs_vii * L_noTF * VIIa  +  koff_viia * VIIa_s  +  kc_xavii * XabVII;

% VIIa_st
 dy(22)  =  +  kon_viia * L_TF * VIIa  -  koff_viia/nbs_vii * VIIa_st  -  ka_s_tf_viia * VIIa_st * TF  +  kd_s_tfviia * TFbVIIa_st  +  kc_s_xavii * Xa_stbVII_st;

% VIIa_s
 dy(23)  =  +  kon_viia * L_noTF * VIIa  -  koff_viia/nbs_vii * VIIa_s  +  kc_xaviii * Xa_sbVII_s;

% VIII
 dy(24)  =  -  kon_viii/nbs_viii * L_TF * VIII  +  koff_viii * VIII_st  -  kon_viii/nbs_viii * L_noTF * VIII  +  koff_viii * VIII_s  -  ka_xa_viii * VIII * Xa  +  kd_xaviii * XabVIII  -  ka_iia_viii * VIII * IIa  +  kd_iiaviii * IIabVIII;

% VIII_st
 dy(25)  =  +  kon_viii * L_TF * VIII  -  koff_viii/nbs_viii * VIII_st  -  ka_s_xa_viii * VIII_st * Xa_st  +  kd_s_xaviii * Xa_stbVIII_st  -  ka_on_s_iia_viii * VIII_st * IIa  +  kd_off_s_iiaviii * IIabVIII_st;

% VIII_s
 dy(26)  =  +  kon_viii * L_noTF * VIII  -  koff_viii/nbs_viii * VIII_s  -  ka_s_xa_viii * VIII_s * Xa_s  +  kd_s_xaviii * Xa_sbVIII_s  -  ka_on_s_iia_viii * VIII_s * IIa  +  kd_off_s_iiaviii * IIabVIII_s;

% VIIIa
 dy(27)  =  -  kon_viiia/nbs_viii * L_TF * VIIIa  +  koff_viiia * VIIIa_st  -  kon_viiia/nbs_viii * L_noTF * VIIIa  +  koff_viiia * VIIIa_s  +  kc_xaviii * XabVIII  +  kc_iiaviii * IIabVIII;

% VIIIa_st
 dy(28)  =  +  kon_viiia * L_TF * VIIIa  -  koff_viiia/nbs_viii * VIIIa_st  +  kc_s_xaviii * Xa_stbVIII_st  -  ka_s_viiia_ixa * VIIIa_st * IXa_st  +  kd_s_viiiaixa * VIIIa_stbIXa_st  +  kc_s_iiaviii * IIabVIII_st;

% VIIIa_s
 dy(29)  =  +  kon_viiia * L_noTF * VIIIa  -  koff_viiia/nbs_viii * VIIIa_s  +  kc_xaviii * Xa_sbVIII_s  -  ka_s_viiia_ixa * VIIIa_s * IXa_s  +  kd_s_viiiaixa * VIIIa_sbIXa_s  +  kc_s_iiaviii * IIabVIII_s;

% IX
 dy(30)  =  -  kon_ix/nbs_ix * L_TF * IX  +  koff_ix * IX_st  -  kon_ix/nbs_ix * L_noTF * IX  +  koff_ix * IX_s  -  ka_xia_ix * IX * XIa  +  kd_xiaix * XIabIX  -  ka_pka_ix * IX * PKa  +  kd_pkaix * PKabIX;

% IX_st
 dy(31)  =  +  kon_ix * L_TF * IX  -  koff_ix/nbs_ix * IX_st  -  ka_s_ix_tfviia * IX_st * TFbVIIa_st  +  kd_s_ixtfviia * IX_stbTFbVIIa_st  -  ka_on_s_xia_ix * IX_st * XIa  +  kd_off_s_xiaix * XIabIX_st;

% IX_s
 dy(32)  =  +  kon_ix * L_noTF * IX  -  koff_ix/nbs_ix * IX_s  -  ka_on_s_xia_ix * IX_s * XIa  +  kd_off_s_xiaix * XIabIX_s  -  ka_s_pka_ix * IX_s * PKa  +  kd_s_pkaix * PKabIX_s;

% IXa
 dy(33)  =  -  kon_ixa/nbs_ixa * L_TF * IXa  +  koff_ixa * IXa_st  -  kon_ixa/nbs_ixa * L_noTF * IXa  +  koff_ixa * IXa_s  +  kc_xiaix * XIabIX  -  ka_ixa_x * IXa * X  +  kd_ixax * IXabX  +  kc_ixax * IXabX  -  ka_ixa_at * IXa * AT  +  kd_ixaat * IXabAT  +  kc_pkaix * PKabIX;

% IXa_st
 dy(34)  =  +  kon_ixa * L_TF * IXa  -  koff_ixa/nbs_ixa * IXa_st  +  kc_s_ixtfviia * IX_stbTFbVIIa_st  -  ka_s_viiia_ixa * VIIIa_st * IXa_st  +  kd_s_viiiaixa * VIIIa_stbIXa_st  +  kc_s_xiaix * XIabIX_st  -  ka_s_ixa_x * IXa_st * X_st  +  kd_s_ixax * IXa_stbX_st  +  kc_ixax * IXa_stbX_st  -  ka_on_s_ixa_at * IXa_st * AT  +  kd_off_s_ixaat * IXa_stbAT;

% IXa_s
 dy(35)  =  +  kon_ixa * L_noTF * IXa  -  koff_ixa/nbs_ixa * IXa_s  -  ka_s_viiia_ixa * VIIIa_s * IXa_s  +  kd_s_viiiaixa * VIIIa_sbIXa_s  +  kc_s_xiaix * XIabIX_s  -  ka_s_ixa_x * IXa_s * X_s  +  kd_s_ixax * IXa_sbX_s  +  kc_s_ixax * IXa_sbX_s  -  ka_on_s_ixa_at * IXa_s * AT  +  kd_off_s_ixaat * IXa_sbAT  +  kc_s_pkaix * PKabIX_s;

% X
 dy(36)  =  -  kon_x/nbs_x * L_TF * X  +  koff_x * X_st  -  kon_x/nbs_x * L_noTF * X  +  koff_x * X_s  -  ka_ixa_x * IXa * X  +  kd_ixax * IXabX;

% X_st
 dy(37)  =  +  kon_x * L_TF * X  -  koff_x/nbs_x * X_st  -  ka_s_x_tfviia * X_st * TFbVIIa_st  +  kd_s_xtfviia * TFbVIIa_stbX_st  -  ka_s_x_viiiaixa * X_st * VIIIa_stbIXa_st  +  kd_s_xviiiaixa * X_stbVIIIa_stbIXa_st  -  ka_s_ixa_x * IXa_st * X_st  +  kd_s_ixax * IXa_stbX_st;

% X_s
 dy(38)  =  +  kon_x * L_noTF * X  -  koff_x/nbs_x * X_s  -  ka_s_x_viiiaixa * X_s * VIIIa_sbIXa_s  +  kd_s_xviiiaixa * X_sbVIIIa_sbIXa_s  -  ka_s_ixa_x * IXa_s * X_s  +  kd_s_ixax * IXa_sbX_s;

% Xa
 dy(39)  =  -  kon_xa/nbs_x * L_TF * Xa  +  koff_xa * Xa_st  -  kon_xa/nbs_x * L_noTF * Xa  +  koff_xa * Xa_s  -  ka_tfpi_xa * Xa * TFPI  +  kd_tfpixa * XabTFPI  -  ka_xa_v * V * Xa  +  kd_xav * XabV  +  kc_xav * XabV  -  ka_xa_vii * VII * Xa  +  kd_xavii * XabVII  +  kc_xavii * XabVII  -  ka_xa_viii * VIII * Xa  +  kd_xaviii * XabVIII  +  kc_xaviii * XabVIII  +  kc_ixax * IXabX  -  ka_tfpi_xa * Xa * TFPIbVh  +  kd_tfpixa * VhbTFPIbXa  -  ka_xa_at * Xa * AT  +  kd_xaat * XabAT;

% Xa_st
 dy(40)  =  +  kon_xa * L_TF * Xa  -  koff_xa/nbs_x * Xa_st  -  ka_s_xa_tfviia * Xa_st * TFbVIIa_st  +  kd_s_xatfviia * TFbVIIa_stbXa_st  -  ka_s_xa_tfvii * Xa_st * TFbVII_st  +  kd_s_xatfvii * TFbVII_stbXa_st  -  ka_on_s_tfpi_xa * Xa_st * TFPI  +  kd_off_s_tfpixa * Xa_stbTFPI  -  ka_s_xa_v * V_st * Xa_st  +  kd_s_xav * Xa_stbV_st  +  kc_s_xav * Xa_stbV_st  -  ka_s_xa_vii * VII_st * Xa_st  +  kd_s_xavii * Xa_stbVII_st  +  kc_s_xavii * Xa_stbVII_st  -  ka_s_xa_viii * VIII_st * Xa_st  +  kd_s_xaviii * Xa_stbVIII_st  +  kc_s_xaviii * Xa_stbVIII_st  +  kc_s_xviiiaixa * X_stbVIIIa_stbIXa_st  -  ka_s_xa_vh * Vh_st * Xa_st  +  kd_s_xavh * Xa_stbVh_st  -  ka_s_xa_va * Va_st * Xa_st  +  kd_s_xava * Xa_stbVa_st  +  kc_ixax * IXa_stbX_st  -  ka_on_s_tfpi_xa * Xa_st * TFPIbVh  +  kd_off_s_tfpixa * VhbTFPIbXa_st  -  ka_s_tfpivh_xa * Xa_st * TFPIbVh_st  +  kd_s_tfpivhxa * Vh_stbTFPIbXa_st  -  ka_on_s_xa_at * Xa_st * AT  +  kd_off_s_xaat * Xa_stbAT;

% Xa_s
 dy(41)  =  +  kon_xa * L_noTF * Xa  -  koff_xa/nbs_x * Xa_s  -  ka_on_s_tfpi_xa * Xa_s * TFPI  +  kd_off_s_tfpixa * Xa_sbTFPI  -  ka_s_xa_v * V_s * Xa_s  +  kd_s_xav * Xa_sbV_s  +  kc_xav * Xa_sbV_s  -  ka_s_xa_vii * VII_s * Xa_s  +  kd_s_xavii * Xa_sbVII_s  +  kc_xaviii * Xa_sbVII_s  -  ka_s_xa_viii * VIII_s * Xa_s  +  kd_s_xaviii * Xa_sbVIII_s  +  kc_xaviii * Xa_sbVIII_s  +  kc_s_xviiiaixa * X_sbVIIIa_sbIXa_s  -  ka_s_xa_vh * Vh_s * Xa_s  +  kd_s_xavh * Xa_sbVh_s  -  ka_s_xa_va * Va_s * Xa_s  +  kd_s_xava * Xa_sbVa_s  +  kc_s_ixax * IXa_sbX_s  -  ka_on_s_tfpi_xa * Xa_s * TFPIbVh  +  kd_off_s_tfpixa * VhbTFPIbXa_s  -  ka_s_tfpivh_xa * Xa_s * TFPIbVh_s  +  kd_s_tfpivhxa * Vh_sbTFPIbXa_s  -  ka_on_s_xa_at * Xa_s * AT  +  kd_off_s_xaat * Xa_sbAT;

% XabTFPI
 dy(42)  =  -  kon_tfpixa/nbs_tfpixa * L_TF * XabTFPI  +  koff_tfpixa * Xa_stbTFPI  -  kon_tfpixa/nbs_tfpixa * L_noTF * XabTFPI  +  koff_tfpixa * Xa_sbTFPI  +  ka_tfpi_xa * Xa * TFPI  -  kd_tfpixa * XabTFPI  -  ka_tfpi_vh * Vh * XabTFPI  +  kd_tfpivh * VhbTFPIbXa  -  ka_on_s_tfpixa_vh * Vh_s * XabTFPI  +  kd_off_s_tfpixavh * Vh_sbTFPIbXa  -  ka_on_s_tfpixa_vh * Vh_st * XabTFPI  +  kd_off_s_tfpixavh * Vh_stbTFPIbXa;

% Xa_stbTFPI
 dy(43)  =  +  kon_tfpixa * L_TF * XabTFPI  -  koff_tfpixa/nbs_tfpixa * Xa_stbTFPI  +  ka_on_s_tfpi_xa * Xa_st * TFPI  -  kd_off_s_tfpixa * Xa_stbTFPI  -  ka_s_tfpixa_tfviia * Xa_stbTFPI * TFbVIIa_st  +  kd_s_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI  -  ka_s_tfpixa_vh * Vh_st * Xa_stbTFPI  +  kd_s_tfpixavh * Vh_stbTFPIbXa_st;

% Xa_sbTFPI
 dy(44)  =  +  kon_tfpixa * L_noTF * XabTFPI  -  koff_tfpixa/nbs_tfpixa * Xa_sbTFPI  +  ka_on_s_tfpi_xa * Xa_s * TFPI  -  kd_off_s_tfpixa * Xa_sbTFPI  -  ka_s_tfpixa_vh * Vh_s * Xa_sbTFPI  +  kd_s_tfpixavh * Vh_sbTFPIbXa_s;

% TFPIbVh
 dy(45)  =  -  kon_tfpivh/nbs_tfpivh * L_TF * TFPIbVh  +  koff_tfpivh * TFPIbVh_st  -  kon_tfpivh/nbs_tfpivh * L_noTF * TFPIbVh  +  koff_tfpivh * TFPIbVh_s  +  ka_tfpi_vh * Vh * TFPI  -  kd_tfpivh * TFPIbVh  -  ka_tfpi_xa * Xa * TFPIbVh  +  kd_tfpixa * VhbTFPIbXa  -  ka_on_s_tfpi_xa * Xa_s * TFPIbVh  +  kd_off_s_tfpixa * VhbTFPIbXa_s  -  ka_on_s_tfpi_xa * Xa_st * TFPIbVh  +  kd_off_s_tfpixa * VhbTFPIbXa_st;

% TFPIbVh_st
 dy(46)  =  +  kon_tfpivh * L_TF * TFPIbVh  -  koff_tfpivh/nbs_tfpivh * TFPIbVh_st  +  ka_on_s_tfpi_vh * Vh_st * TFPI  -  kd_off_s_tfpivh * TFPIbVh_st  -  ka_s_tfpivh_xa * Xa_st * TFPIbVh_st  +  kd_s_tfpivhxa * Vh_stbTFPIbXa_st;

% TFPIbVh_s
 dy(47)  =  +  kon_tfpivh * L_noTF * TFPIbVh  -  koff_tfpivh/nbs_tfpivh * TFPIbVh_s  +  ka_on_s_tfpi_vh * Vh_s * TFPI  -  kd_off_s_tfpivh * TFPIbVh_s  -  ka_s_tfpivh_xa * Xa_s * TFPIbVh_s  +  kd_s_tfpivhxa * Vh_sbTFPIbXa_s;

% IXabAT
 dy(48)  =  -  kon_ixat/nbs_ixaAT * L_noTF * IXabAT  +  koff_ixat * IXa_sbAT  -  kon_ixat/nbs_ixaAT * L_TF * IXabAT  +  koff_ixat * IXa_stbAT  +  ka_ixa_at * IXa * AT  -  kd_ixaat * IXabAT;

% IXa_sbAT
 dy(49)  =  +  kon_ixat * L_noTF * IXabAT  -  koff_ixat/nbs_ixaAT * IXa_sbAT  +  ka_on_s_ixa_at * IXa_s * AT  -  kd_off_s_ixaat * IXa_sbAT;

% IXa_stbAT
 dy(50)  =  +  kon_ixat * L_TF * IXabAT  -  koff_ixat/nbs_ixaAT * IXa_stbAT  +  ka_on_s_ixa_at * IXa_st * AT  -  kd_off_s_ixaat * IXa_stbAT;

% XabAT
 dy(51)  =  -  kon_xaat/nbs_xaAT * L_noTF * XabAT  +  koff_xaat * Xa_sbAT  -  kon_xaat/nbs_xaAT * L_TF * XabAT  +  koff_xaat * Xa_stbAT  +  ka_xa_at * Xa * AT  -  kd_xaat * XabAT;

% Xa_sbAT
 dy(52)  =  +  kon_xaat * L_noTF * XabAT  -  koff_xaat/nbs_xaAT * Xa_sbAT  +  ka_on_s_xa_at * Xa_s * AT  -  kd_off_s_xaat * Xa_sbAT;

% Xa_stbAT
 dy(53)  =  +  kon_xaat * L_TF * XabAT  -  koff_xaat/nbs_xaAT * Xa_stbAT  +  ka_on_s_xa_at * Xa_st * AT  -  kd_off_s_xaat * Xa_stbAT;

% PS
 dy(54)  =  -  kon_ps/nbs_PS * L_TF * PS  +  koff_ps * PS_st  -  kon_ps/nbs_PS * L_noTF * PS  +  koff_ps * PS_s;

% PS_st
 dy(55)  =  +  kon_ps * L_TF * PS  -  koff_ps/nbs_PS * PS_st;

% PS_s
 dy(56)  =  +  kon_ps * L_noTF * PS  -  koff_ps/nbs_PS * PS_s;

% TF
 dy(57)  =  -  ka_s_tf_vii * VII_st * TF  +  kd_s_tfvii * TFbVII_st  -  ka_s_tf_viia * VIIa_st * TF  +  kd_s_tfviia * TFbVIIa_st;

% TFbVII_st
 dy(58)  =  +  ka_s_tf_vii * VII_st * TF  -  kd_s_tfvii * TFbVII_st  -  ka_s_xa_tfvii * Xa_st * TFbVII_st  +  kd_s_xatfvii * TFbVII_stbXa_st;

% TFbVIIa_st
 dy(59)  =  +  ka_s_tf_viia * VIIa_st * TF  -  kd_s_tfviia * TFbVIIa_st  -  ka_s_ix_tfviia * IX_st * TFbVIIa_st  +  kd_s_ixtfviia * IX_stbTFbVIIa_st  +  kc_s_ixtfviia * IX_stbTFbVIIa_st  -  ka_s_x_tfviia * X_st * TFbVIIa_st  +  kd_s_xtfviia * TFbVIIa_stbX_st  -  ka_s_xa_tfviia * Xa_st * TFbVIIa_st  +  kd_s_xatfviia * TFbVIIa_stbXa_st  -  ka_s_tfpixa_tfviia * Xa_stbTFPI * TFbVIIa_st  +  kd_s_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI;

% IX_stbTFbVIIa_st
 dy(60)  =  +  ka_s_ix_tfviia * IX_st * TFbVIIa_st  -  kd_s_ixtfviia * IX_stbTFbVIIa_st  -  kc_s_ixtfviia * IX_stbTFbVIIa_st;

% TFbVIIa_stbX_st
 dy(61)  =  +  ka_s_x_tfviia * X_st * TFbVIIa_st  -  kd_s_xtfviia * TFbVIIa_stbX_st  -  kc_s_xtfviia * TFbVIIa_stbX_st;

% TFbVIIa_stbXa_st
 dy(62)  =  +  kc_s_xtfviia * TFbVIIa_stbX_st  +  ka_s_xa_tfviia * Xa_st * TFbVIIa_st  -  kd_s_xatfviia * TFbVIIa_stbXa_st  +  kc_s_xatfvii * TFbVII_stbXa_st  -  ka_on_s_tfpi_tfviiaxa * TFbVIIa_stbXa_st * TFPI  +  kd_off_s_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI;

% TFbVII_stbXa_st
 dy(63)  =  +  ka_s_xa_tfvii * Xa_st * TFbVII_st  -  kd_s_xatfvii * TFbVII_stbXa_st  -  kc_s_xatfvii * TFbVII_stbXa_st;

% TFPI
 dy(64)  =  -  ka_on_s_tfpi_xa * Xa_st * TFPI  +  kd_off_s_tfpixa * Xa_stbTFPI  -  ka_on_s_tfpi_xa * Xa_s * TFPI  +  kd_off_s_tfpixa * Xa_sbTFPI  -  ka_tfpi_xa * Xa * TFPI  +  kd_tfpixa * XabTFPI  -  ka_on_s_tfpi_tfviiaxa * TFbVIIa_stbXa_st * TFPI  +  kd_off_s_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI  -  ka_tfpi_vh * Vh * TFPI  +  kd_tfpivh * TFPIbVh  -  ka_on_s_tfpi_vh * Vh_s * TFPI  +  kd_off_s_tfpivh * TFPIbVh_s  -  ka_on_s_tfpi_vh * Vh_st * TFPI  +  kd_off_s_tfpivh * TFPIbVh_st  -  ka_on_s_tfpi_xa * TFPI * Xa_sbVh_s  +  kd_off_s_tfpixa * TFPIbXa_sbVh_s  -  ka_on_s_tfpi_xa * TFPI * Xa_stbVh_st  +  kd_off_s_tfpixa * TFPIbXa_stbVh_st  -  ka_on_s_tfpi_vh * TFPI * Xa_sbVh_s  +  kd_off_s_tfpivh * Xa_sbVh_sbTFPI  -  ka_on_s_tfpi_vh * TFPI * Xa_stbVh_st  +  kd_off_s_tfpivh * Xa_stbVh_stbTFPI;

% TFbVIIa_stbXa_stbTFPI
 dy(65)  =  +  ka_s_tfpixa_tfviia * Xa_stbTFPI * TFbVIIa_st  -  kd_s_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI  -  k_p1_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI  +  k_m1_xatfpitfviia * Xa_stbTFPIbTFbVIIa_st  +  ka_on_s_tfpi_tfviiaxa * TFbVIIa_stbXa_st * TFPI  -  kd_off_s_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI;

% Xa_stbTFPIbTFbVIIa_st
 dy(66)  =  +  k_p1_tfviiaxatfpi * TFbVIIa_stbXa_stbTFPI  -  k_m1_xatfpitfviia * Xa_stbTFPIbTFbVIIa_st;

% Xa_stbV_st
 dy(67)  =  +  ka_s_xa_v * V_st * Xa_st  -  kd_s_xav * Xa_stbV_st  -  kc_s_xav * Xa_stbV_st;

% Xa_stbVII_st
 dy(68)  =  +  ka_s_xa_vii * VII_st * Xa_st  -  kd_s_xavii * Xa_stbVII_st  -  kc_s_xavii * Xa_stbVII_st;

% Xa_stbVIII_st
 dy(69)  =  +  ka_s_xa_viii * VIII_st * Xa_st  -  kd_s_xaviii * Xa_stbVIII_st  -  kc_s_xaviii * Xa_stbVIII_st;

% VIIIa_stbIXa_st
 dy(70)  =  +  ka_s_viiia_ixa * VIIIa_st * IXa_st  -  kd_s_viiiaixa * VIIIa_stbIXa_st  -  ka_s_x_viiiaixa * X_st * VIIIa_stbIXa_st  +  kd_s_xviiiaixa * X_stbVIIIa_stbIXa_st  +  kc_s_xviiiaixa * X_stbVIIIa_stbIXa_st;

% X_stbVIIIa_stbIXa_st
 dy(71)  =  +  ka_s_x_viiiaixa * X_st * VIIIa_stbIXa_st  -  kd_s_xviiiaixa * X_stbVIIIa_stbIXa_st  -  kc_s_xviiiaixa * X_stbVIIIa_stbIXa_st;

% Xa_stbVh_st
 dy(72)  =  +  ka_s_xa_vh * Vh_st * Xa_st  -  kd_s_xavh * Xa_stbVh_st  -  ka_s_ii_xavh * II_st * Xa_stbVh_st  +  kd_s_iixavh * II_stbXa_stbVh_st  +  kc_s_iixavh * II_stbXa_stbVh_st  -  ka_on_s_tfpi_xa * TFPI * Xa_stbVh_st  +  kd_off_s_tfpixa * TFPIbXa_stbVh_st  -  ka_on_s_tfpi_vh * TFPI * Xa_stbVh_st  +  kd_off_s_tfpivh * Xa_stbVh_stbTFPI;

% II_stbXa_stbVh_st
 dy(73)  =  +  ka_s_ii_xavh * II_st * Xa_stbVh_st  -  kd_s_iixavh * II_stbXa_stbVh_st  -  kc_s_iixavh * II_stbXa_stbVh_st;

% IIa
 dy(74)  =  +  kc_s_iixavh * II_stbXa_stbVh_st  -  ka_on_s_iia_v * V_st * IIa  +  kd_off_s_iiav * IIabV_st  +  kc_s_iiav * IIabV_st  -  ka_on_s_iia_vh * Vh_st * IIa  +  kd_off_s_iiavh * IIabVh_st  +  kc_s_iiavh * IIabVh_st  -  ka_on_s_iia_viii * VIII_st * IIa  +  kd_off_s_iiaviii * IIabVIII_st  +  kc_s_iiaviii * IIabVIII_st  +  kc_s_iixava * II_stbXa_stbVa_st  +  kc_s_iixavh * II_sbXa_sbVh_s  -  ka_on_s_iia_v * V_s * IIa  +  kd_off_s_iiav * IIabV_s  +  kc_s_iiav * IIabV_s  -  ka_on_s_iia_vh * Vh_s * IIa  +  kd_off_s_iiavh * IIabVh_s  +  kc_s_iiavh * IIabVh_s  -  ka_on_s_iia_viii * VIII_s * IIa  +  kd_off_s_iiaviii * IIabVIII_s  +  kc_s_iiaviii * IIabVIII_s  +  kc_s_iixava * II_sbXa_sbVa_s  -  ka_iia_v * V * IIa  +  kd_iiav * IIabV  +  kc_iiav * IIabV  -  ka_iia_vh * Vh * IIa  +  kd_iiavh * IIabVh  +  kc_iiavh * IIabVh  -  ka_iia_viii * VIII * IIa  +  kd_iiaviii * IIabVIII  +  kc_iiaviii * IIabVIII  -  ka_iia_xi * IIa * XI  +  kd_iiaxi * IIabXI  +  kc_iiaxi * IIabXI  -  ka_iia_at * IIa * AT  +  kd_iiaat * IIabAT;

% IIabV_st
 dy(75)  =  +  ka_on_s_iia_v * V_st * IIa  -  kd_off_s_iiav * IIabV_st  -  kc_s_iiav * IIabV_st;

% IIabVh_st
 dy(76)  =  +  ka_on_s_iia_vh * Vh_st * IIa  -  kd_off_s_iiavh * IIabVh_st  -  kc_s_iiavh * IIabVh_st;

% IIabVIII_st
 dy(77)  =  +  ka_on_s_iia_viii * VIII_st * IIa  -  kd_off_s_iiaviii * IIabVIII_st  -  kc_s_iiaviii * IIabVIII_st;

% XIa
 dy(78)  =  -  ka_on_s_xia_ix * IX_st * XIa  +  kd_off_s_xiaix * XIabIX_st  +  kc_s_xiaix * XIabIX_st  -  ka_on_s_xia_ix * IX_s * XIa  +  kd_off_s_xiaix * XIabIX_s  +  kc_s_xiaix * XIabIX_s  -  ka_xia_ix * IX * XIa  +  kd_xiaix * XIabIX  +  kc_xiaix * XIabIX  +  kc_iiaxi * IIabXI  -  ka_xia_at * XIa * AT  +  kd_xiaat * XIabAT  +  kc_xiiaxi * XIIabXI  -  ka_xia_c1inh * XIa * C1INH  +  koff_sn_xiahk * XIabHK_sn;

% XIabIX_st
 dy(79)  =  +  ka_on_s_xia_ix * IX_st * XIa  -  kd_off_s_xiaix * XIabIX_st  -  kc_s_xiaix * XIabIX_st;

% Xa_stbVa_st
 dy(80)  =  +  ka_s_xa_va * Va_st * Xa_st  -  kd_s_xava * Xa_stbVa_st  -  ka_s_ii_xava * II_st * Xa_stbVa_st  +  kd_s_iixava * II_stbXa_stbVa_st  +  kc_s_iixava * II_stbXa_stbVa_st;

% II_stbXa_stbVa_st
 dy(81)  =  +  ka_s_ii_xava * II_st * Xa_stbVa_st  -  kd_s_iixava * II_stbXa_stbVa_st  -  kc_s_iixava * II_stbXa_stbVa_st;

% IXa_stbX_st
 dy(82)  =  +  ka_s_ixa_x * IXa_st * X_st  -  kd_s_ixax * IXa_stbX_st  -  kc_ixax * IXa_stbX_st;

% Xa_sbV_s
 dy(83)  =  +  ka_s_xa_v * V_s * Xa_s  -  kd_s_xav * Xa_sbV_s  -  kc_xav * Xa_sbV_s;

% Xa_sbVII_s
 dy(84)  =  +  ka_s_xa_vii * VII_s * Xa_s  -  kd_s_xavii * Xa_sbVII_s  -  kc_xaviii * Xa_sbVII_s;

% Xa_sbVIII_s
 dy(85)  =  +  ka_s_xa_viii * VIII_s * Xa_s  -  kd_s_xaviii * Xa_sbVIII_s  -  kc_xaviii * Xa_sbVIII_s;

% VIIIa_sbIXa_s
 dy(86)  =  +  ka_s_viiia_ixa * VIIIa_s * IXa_s  -  kd_s_viiiaixa * VIIIa_sbIXa_s  -  ka_s_x_viiiaixa * X_s * VIIIa_sbIXa_s  +  kd_s_xviiiaixa * X_sbVIIIa_sbIXa_s  +  kc_s_xviiiaixa * X_sbVIIIa_sbIXa_s;

% X_sbVIIIa_sbIXa_s
 dy(87)  =  +  ka_s_x_viiiaixa * X_s * VIIIa_sbIXa_s  -  kd_s_xviiiaixa * X_sbVIIIa_sbIXa_s  -  kc_s_xviiiaixa * X_sbVIIIa_sbIXa_s;

% Xa_sbVh_s
 dy(88)  =  +  ka_s_xa_vh * Vh_s * Xa_s  -  kd_s_xavh * Xa_sbVh_s  -  ka_s_ii_xavh * II_s * Xa_sbVh_s  +  kd_s_iixavh * II_sbXa_sbVh_s  +  kc_s_iixavh * II_sbXa_sbVh_s  -  ka_on_s_tfpi_xa * TFPI * Xa_sbVh_s  +  kd_off_s_tfpixa * TFPIbXa_sbVh_s  -  ka_on_s_tfpi_vh * TFPI * Xa_sbVh_s  +  kd_off_s_tfpivh * Xa_sbVh_sbTFPI;

% II_sbXa_sbVh_s
 dy(89)  =  +  ka_s_ii_xavh * II_s * Xa_sbVh_s  -  kd_s_iixavh * II_sbXa_sbVh_s  -  kc_s_iixavh * II_sbXa_sbVh_s;

% IIabV_s
 dy(90)  =  +  ka_on_s_iia_v * V_s * IIa  -  kd_off_s_iiav * IIabV_s  -  kc_s_iiav * IIabV_s;

% IIabVh_s
 dy(91)  =  +  ka_on_s_iia_vh * Vh_s * IIa  -  kd_off_s_iiavh * IIabVh_s  -  kc_s_iiavh * IIabVh_s;

% IIabVIII_s
 dy(92)  =  +  ka_on_s_iia_viii * VIII_s * IIa  -  kd_off_s_iiaviii * IIabVIII_s  -  kc_s_iiaviii * IIabVIII_s;

% XIabIX_s
 dy(93)  =  +  ka_on_s_xia_ix * IX_s * XIa  -  kd_off_s_xiaix * XIabIX_s  -  kc_s_xiaix * XIabIX_s;

% Xa_sbVa_s
 dy(94)  =  +  ka_s_xa_va * Va_s * Xa_s  -  kd_s_xava * Xa_sbVa_s  -  ka_s_ii_xava * II_s * Xa_sbVa_s  +  kd_s_iixava * II_sbXa_sbVa_s  +  kc_s_iixava * II_sbXa_sbVa_s;

% II_sbXa_sbVa_s
 dy(95)  =  +  ka_s_ii_xava * II_s * Xa_sbVa_s  -  kd_s_iixava * II_sbXa_sbVa_s  -  kc_s_iixava * II_sbXa_sbVa_s;

% IXa_sbX_s
 dy(96)  =  +  ka_s_ixa_x * IXa_s * X_s  -  kd_s_ixax * IXa_sbX_s  -  kc_s_ixax * IXa_sbX_s;

% XabV
 dy(97)  =  +  ka_xa_v * V * Xa  -  kd_xav * XabV  -  kc_xav * XabV;

% XabVII
 dy(98)  =  +  ka_xa_vii * VII * Xa  -  kd_xavii * XabVII  -  kc_xavii * XabVII;

% XabVIII
 dy(99)  =  +  ka_xa_viii * VIII * Xa  -  kd_xaviii * XabVIII  -  kc_xaviii * XabVIII;

% IIabV
 dy(100)  =  +  ka_iia_v * V * IIa  -  kd_iiav * IIabV  -  kc_iiav * IIabV;

% IIabVh
 dy(101)  =  +  ka_iia_vh * Vh * IIa  -  kd_iiavh * IIabVh  -  kc_iiavh * IIabVh;

% IIabVIII
 dy(102)  =  +  ka_iia_viii * VIII * IIa  -  kd_iiaviii * IIabVIII  -  kc_iiaviii * IIabVIII;

% XIabIX
 dy(103)  =  +  ka_xia_ix * IX * XIa  -  kd_xiaix * XIabIX  -  kc_xiaix * XIabIX;

% XI
 dy(104)  =  -  ka_iia_xi * IIa * XI  +  kd_iiaxi * IIabXI  -  ka_xiia_xi * XI * XIIa  +  kd_xiiaxi * XIIabXI  -  ka_xi_hk * XI * HK  +  kd_xihk * XIbHK;

% IIabXI
 dy(105)  =  +  ka_iia_xi * IIa * XI  -  kd_iiaxi * IIabXI  -  kc_iiaxi * IIabXI;

% IXabX
 dy(106)  =  +  ka_ixa_x * IXa * X  -  kd_ixax * IXabX  -  kc_ixax * IXabX;

% TFPIbXa_sbVh_s
 dy(107)  =  +  ka_on_s_tfpi_xa * TFPI * Xa_sbVh_s  -  kd_off_s_tfpixa * TFPIbXa_sbVh_s;

% TFPIbXa_stbVh_st
 dy(108)  =  +  ka_on_s_tfpi_xa * TFPI * Xa_stbVh_st  -  kd_off_s_tfpixa * TFPIbXa_stbVh_st;

% Xa_sbVh_sbTFPI
 dy(109)  =  +  ka_on_s_tfpi_vh * TFPI * Xa_sbVh_s  -  kd_off_s_tfpivh * Xa_sbVh_sbTFPI;

% Xa_stbVh_stbTFPI
 dy(110)  =  +  ka_on_s_tfpi_vh * TFPI * Xa_stbVh_st  -  kd_off_s_tfpivh * Xa_stbVh_stbTFPI;

% VhbTFPIbXa
 dy(111)  =  +  ka_tfpi_vh * Vh * XabTFPI  -  kd_tfpivh * VhbTFPIbXa  +  ka_tfpi_xa * Xa * TFPIbVh  -  kd_tfpixa * VhbTFPIbXa;

% Vh_sbTFPIbXa
 dy(112)  =  +  ka_on_s_tfpixa_vh * Vh_s * XabTFPI  -  kd_off_s_tfpixavh * Vh_sbTFPIbXa;

% Vh_stbTFPIbXa
 dy(113)  =  +  ka_on_s_tfpixa_vh * Vh_st * XabTFPI  -  kd_off_s_tfpixavh * Vh_stbTFPIbXa;

% Vh_sbTFPIbXa_s
 dy(114)  =  +  ka_s_tfpixa_vh * Vh_s * Xa_sbTFPI  -  kd_s_tfpixavh * Vh_sbTFPIbXa_s  +  ka_s_tfpivh_xa * Xa_s * TFPIbVh_s  -  kd_s_tfpivhxa * Vh_sbTFPIbXa_s;

% Vh_stbTFPIbXa_st
 dy(115)  =  +  ka_s_tfpixa_vh * Vh_st * Xa_stbTFPI  -  kd_s_tfpixavh * Vh_stbTFPIbXa_st  +  ka_s_tfpivh_xa * Xa_st * TFPIbVh_st  -  kd_s_tfpivhxa * Vh_stbTFPIbXa_st;

% VhbTFPIbXa_s
 dy(116)  =  +  ka_on_s_tfpi_xa * Xa_s * TFPIbVh  -  kd_off_s_tfpixa * VhbTFPIbXa_s;

% VhbTFPIbXa_st
 dy(117)  =  +  ka_on_s_tfpi_xa * Xa_st * TFPIbVh  -  kd_off_s_tfpixa * VhbTFPIbXa_st;

% AT
 dy(118)  =  -  ka_iia_at * IIa * AT  +  kd_iiaat * IIabAT  -  ka_ixa_at * IXa * AT  +  kd_ixaat * IXabAT  -  ka_on_s_ixa_at * IXa_s * AT  +  kd_off_s_ixaat * IXa_sbAT  -  ka_on_s_ixa_at * IXa_st * AT  +  kd_off_s_ixaat * IXa_stbAT  -  ka_xa_at * Xa * AT  +  kd_xaat * XabAT  -  ka_on_s_xa_at * Xa_s * AT  +  kd_off_s_xaat * Xa_sbAT  -  ka_on_s_xa_at * Xa_st * AT  +  kd_off_s_xaat * Xa_stbAT  -  ka_xia_at * XIa * AT  +  kd_xiaat * XIabAT  -  ka_sn_xiia_at * AT * XIIa_sn  -  ka_xiia_at * AT * XIIa;

% IIabAT
 dy(119)  =  +  ka_iia_at * IIa * AT  -  kd_iiaat * IIabAT;

% XIabAT
 dy(120)  =  +  ka_xia_at * XIa * AT  -  kd_xiaat * XIabAT;

% XII
 dy(121)  =  -  kon_sn_xii * XII * SN  +  koff_sn_xii * XII_sn  -  ka_pka_xii * XII * PKa  +  kd_pkaxii * PKabXII;

% SN
 dy(122)  =  -  kon_sn_xii * XII * SN  +  koff_sn_xii * XII_sn  -  kon_sn_xiia * SN * XIIa  +  koff_sn_xiia * XIIa_sn  -  kon_sn_xihk * SN * XIbHK  +  koff_sn_xihk * XIbHK_sn  +  koff_sn_xiahk * XIabHK_sn  -  kon_sn_pkhk * SN * PKbHK  +  koff_sn_pkhk * PKbHK_sn  +  koff_sn_pkahk_split * PKabHK_sn  -  kon_sn_pkahk * SN * PKabHK  +  koff_sn_pkahk * PKabHK_sn  +  kc_sn_pkahkxii * PKabHK_snbXII_sn;

% XII_sn
 dy(123)  =  +  kon_sn_xii * XII * SN  -  koff_sn_xii * XII_sn  -  ka_sn_xiia_xii * XII_sn * XIIa_sn  +  kd_sn_xiiaxii * XIIa_snbXII_sn  -  kc_sn_xii * XII_sn  -  ka_sn_pkahk_xii * XII_sn * PKabHK_sn  +  kd_sn_pkahkxii * PKabHK_snbXII_sn;

% XIIa
 dy(124)  =  -  kon_sn_xiia * SN * XIIa  +  koff_sn_xiia * XIIa_sn  -  ka_xiia_pk * XIIa * PK  +  kd_xiiapk * XIIabPK  +  kc_xiiapk * XIIabPK  +  kc_pkaxii * PKabXII  -  ka_xiia_xi * XI * XIIa  +  kd_xiiaxi * XIIabXI  +  kc_xiiaxi * XIIabXI  -  ka_xiia_at * AT * XIIa  -  ka_xiia_c1inh * XIIa * C1INH;

% XIIa_sn
 dy(125)  =  +  kon_sn_xiia * SN * XIIa  -  koff_sn_xiia * XIIa_sn  -  ka_sn_xiia_xii * XII_sn * XIIa_sn  +  kd_sn_xiiaxii * XIIa_snbXII_sn  +  kc_sn_xiiaxii * 2 * XIIa_snbXII_sn  +  kc_sn_xii * XII_sn  -  ka_sn_xiia_at * AT * XIIa_sn  -  ka_sn_xiia_c1inh * XIIa_sn * C1INH  -  ka_sn_xihk_xiia * XIIa_sn * XIbHK_sn  +  kd_sn_xihkxiia * XIbHK_snbXIIa_sn  +  kc_sn_xihkxiia * XIbHK_snbXIIa_sn  -  ka_sn_pkhk_xiia * XIIa_sn * PKbHK_sn  +  kd_sn_pkhkxiia * PKbHK_snbXIIa_sn  +  kc_sn_pkhkxiia * PKbHK_snbXIIa_sn  +  kc_sn_pkahkxii * PKabHK_snbXII_sn;

% XIIa_snbXII_sn
 dy(126)  =  +  ka_sn_xiia_xii * XII_sn * XIIa_sn  -  kd_sn_xiiaxii * XIIa_snbXII_sn  -  kc_sn_xiiaxii * XIIa_snbXII_sn;

% PK
 dy(127)  =  -  ka_xiia_pk * XIIa * PK  +  kd_xiiapk * XIIabPK  -  ka_pk_hk * PK * HK  +  kd_pkhk * PKbHK;

% XIIabPK
 dy(128)  =  +  ka_xiia_pk * XIIa * PK  -  kd_xiiapk * XIIabPK  -  kc_xiiapk * XIIabPK;

% PKa
 dy(129)  =  +  kc_xiiapk * XIIabPK  -  ka_pka_xii * XII * PKa  +  kd_pkaxii * PKabXII  +  kc_pkaxii * PKabXII  -  ka_pka_ix * IX * PKa  +  kd_pkaix * PKabIX  +  kc_pkaix * PKabIX  -  ka_pka_c1inh * PKa * C1INH  -  ka_s_pka_ix * IX_s * PKa  +  kd_s_pkaix * PKabIX_s  +  kc_s_pkaix * PKabIX_s  +  koff_sn_pkahk_split * PKabHK_sn  -  ka_pka_hk * PKa * HK  +  kd_pkahk * PKabHK  +  kc_sn_pkahkxii * PKabHK_snbXII_sn;

% PKabXII
 dy(130)  =  +  ka_pka_xii * XII * PKa  -  kd_pkaxii * PKabXII  -  kc_pkaxii * PKabXII;

% XIIabXI
 dy(131)  =  +  ka_xiia_xi * XI * XIIa  -  kd_xiiaxi * XIIabXI  -  kc_xiiaxi * XIIabXI;

% PKabIX
 dy(132)  =  +  ka_pka_ix * IX * PKa  -  kd_pkaix * PKabIX  -  kc_pkaix * PKabIX;

% XIIa_snbAT
 dy(133)  =  +  ka_sn_xiia_at * AT * XIIa_sn;

% XIIabAT
 dy(134)  =  +  ka_xiia_at * AT * XIIa;

% C1INH
 dy(135)  =  -  ka_xia_c1inh * XIa * C1INH  -  ka_xiia_c1inh * XIIa * C1INH  -  ka_sn_xiia_c1inh * XIIa_sn * C1INH  -  ka_pka_c1inh * PKa * C1INH;

% XIabC1INH
 dy(136)  =  +  ka_xia_c1inh * XIa * C1INH;

% XIIabC1INH
 dy(137)  =  +  ka_xiia_c1inh * XIIa * C1INH;

% XIIa_snbC1INH
 dy(138)  =  +  ka_sn_xiia_c1inh * XIIa_sn * C1INH;

% PKabC1INH
 dy(139)  =  +  ka_pka_c1inh * PKa * C1INH;

% PKabIX_s
 dy(140)  =  +  ka_s_pka_ix * IX_s * PKa  -  kd_s_pkaix * PKabIX_s  -  kc_s_pkaix * PKabIX_s;

% HK
 dy(141)  =  -  ka_xi_hk * XI * HK  +  kd_xihk * XIbHK  +  koff_sn_xiahk * XIabHK_sn  -  ka_pk_hk * PK * HK  +  kd_pkhk * PKbHK  -  ka_pka_hk * PKa * HK  +  kd_pkahk * PKabHK  +  kc_sn_pkahkxii * PKabHK_snbXII_sn;

% XIbHK
 dy(142)  =  +  ka_xi_hk * XI * HK  -  kd_xihk * XIbHK  -  kon_sn_xihk * SN * XIbHK  +  koff_sn_xihk * XIbHK_sn;

% XIbHK_sn
 dy(143)  =  +  kon_sn_xihk * SN * XIbHK  -  koff_sn_xihk * XIbHK_sn  -  ka_sn_xihk_xiia * XIIa_sn * XIbHK_sn  +  kd_sn_xihkxiia * XIbHK_snbXIIa_sn;

% XIbHK_snbXIIa_sn
 dy(144)  =  +  ka_sn_xihk_xiia * XIIa_sn * XIbHK_sn  -  kd_sn_xihkxiia * XIbHK_snbXIIa_sn  -  kc_sn_xihkxiia * XIbHK_snbXIIa_sn;

% XIabHK_sn
 dy(145)  =  +  kc_sn_xihkxiia * XIbHK_snbXIIa_sn  -  koff_sn_xiahk * XIabHK_sn;

% PKbHK
 dy(146)  =  +  ka_pk_hk * PK * HK  -  kd_pkhk * PKbHK  -  kon_sn_pkhk * SN * PKbHK  +  koff_sn_pkhk * PKbHK_sn;

% PKbHK_sn
 dy(147)  =  +  kon_sn_pkhk * SN * PKbHK  -  koff_sn_pkhk * PKbHK_sn  -  ka_sn_pkhk_xiia * XIIa_sn * PKbHK_sn  +  kd_sn_pkhkxiia * PKbHK_snbXIIa_sn;

% PKbHK_snbXIIa_sn
 dy(148)  =  +  ka_sn_pkhk_xiia * XIIa_sn * PKbHK_sn  -  kd_sn_pkhkxiia * PKbHK_snbXIIa_sn  -  kc_sn_pkhkxiia * PKbHK_snbXIIa_sn;

% PKabHK_sn
 dy(149)  =  +  kc_sn_pkhkxiia * PKbHK_snbXIIa_sn  -  koff_sn_pkahk_split * PKabHK_sn  +  kon_sn_pkahk * SN * PKabHK  -  koff_sn_pkahk * PKabHK_sn  -  ka_sn_pkahk_xii * XII_sn * PKabHK_sn  +  kd_sn_pkahkxii * PKabHK_snbXII_sn;

% HKa
 dy(150)  =  +  koff_sn_pkahk_split * PKabHK_sn;

% BK
 dy(151)  =  +  koff_sn_pkahk_split * PKabHK_sn;

% PKabHK
 dy(152)  =  +  ka_pka_hk * PKa * HK  -  kd_pkahk * PKabHK  -  kon_sn_pkahk * SN * PKabHK  +  koff_sn_pkahk * PKabHK_sn;

% PKabHK_snbXII_sn
 dy(153)  =  +  ka_sn_pkahk_xii * XII_sn * PKabHK_sn  -  kd_sn_pkahkxii * PKabHK_snbXII_sn  -  kc_sn_pkahkxii * PKabHK_snbXII_sn;





end