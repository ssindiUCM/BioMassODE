% Initial Conditions 
L_TF_IC = 0; 
X_IC = 0; 
X_st_IC = 0; 
IIa_IC = 5.0; 
V_IC = 0.0; 
PL_IC = 0; 
P_SUB_IC = 10.0; 
PL_S_IC = 0; 
PL_V_IC = 0; 
p2_IC = 0; 
IIa_p_IC = 0; 
p5_IC = 0; 
V_p_IC = 0; 
Vh_IC = 0; 

init_cond = [ L_TF_IC, X_IC, X_st_IC, IIa_IC, V_IC, PL_IC, P_SUB_IC, PL_S_IC, PL_V_IC, p2_IC, IIa_p_IC, p5_IC, V_p_IC, Vh_IC ];


% Flow Rate Parameters 
IIa_up = IIa_IC; 
V_up = V_IC; 
PL_up = 1.0; 

flowUp = [ IIa_up, V_up, PL_up ];


