% Kinetic Parameters 
kon_x = 1; 
koff_x = 1; 
kflow = 1; 
k_pla_plus = 0.5; 
k_pla_minus = 10.0; 
k_pla_act = 1; 
kact_e2 = 1; 
kon_IIa_sp = 1; 
koff_IIa_sp = 1; 
kon_v_sp = 1; 
koff_v_sp = 1; 
kon_se = 1; 

p = [ kon_x, koff_x, kflow, k_pla_plus, k_pla_minus, k_pla_act, kact_e2, kon_IIa_sp, koff_IIa_sp, kon_v_sp, koff_v_sp, kon_se ];


% Binding Site Parameters 
nbs_x = 1; 
np2 = 200.0; 
np5 = 500.0; 
nv = 1; 
nhv = 1; 

nbs = [ nbs_x, np2, np5, nv, nhv ];


% Function Arguments 
e2P = 0.001; 
VolP = 10.0; 

otherArgs = [ e2P, VolP ];


