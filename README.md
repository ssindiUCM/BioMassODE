# BioMassODE
    Python script to generate MATLAB code for Coagulation Biochemical Reactions.

    Usage:
        python3 createCoagModel.py StaticCoag.txt

    Input File: StaticCoag.txt (List of Biochemical Reactions)
    
    Output File(s): Separates outputs for initial conditions, parameters, and code.
        (Assumes prefix based on input file)
        - StaticCoagMatlab.m (assumes prefix).
        - StaticCoagIC.m
        - StaticCoagParams.m
        - StaticCoagRename.m

    -*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-
    
    Input File Specifications:
    ----------------------------
    There are 3 classes of objects that can be specified:
    (1) [OPTIONAL] Variable Specifications: Initial Conditions, Parameter Values, Special Species
    (2) [OPTIONAL] Functions (to be used later for dilution/non-mass action biochemical kinetics)
    (3) [REQUIRED] Biochemical Reactions
    
    ----------------------------
    (1) Variable Specification:
    ----------------------------
    - Comments/Whitespace: Anything following "#" is ignored; Whitespace is ignored.
    - Parameter/Species names must start with a letter (no leading numbers).
    - Parameter/Species names may contain: [0-9a-zA-Z_:]

    - Users can define:
        * Initial Conditions (must be a non-negative real number; DEFAULT = 0)
            Ex:
                IIa_IC  = 5.0;  #mu M
                V_IC    = 0;
        * Parameter Values (real values of Initial Conditions; DEFAULT = 1)
            Ex:
                IIa_up  = IIa_IC #mu M; Let's see if this works!
                V_up    = V_IC
                PL_up   = 1.0    #Set to be a random value.
        * Special classes of species: LIPID, PLATELET, PLATELET_SITE
            Ex:
                PL      = PLATELET          #Platlet in Solution
                PL_S    = PROCOAG_PLATELET  #Platelet in Subendothelium
                PL_V    = PROCOAG_PLATELET  #Platelet in Volume
                p2      = PLATELET_SITE
                p5      = PLATELET_SITE
                
    ---------------
    (2) Functions:
    ---------------
    - Functions are defined as a name, list of arguments (comma separated) and body
    - Function arguments can be a parameter, species name or a dummy variable (dummy:x)
    - Examples:
        FUNCTION A(IIa,e2P) = IIa/(e2P + IIa) #1 nM = 0.001 mu M
        FUNCTION B(dummy:x, e2P) = x/(e2P + x)
    
    --------------
    (3) Diultion:
    --------------
    - A flag that turns on a diultion equation for every species.
    - The diultion function (and any dependent reactions) must already be defined.
    
    Ex:
    DILUTION = Dilution(VolP,PL, PL_S, PL_V,IIa,k_pla_act,k_pla_plus,kact_e2,e2P)

    FUNCTION Dilution(VolP,PL, PL_S, PL_V,IIa,k_pla_act,k_pla_plus,kact_e2,e2P) = (VolP)/((1-VolP)*(PL_S + PL_V))*dPdt(PL, PL_S, PL_V,IIa,k_pla_act,k_pla_plus,kact_e2,e2P)

    FUNCTION dPdt(PL, PL_S, PL_V,IIa,k_pla_act,k_pla_plus,kact_e2,e2P) = +  k_pla_act * PL * PL_S  + k_pla_act * PL * PL_V  +  kact_e2 * A(IIa,e2P) * PL +  k_pla_plus * PL * P_SUB

    ----------------------------
    (4) Biochemical Equations:
    ----------------------------
    
    - Assumes biochemical reactions (one per line) of the form:
            LHS <-> RHS, Rates, TYPE
    - Supports forward and reversible reactions.
    - Reaction Operators allowed: '*', '+', '<->', '->'
    - 5 Reaction Types:
        * MASS_ACTION: Default (if no type given)
        * FLOW:        species entering or exiting reaction zone.
        * LIPID:       Binding on/off lipid (competition)
        * FUNCTION:    Concentration changing due to reaction zone

    Example Reactions:
        A + 2 * B -> C , k_1 #Forward
        A + 2 * B <-> C , k_1 , k_2 #Reversible
        
    Example Initial Conditions/Parameter Setting:
        A_IC = 10.0 #\\mu M
        B_IC = 1.0  #\\mu M

    Lipid Binding Support:
        - Allows non-mass action terms for lipid binding:
            L_TF + II <-> II_st, kon_ii, koff_ii, nbs_ii, LIPID
        - kon_ii units: 1/(concentration * time * binding sites)

    Flow Reactant Support:
        - Allows specification of flow species as: list, reactions
        - Requires that the upstream species for S has the name S_up
   
    Example Flow List:
        FLOW, kflow, IIa
    
    Example Flow Reactions:
          -> K, kflow, K_up, FLOW #Species Flowing In
        K ->  , kflow, FLOW       #Flowing Flowing Out
        
    Supported Features:
    -----------------------------------
      - Input file lines can be in any order. (For cases w/ duplicates first entry read is retained)
      - Support for non-mass action lipid binding binding.
      - Outputs lipid/platlet binding sites as a separate parameter vector (nbs)
      - Can handle non-constant coefficients on the RHS for platelet sites and stores.
      - Outputs Species and rates output in input order.
      - Supports pure synthesis/degradation/in-out flow (e.g., "-> A", "B ->").
      - Consolidates duplicate kinetic rates.
      - Splits stoichiometric matrix for A + B -> A + C reactions.
      - Removes duplicate reactions (even one side of bidirectional ones).
      - Checks reaction rate dimensions for consistency.
      - Initial conditions can be set in the input file.

    In-Progress Features (Still Working On):
    -----------------------------------
      - Support for inline kinetic rate values:
        Example:
          A + B -> C, k1=0.1
          A + B <-> C, kon, koff=100

    Feature to Consider Developing:
    -----------------------------------
    - Allow setting rates and initial conditions seprate external file.
    - Improve MATLAB text wrapping for long lines.
    - Add back support for Python code.

    Current Version:
        Suzanne Sindi, 05/06/2025

Usage: python3 createMatlabFile.py StaticCoag.txt
