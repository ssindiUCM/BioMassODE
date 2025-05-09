#!/usr/bin/env python3

import sys, re, io, os, textwrap, math, csv
import numpy as np
from collections import defaultdict
from collections import Counter

###############################
# Error and Usage Information #
###############################

# Check if the correct number of arguments is provided
if len(sys.argv) < 2:
    print(textwrap.dedent("""

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
                PL      = PLATELET          #Platelet in Solution
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
    (3) Dilution:
    --------------
    - A flag that turns on a dilution equation for every species.
    - The dilution function (and any dependent reactions) must already be defined.
    
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
        - Allows specification of flow species in two ways: list, reactions
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
      - Outputs lipid/platelet binding sites as a separate parameter vector (nbs)
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
      - Allow setting rates and initial conditions separate external file.
      - Improve MATLAB text wrapping for long lines.
      - Add back support for Python code.
    
    Current Version:
        Suzanne Sindi, 05/06/2025

    """))
    sys.exit("Usage: python3 createCoagModel.py StaticCoag.txt")

########################
# Supporting Functions #
########################
def getSubset(species_dict, *types):
    """
    Returns a subset of species from species_dict that match one or more specified types.
    
    Parameters:
    - species_dict: The dictionary of species (e.g., specialSpecies).
    - *types: One or more species types to filter by.
    
    Returns:
    - A dictionary containing species matching the specified types.
    """
    # Ensure types is a tuple, even if only one type is passed
    types_set = set(types)

    # Filter the species_dict for entries matching the specified types
    subset = {
        name: info
        for name, info in species_dict.items()
        if info["type"] in types_set
    }

    return subset


def appendToSpecialSpecies(special_species_dict, species_list, species_type):
    for species in species_list:
        if species in special_species_dict:
            existing_type = special_species_dict[species]
            print(f"Error: Species '{species}' already assigned as type '{existing_type}', cannot reassign as '{species_type}'.")
        else:
            special_species_dict[species] = species_type

#Matlab can not handle ":"'s in variable names!
def transform_string(s):
    return s.replace(':', 'b')

#Get Args: A(IIa,e2P) ->
def getArgs(function_expr):
    match = re.search(r'\((.*?)\)', function_expr)
    if match:
        args_str = match.group(1)  # everything inside the parentheses
        return [arg.strip() for arg in args_str.split(',') if arg.strip()]
    return []

#Converts List to a String;
def list_to_string(lst):
    """Convert a list into a comma-separated string."""
    return ', '.join(map(str, lst))

# Flatten a list of lists into a single list
def flatten_list(lst_of_lists):
    """Flatten a list of lists into a single list."""
    return [item for sublist in lst_of_lists for item in sublist]

#Keep only unique entries of a list
def unique_entries_only(lst):
    """Return only unique entries from a list."""
    return list(set(lst))
    
#Keep only unique entries of a list IN ORDER
def unique_entries_in_order(lst):
    """Return only unique entries from a list, preserving their order."""
    seen = set()
    return [x for x in lst if not (x in seen or seen.add(x))]

def create_reaction(reactants, reactant_coeffs, products, product_coeffs, names, values, reaction_type, reaction_modifiers, reverse=False):
    if reverse:
        reactants, products = products, reactants
        reactant_coeffs, product_coeffs = product_coeffs, reactant_coeffs
        name = names[1]
        value = values[1]
    else:
        name = names[0]
        value = values[0]
    return Reaction(
        equation = formatFwdReaction(reactants, reactant_coeffs, products, product_coeffs),
        rateName = name,
        rateValue = value,
        reactants = reactants,
        reactant_coeffs = reactant_coeffs,
        products = products,
        product_coeffs = product_coeffs,
        reactionType = reaction_type,
        reactionModifiers = reaction_modifiers
    )

def formatFwdReaction(reactants, reactant_coeffs, products, product_coeffs):
    # Join reactants with their coefficients
    reactant_str = " + ".join(
        f"{coeff} * {reactant}" if coeff > 1 else reactant
        for reactant, coeff in zip(reactants, reactant_coeffs)
    )

    # Join products with their coefficients
    product_str = " + ".join(
        f"{coeff} * {product}" if coeff > 1 else product
        for product, coeff in zip(products, product_coeffs)
    )

    # Combine reactants and products into the final reaction string
    final_string = f"{reactant_str} -> {product_str}"

    return final_string


def parseFlowReactions(flowList, verbose=False):
    
    formattedReactions = []

    if verbose:
        print(f"\n{'=' * 50}")
        print(f"Step 1(d): Convert FLOW reactions to biochemical format")
        print(f"{'=' * 50}")

    for i, line in enumerate(flowList):
        parts = [p.strip() for p in line.split(',')]

        if len(parts) < 3:
            print(f"\tLine {i+1}: Invalid format (too few parts) -> {line}")
            continue

        keyword, rate, *species = parts

        if keyword != "FLOW":
            print(f"\tLine {i+1}: Invalid keyword '{keyword}', expected 'FLOW'")
            continue

        for sp in species:
            # First reaction: species exits the system
            outRxn = f"{sp} -> , {rate}, FLOW"
            # Second reaction: species appears along with its "up" version
            sp_up = f"{sp}_up"
            inRxn = f"-> {sp}, {rate}, {sp_up}, FLOW"

            formattedReactions.append(outRxn)
            formattedReactions.append(inRxn)

            if verbose:
                print(f"\tLine {i+1}: {sp}")
                print(f"\t\tGenerated outRxn: {outRxn}")
                print(f"\t\tGenerated inRxn : {inRxn}")

    if verbose:
        print(f"{'=' * 50}")

    return formattedReactions


def parseReactions(reactions, plateletSites, verbose, verboseDetailed):

    if verbose:
        print(f"Step 2: Parsing the Biochemical Reacitons")

    if verbose:
        print(f"\tNumber of biochemical reactions: {len(biochemicalReactions)}")
        print(f"\tSome of the Biochemical Reactions:")
        num_reactions_to_print = min(previewVal, len(reactions))
        for i in range(num_reactions_to_print):
            print(f"\t\t" + biochemicalReactions[i])
        print('-' * 50)

    parsed_reactions  = []
    newRateParameters = []

    for reaction in reactions:
        if verbose: print(f"Processing Reaction: '{reaction}'")
        
        # Parse the equation
        result = parseEquation(reaction,plateletSites,verboseDetailed)
        
        # Ensure result is valid before unpacking
        if not result or all(val is None for val in result):
            print(f"\tError: Invalid equation format in reaction: {reaction}")
            continue  # Skip this reaction
            
        (reactants,reactantCoeffs,products,productCoeffs,rates,reaction_type,reaction_modifiers)=result
        if verboseDetailed:
            print(f"\tReactants = {reactants})")
            print(f"\tCoeffs = {reactantCoeffs})")
            print(f"\tProducts = {products}")
            print(f"\tProducts Coeffs = {productCoeffs}")
            print(f"\tRates = {rates}")
            print(f"\tReaction Type = {reaction_type}")
            print(f"\tReaction Modifiers = {reaction_modifiers}")
            print(f"**********")

        # Store (if needed) names and values for each rate
        names  = []
        values = []
        for key, rate in rates.items():
            if '=' in rate:
                name, value = rate.split('=')
                names.append(name.strip())  # Store the name part
                values.append(float(value.strip()))  # Convert the value to a float
                newRateParameters.append(Parameter(name, value))
            else:
                names.append(rate.strip())  # Store the rate if no '='
                values.append(-1)  # Assign -1 as a placeholder for missing values

        reactionCount = len(rates);

        if reactionCount == 1:  # We add only 1 case, easy
            parsed_reactions.append(create_reaction(reactants, reactantCoeffs, products, productCoeffs,names,values,reaction_type, reaction_modifiers) )
        elif reactionCount == 2:  # We add 2 objects for bi-directional reactions
            parsed_reactions.append(create_reaction(reactants, reactantCoeffs, products, productCoeffs,names,values,reaction_type, reaction_modifiers,reverse=False) )
            parsed_reactions.append(create_reaction(reactants, reactantCoeffs, products, productCoeffs,names,values,reaction_type, reaction_modifiers,reverse=True) )
        else:
            print(f"\tError: Invalid reaction count in: {reaction}")
            continue
    
    return parsed_reactions, newRateParameters

def parse_biochemical_equation(line):
    parts = [p.strip() for p in line.split(",")]

    if len(parts) < 2:
        print(f"Error: Invalid format {line} (must have at least equation and rate)")
        return None, None, None, None, None

    equation_part = parts[0]

    # Determine arrow type and split equation
    if '<->' in equation_part:
        arrow = '<->'
        if len(parts) < 3:
            print(f"Error: Bidirectional reaction must include both FWD and REV rates: {line}")
            return None, None, None, None, None
        lhs, rhs = [s.strip() for s in equation_part.split('<->')]
        rates = {
            "FWD_RATE": parts[1],
            "REV_RATE": parts[2]
        }
        other = parts[3:]

    elif '->' in equation_part:
        arrow = '->'
        lhs, rhs = [s.strip() for s in equation_part.split('->')]
        rates = {
            "RATE": parts[1]
        }
        other = parts[2:]

    else:
        print(f"Error: Equation must contain '->' or '<->': {line}")
        return None, None, None, None, None

    return lhs, rhs, arrow, rates, other

def parse_reaction_type(other, reactants, products, coefficients):
    known_types = {"LIPID", "FLOW", "MASS_ACTION", "FUNCTION"}
    reactionModifiers = {}

    if len(other) == 0:
        maybe_type = "MASS_ACTION"
    else:
        maybe_type = other[-1]

    if maybe_type not in known_types:
        print(f"Error: Unknown reaction type '{maybe_type}'")
        return None, None
    
    # Handle based on type
    if maybe_type == "LIPID":
        if len(other) != 2:
            print(f"Error: LIPID reactions must have exactly one modifier before the type — got {other[:-1]}")
            return None, None
        reactionModifiers["bindingSites"] = other[0]
    
    elif maybe_type == "FLOW":
        if len(other) == 2 and len(reactants) == 0:  #In: ->K, kflow, Other={K_up, FLOW
            reactionModifiers["upstream"] = other[0]
        elif len(other) == 1 and len(reactants) == 1: #Out: K->, kflow, Other={FLOW}
            reactionModifiers["downstream"] = "outflow"
        else:
            print(f"Error: FLOW reaction doesn't have expected in/out flow structure")
            return None, None
        
    elif maybe_type == "FUNCTION":
        if len(other) > 1:
            function_expr = ",".join(other[:-1])  # Join all but the last element
            reactionModifiers["function"] = function_expr
            reactionModifiers["args"] = getArgs(function_expr)
            
        else:
            print(f"Error: FUNCTION reactions must include a function expression — got none")
            return None, None

    # MASS_ACTION doesn't need anything extra
    
    #Add if we have any values in products or coefficients
    for key, value in zip(products, coefficients):
        reactionModifiers[key] = value

    return maybe_type, reactionModifiers


def parseEquation(equationString, plateletSites, verbose=False):

    #Remove trailing whitespace
    equationString = equationString.split("#")[0].strip()
    
    #Pars into components: lhs arrow rhs, rates, other
    lhs, rhs, arrow, rates, other = parse_biochemical_equation(equationString)
    if verbose:
        print(f"NEW: Equation = {equationString}")
        print(f"NEW: LHS = {lhs}")
        print(f"NEW: RHS = {rhs}")
        print(f"NEW: Arrow = {arrow}")
        print(f"NEW: Rates = {rates}")
        print(f"NEW: Other = {other}")

    #(1) Process the LHS and RHS of the biochemical equation:
    reactants = [r.strip() for r in lhs.split('+')] if lhs.strip() else []
    products = [p.strip() for p in rhs.split('+')] if rhs.strip() else []

    # Extract_coefficients to process reactants and products
    # Note: We can only handle non integer coefficients on the RHS. Will throw an error
    reactantCoeffs, reactants, plateletSiteReactants, plateletSiteReactantCoeffs = extract_coefficients(reactants,plateletSites,"LHS")
    productCoeffs, products, plateletSiteProducts, plateletSiteProductCoeffs   = extract_coefficients(products,plateletSites,"RHS")
        
    # Check LHS
    if any(x is None for x in [
        reactantCoeffs, reactants, plateletSiteReactants, plateletSiteReactantCoeffs
    ]):
        print(f"Error: Malformed equation on LHS - {equation}")
        return None, None, None, None, None, None, None

    # Check RHS
    if any(x is None for x in [
        productCoeffs, products, plateletSiteProducts, plateletSiteProductCoeffs
    ]):
        print(f"Error: Malformed equation on RHS - {equation}")
        return None, None, None, None, None, None, None

    #(2) Parse OTHER:
    reactionType, reactionModifiers  = parse_reaction_type(other,reactants, plateletSiteProducts, plateletSiteProductCoeffs)

    if verbose:
        print(f"reactantCoeffs: {reactantCoeffs}")
        print(f"reactants: {reactants}")
        print(f"plateletSiteReactants: {plateletSiteReactants}")
        print(f"plateletSiteReactantCoeffs: {plateletSiteReactantCoeffs}")
        print(f"productCoeffs: {productCoeffs}")
        print(f"products: {products}")
        print(f"plateletSiteProducts: {plateletSiteProducts}")
        print(f"plateletSiteProductCoeffs: {plateletSiteProductCoeffs}")
        print(f"reactionType: {reactionType}")
        print(f"reactionModifiers: {reactionModifiers}")
    
    #(3) Check if we have the right number of rates.
    # Rate is a dictionary of either "RATE" or "FWD_RATE" "REV_RATE"
    reaction_count = len(rates)
    if reaction_count > 1 and reactionType not in {"MASS_ACTION", "LIPID"}:
        print(f"Error in Equation = {equation}: Only MASS_ACTION and LIPID reactionTypes can be bidirectional")
        return None, None, None, None, None, None, None

    if verbose: print(f"**************")

    return reactants, reactantCoeffs, products, productCoeffs, rates, reactionType, reactionModifiers

# Extract the Coefficients for a Set of Reactants
# Added: plateletSiteBool to check if there's plateletSites in the equation!
# This should only happen with type PLATELET_ACTIVATION
def extract_coefficients(terms, plateletSites, side=None):
    plateletSitesFound = []  # Stores species found in plateletSites
    plateletSiteCoeffs = []  # Stores corresponding coefficients

    # Check if the input contains only an empty string
    if len(terms) == 1 and terms[0] == '':
        return [], [], plateletSitesFound, plateletSiteCoeffs

    coeffs = []
    species = []

    for term in terms:
        term = term.strip()
        found_platelet_species = None
        platelet_coeff = "1"  # Default coefficient is 1

        if '*' in term:
            parts = term.split('*')
            for part in parts:
                stripped_part = part.strip()
                if stripped_part in plateletSites:
                    found_platelet_species = stripped_part
                else:
                    platelet_coeff = stripped_part
        else:
            if term in plateletSites:
                found_platelet_species = term

        if found_platelet_species:
            # Error check if we're on the LHS
            if side == "LHS":
                try:
                    int_val = int(platelet_coeff)  # Check if it's numeric
                except ValueError:
                    print(f"Error: Non-numeric coefficient '{platelet_coeff}' for platelet site '{found_platelet_species}' on LHS of equation {terms}.")
                    return None, None, None, None
                else:
                    coeffs.append(int_val)  # Store numeric version
            else:
                coeffs.append(1)  # On RHS or general case, use default 1

            species.append(found_platelet_species)
            plateletSitesFound.append(found_platelet_species)
            plateletSiteCoeffs.append(platelet_coeff)  # Store original value as string

        else:
            # Match terms with coefficients followed by '*'
            match = re.match(r'(\d*)\s*\*\s*(\S+)', term)
            if match:
                coeff_str, species_name = match.groups()
                coeff = int(coeff_str) if coeff_str else 1
            else:
                # Match terms with optional coefficients without '*'
                match = re.match(r'(\d*)\s*(\S+)', term)
                if match:
                    coeff_str, species_name = match.groups()
                    coeff = int(coeff_str) if coeff_str else 1
                else:
                    print(f"Error: Invalid term format: {term}")
                    coeffs.append(1)
                    species.append(term.strip())
                    continue
            
            coeffs.append(coeff)
            species.append(species_name.strip())

    return coeffs, species, plateletSitesFound, plateletSiteCoeffs

def parseInputFile(verbose=False):
    # Initialize biochemical arrays and counters
    biochemicalReactions = []
    initialConditions = []      # User-specified initial conditions
    specialSpecies = []         # Special Species: Lipid, Platelet, Sites, Stores
    parameterCondition = []     # User-specified parameters (not in-line)
    functionCondition = []      # Dilution and platelet activation
    flowList = []
    dilutionCondition = []      # Where we store dilution function

    numReactionsReadIn = 0
    numInitialConditionsReadIn = 0
    numLipidSpecies = 0
    numPlateletSpecies = 0
    numPlateletSites = 0
    numPlateletStores = 0
    numParametersReadIn = 0
    numFunctionsDefined = 0
    numFlowListDefined = 0
    numDilutionCondition = 0
    
    DILUTION = False #Flag for if we turn on dilution or not.

    if verbose:
        print(f"{'-' * 50}")
        print(f"Step 0 Preprocessing of {sys.argv[1]}")

    try:
        with open(sys.argv[1], 'r') as file:
            for line in file:
                line = line.rstrip()  # Remove trailing newline characters
                
                # Ignore blank lines or comments
                if not line or line.startswith('#'):
                    continue
                
                # Remove inline comments after '#'
                line = line.split('#')[0].rstrip()

                # Split the line by commas
                fields = line.split(',')
                
                ##(Q) Is this a special species?
                if "= LIPID" in line:
                    specialSpecies.append({"type": "LIPID", "line": line.replace("= LIPID", "").strip()})
                    numLipidSpecies += 1
                    if verbose: print(f"\t\tLipid Species: {line}")
                elif "= PLATELET_SITE" in line:
                    specialSpecies.append({"type": "PLATELET_SITE", "line": line.replace("= PLATELET_SITE", "").strip()})
                    numPlateletSites += 1
                    if verbose: print(f"\t\tPlatelet Site: {line}")
                elif "= PLATELET_STORE" in line:
                    specialSpecies.append({"type": "PLATELET_STORE", "line": line.replace("= PLATELET_STORE", "").strip()})
                    numPlateletStores += 1
                    if verbose: print(f"\t\tPlatelet Store: {line}")
                elif "= PLATELET" in line:
                    specialSpecies.append({"type": "PLATELET", "line": line.replace("= PLATELET", "").strip()})
                    numPlateletSpecies += 1
                    if verbose: print(f"\t\tPlatelet Species: {line}")
                elif "= PROCOAG_PLATELET" in line:
                    specialSpecies.append({"type": "PROCOAG_PLATELET", "line": line.replace("= PROCOAG_PLATELET", "").strip()})
                    numPlateletSpecies += 1
                    if verbose: print(f"\t\tPlatelet Species: {line}")

                #(Q2): Do we start with function? Then it's a function!
                elif line.lstrip().startswith("FUNCTION"):
                    functionCondition.append(line)
                    numFunctionsDefined += 1
                    if verbose: print(f"\t\tFunctions: {line}")
                    
                #(Q3): Do we start the line with FLOW? Then it's a flow LIST!
                elif line.lstrip().startswith("FLOW"):
                    flowList.append(line);
                    numFlowListDefined += 1
                    if verbose: print(f"\t\tFlow List: {line}")

                #(Q4): Do we start with DIULTION? Then DILUTION = True and we go!
                elif line.lstrip().startswith("DILUTION"):
                    DILUTION = True;
                    dilutionCondition.append(line)
                    if verbose: print(f"\t\tDilution is True: {line}")
            
                #(Q5): Is this an intial condition?
                elif "_IC" in line.split('=')[0].strip():  # Check if it's an initial condition
                    initialConditions.append(line)
                    numInitialConditionsReadIn += 1
                    if verbose: print(f"\t\tInitial Condition: {line}")
                
                #(Q6): Does it contain a comma? Then Try for biochemical equation
                elif len(fields) > 1:
                    biochemicalReactions.append(line)
                    numReactionsReadIn += 1
                    if verbose: print(f"\t\tBiochemical Equation: {line}")
                else:
                   # Assume it's a parameter if it's neither a species nor an initial condition
                    parameterCondition.append(line)
                    numParametersReadIn += 1
                    if verbose: print(f"\t\tParameter: {line}")

    except IOError:
        sys.exit(f"Couldn't open {sys.argv[1]}")

    # Print output summary
    if verbose:
        print(f"DONE: Step 0 Preprocessing of {sys.argv[1]}")
        print(f"{'-' * 50}")

        # Define headers
        headers = ["Category", "Count"]
        data = [
            ("Biochemical Reactions", numReactionsReadIn),
            ("Initial Conditions", numInitialConditionsReadIn),
            ("Lipid Species", numLipidSpecies),
            ("Platelet Binding Sites",numPlateletSites),
            ("Platelet Stores", numPlateletStores),
            ("Platelet Species", numPlateletSpecies),
            ("Parameters", numParametersReadIn),
            ("Functions", numFunctionsDefined),
            ("Flow List", numFlowListDefined),
            ("Dilution Status", DILUTION )
        ]

        # Print table header
        print(f"\n{'=' * 50}")
        print(f"{'Step 0: Processed Input File':^30}")
        print(f"{'=' * 50}")
        print(f"{headers[0]:<30} | {headers[1]:>10}")
        print("-" * 50)

        # Later, during printing:
        for category, count in data:
            display_value = str(count) if isinstance(count, bool) else count
            if count>0:
                print(f"{category:<30} | {display_value:>10}")

        # Print done message and finish with a line
        print(f"{'=' * 50}\n")

    return {
        "biochemicalReactions": biochemicalReactions,
        "initialConditions": initialConditions,
        "specialSpecies": specialSpecies,
        "parameterCondition": parameterCondition,
        "functionCondition": functionCondition,
        "flowList": flowList,
        "DILUTION" : DILUTION,
        "dilutionCondition": dilutionCondition,
        "counts": {
            "numReactionsReadIn": numReactionsReadIn,
            "numInitialConditionsReadIn": numInitialConditionsReadIn,
            "numLipidSpecies": numLipidSpecies,
            "numPlateletSpecies": numPlateletSpecies,
            "numPlateletSites": numPlateletSites,
            "numParametersReadIn": numParametersReadIn,
            "numFunctionsDefined": numFunctionsDefined,
            "numFlowListDefined": numFlowListDefined,
        }
    }


def parseSpecies(lines, verbose=False):
    allowed_types = { "PLATELET", "PROCOAG_PLATELET", "LIPID",
        "PLATELET_SITE", "PLATELET_STORE"}

    species_dict = {}
    error_flag = False

    for entry in lines:
        species_type = entry.get('type')
        line = entry.get('line', '').strip()

        if species_type not in allowed_types:
            print(f"Error: Invalid species type '{species_type}' in entry: {entry}")
            error_flag = True
            continue

        parts = [p.strip() for p in line.split(',')]

        if not parts or not parts[0]:
            print(f"Error: Empty or invalid species name in entry: {entry}")
            error_flag = True
            continue

        name = parts[0]
        modifier = parts[1] if len(parts) > 1 else ""

        if name in species_dict:
            print(f"Error: Duplicate species name '{name}' detected.")
            error_flag = True
            continue

        if species_type in {"PLATELET_STORE", "PLATELET_SITE"} and modifier == "":
            print(f"Error: Species '{name}' of type '{species_type}' must have a modifier.")
            error_flag = True
            continue

        species_dict[name] = {
            "type": species_type,
            "modifier": modifier
        }

    # Verbose: Print nicely formatted table
    if verbose:
        print(f"{'=' * 50}")
        print(f"{'Step 1(b): Special Species':^30}")
        print(f"{'=' * 50}")
        print(f"{'Name':<12}| {'Type':<20} | Modifier")
        print("-" * 50)

        # Sort by type, then by name
        sorted_items = sorted(species_dict.items(), key=lambda x: (x[1]["type"], x[0]))

        for name, info in sorted_items:
            type_str = info["type"]
            mod_str = info["modifier"] if info["modifier"] else "[None]"
            print(f"{name:<12}| {type_str:<20} | {mod_str}")

        print("=" * 50)

    return species_dict, error_flag


def parseInitialConditions(initialConditions, verbose=False):
    
    parsed_ICs = []
    seen_names = set()  # A set to track names that have been parsed
    duplicates = []     # List to store any duplicates found

    for ic_str in initialConditions:
        ic_str = ic_str.strip()  # Remove leading/trailing spaces
        
        # Remove comments (anything after #)
        ic_str = ic_str.split("#", 1)[0].strip()
        
        if not ic_str:
            continue  # Skip empty or fully commented-out lines

        # Ensure there is exactly one '='
        if ic_str.count("=") != 1:
            print(f"\t⚠️ Error: Invalid format in '{ic_str}'. Must contain exactly one '='.")
            continue

        try:
            name, value_str = ic_str.split("=")  # Guaranteed to have one '=' due to check above
            name = name.strip()
            value_str = value_str.strip()

            # Ensure the name matches the correct format (String_IC)
            if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*_IC$", name):
                corrected_name = re.sub(r"(_[iI][cC])$", "_IC", name)  # Fix capitalization
                if corrected_name != name:
                    print(f"\t🔄 Auto-corrected '{name}' to '{corrected_name}'")
                    name = corrected_name
                else:
                    print(f"\t⚠️ Error: '{name}' is invalid. Should match 'String_IC'. Did you mean '{name}_IC'?")
                    continue  # Skip invalid entries

            # Clean the value string: Remove any extra characters or units (e.g., semicolon)
            value_str = ''.join(filter(lambda x: x.isdigit() or x == '.', value_str))  # Keep digits and decimal points
            
            # Convert value to float and check if non-negative
            value = float(value_str)
            if value < 0:
                print(f"\t⚠️ Error: Value for '{name}' must be non-negative. Found {value}.")
                continue

            # Check for duplicates
            if name in seen_names:
                duplicates.append(name)  # Add to duplicates list
                continue  # Skip adding this duplicate to the final list
            else:
                seen_names.add(name)  # Add the name to the set of seen names

            # Create InitialCondition object and append to the result list
            parsed_ICs.append(InitialCondition(name, value))

        except ValueError:
            print(f"\t⚠️ Error: Could not parse value in '{ic_str}'. Ensure RHS is a non-negative number.")
        except Exception as e:
            print(f"\t❌ Unexpected error while parsing '{ic_str}': {e}")

    # Report duplicates if any
    if duplicates:
        print(f"\t⚠️ Error: Found duplicates in initial conditions: {', '.join(duplicates)} (retaining only 1st value")


    if verbose:
        
        # Define headers for the initial conditions table
        ic_headers = ["Initial Condition Name", "Value"]
    
        # Create a data list for parsed initial conditions
        ic_data = [(ic.name, ic.value) for ic in parsed_ICs]  # Assuming parsed_ICs is a list of InitialCondition objects

        # Print table header for initial conditions
        print(f"{'=' * 50}")
        print(f"Step 1(a): Specified Intial Conditions")
        print(f"{'=' * 50}")
        print(f"{ic_headers[0]:<30} | {ic_headers[1]:>10}")
        print("-" * 50)

        # Print each parsed initial condition in the data list
        for name, value in ic_data:
            print(f"{name:<30} | {value:>10.3f}")  # Format value to 2 decimal places
        # Print done message and finish with a line
        print(f"{'=' * 50}\n")

    return parsed_ICs

def parseFunction(functionsDefined,verbose=False):
    parsed_functions = []
    keyword = "FUNCTION"
    function_pattern = re.compile(
        rf"^{keyword}\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(([^)]*)\)\s*=\s*(.+)$"
    )

    for func in functionsDefined:
        match = function_pattern.match(func.strip())

        if not match:
            if verbose:
                print(f"\tInvalid function format: {func}")
            continue  # Skip invalid functions

        name, args, body = match.groups()

        # Validate function name (MATLAB/Python compatibility)
        if not name.isidentifier():
            if verbose:
                print(f"\tInvalid function name: {name}")
            continue

        # Validate arguments (should be valid variable names, comma-separated)
        args_list = []
        dummy_args = []

        for raw_arg in args.split(","):
            arg = raw_arg.strip()
            if not arg:
                continue
            if arg.startswith("dummy:"):
                arg_name = arg[len("dummy:"):]
                dummy_args.append(arg_name)
                args_list.append(arg_name)
            else:
                args_list.append(arg)

        # Validate argument names
        if not all(arg.isidentifier() for arg in args_list):
            if verbose:
                print(f"\tInvalid function arguments: {args}")
            continue

        # Validate expression (basic check: should not be empty)
        body = body.strip()
        if not body:
            if verbose:
                print(f"\tInvalid function body: {func}")
            continue

        # Store valid function details
        parsed_functions.append({"name": name, "args": args_list, "dummy_args": dummy_args,"body": body})

        #if verbose:
            #print(f"\tValid function parsed: {name}({', '.join(args_list)}) = {body}")

    if verbose:
        # Create a data list for parsed initial conditions
        print(f"\n{'=' * 50}")
        print(f"Step 1(c): Parsed Functions")
        print(f"{'=' * 50}")
        for fVal in parsed_functions:  # Assuming parsed_ICs is a list of InitialCondition objects
            # Print table header for initial conditions
            print(f"{fVal}")
            
    return parsed_functions
    
def parseDilution(dilutionCondition,verbose=False):
    #Should look like this:
    #DILUTION = FUNCTION Dilution(Vp,PL, PL_S, PL_V,IIa,k_pla_act,k_pla_plus,kact_e2,e2P)
    parsed_dilution = []

    # Pattern: DILUTION = FunctionName(arg1, arg2, ...)
    dilution_pattern = re.compile(
        r"^\s*DILUTION\s*=\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\(([^)]*)\)\s*$"
    )

    for func in dilutionCondition:
        match = dilution_pattern.match(func.strip())

        if not match:
            if verbose:
                print(f"\tInvalid dilution format: {func}")
            continue  # Skip invalid functions

        name, args = match.groups()

        # Validate function name (MATLAB/Python compatibility)
        if not name.isidentifier():
            if verbose:
                print(f"\tInvalid dilution function name: {name}")
            continue

        # Validate arguments (should be valid variable names, comma-separated)
        args_list = [arg.strip() for arg in args.split(",") if arg.strip()]
        if not all(arg.isidentifier() for arg in args_list):
            if verbose:
                print(f"\tInvalid dilution function arguments: {args}")
            continue

        # Store valid function details
        parsed_dilution.append({"name": name, "args": args_list})

    if len(parsed_dilution)>1:
        print(f"Warning: More than 1 diultion function given. Will only use the first.")

    if verbose:
        # Create a data list for parsed initial conditions
        print(f"{'=' * 50}")
        print(f"Step 1(g): Dilution Function")
        print(f"{'=' * 50}")
        for fVal in parsed_dilution:  # Assuming parsed_ICs is a list of InitialCondition objects
            # Print table header for initial conditions
            print(f"{fVal}")
        
    return parsed_dilution


def parseParameters(parameterLines, verbose=False):

    if verbose:
        print("-" * 50)
        print(f"Step 1(f): Parse Specified Parameters")
    
    parsed_params = []
    seen_names = set()  # Track names to detect duplicates
    duplicates = []  # Store duplicate names

    for param_str in parameterLines:
        param_str = param_str.strip()  # Remove leading/trailing spaces

        # Remove comments (anything after #)
        param_str = param_str.split("#", 1)[0].strip()
        
        if not param_str:
            continue  # Skip empty or fully commented-out lines

        # Ensure there is exactly one '='
        if param_str.count("=") != 1:
            print(f"\t⚠️ Error: Invalid format in '{param_str}'. Must contain exactly one '='.")
            continue

        try:
            name, value_str = param_str.split("=")  # Guaranteed to have one '='
            name = name.strip()
            value_str = value_str.strip()

            # Ensure the name is a valid identifier (letters, numbers, underscores, but cannot start with a number)
            if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", name):
                print(f"\t⚠️ Error: Invalid parameter name '{name}'. Must start with a letter or underscore.")
                continue

            # Attempt to parse value as a number
            if re.match(r"^[-+]?(?:\d*\.\d+|\d+)$", value_str):  # Matches integers, decimals, and .5 style numbers

                value = float(value_str)
                if value < 0:
                    print(f"\t⚠️ Error: Parameter '{name}' must be non-negative. Found {value}.")
                    continue
            else:
                value = value_str  # Keep as a string if it's not a valid number

            # Check for duplicates
            if name in seen_names:
                duplicates.append(name)
                print(f"\t⚠️ Warning: Duplicate parameter '{name}' found. Keeping first occurrence.")
                continue
            else:
                seen_names.add(name)  # Mark as seen

            # Create Parameter object and append to list
            parsed_params.append(Parameter(name, value))

        except ValueError:
            print(f"\t⚠️ Error: Could not parse value in '{param_str}'. Ensure RHS is a valid number or string.")
        except Exception as e:
            print(f"\t❌ Unexpected error while parsing '{param_str}': {e}")

    # Report duplicates if any
    if duplicates:
        print(f"\t⚠️ Warning: Found duplicate parameters: {', '.join(duplicates)} (retaining only 1st value).")
        
    if verbose:
        print(f"DONE: 1(f): Parse Specified Parameters")
        print(f"{'-' * 50}")
        
        # Define headers for the parameters table
        p_headers = ["Parameter Name", "Value", "Type"]
    
        # Create a data list for parsed parameters, including type
        p_data = [(p.name, p.value, "Float" if isinstance(p.value, float) else "String") for p in parsed_params]

        # Print table header
        print(f"\n{'=' * 65}")
        print(f"{'Step 1(f): Parse Specified Parameters':^50}")
        print(f"{'=' * 65}")
        print(f"{p_headers[0]:<30} | {p_headers[1]:>10} | {p_headers[2]:>10}")
        print("-" * 65)

        for name, value, value_type in p_data:
            if isinstance(value, float):
                if value.is_integer():
                    print(f"{name:<30} | {value:>10.1f} | {value_type:>10}")
                else:
                    print(f"{name:<30} | {value:>10.1e} | {value_type:>10}")
            else:
                print(f"{name:<30} | {value:>10} | {value_type:>10}")

        # Print done message and finish with a line
        print(f"{'=' * 65}\n")
        
    return parsed_params

def add_extra_parameters(parsed_params, extra_params, verbose = False):
    existing_names = {param.name for param in parsed_params}
    
    for extra in extra_params:
        if extra.name in existing_names:
            print(f"Already exists in parsedParameters - will retain the value already stored: {extra.name}")
        else:
            try:
                # Attempt to convert the value to float
                numeric_value = float(extra.value)
                parsed_params.append(Parameter(extra.name, numeric_value))
                if verbose:
                    print(f"- Adding value for in-line defined parameter: {Parameter(extra.name, numeric_value)}")
            except (ValueError, TypeError):
                print(f"Warning: Parameter '{extra.name}' has a non-numeric value: {extra.value} (not added)")


def validate_parameters(parsedParameters, uniqueSpecies, verbose = False):
    numErrors = 0
    if verbose:
        print(f"(Q): Do all parameters defined as a string occur as species_IC?")
        
    for param in parsedParameters:
        if verbose:
            print(f"Parameter name: {param.name}\tParameter Value: {param.value}")
        if isinstance(param.value, str):  # Check if value is a string
            match = re.fullmatch(r"([A-Za-z0-9_]+)_IC", param.value)
            if match:
                species_name = match.group(1)  # Extract the species part
                if species_name not in uniqueSpecies:
                    print(f"\tWarning: Parameter '{param.name}' refers to an unknown species '{species_name}' in value '{param.value}'.")
                    numErrors = numErrors + 1
            else:
                print(f"\tWarning: Parameter '{param.name}' has an invalid string value '{param.value}'. Expected format: 'STRING_IC'.")
                numErrors = numErrors + 1

        else:  # Check if value is a float or convertible to float
            try:
                float(param.value)  # Attempt conversion
            except ValueError:
                print(f"\tWarning: Parameter '{param.name}' has an invalid numeric value '{param.value}'.")
                numErrors = numErrors + 1
        
    if verbose:
        if numErrors == 0:
            print(f"\tYes...Finished parsing with no warnings!")
        else:
            print(f"\tProceeded with caution (reported {numErrors} warnings.")

def createBiochemicalMatrices(unique_species, parsed_reactions, flowRate, prefix, verbose=False):
    # Define the output file names
    stoich_output_file   = f"{prefix}Stoich.csv"
    reactant_output_file = f"{prefix}Reactants.csv"
    product_output_file  = f"{prefix}Products.csv"

    if verbose:
        print("-" * 50)
        print(f"Step 2: Creating Stoichiometry Matrix")

    s = Stoich(unique_species,parsed_reactions)
    
    bad_rates = s.check_rates()
    
    if len(set(flowRate)) > 1:
       if verbose: print(f"Error: Flow list contains multiple unique values: {flowRate}.")
       exit(-1)
    elif len(set(flowRate)) > 0:
        uniqueFlowRate = flowRate[0];
        if verbose: print(f"Processing Unique Flow Rate: {uniqueFlowRate}")
    else:
        if verbose: print(f"Model has NO Flow Rates")

    #Do we have any bad reaction rates?
    filtered_bad_rates = bad_rates.copy()
    #filtered_bad_rates.pop(uniqueFlowRate, None)  # Safely remove if it exists
    try:
        filtered_bad_rates.pop(uniqueFlowRate, None)
    except NameError:
        pass  # uniqueFlowRate was never defined, so do nothing
    
    if len(filtered_bad_rates) > 0:
        print(f"Error: We found some biochemical rates with different dimensions.")
        for rate, column_sum in filtered_bad_rates.items():
            print(f"Rate {rate}: Sum of columns = {column_sum}")
    else:
        if verbose: print(f"All rates are correct dimension")

    if verbose:
        s.to_csv(stoich_output_file,"S")
        s.to_csv(reactant_output_file,"R")
        s.to_csv(product_output_file,"P")
        print(f"DONE! Successfully created Stoichiometry Matrix (output file {stoich_output_file})")
        print('-' * 50)
    
    return s;

#################
# Class Objects #
#################

class InitialCondition:
    def __init__(self, name, value):
        self.name  = name
        self.value = value

    def __repr__(self):
        return f"InitialCondition(name={self.name}, value={self.value})"

class Parameter:
    def __init__(self, name, value):
        self.name  = name
        self.value = value
    
    def __repr__(self):
        return f"Parameter(name={self.name}, value={self.value})"

class Reaction:
    index_counter = 0  # Class variable to keep track of the index

    def __init__(self, equation, rateName, rateValue, reactants, reactant_coeffs, products, product_coeffs,reactionType, reactionModifiers): #, rateModifier,rateModifierValue):
        self.equation = equation
        self.rateName  = rateName
        self.rateValue = rateValue
        self.reactants = reactants
        self.reactant_coeffs = reactant_coeffs
        self.products = products
        self.product_coeffs = product_coeffs
        self.reactionType = reactionType
        self.reactionModifiers = reactionModifiers
        #self.rateModifier = rateModifier
        #self.rateModifierValue = rateModifierValue
        self.index = Reaction.index_counter  # Assign the current index
        Reaction.index_counter += 1  # Increment the index for the next reaction

    def __repr__(self):
        return (f"Reaction(index={self.index}, equation={self.equation}, rateName={self.rateName}, "
                f"rateValue={self.rateValue}, reactants={self.reactants}, "
                f"reactant_coeffs={self.reactant_coeffs}, products={self.products}, "
                f"product_coeffs={self.product_coeffs}, rxn_type={self.reactionType}, "
                f"reactionModifiers={self.reactionModifiers}")
                #f"rate_modifier = {self.rateModifier}, rate_modifierValue = {self.rateModifierValue}")

    def __eq__(self, other):
        if isinstance(other, Reaction):
            sorted_self_reactants = sorted(zip(self.reactants, self.reactant_coeffs))
            sorted_other_reactants = sorted(zip(other.reactants, other.reactant_coeffs))
            sorted_self_products = sorted(zip(self.products, self.product_coeffs))
            sorted_other_products = sorted(zip(other.products, other.product_coeffs))

            # Compare only the content (ignoring index)
            return (self.rateName == other.rateName and
                    self.rateValue == other.rateValue and
                    sorted_self_reactants == sorted_other_reactants and
                    sorted_self_products == sorted_other_products)
        return False

    def __hash__(self):
        sorted_reactants = tuple(sorted(zip(self.reactants, self.reactant_coeffs)))
        sorted_products = tuple(sorted(zip(self.products, self.product_coeffs)))
        
        # Hash based on content (ignoring index)
        return hash((self.rateName, self.rateValue, sorted_reactants, sorted_products))


class Stoich:

    def __init__(self, uniqueSpecies: list, parsedReactions: list) -> None:
        self.species = dict()
        c = 0
        for s in uniqueSpecies:
            if s not in self.species:
                self.species[s] = c
                c += 1
    
        self.rates = []
        for r in parsedReactions:
            self.rates.append(r.rateName)
    
        self.N = len(self.species)  #Species
        self.M = len(self.rates)    #Reactions
        
        self.stoich    = [[math.nan for i in range(self.M)] for j in range(self.N)]
        self.reactants = [[math.nan for i in range(self.M)] for j in range(self.N)]
        self.products  = [[math.nan for i in range(self.M)] for j in range(self.N)]

        for r, reaction in enumerate(parsedReactions):

            #Process Reactants
            for i,spec in enumerate(reaction.reactants):
                spec_ind = self.species[spec]
            
                # Subtract species multiplicities from stoich matrix
                if math.isnan(self.stoich[spec_ind][r]): self.stoich[spec_ind][r] = 0
                self.stoich[spec_ind][r]    += - reaction.reactant_coeffs[i]

                if math.isnan(self.reactants[spec_ind][r]): self.reactants[spec_ind][r] = 0
                self.reactants[spec_ind][r] += - reaction.reactant_coeffs[i]

            # Process Products
            for i, spec in enumerate(reaction.products):
                spec_ind = self.species[spec]

                # Add species multiplicities from stoich matrix
                if math.isnan(self.stoich[spec_ind][r]): self.stoich[spec_ind][r] = 0
                self.stoich[spec_ind][r] += reaction.product_coeffs[i]
                
                if math.isnan(self.products[spec_ind][r]): self.products[spec_ind][r] = 0
                self.products[spec_ind][r] += reaction.product_coeffs[i]
    
    # Define how to print stoich objects
    def __str__(self):
        return f"Stoich Info:\nSpecies: {list(self.species)}\nRates: {list(self.rates)}\nMatrix: {self.stoich}\nMatrix: {self.reactants}\nMatrix: {self.products}"
    
    def __repr__(self):
        return self.__str__()

    def check_rates(self):
        # Create a dictionary to map rates to their corresponding columns
        rate_to_columns = defaultdict(list)
        for i, rate in enumerate(self.rates):
            rate_to_columns[rate].append(i)

        # Calculate sum for each column where rates are the same
        column_sums = {}
        for rate, columns in rate_to_columns.items():
            # Initialize a list for storing sums, with length equal to the number of rows
            column_sum = [0] * len(columns)
            
            # Iterate over each column index for the current rate
            for j,col in enumerate(columns):
                for row in range(self.N):
                    # Sum the values in the column only if they are not NaN
                    if not math.isnan(self.reactants[row][col]):
                        column_sum[j] += self.reactants[row][col]
            # Store the sum for the current rate
            column_sums[rate] = column_sum
    
        # Identify bad rates
        bad_rates = {}
        for rate, column_sum in column_sums.items():
            unique_sums = set(column_sum)
            if len(unique_sums) > 1:
                bad_rates[rate] = column_sum

        return bad_rates
            
    def to_csv(self, filename: str, option: str):
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write the header row
            writer.writerow([''] + list(self.rates))  # Header with rate names
            
            if option == "S":
                # Write each row for species
                for i, species in enumerate(self.species):
                    row = [species] + [self.stoich[i][j] if not math.isnan(self.stoich[i][j]) else '' for j in range(self.M)]
                    writer.writerow(row)
            elif option == "R":
                # Write each row for species
                for i, species in enumerate(self.species):
                    row = [species] + [self.reactants[i][j] if not math.isnan(self.reactants[i][j]) else '' for j in range(self.M)]
                    writer.writerow(row)
            
            elif option == "P":
                # Write each row for species
                for i, species in enumerate(self.species):
                    row = [species] + [self.products[i][j] if not math.isnan(self.products[i][j]) else '' for j in range(self.M)]
                    writer.writerow(row)
            else:
                print(f"Error: Invalid Option Given: {option}")


#######################
# Create Output Files #
#######################

def write_species_ode(f, i, species, s, rates, parsed_reactions, specialSpecies, DILUTION, dil_string=None, verbose=False):
    Ns = len(s.species)
    Nr = len(s.rates)
    
    speciesName = species[i]
    species_info = specialSpecies.get(speciesName, {})
    species_type = species_info.get("type")
    
    isLipidSpecies = species_type == "LIPID"
    isPlateletSite = species_type == "PLATELET_SITE"
    isPlateletStore = species_type == "PLATELET_STORE"

    isOnEmptyLipid = "_s" in speciesName and not "_st" in speciesName and not "_sn" in speciesName
    isOnTFLipid = "_st" in speciesName
    isOnSilica = "_sn" in speciesName
    
    f.write("% ")
    f.write(str(transform_string(species[i])))
    f.write(f"\n d{transform_string(species[i])} = ")
    
    for j in range(Nr):  # Assuming Nr = len(rates)
        isLipidReaction    = (parsed_reactions[j].reactionType=="LIPID")
        isFlowReaction     = (parsed_reactions[j].reactionType=="FLOW")
        isFunctionReaction = (parsed_reactions[j].reactionType=="FUNCTION")
            
        #################
        # Case 1: Reaction j DECREASES the amount of species i;
        #################
        if (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) < 0):
            ##################
            #Part 1.1: Sign, Rate Constant & Stochiometric Change
            ##################
            f.write("  -  ")
            if abs(int(s.stoich[i][j]))>1:
                f.write(str(abs(int(s.stoich[i][j])))) #Species i change;
                f.write(" * ")
            f.write(str(rates[j])) #Reaction Rate
                
            ##################
            #Part 1.2: Modifying the Reaction Rate if Needed for Each Type
            ##################
            if isLipidReaction and not (isLipidSpecies or isOnEmptyLipid or isOnTFLipid): #V_s
                f.write("/")
                f.write(parsed_reactions[j].reactionModifiers["bindingSites"])
                 
            if isFunctionReaction:
                f.write(" * ")
                f.write(parsed_reactions[j].reactionModifiers["function"])
            
            ##################
            #Part 1.3: Calculate the forward rate from reactant concentration
            ##################
            for k in range(Ns): #Species k
                if (not math.isnan(s.reactants[k][j])) and (int(s.reactants[k][j]) <= 0):
                    f.write(" * ")
                    f.write(transform_string(species[k]))
                    if (abs(int(s.reactants[k][j])) != 1) and (s.reactants[k][j] != 0):
                        f.write("^")
                        f.write(str(abs(int(s.reactants[k][j]))))
        ##################
        # Case 2: Reaction j INCREASES the amount of species i;
        #################
        elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) > 0):
            ##################
            #Part 2.1: Sign, Rate Constant & Stochiometric Change
            ##################
            f.write("  +  ")
            if abs(int(s.stoich[i][j]))>1:
                f.write(str(abs(int(s.stoich[i][j])))) #Species i change;
                f.write(" * ")
            f.write(str(rates[j])) #Biochemical Reaction Rate
                
            ##################
            #Part 2.2: Modifying the Reaction Rate if Needed for Each Type
            ##################
            #Writing the Reaction Rate isLipidSpecies (true/false); isLipidReaction (true/false)
            if isLipidSpecies: #This must be a koff
                f.write("*")
                f.write(parsed_reactions[j].reactionModifiers["bindingSites"])
                
            if isLipidReaction and ( not isLipidSpecies ) and (isOnEmptyLipid or isOnTFLipid): #V_s
                f.write("/")
                f.write(parsed_reactions[j].reactionModifiers["bindingSites"])
                
            if isFlowReaction:
                f.write(" * ")
                f.write(parsed_reactions[j].reactionModifiers["upstream"])
            
            if isFunctionReaction:
                f.write(" * ")
                f.write(parsed_reactions[j].reactionModifiers["function"])

            if isPlateletSite and speciesName in parsed_reactions[j].reactionModifiers:
                #print(f"We are processing a species {speciesName} for a reaction {parsed_reactions[j]}")
                value = parsed_reactions[j].reactionModifiers[speciesName]
                if value:  # make sure it's not empty/None
                    f.write(" * ")
                    f.write(value)

            ##################
            #Part 2.3: Calculate the forward rate from reactant concentration
            ##################
            for k in range(Ns):
                ##Double check that this makes sense; perhaps should be reactants.
                if (not math.isnan(s.reactants[k][j])) and (int(s.reactants[k][j]) <= 0):
                    f.write(" * ")
                    f.write(transform_string(species[k]))
                    if (abs(int(s.reactants[k][j])) != 1) and (s.reactants[k][j] != 0):
                        f.write("^")
                        f.write(str(abs(int(s.reactants[k][j]))))
            
        #################
        # Case 3: Reaction j DEPENDS on Species i; but doesn't change it's concentration
        #################
        elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) == 0):
            f.write("  +  ")
            f.write(" 0 ")

    #End Looping over all reactions for species [i];
    
    ############################################
    # Case 4: Platelet Sites & Stores Creation
    #############################################
    if ( isPlateletSite or isPlateletStore ):
        #print(f"We are here with special species: {specialSpecies[speciesName]}")
        infoSite = specialSpecies[speciesName]
        for j in range(Ns):
            if species[j] in specialSpecies:
                info = specialSpecies[species[j]]
                if info['type'] == "PROCOAG_PLATELET":
                    f.write(f" + {infoSite['modifier']} * d{transform_string(species[j])} ")

    ############################################
    # Case 5: Adding in the Diultion
    #############################################
    if DILUTION:
        f.write(" - ")
        f.write(transform_string(species[i]))
        f.write(" * ")
        f.write(dil_string)

    if verbose:
        print(f"\t\tSpecies {species[i]} ODE Complete")

    f.write(";\n\n")

def create_matlab_multipleFileOutput(input_file: str, outputPrefix: str,s: Stoich, parsed_reactions: list, species: list, rates: list, uniqueRates: list, uniqueNBS: list, uniqueFlow: list, uniqueftnArgs: list, parsed_ICs: list, parsed_params: list, specialSpecies: list, function_dicts: list, dilution = False, dilution_list = [], verbose: bool = False):
    
    if verbose:
        print('-' * 50)
        print(f"Step 3: Creating Matlab Output Files")
    
    Ns = len(s.species)
    Nr = len(s.rates)
    
    NumReactions = len(parsed_reactions)
    if( not NumReactions == Nr):
        print(f"Error in Creating Matlab File: Total Rates do not Match Reactions.\n")
        exit(-1);
        
    # Dictionaries for specified IC's and params
    ic_dict = {ic.name: ic.value for ic in parsed_ICs}                # Convert IC to dict
    param_dict = {param.name: param.value for param in parsed_params} #Convert Params to dict
        
    matlabFilePrefix = outputPrefix + "Matlab"
    ICFilePrefix     = outputPrefix + "IC"
    ParamPrefix      = outputPrefix + "Params"
    RenamePrefix     = outputPrefix + "Rename"
    
    fIC     = open(ICFilePrefix + ".m", 'w')
    fParam  = open(ParamPrefix + ".m", 'w')
    fRename = open(RenamePrefix + ".m", 'w')
    f       = open(matlabFilePrefix + ".m", 'w')

    f.write("function [time,y] = ")
    f.write(matlabFilePrefix)
    f.write("(t_start,t_final)\n")
    f.write("% Solves a system of ODEs from t=t_start to t=t_final \n")
    f.write("% If no start time is given, then t_start = 0 \n")
    f.write("% If no start or final time is given, then t_start = 0, t_final = 30*60 \n")
    f.write("%\n")
    f.write("%\n")
    f.write("% This file was created by issuing command: \n")
    f.write("%     python createMatlabFile.py ")
    f.write(input_file)
    f.write("\n")
    f.write("%\n")
    f.write("\n")
    f.write("if nargin == 0\n")
    f.write("     t_start = 0;     % Default start time is 0\n")
    f.write("     t_final = 30*60; % Default final time is 30*60\n")
    f.write("elseif nargin~=2\n")
    f.write("   disp('Need to Specify t_start, t_end')\n")
    f.write("   return\n")
    f.write("end\n\n\n")

    fParam.write("% Kinetic Parameters \n")
    for i in range(len(uniqueRates)):
        rate_name = uniqueRates[i]
        value = param_dict.get(rate_name, 1)  # Default to 1 if not explicitly provided
        fParam.write(f"{rate_name} = {value}; \n")

    fParam.write("\np = [ ")
    fParam.write(uniqueRates[0])
    for i in range(1, len(uniqueRates)):
        fParam.write(", ")
        fParam.write(uniqueRates[i])
    fParam.write(" ];\n\n\n")
    
    if len(uniqueNBS) >=1:
        fParam.write("% Binding Site Parameters \n")
        for i in range(len(uniqueNBS)):
            species_name = uniqueNBS[i]
            value = param_dict.get(species_name, 1)  # Default to 1 if not explicitly provided
            fParam.write(f"{species_name} = {value}; \n")
    
        fParam.write("\n")
        fParam.write("nbs = [ ")
        fParam.write(uniqueNBS[0])
        for i in range(1, len(uniqueNBS)):
            fParam.write(", ")
            fParam.write(uniqueNBS[i])
        fParam.write(" ];\n\n\n")

    if len(uniqueftnArgs) >=1:
        fParam.write("% Function Arguments \n")
        for i in range(len(uniqueftnArgs)):
            species_name = uniqueftnArgs[i]
            value = param_dict.get(species_name, 1)  # Default to 1 if not explicitly provided
            fParam.write(f"{species_name} = {value}; \n")
        
        fParam.write("\n")
        fParam.write("otherArgs = [ ")
        fParam.write(uniqueftnArgs[0])
        for i in range(1, len(uniqueftnArgs)):
            fParam.write(", ")
            fParam.write(uniqueftnArgs[i])
        fParam.write(" ];\n\n\n")

    f.write("% Set the Kinetic Parameters\n")
    f.write(f"{ParamPrefix}\n\n")

    fIC.write("% Initial Conditions \n")
    for i in range(Ns):
        species_name = transform_string(species[i]) + "_IC"
        value = ic_dict.get(species_name, 0)  # Default to 0 if not explicitly provided
        fIC.write(f"{species_name} = {value}; \n")

    ##Line that's too long;
    fIC.write("\ninit_cond = [ ")
    fIC.write(transform_string(species[0]))
    fIC.write("_IC")
    for i in range(1,Ns):
        fIC.write(", ")
        fIC.write(transform_string(species[i]))
        fIC.write("_IC")
    fIC.write(" ];\n\n\n")
    
    if len(uniqueFlow) >= 1:
        fIC.write("% Flow Rate Parameters \n")
        for name in uniqueFlow:
            value = param_dict.get(name, 1)     # Get value if exists, otherwise default to 1
            fIC.write(f"{name} = {value}; \n")  # Use the actual value if available

        fIC.write("\n")
        
        fIC.write("flowUp = [ ")
        fIC.write(uniqueFlow[0])
        for i in range(1, len(uniqueFlow)):
            fIC.write(", ")
            fIC.write(uniqueFlow[i])
        fIC.write(" ];\n\n\n")
    
    f.write("% Set the Initial Conditions\n")
    f.write(f"{ICFilePrefix}\n\n")

    ######################################
    # Set Function Mathematical Function #
    ######################################
    args = ['p']
    if len(uniqueNBS) >= 1:
        args.append('nbs')
    if len(uniqueFlow) >= 1:
        args.append('flowUp')
    if len(uniqueftnArgs) >=1:
        args.append('otherArgs')

    f.write("options = odeset('RelTol',1e-12,'AbsTol',1e-23);\n\n\n")

    f.write("%------------------------- Main Solve ----------------------%\n")
    f.write(f"[time,y] = ode15s(@(t,y)RHS(t,y,{','.join(args)}), t_start:1:t_final, init_cond, options);\n")
    f.write("%-----------------------------------------------------------%\n\n\n")

    fRename.write("% Rename solution components\n") ##Modify
    for i in range(Ns):
        fRename.write(transform_string(species[i]))
        fRename.write(" = y(:,")
        fRename.write(str(i+1))
        fRename.write("); \n")

    f.write("% Rename solution components\n") ##Modify
    f.write(f"{RenamePrefix}\n") ##Modify

    f.write("%  \n")
    f.write("% Place plots or other calculations here\n")
    f.write("%   \n")
    f.write("% Example: \n")
    f.write("% plot(time, ")
    f.write(str(transform_string(species[0])))
    f.write(", 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('")
    f.write(str(species[0]))
    f.write("');\n\n\n")

    f.write("end\n\n\n\n")

    f.write("%-----------------------------------------------------%\n")
    f.write("%-------------------- RHS Function -------------------%\n")
    f.write("%-----------------------------------------------------%\n\n")

    f.write(f"function dy = RHS(t,y,{','.join(args)})\n\n")
    
    f.write("dy = zeros(")
    f.write(str(Ns))
    f.write(",1);\n")
    f.write("\n\n")

    f.write("% Rename Variables \n\n") ## Modify

    for i in range(Ns):
        f.write(str(transform_string(species[i])))
        f.write("   = y(")
        f.write(str(i+1))
        f.write("); \n")
    f.write("\n\n")

    f.write("% Rename Kinetic Parameters \n")
    for i in range(len(uniqueRates)):
        f.write(str(uniqueRates[i]))
        f.write(" = p(")
        f.write(str(i+1))
        f.write(");  \n")
    f.write("\n\n")

    if len(uniqueNBS) >=1:
        f.write("% Rename Binding Site \n")
        for i in range(len(uniqueNBS)):
            f.write(str(uniqueNBS[i]))
            f.write(" = nbs(")
            f.write(str(i+1))
            f.write(");  \n")
        f.write("\n\n")

    if len(uniqueftnArgs) >=1:
        f.write("% Rename Function Arguments Site \n")
        for i in range(len(uniqueftnArgs)):
            f.write(str(uniqueftnArgs[i]))
            f.write(" = otherArgs(")
            f.write(str(i+1))
            f.write(");  \n")
        f.write("\n\n")

    if len(uniqueFlow) >=1:
        f.write("% Rename Flow Up \n")
        for i in range(len(uniqueFlow)):
            f.write(str(uniqueFlow[i]))
            f.write(" = flowUp(")
            f.write(str(i+1))
            f.write(");  \n")
        f.write("\n\n")

    f.write("% ODEs from reaction equations \n\n")

    if verbose: print("\tWriting ODEs now....")

    #Makes Sure if we transformed a species name
    if DILUTION:
        first_func = dilution_list[0]
        transformed_args = [transform_string(arg) for arg in first_func['args']]
        dil_string = f"{first_func['name']}({', '.join(transformed_args)})"
    else:
        dil_string = None
        
    # First pass: all species except platelet sites & platelet stores
    for i in range(Ns):
        speciesName = species[i]
        species_info = specialSpecies.get(speciesName, {})
        species_type = species_info.get("type")

        if species_type in ["PLATELET_SITE", "PLATELET_STORE"]:
            continue
        write_species_ode(f, i, species, s, rates, parsed_reactions, specialSpecies, DILUTION, dil_string, verbose)

    # Second pass: platelet sites & platelet stores
    for i in range(Ns):
        speciesName = species[i]
        species_info = specialSpecies.get(speciesName, {})
        species_type = species_info.get("type")

        if species_type not in ["PLATELET_SITE", "PLATELET_STORE"]:
            continue
        write_species_ode(f, i, species, s, rates, parsed_reactions, specialSpecies, DILUTION, dil_string, verbose)
    
    
    f.write(f" dy = [ d{transform_string(species[0])}")
    for i in range(1,Ns):
        f.write(f", d{transform_string(species[i])}")
    f.write(" ]';\n\n\n")
    
    f.write("end")   ###END RHS
    f.write("\n\n")

    ##Output Helper Functions:
    f.write("%Beginning of Helper Functions\n")

    for function_dict in function_dicts:
        # Generate the MATLAB code from the function dictionary
        matlab_code = create_matlab_function(function_dict)
    
        # Write the generated function code to the file
        f.write(matlab_code)
        f.write("\n\n")  # Add some spacing between functions if desired

    f.write("%End of Helper Functions\n")
    
    if verbose:
        print(f"DONE! Successfully created Matlab Files")
        print('-' * 50)
    
def create_matlab_multipleFileOutputOLD(input_file: str, outputPrefix: str,s: Stoich, parsed_reactions: list, species: list, rates: list, uniqueRates: list, uniqueNBS: list, uniqueFlow: list, uniqueftnArgs: list, parsed_ICs: list, parsed_params: list, specialSpecies: list, function_dicts: list, dilution = False, dilution_list = [], verbose: bool = False):
    
    if verbose:
        print('-' * 50)
        print(f"Step 3: Creating Matlab Output Files")
    
    Ns = len(s.species)
    Nr = len(s.rates)
    
    NumReactions = len(parsed_reactions)
    if( not NumReactions == Nr):
        print(f"Error in Creating Matlab File: Total Rates do not Match Reactions.\n")
        exit(-1);
        
    matlabFilePrefix = outputPrefix + "Matlab"
    ICFilePrefix     = outputPrefix + "IC"
    ParamPrefix      = outputPrefix + "Params"
    RenamePrefix     = outputPrefix + "Rename"
    
    fIC     = open(ICFilePrefix + ".m", 'w')
    fParam  = open(ParamPrefix + ".m", 'w')
    fRename = open(RenamePrefix + ".m", 'w')
    f       = open(matlabFilePrefix + ".m", 'w')

    f.write("function [time,y] = ")
    f.write(matlabFilePrefix)
    f.write("(t_start,t_final)\n")
    f.write("% Solves a system of ODEs from t=t_start to t=t_final \n")
    f.write("% If no start time is given, then t_start = 0 \n")
    f.write("% If no start or final time is given, then t_start = 0, t_final = 30*60 \n")
    f.write("%\n")
    f.write("%\n")
    f.write("% This file was created by issuing command: \n")
    f.write("%     python createMatlabFile.py ")
    f.write(input_file)
    f.write("\n")
    f.write("%\n")
    f.write("\n")
    #f.write("\nif nargin == 1\n")
    #f.write("     t_start = 0;  % Default start time is 0 \n")
    f.write("if nargin == 0\n")
    f.write("     t_start = 0;     % Default start time is 0\n")
    f.write("     t_final = 30*60; % Default final time is 30*60\n")
    f.write("elseif nargin~=2\n")
    f.write("   disp('Need to Specify t_start, t_end')\n")
    f.write("   return\n")
    f.write("end\n\n\n")

    fParam.write("% Kinetic Parameters \n")
    for i in range(len(uniqueRates)):
        fParam.write(uniqueRates[i])
        fParam.write(" = 1; \n")                            #<- Should be able to set initial parameter

    fParam.write("\np = [ ")
    fParam.write(uniqueRates[0])
    for i in range(1, len(uniqueRates)):
        fParam.write(", ")
        fParam.write(uniqueRates[i])
    fParam.write(" ];\n\n\n")
    
    if len(uniqueNBS) >=1:
        fParam.write("% Binding Site Parameters \n")
        for i in range(len(uniqueNBS)):
            fParam.write(uniqueNBS[i])
            fParam.write(" = 1; \n")                          #<- Should be able to set initial parameter
        
        fParam.write("\n")
        fParam.write("nbs = [ ")
        fParam.write(uniqueNBS[0])
        for i in range(1, len(uniqueNBS)):
            fParam.write(", ")
            fParam.write(uniqueNBS[i])
        fParam.write(" ];\n\n\n")

    if len(uniqueftnArgs) >=1:
        fParam.write("% Function Arguments \n")
        for i in range(len(uniqueftnArgs)):
            fParam.write(uniqueftnArgs[i])
            fParam.write(" = 1; \n")                          #<- Should be able to set initial parameter
        
        fParam.write("\n")
        fParam.write("otherArgs = [ ")
        fParam.write(uniqueftnArgs[0])
        for i in range(1, len(uniqueftnArgs)):
            fParam.write(", ")
            fParam.write(uniqueftnArgs[i])
        fParam.write(" ];\n\n\n")

    f.write("% Set the Kinetic Parameters\n")
    f.write(f"{ParamPrefix}\n\n")

    ########################################################################################
    # Set Initial Conditions and Upstream Flow Values: Default value of 0 if not specified #
    ########################################################################################
    
    ic_dict = {ic.name: ic.value for ic in parsed_ICs}                # Convert IC to dict
    param_dict = {param.name: param.value for param in parsed_params} #Convert Params to dict

    fIC.write("% Initial Conditions \n")
    for i in range(Ns):
        species_name = transform_string(species[i]) + "_IC"
        value = ic_dict.get(species_name, 0)  # Default to 0 if not explicitly provided
        fIC.write(f"{species_name} = {value}; \n")

    ##Line that's too long;
    fIC.write("\ninit_cond = [ ")
    fIC.write(transform_string(species[0]))
    fIC.write("_IC")
    for i in range(1,Ns):
        fIC.write(", ")
        fIC.write(transform_string(species[i]))
        fIC.write("_IC")
    fIC.write(" ];\n\n\n")
    
    if len(uniqueFlow) >= 1:
        fIC.write("% Flow Rate Parameters \n")
        for name in uniqueFlow:
            value = param_dict.get(name, 1)     # Get value if exists, otherwise default to 1
            fIC.write(f"{name} = {value}; \n")  # Use the actual value if available

        fIC.write("\n")
        
        fIC.write("flowUp = [ ")
        fIC.write(uniqueFlow[0])
        for i in range(1, len(uniqueFlow)):
            fIC.write(", ")
            fIC.write(uniqueFlow[i])
        fIC.write(" ];\n\n\n")
    
    f.write("% Set the Initial Conditions\n")
    f.write(f"{ICFilePrefix}\n\n")

    ######################################
    # Set Function Mathematical Function #
    ######################################
    args = ['p']
    if len(uniqueNBS) >= 1:
        args.append('nbs')
    if len(uniqueFlow) >= 1:
        args.append('flowUp')
    if len(uniqueftnArgs) >=1:
        args.append('otherArgs')

    f.write("options = odeset('RelTol',1e-12,'AbsTol',1e-23);\n\n\n")

    f.write("%------------------------- Main Solve ----------------------%\n")
    f.write(f"[time,y] = ode15s(@(t,y)RHS(t,y,{','.join(args)}), t_start:1:t_final, init_cond, options);\n")
    f.write("%-----------------------------------------------------------%\n\n\n")

    fRename.write("% Rename solution components\n") ##Modify
    for i in range(Ns):
        fRename.write(transform_string(species[i]))
        fRename.write(" = y(:,")
        fRename.write(str(i+1))
        fRename.write("); \n")

    f.write("% Rename solution components\n") ##Modify
    f.write(f"{RenamePrefix}\n") ##Modify

    f.write("%  \n")
    f.write("% Place plots or other calculations here\n")
    f.write("%   \n")
    f.write("% Example: \n")
    f.write("% plot(time, ")
    f.write(str(transform_string(species[0])))
    f.write(", 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('")
    f.write(str(species[0]))
    f.write("');\n\n\n")

    f.write("end\n\n\n\n")

    f.write("%-----------------------------------------------------%\n")
    f.write("%-------------------- RHS Function -------------------%\n")
    f.write("%-----------------------------------------------------%\n\n")

    f.write(f"function dy = RHS(t,y,{','.join(args)})\n\n")
    
    f.write("dy = zeros(")
    f.write(str(Ns))
    f.write(",1);\n")
    f.write("\n\n")

    f.write("% Rename Variables \n\n") ## Modify

    for i in range(Ns):
        f.write(str(transform_string(species[i])))
        f.write("   = y(")
        f.write(str(i+1))
        f.write("); \n")
    f.write("\n\n")

    f.write("% Rename Kinetic Parameters \n")
    for i in range(len(uniqueRates)):
        f.write(str(uniqueRates[i]))
        f.write(" = p(")
        f.write(str(i+1))
        f.write(");  \n")
    f.write("\n\n")

    if len(uniqueNBS) >=1:
        f.write("% Rename Binding Site \n")
        for i in range(len(uniqueNBS)):
            f.write(str(uniqueNBS[i]))
            f.write(" = nbs(")
            f.write(str(i+1))
            f.write(");  \n")
        f.write("\n\n")

    if len(uniqueftnArgs) >=1:
        f.write("% Rename Function Arguments Site \n")
        for i in range(len(uniqueftnArgs)):
            f.write(str(uniqueftnArgs[i]))
            f.write(" = otherArgs(")
            f.write(str(i+1))
            f.write(");  \n")
        f.write("\n\n")

    if len(uniqueFlow) >=1:
        f.write("% Rename Flow Up \n")
        for i in range(len(uniqueFlow)):
            f.write(str(uniqueFlow[i]))
            f.write(" = flowUp(")
            f.write(str(i+1))
            f.write(");  \n")
        f.write("\n\n")

    f.write("% ODEs from reaction equations \n\n")

    if verbose: print("\tWriting ODEs now....")

    ####################
    # Process Diultion #
    ####################
    
    if DILUTION:
        first_func = dilution_list[0]  # Access the first element
        dil_string = f"{first_func['name']}({', '.join(first_func['args'])})"
    
    for i in range(Ns): #For each Species
        speciesName = species[i]
        isLipidSpecies = specialSpecies.get(speciesName) == "LIPID"
        isPlateletSite = specialSpecies.get(speciesName) == "PLATELET_SITE"
        isOnEmptyLipid = "_s" in speciesName and (not "_st" in speciesName) and (not "_sn" in speciesName)
        isOnTFLipid    = "_st" in speciesName
        isOnSilica     = "_sn" in speciesName
        f.write("% ")
        f.write(str(transform_string(species[i])))
        f.write("\n dy(")
        f.write(str(i+1))
        f.write(")  =")
        for j in range(Nr): #For Each Reaction Rate....
            isLipidReaction    = (parsed_reactions[j].reactionType=="LIPID")
            isFlowReaction     = (parsed_reactions[j].reactionType=="FLOW")
            isFunctionReaction = (parsed_reactions[j].reactionType=="FUNCTION")
            
            #################
            # Case 1: Reaction j DECREASES the amount of species i;
            #################
            if (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) < 0):
                ##################
                #Part 1: Sign, Rate Constant & Stochiometric Change
                ##################
                f.write("  -  ")
                if abs(int(s.stoich[i][j]))>1:
                    f.write(str(abs(int(s.stoich[i][j])))) #Species i change;
                    f.write(" * ")
                f.write(str(rates[j])) #Reaction Rate
                
                ##################
                #Part 2: Modifying the Reaction Rate if Needed for Each Type
                #################
                if isLipidReaction and not (isLipidSpecies or isOnEmptyLipid or isOnTFLipid): #V_s
                    f.write("/")
                    f.write(parsed_reactions[j].reactionModifiers["bindingSites"])
                 
                if isFunctionReaction:
                    f.write(" * ")
                    f.write(parsed_reactions[j].reactionModifiers["function"])
                
                ##################
                #Part 3: Calculate the forward rate from reactant concentration
                ##################
                for k in range(Ns): #Species k
                    if (not math.isnan(s.reactants[k][j])) and (int(s.reactants[k][j]) <= 0):
                        f.write(" * ")
                        f.write(transform_string(species[k]))
                        if (abs(int(s.reactants[k][j])) != 1) and (s.reactants[k][j] != 0):
                            f.write("^")
                            f.write(str(abs(int(s.reactants[k][j]))))

            #################
            # Case 2: Reaction j INCREASES the amount of species i;
            #################
            elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) > 0):
                ##################
                #Part 1: Sign, Rate Constant & Stochiometric Change
                ##################
                f.write("  +  ")
                if abs(int(s.stoich[i][j]))>1:
                    f.write(str(abs(int(s.stoich[i][j])))) #Species i change;
                    f.write(" * ")
                f.write(str(rates[j])) #Biochemical Reaction Rate
                
                ##################
                #Part 2: Modifying the Reaction Rate if Needed for Each Type
                ##################
                #Writing the Reaction Rate isLipidSpecies (true/false); isLipidReaction (true/false)
                if isLipidSpecies: #This must be a koff
                    f.write("*")
                    f.write(parsed_reactions[j].reactionModifiers["bindingSites"])
                    
                if isLipidReaction and ( not isLipidSpecies ) and (isOnEmptyLipid or isOnTFLipid): #V_s
                    f.write("/")
                    f.write(parsed_reactions[j].reactionModifiers["bindingSites"])
                    
                if isFlowReaction:
                    f.write(" * ")
                    f.write(parsed_reactions[j].reactionModifiers["upstream"])

                if isFunctionReaction:
                    f.write(" * ")
                    f.write(parsed_reactions[j].reactionModifiers["function"])

                if isPlateletSite and speciesName in parsed_reactions[j].reactionModifiers:
                    #print(f"We are processing a species {speciesName} for a reaction {parsed_reactions[j]}")
                    value = parsed_reactions[j].reactionModifiers[speciesName]
                    if value:  # make sure it's not empty/None
                        f.write(" * ")
                        f.write(value)

                ##################
                #Part 3: Calculate the forward rate from reactant concentration
                ##################
                for k in range(Ns):
                    ##Double check that this makes sense; perhaps should be reactants.
                    if (not math.isnan(s.reactants[k][j])) and (int(s.reactants[k][j]) <= 0):
                        f.write(" * ")
                        f.write(transform_string(species[k]))
                        if (abs(int(s.reactants[k][j])) != 1) and (s.reactants[k][j] != 0):
                            f.write("^")
                            f.write(str(abs(int(s.reactants[k][j]))))
            
            #################
            # Case 3: Reaction j DEPENDS on Species i; but doesn't change it's concentration
            #################
            elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) == 0):
                f.write("  +  ")
                f.write(" 0 ")
                
        #####################
        # Consider Dilution #
        #####################
        if DILUTION:
            f.write(" - ")
            f.write(transform_string(species[i]))
            f.write(" * ")
            f.write(dil_string)
        
        if verbose: print(f"\t\tSpecies {species[i]} ODE Complete")
        f.write(";\n\n")

    f.write("\n\n\n\n")
    f.write("end")
    f.write("\n\n")

    ##Output Helper Functions:
    f.write("%Beginning of Helper Functions\n")

    for function_dict in function_dicts:
        # Generate the MATLAB code from the function dictionary
        matlab_code = create_matlab_function(function_dict)
    
        # Write the generated function code to the file
        f.write(matlab_code)
        f.write("\n\n")  # Add some spacing between functions if desired

    f.write("%End of Helper Functions\n")
    
    if verbose:
        print(f"DONE! Successfully created Matlab Files")
        print('-' * 50)

#To include to text-wrap in Matlab file when output. Not being used yet.
def add_line_continuations(code: str, max_line_length: int = 80) -> str:
    lines = code.split('\n')
    new_lines = []

    for line in lines:
        while len(line) > max_line_length:
            split_index = line.rfind(' ', 0, max_line_length)
            if split_index == -1:
                split_index = max_line_length
            new_lines.append(line[:split_index] + ' ...')
            line = line[split_index:].strip()
        new_lines.append(line)
    
    return '\n'.join(new_lines)

def create_matlab_function(function_dict):
    function_name = function_dict['name']
    args = ', '.join(function_dict['args'])
    body = function_dict['body']

    # Generate the MATLAB function
    matlab_function = f"function output = {function_name}({args})\n"
    matlab_function += f"    % Function: {function_name}\n"
    matlab_function += f"    % Arguments: {', '.join(function_dict['args'])}\n"
    matlab_function += f"    % Body: {body}\n"
    matlab_function += f"    output = {body};\n"
    matlab_function += "end\n"
    
    return matlab_function

#################################################################################
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* Main Code *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* #
#################################################################################

verbose                           = True
verboseReaction                   = False
determineLipidConservedQuantities = False
previewVal                        = 10

######################################
# Step 0: Preprocessing              #
#   - Read in Txt File               #
#   - Remove white space             #
#   - Store the each type of line    #
######################################

parsed_data = parseInputFile(verbose)

biochemicalReactions = parsed_data["biochemicalReactions"]
initialConditions    = parsed_data["initialConditions"]
speciesDefined       = parsed_data["specialSpecies"] #Lipid (Active/Inactive) Platelet (Procoag), Sites & Stores
functionsDefined     = parsed_data["functionCondition"]
parameters           = parsed_data["parameterCondition"]
flowList             = parsed_data["flowList"]
dilutionCondition    = parsed_data["dilutionCondition"]
DILUTION             = parsed_data["DILUTION"]

#########################################################
# Step 1: Parse and Check for Consistency in all Inputs #
#########################################################

parsed_ICs                   = parseInitialConditions(initialConditions,verbose)
specialSpecies, errorFlag    = parseSpecies(speciesDefined,verbose)
parsedFunctions              = parseFunction(functionsDefined,verbose)
parsedFlowReactions          = parseFlowReactions(flowList,verbose)
parsedParameters             = parseParameters(parameters, verbose)
parsedDilution               = parseDilution(dilutionCondition, verbose)

#######################################################
# Step 2: Parse and Error Check Biochemical Reactions #
#######################################################

#Add the flow reactions to the full set of reactions:
biochemicalReactions.extend(parsedFlowReactions)

parsed_reactions, extraParsedParameters = parseReactions(biochemicalReactions, getSubset(specialSpecies,"PLATELET_SITE"), verbose, verboseReaction)

#Add any in-line parsed parameters to the existing list of parameters
add_extra_parameters(parsedParameters, extraParsedParameters, verbose)

# (2b) Check for & Remove Duplicates in reaction list
unique_reactions = set()
duplicates = [reaction for reaction in parsed_reactions if reaction in unique_reactions or unique_reactions.add(reaction)]

if verbose:
    print("-" * 50)
    print(f"Begin: Consistency Checking Biochemical Reactions:")
    print("\t(1) Sanity Check: Are there duplicate Reactions?")

    if duplicates:
        print("\t\t---> Suspicious Duplicate Reactions Found. Removing them:")
        print("\n".join(r.equation for r in duplicates))
    else:
        print("\t\t\t--->Sanity Check PASSED: No Duplicate Reactions Found.")

# (2c) Determine the Number of reaction types, List of Unique Species, Unique Reaction Rates
species      = []
rates        = []
nbs          = []
flow         = []
flowRate     = []
functionArgs = []

# (2d) Dictionary to track the flow reactions for each species
flow_species_count  = {}
badFlowReactions    = []

# Dictionary counting occurrences of each reaction type
reaction_type_count = {"LIPID": 0, "PLATELET": 0, "PLATELET_ACTIVATION":0, "MASS_ACTION": 0, "FLOW": 0, "FUNCTION": 0}

# Iterate through each reaction (now unidirectional after parsing)
for reaction in parsed_reactions:
    rates.append(reaction.rateName)
    
    # Increment the count for the reaction type
    if reaction.reactionType in reaction_type_count:
        reaction_type_count[reaction.reactionType] += 1
     
    if(reaction.reactionType=="LIPID" or reaction.reactionType=="PLATELET"):
        nbs.append(reaction.reactionModifiers["bindingSites"])
        
    if(reaction.reactionType == "MASS_ACTION"):
        # Check for plateletSites in reactionModifiers
        for site in getSubset(specialSpecies,"PLATELET_SITE"):
            if site in reaction.reactionModifiers:
                val = reaction.reactionModifiers[site]
                # Check if val is NOT an integer (could be symbolic expression or parameter name)
                try:
                    int(val)
                except ValueError:
                    nbs.append(val)  # Add non-integer-valued platelet site modifier
    
    if(reaction.reactionType == "FUNCTION"):
        functionArgs.extend(reaction.reactionModifiers["args"])

    if(reaction.reactionType=="FLOW"):
        if "upstream" in reaction.reactionModifiers:
            flow.append(reaction.reactionModifiers["upstream"])
        flowRate.append(reaction.rateName)
        if not (len(reaction.reactants) == 1 and len(reaction.products) == 0) and \
           not (len(reaction.reactants) == 0 and len(reaction.products) == 1):
            badFlowReactions.append(reaction)
    
    # Append Reactants and Products for all reactions
    if reaction.reactants:
        species.extend(reaction.reactants)
        if reaction.reactionType == "FLOW":
            for reactant in reaction.reactants:
                if reactant in flow_species_count:
                    flow_species_count[reactant] += 1
                else:
                    flow_species_count[reactant] = 1

    if reaction.products:
        species.extend(reaction.products)
        if reaction.reactionType == "FLOW":
            for product in reaction.products:
                if product in flow_species_count:
                    flow_species_count[product] += 1
                else:
                    flow_species_count[product] = 1

##Now we should also add for all the plateletSites and plateletStores the information.
for name, info in specialSpecies.items():
    if info["type"] in ["PLATELET_SITE", "PLATELET_STORE"]:
        species.append(name)
        nbs.append(info["modifier"])

unique_species = unique_entries_in_order(species)
unique_rates   = unique_entries_in_order(rates)
unique_nbs     = unique_entries_in_order(nbs)
unique_flow    = unique_entries_in_order(flow)
unique_ftnArgs = unique_entries_in_order(functionArgs)

#############################################
# Error Check: Am I using all unique names? #
#############################################

if verbose:
    print(f"\t(2) Are all species and parameter names unique?")
    
all_variable_names = unique_species + unique_rates + unique_nbs + unique_flow
name_counts        = Counter(all_variable_names)
duplicates         = [name for name, count in name_counts.items() if count > 1]

if verbose:
    if duplicates:
        print(f"\t\t - Error: Duplicate variable names found across categories:", duplicates)
        print(f"\t\t Exiting due to duplicate variable error.")
        exit(-1)
    else:
        print(f"\t\t - All variable names are unique across categories.")

############################################
# Do we have Extra Parameters in Functions #
############################################

##Add in ANY extra arguments from any function in the list
for func in parsedFunctions:
    # Iterate through each argument in the function's args list
    for arg in func['args']:
        # Add to unique_ftnArgs if it's (1) NOT a dummy and (2) NOT already in unique_ftnArgs;
        if arg not in func['dummy_args'] and arg not in unique_ftnArgs:
        #if arg not in unique_ftnArgs:
            unique_ftnArgs.append(arg)
 
#Update: Retain only Unique_FtnArgs that are NOT another type of variable.
if verbose:
    print(f"\t(3) Do we need to add any extra parameters?")
    #print(f"\t\t - All Function Arguments: {unique_ftnArgs}")

unique_ftnArgs = [arg for arg in unique_ftnArgs if arg not in all_variable_names]

if verbose:
    print(f"\t\t - Function Arguments that are New Parameters: {unique_ftnArgs}")


# Error Check on Flow Reactions:
if reaction_type_count["FLOW"] > 0:
    if verbose: print("\t(4) Error Check: If FLOW reactions are included do they make sense?")

    #(1) Check that all flow reactions have the same term
    if verbose: print(f"\t\t(Q): Is there only 1 rate for Flow reactions?")
    if len(set(flowRate)) > 1:
        print("\t\tError: Flow list contains multiple unique values.")
        exit(-1)
    if verbose: print(f"\t\t\tYES.")

    #(2) Check that each species involved in flow reactions has exactly 2 reactions
    if verbose: print(f"\t\t(Q): Do all flow species have 2 rates (in and out)?")
    for species_name, count in flow_species_count.items():
        if count != 2:
            print(f"\t\tError: Species '{species_name}' is involved in {count} flow reactions, but exactly 2 are expected.")
            exit(-1)
    if verbose: print(f"\t\t\tYES.")
        
    #(3) Check that all flow reactions are of the form; "A -> ", or " -> A"
    if verbose: print(f"\t\t(Q): Are all flow reactions of the form A->0 or 0->A?")
    if len(badFlowReactions)>0:
        print(f"\t\t\tExpected (1 reactant, 0 products) or (0 reactants, 1 product).")
        for badRxn in badFlowReactions:
            print(f"\t\t\tError: FLOW reaction format error: {badRxn.equation}")
        exit(-1)
    else:
        if verbose: print(f"\t\t\tYES.")
    
    #########################################
    # What other Error Checks to Include??  #
    #########################################
    #(4) Are all flow reaction modifiers valid (non-negative real numbers or Species_IN)
    ##if verbose: print(f"\t\t(Q): Are all flow reaction modifiers valid?")
        
    if verbose: print("\t\t--->Error Check PASSED")

print("=" * 50)
print(f"{'Reaction Type Counts':^50}")
print("=" * 50)
print(f"{'Reaction Type':<25} | {'Count':>10}")
print("-" * 50)

for reaction_type, count in reaction_type_count.items():
    if count>0:
        print(f"{reaction_type:<25} | {count:>10} reactions")

if DILUTION:
    print(f"{'DILUTION':<25} | {'True':>10}")
print("=" * 50)

####################
# Things to Modify #
####################
# (1) Error Check to Add: We allow upstream in values to be set to other parameters
#     (Example: Cup = C_IN) symbolically. We do NOT check that this exists.
#     Add check
# (2) Pull all the in-line parameters from the equations (if the values do NOT equal to the numerical value)
#(6) Make it possible to be able to in-line declare a biochemical rate:
#    - Check that any rate for a parameter is the same (error check)
#(3) All specified parameters should occur somewhere! (Rxn rates, in a function).

##########################################
# Step 3: Create Stoichiometric Matrices #
#   - Columns: Species                   #
#   - Rows:    Reactions                 #
##########################################

# Extract the PREFIX from the input file name
prefix, _  = os.path.splitext(sys.argv[1])
stoich = createBiochemicalMatrices(unique_species,parsed_reactions,flowRate,prefix,verbose)

print("=" * 50)
print(f"System Details:")

if reaction_type_count["FLOW"] > 0:
    print(f"(*) Flow Details:")
    print(f"\t- Flow Reactions occur at rate: {flowRate[0]}.")
    flow_species_list = sorted(flow_species_count.keys())
    print(f"\t- Species that are involved in Flow: {', '.join(flow_species_list)}")
    print(f"\t- Upstream Flow Values: {', '.join(unique_flow)}")

print(f"(*) Species Details:")
print(f"\t- All Species: {', '.join(unique_species)}")
lipidSpeciesList    = getSubset(specialSpecies,"LIPID").keys()
plateletSpeciesList = getSubset(specialSpecies,"PLATELET").keys()
procoagPlateletSpeciesList = getSubset(specialSpecies,"PROCOAG_PLATELET").keys()
plateletSitesList   = getSubset(specialSpecies,"PLATELET_SITE").keys()
plateletStoreList   = getSubset(specialSpecies,"PLATELET_STORE").keys()
if len(lipidSpeciesList)>0:
    print(f"\t- Lipid Species: {', '.join(lipidSpeciesList)}")
if len(plateletSpeciesList)>0:
    print(f"\t- Platelet Species: {', '.join(plateletSpeciesList)}")
if len(procoagPlateletSpeciesList)>0:
    print(f"\t- Procoag-Platelet Species: {', '.join(procoagPlateletSpeciesList)}")
if len(plateletSitesList)>0:
    print(f"\t- Platelet Sites: {', '.join(plateletSitesList)}")
if len(plateletStoreList)>0:
    print(f"\t- Platelet Stores: {', '.join(plateletStoreList)}")

print(f"(*) Parameter Details:")
if reaction_type_count["LIPID"] > 0 or reaction_type_count["PLATELET"] > 0:
    print(f"\t- Lipid & Platelet Binding Parameters: {', '.join(unique_nbs)}")

if len(parsedParameters)>0:
    param_list = [f"{p.name}={p.value}" for p in parsedParameters]
    print("\t- User Defined Parameters: " + ", ".join(param_list))

if len(parsedFunctions)>0:
    print(f"(*) Function Details:")
    print(f"\t- Parameters Added by Functions: {', '.join(unique_ftnArgs)}")
    func_strings = [f"{func['name']}({', '.join(func['args'])})" for func in parsedFunctions]
    print(f"\t- Functions Defined: {func_strings}")
    
if DILUTION:
    first_func = parsedDilution[0]  # Access the first element
    dil_string = f"{first_func['name']}({', '.join(first_func['args'])})"
    print(f"\t- Dilution Rate: {dil_string}")

print("=" * 50)

if verbose:
    print(f"(Q): Do all Specified Species Occur in the Biochemical Reactions?")

    # Check Lipid Species
    for lipid in lipidSpeciesList:
        if lipid not in unique_species:
            print(f"\tWarning: Lipid species '{lipid}' not found in biochemical reactions.")

    # Check Platelet Species
    for platelet in plateletSpeciesList:
        if platelet not in unique_species:
            print(f"\tWarning: Platelet species '{platelet}' not found in biochemical reactions.")

    # Check Platelet Sites
    for site in plateletSitesList:
        if site not in unique_species:
            print(f"\tWarning: Platelet site '{site}' not found in biochemical reactions.")

##What does Validate Parameters do???
validate_parameters(parsedParameters,unique_species,verbose)

#####################################
# Step 4: Create Matlab File Output #
#####################################

create_matlab_multipleFileOutput(sys.argv[1], prefix, stoich, parsed_reactions, unique_species, rates, unique_rates, unique_nbs, unique_flow, unique_ftnArgs, parsed_ICs, parsedParameters, specialSpecies, parsedFunctions, DILUTION, parsedDilution, verbose)

print(f"FINISHED SUCCESSFULLY!")

#######################################################
## Step 5: Determine Conserved Quantities (If Needed) #
#######################################################

if determineLipidConservedQuantities:
    s_count = 0
    st_count = 0
    sn_count = 0;

    sLipidSites  = "Ls_Sites = ";
    stLipidSites = "Lst_Sites = ";
    snLipidSites = "Sn_Sites = ";

    #Goal:
    # For each species bound to lipid; determine the number of bindings sites.
    # The binding site sizes we have are given in this list:
    # unique_nbs

    for spec in unique_species:
        #Determine which strings in this list we contain a match to;
        #matches = [nbs for nbs in unique_nbs if nbs in spec]
        #print(f"Species {spec} and matches = {matches}")

        # Count occurrences of "_st" first
        local_st_count = spec.count("_st")
    
        # Replace "_st" with a placeholder to avoid counting it as "_s"
        temp_spec = spec.replace("_st", "")
    
        # Count occurrences of "_st" first
        local_sn_count = temp_spec.count("_sn")
    
        # Replace "_sn" with a placeholder to avoid counting it as "_s"
        temp_spec = temp_spec.replace("_sn", "")
    
        # Count occurrences of "_s" in the modified string
        local_s_count = temp_spec.count("_s")
    
        if local_st_count > 0:
            if local_st_count > 1:
                stLipidSites += f"+ {local_st_count} * {transform_string(spec)} "
            else:
                stLipidSites += f"+ {transform_string(spec)} "
        if local_sn_count > 0:
            if local_sn_count > 1:
                snLipidSites += f"+ {local_sn_count} * {transform_string(spec)} "
            else:
                snLipidSites += f"+ {transform_string(spec)} "
        if local_s_count > 0:
            if local_s_count > 1:
                sLipidSites += f"+ {local_s_count} * {transform_string(spec)} "
            else:
                sLipidSites += f"+ {transform_string(spec)} "
    
        st_count += local_st_count
        sn_count += local_sn_count
        s_count += local_s_count
        
    
        if verbose: print(f"{spec}\t{local_st_count}\t{local_sn_count}\t{local_s_count}")

        if verbose:
            print(f'The string "_s" occurs {s_count} times (excluding "_st" and "_sn").')
            print(f'The string "_st" occurs {st_count} times.')
            print(f'The string "_sn" occurs {sn_count} times.')

            print(f'The total bound lipid: {stLipidSites}')
            print(f'The total bound ubBoundlipid: {sLipidSites}')
            print(f'The total bound Sn Sites: {snLipidSites}')



