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
    - Species names may contain: [0-9a-zA-Z_:]
    - Species names must start with a letter (no leading numbers).

    - Users can define:
        * Initial Conditions (must be a non-negative real number; DEFAULT = 0)
            Ex:
                IIa_IC  = 5.0;  #mu M
                V_IC    = 0;
        * Parameter Values (real valued or as previous initial condition; DEFAULT = 1)
            Ex:
                IIa_up  = IIa_IC #mu M; Let's see if this works!
                V_up    = V_IC
                PL_up   = 1.0    #Set to be a random value.
        * Special classes of species: LIPID, PLATELET, PLATELET_SITE
            Ex:
                PL      = PLATELET  #Platlet in Solution
                PL_S    = PLATELET  #Platelet in Subendothelium
                PL_V    = PLATELET  #Platelet in Volume
                p2avail = PLATELET_SITE
                p5avail = PLATELET_SITE
                
    ---------------
    (2) Functions:  NOT FULLY SUPPORTED
    ---------------
    - Functions are defined as a name, list of arguments (comma separated) and body:
    - Examples:
        FUNCTION A(x,e2P) = x/(e2P + x) #1 nM = 0.001 mu M
    - Function arguments can be a parameter, species name, etc

    ----------------------------
    (3) Biochemical Equations:
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
        * [?not needed?] PLATELET:    Binding on/off platelet (no competition)

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

    Example Flow Reactions:
          -> K, kflow, K_up, FLOW #Species Flowing In
        K ->  , kflow, FLOW       #Flowing Flowing Out
        
    Supported Features:
    -----------------------------------
      - Support for non-mass action lipid/platlet binding.
      - Outputs lipid/platlet binding sites as a separate parameter vector (nbs)
      - Can handle non-constant coefficients for platelet site activation.
      - Outputs Species and rates output in input order.
      - Supports pure synthesis/degradation/in-out flow (e.g., "-> A", "B ->").
      - Consolidates duplicate kinetic rates.
      - Splits stoichiometric matrix for A + B -> A + C reactions.
      - Removes duplicate reactions (even one side of bidirectional ones).
      - Checks reaction rate dimensions for consistency.
      - Initial conditions can be set in the input file.

    In-Progress Features (Should Check):
    -----------------------------------
      - Parameter values set set by constants.
      - Support for inline kinetic rate values:
        Example:
          A + B -> C, k1=0.1
          A + B <-> C, kon, koff=100

    Feature to Consider Developing:
    -----------------------------------
    - Allow setting rates and initial conditions from text file or external file.
    - Investigate prior use of "=" operator in reactions.
    - Improve MATLAB text wrapping for long lines.
    - Add back support for Python code.
    - Consider Flow reactions of the form: C <-> C, kflow, C_up, FLOW

    Problems to Resolve:
    -----------------------------------
    - Check for valid reaction types (i.e., flow types must all have the same FLOW rate)
        * Really there are 2 flow rates for biochemical species and platelets.
    
    Current Version:
        Suzanne Sindi, 04/08/2025

    """))
    sys.exit("Usage: python3 createCoagModel.py StaticCoag.txt")

########################
# Supporting Functions #
########################

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

def parseReactions(reactions, plateletSites, verbose):

    if verbose:
        print(f"Step 2: Parsing the Biochemical Reacitons")

    if verbose:
        print(f"\tNumber of biochemical reactions: {len(biochemicalReactions)}")
        print(f"\tSome of the Biochemical Reactions:")
        num_reactions_to_print = min(previewVal, len(reactions))
        for i in range(num_reactions_to_print):
            print(f"\t\t" + biochemicalReactions[i])
        print('-' * 50)

    parsed_reactions = []

    for reaction in reactions:
        if verbose: print(f"Processing Reaction: '{reaction}'")
        
        # Parse the equation
        result = parseEquation(reaction,plateletSites,verbose)
        
        # Ensure result is valid before unpacking
        if not result or all(val is None for val in result):
            print(f"\tError: Invalid equation format in reaction: {reaction}")
            continue  # Skip this reaction
            
        (reactants,reactantCoeffs,products,productCoeffs,rates,reaction_type,reaction_modifiers)=result
        if verbose:
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
            else:
                names.append(rate.strip())  # Store the rate if no '='
                values.append(-1)  # Assign -1 as a placeholder for missing values

        reactionCount = len(rates);

        if reactionCount == 1:  # We add only 1 case, easy
            parsed_reactions.append(create_reaction(reactants, reactantCoeffs, products, productCoeffs,names,values,reaction_type, reaction_modifiers) )
            # Create a Reaction object and add to the list
            #reaction_obj = Reaction(equation = formatFwdReaction(reactants, reactantCoeffs, products, productCoeffs),rateName=names[0],rateValue=values[0],reactants=reactants,reactant_coeffs=reactantCoeffs,products=products, product_coeffs=productCoeffs, reactionType=reaction_type, reactionModifiers = reaction_modifiers)
            #parsed_reactions.append(reaction_obj)
        
        elif reactionCount == 2:  # We add 2 objects for bi-directional reactions
            parsed_reactions.append(create_reaction(reactants, reactantCoeffs, products, productCoeffs,names,values,reaction_type, reaction_modifiers,reverse=False) )
            parsed_reactions.append(create_reaction(reactants, reactantCoeffs, products, productCoeffs,names,values,reaction_type, reaction_modifiers,reverse=True) )

            #reaction_fwd = Reaction( equation=formatFwdReaction(reactants, reactantCoeffs, products, productCoeffs), rateName=names[0], rateValue=values[0], reactants=reactants, reactant_coeffs=reactantCoeffs, products=products, product_coeffs=productCoeffs, reactionType=reaction_type,                 reactionModifiers=reaction_modifiers)
            
            #reaction_rev = Reaction( equation=formatFwdReaction(products, productCoeffs, reactants, reactantCoeffs), rateName=names[1], rateValue=values[1], reactants=products, reactant_coeffs=productCoeffs, products=reactants, product_coeffs=reactantCoeffs, reactionType=reaction_type, reactionModifiers=reaction_modifiers)
            #parsed_reactions.append(reaction_fwd)
            #parsed_reactions.append(reaction_rev)
        
        else:
            print(f"\tError: Invalid reaction count in: {reaction}")
            continue
    
    return parsed_reactions

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
        print(f"NEW reactantCoeffs: {reactantCoeffs}")
        print(f"NEW reactants: {reactants}")
        print(f"NEW plateletSiteReactants: {plateletSiteReactants}")
        print(f"NEW plateletSiteReactantCoeffs: {plateletSiteReactantCoeffs}")
        print(f"NEW productCoeffs: {productCoeffs}")
        print(f"NEW products: {products}")
        print(f"NEW plateletSiteProducts: {plateletSiteProducts}")
        print(f"NEW plateletSiteProductCoeffs: {plateletSiteProductCoeffs}")
        print(f"NEW reactionType: {reactionType}")
        print(f"NEW reactionModifiers: {reactionModifiers}")
    
    #(3) Check if we have the right number of rates.
    # Rate is a dictionary of either "RATE" or "FWD_RATE" "REV_RATE"
    reaction_count = len(rates)
    if reaction_count > 1 and reactionType not in {"MASS_ACTION", "LIPID"}:
        print(f"Error in Equation = {equation}: Only MASS_ACTION and LIPID reactionTypes can be bidirectional")
        return None, None, None, None, None, None, None

    if verbose: print(f"NEW **************")

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
    initialConditions = []    # User-specified initial conditions
    lipidSpecies = []         # Lipids must be handled differently
    plateletSpecies = []      # Platelets must be handled differently
    plateletSites = []        # Platelet Sites (not the species)
    parameterCondition = []   # User-specified parameters (not in-line)
    functionCondition = []    # Dilution and platelet activation

    numReactionsReadIn = 0
    numInitialConditionsReadIn = 0
    numLipidSpecies = 0
    numPlateletSpecies = 0
    numPlateletSites = 0
    numParametersReadIn = 0
    numFunctionsDefined = 0

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
                
                if len(fields) > 1 and not line.lstrip().startswith("FUNCTION"):
                    # If a comma and NOT the word function it's a biochemical equation
                    biochemicalReactions.append(line)
                    numReactionsReadIn += 1
                    if verbose: print(f"\t\tBiochemical Equation: {line}")
                else:
                    # Process single-field lines
                    if "= LIPID" in line:
                        lipidSpecies.append(line)
                        numLipidSpecies += 1
                        if verbose: print(f"\t\tLipid Species: {line}")

                    # Check for platelet sites before checking for platelet species
                    elif "= PLATELET_SITE" in line:
                        plateletSites.append(line)
                        numPlateletSites += 1
                        if verbose:
                            print(f"\t\tPlatelet Sites: {line}")

                    elif "= PLATELET" in line:
                        plateletSpecies.append(line)
                        numPlateletSpecies += 1
                        if verbose: print(f"\t\tPlatelet Species: {line}")

                    elif "_IC" in line.split('=')[0].strip():  # Check if it's an initial condition
                        initialConditions.append(line)
                        numInitialConditionsReadIn += 1
                        if verbose: print(f"\t\tInitial Condition: {line}")
                    
                    # Do we start with function? Then it's a function!
                    elif line.lstrip().startswith("FUNCTION"):
                        functionCondition.append(line)
                        numFunctionsDefined += 1
                        if verbose:
                            print(f"\t\tFunctions: {line}")

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
            ("Specified Initial Conditions", numInitialConditionsReadIn),
            ("Specified Lipid Species", numLipidSpecies),
            ("Specified Platelet Species", numPlateletSpecies),
            ("Specified Parameters", numParametersReadIn),
            ("Functions Specified", numFunctionsDefined)
        ]

        # Print table header
        print(f"\n{'=' * 50}")
        print(f"{'Step 0: Processed Input File':^30}")
        print(f"{'=' * 50}")
        print(f"{headers[0]:<30} | {headers[1]:>10}")
        print("-" * 50)

        # Print each row in the data
        for category, count in data:
            if count > 0:  # Only print if count is greater than zero
                print(f"{category:<30} | {count:>10}")

        # Print done message and finish with a line
        print(f"{'=' * 50}\n")

    return {
        "biochemicalReactions": biochemicalReactions,
        "initialConditions": initialConditions,
        "lipidSpecies": lipidSpecies,
        "plateletSpecies": plateletSpecies,
        "plateletSites": plateletSites,
        "parameterCondition": parameterCondition,
        "functionCondition": functionCondition,
        "counts": {
            "numReactionsReadIn": numReactionsReadIn,
            "numInitialConditionsReadIn": numInitialConditionsReadIn,
            "numLipidSpecies": numLipidSpecies,
            "numPlateletSpecies": numPlateletSpecies,
            "numPlateletSites": numPlateletSites,
            "numParametersReadIn": numParametersReadIn,
            "numFunctionsDefined": numFunctionsDefined,
        }
    }


def parseList(lines, mode=None, verbose=False):
    unique_items = set()  # Track unique items for quick lookup
    ordered_list = []  # Preserve insertion order

    step_details = {
        "LIPID": {
            "header": "Step 1(b): Parse Specified Lipid Species",
            "done": "DONE: 1(b): Parse Specified Lipid Species",
            "header_display": "Step 1(b): Specified Lipid Species"
        },
        "PLATELET": {
            "header": "Step 1(c): Parse Specified Platelet Species",
            "done": "DONE: 1(c): Parse Specified Platelet Species",
            "header_display": "Step 1(c): Specified Platelet Species"
        },
        "PLATELET_SITE": {
            "header": "Step 1(d): Parse Specified Platelet Binding Sites",
            "done": "DONE: 1(d): Parse Specified Platelet Binding Sites",
            "header_display": "Step 1(d): Specified Platelet Binding Sites"
        }
        
    }

    if verbose and mode in step_details:
        print("-" * 50)
        print(step_details[mode]["header"])

    for line in lines:
        line = line.strip()

        # Determine parsing behavior based on mode
        if mode in ["LIPID", "PLATELET","PLATELET_SITE"]:
            keyword = f"= {mode}"
            if keyword in line:
                item_name = line.split('=')[0].strip()  # Get the part before `=`
            else:
                return f"⚠️ Error: Expected '{keyword}' but not found in '{line}'."
        else:
            item_name = line  # Treat as a direct list of items

        # Ensure valid name (cannot start with a number)
        if re.match(r"^[^\d].*", item_name):
            if item_name not in unique_items:  # Only add if unique
                unique_items.add(item_name)
                ordered_list.append(item_name)
                if verbose:
                    print(f"\t✔ Added: {item_name}")
        else:
            print(f"\t⚠️ Error: Invalid name '{item_name}' (cannot start with a number).")

    if verbose:
        print(f"\t✅ Total Unique Items: {len(ordered_list)}")
        print(f"\t🔹 Unique List (Order Preserved): {ordered_list}")

    if verbose and mode in step_details:
        print(step_details[mode]["done"])
        print("-" * 50)
    
        print(f"\n{'=' * 50}")
        print(f"{step_details[mode]['header_display']:^30}")
        print(f"{'=' * 50}")
        print(f"{ordered_list}")
        print(f"{'=' * 50}\n")

    return ordered_list
    
def parseInitialConditions(initialConditions, verbose=False):

    if verbose:
        print("-" * 50)
        print(f"Step 1(a): Parse Intial Conditions")
    
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
        print(f"DONE: Step 1(a): Parse Initial Conditions")
        print('-' * 50)

        # Define headers for the initial conditions table
        ic_headers = ["Initial Condition Name", "Value"]
    
        # Create a data list for parsed initial conditions
        ic_data = [(ic.name, ic.value) for ic in parsed_ICs]  # Assuming parsed_ICs is a list of InitialCondition objects

        # Print table header for initial conditions
        print(f"\n{'=' * 50}")
        print(f"{'Step 1(a): Parse Intial Conditions':^30}")
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

    if verbose:
        print('-' * 50)
        print(f"Step 1(e): Parse Specified Functions")

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
        args_list = [arg.strip() for arg in args.split(",") if arg.strip()]
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
        parsed_functions.append({"name": name, "args": args_list, "body": body})

        if verbose:
            print(f"\tValid function parsed: {name}({', '.join(args_list)}) = {body}")

    if verbose:
        # Create a data list for parsed initial conditions
        print(f"\n{'=' * 50}")
        print(f"{'Step 1(e): Parse Function Conditions':^30}")
        print(f"{'=' * 50}")
        for fVal in parsed_functions:  # Assuming parsed_ICs is a list of InitialCondition objects
            # Print table header for initial conditions
            print(f"{fVal}")
        
        print(f"DONE: Step 1(e): Parse Specified Functions")
        print('-' * 50)

    return parsed_functions


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
            if re.match(r"^\d+(\.\d+)?$", value_str):  # Matches integers and decimals
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

        # Print each parsed parameter
        for name, value, value_type in p_data:
            if isinstance(value, float):  # Format floats to 2 decimal places
                print(f"{name:<30} | {value:>10.2f} | {value_type:>10}")
            else:  # Print strings without decimal formatting
                print(f"{name:<30} | {value:>10} | {value_type:>10}")
        
        # Print done message and finish with a line
        print(f"{'=' * 65}\n")
        
    return parsed_params

def validate_parameters(parsedParameters, uniqueSpecies, verbose = False):
    numErrors = 0
    if verbose:
        print(f"(Q): Do all parameters defined as a string occur as species_IC?")
        
    for param in parsedParameters:
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


    if len(bad_rates) > 0:
        print(f"Error: We found some biochemical rates with different dimensions.")
        for rate, column_sum in bad_rates.items():
            #if(not rate == uniqueFlowRate):
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

#def create_matlab_multipleFileOutput(input_file: str, parsed_reactions: list, outputPrefix: str, s: Stoich, species: list, rates: list, uniqueRates: list, uniqueNBS: list, v: bool = False):

def create_matlab_multipleFileOutput(input_file: str, outputPrefix: str,s: Stoich, parsed_reactions: list, species: list, rates: list, uniqueRates: list, uniqueNBS: list, uniqueFlow: list, uniqueftnArgs: list, parsed_ICs: list, parsed_params: list, specialSpecies: list, function_dicts: list, verbose: bool = False):
    
    if verbose:
        print('-' * 50)
        print(f"Step 3: Creating Matlab Output Files")
    
    Ns = len(s.species)
    Nr = len(s.rates)
    
    NumReactions = len(parsed_reactions)
    if( not NumReactions == Nr):
        print(f"Error in Creating Matlab File: Total Rates do not Match Reactions.\n")
        exit(-1);

    #species = [i.name for i in s.species.keys()]
    #rates = [i for i in s.rates.keys()]
    #if v: print('\nOutput File: \n', output_file)
    #if v: print("\n")
        
    matlabFilePrefix = outputPrefix + "Matlab"
    ICFilePrefix     = outputPrefix + "IC"
    ParamPrefix      = outputPrefix + "Params"
    RenamePrefix     = outputPrefix + "Rename"
    
    fIC     = open(ICFilePrefix + ".m", 'w')
    fParam  = open(ParamPrefix + ".m", 'w')
    fRename = open(RenamePrefix + ".m", 'w')
    f       = open(matlabFilePrefix + ".m", 'w')

    #if v: print('\nOutput File: \n', output_file)
    #if v: print("\n")

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
        #fIC.write(transform_string(species[i]))
        #fIC.write("_IC")
        #fIC.write(" = 0; \n")
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

    ################################
    # NEW: Working on the Dilution #
    ################################
    
    DILUTION = False #True
    ####Pre-Define the Platelet Species for the Dilution###
    if DILUTION:
        platelet_indices = [i for i in range(Ns) if species[i] in parsedPlateletSpecies]
        for i in platelet_indices:
            print(f"\tFound Platelet Species {species[i]} at index {i}....")

    #exit(-1)


    for i in range(Ns): #For each Species
        speciesName = species[i]
        isLipidSpecies = specialSpecies.get(speciesName) == "LIPID"
        isPlateletSite = specialSpecies.get(speciesName) == "PLATELET_SITE"
        #isLipidSpecies = "L_noTF" in speciesName or "L_TF" in speciesName
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
                ##################
                #Writing the Reaction Rate isLipidSpecies (true/false); isLipidReaction (true/false)
                #if isLipidReaction: #This must be a koff
                # V_s -> V + L; L is correct
                # (V_s)_dt = -rate
                # (V)_dt   = +rate
                # L_dt     = +rate*nbs
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
    """
    This function parses the dictionary and generates a corresponding MATLAB function.
    
    Parameters:
    - function_dict (dict): A dictionary with 'name', 'args', and 'body' keys.

    Returns:
    - str: A string that represents the MATLAB function code.
    """
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

verbose            = True
determineConserved = False
previewVal         = 10

######################################
# Step 0: Preprocessing              #
#   - Read in Txt File               #
#   - Remove white space             #
#   - Store the each type of line    #
######################################

parsed_data = parseInputFile(verbose)

biochemicalReactions = parsed_data["biochemicalReactions"]
initialConditions    = parsed_data["initialConditions"]
lipidSpecies         = parsed_data["lipidSpecies"]
plateletSpecies      = parsed_data["plateletSpecies"]
plateletSites        = parsed_data["plateletSites"]
functionsDefined     = parsed_data["functionCondition"]
parameters           = parsed_data["parameterCondition"]

#########################################################
# Step 1: Parse and Check for Consistency in all Inputs #
#########################################################

parsed_ICs            = parseInitialConditions(initialConditions,verbose)
parsedLipidSpecies    = parseList(lipidSpecies,"LIPID", verbose)
parsedPlateletSpecies = parseList(plateletSpecies,"PLATELET", verbose)
parsedPlateletSites   = parseList(plateletSites,"PLATELET_SITE", verbose)
parsedFunctions       = parseFunction(functionsDefined,verbose)
parsedParameters      = parseParameters(parameters, verbose)

#Store all special named species
specialSpecies = {} #Dictionary for special species
appendToSpecialSpecies(specialSpecies,parsedLipidSpecies,"LIPID")
appendToSpecialSpecies(specialSpecies,parsedPlateletSpecies,"PLATELET")
appendToSpecialSpecies(specialSpecies,parsedPlateletSites,"PLATELET_SITE")

#######################################################
# Step 2: Parse and Error Check Biochemical Reactions #
#######################################################

# (2a) Do initial parsing of reactions
parsed_reactions = parseReactions(biochemicalReactions, parsedPlateletSites, True)

# (2b) Check for & Remove Duplicates in reaction list
unique_reactions = set()
duplicates = [reaction for reaction in parsed_reactions if reaction in unique_reactions or unique_reactions.add(reaction)]

print("-" * 50)
print(f"Begin: Consistency Checking Biochemical Reactions:")

print("\t(1) Sanity Check: Are there duplicate Reactions?")
if duplicates:
    print("\t\t---> Suspicious Duplicate Reactions Found. Removing them:")
    #print("\n".join(duplicates))
    print("\n".join(r.equation for r in duplicates))
else:
    print("\t--->Sanity Check PASSED: No Duplicate Reactions Found.")

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
        for site in parsedPlateletSites:
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

unique_species = unique_entries_in_order(species)
unique_rates   = unique_entries_in_order(rates)
unique_nbs     = unique_entries_in_order(nbs)
unique_flow    = unique_entries_in_order(flow)
unique_ftnArgs = unique_entries_in_order(functionArgs)

#Error Check: Am I using all unique names?
all_variable_names = unique_species + unique_rates + unique_nbs + unique_flow
name_counts        = Counter(all_variable_names)
duplicates         = [name for name, count in name_counts.items() if count > 1]

if duplicates:
    print(f"Error: Duplicate variable names found across categories:", duplicates)
else:
    print(f"All variable names are unique across categories.")

#Update: Retain only Unique_FtnArgs that do NOT occur elsewhere
print(f"Before Filtering: {unique_ftnArgs}")
unique_ftnArgs = [arg for arg in unique_ftnArgs if arg not in all_variable_names]
print(f"After Filtering: {unique_ftnArgs}")

# Error Check on Flow Reactions:
if reaction_type_count["FLOW"] > 0:
    print("\t(2) Error Check: If FLOW reactions are included do they make sense?")

    #(1) Check that all flow reactions have the same term
    print(f"\t\t(Q): Is there only 1 rate for Flow reactions?")
    if len(set(flowRate)) > 1:
        print("\t\tError: Flow list contains multiple unique values.")
        exit(-1)
    print(f"\t\t\tYES.")

    #(2) Check that each species involved in flow reactions has exactly 2 reactions
    print(f"\t\t(Q): Do all flow species have 2 rates (in and out)?")
    for species_name, count in flow_species_count.items():
        if count != 2:
            print(f"\t\tError: Species '{species_name}' is involved in {count} flow reactions, but exactly 2 are expected.")
            exit(-1)
    print(f"\t\t\tYES.")
        
    #(3) Check that all flow reactions are of the form; "A -> ", or " -> A"
    print(f"\t\t(Q): Are all flow reactions of the form A->0 or 0->A?")
    if len(badFlowReactions)>0:
        print(f"\t\t\tExpected (1 reactant, 0 products) or (0 reactants, 1 product).")
        for badRxn in badFlowReactions:
            print(f"\t\t\tError: FLOW reaction format error: {badRxn.equation}")
        exit(-1)
    else:
        print(f"\t\t\tYES.")
        
    #(4) Are all flow reaction modifiers valid (non-negative real numbers or Species_IN)
    print(f"\t\t(Q): Are all flow reaction modifiers valid?")

    print("\t--->Error Check PASSED")

print("=" * 50)
print(f"{'Reaction Type Counts':^50}")
print("=" * 50)
print(f"{'Reaction Type':<25} | {'Count':>10}")
print("-" * 50)

for reaction_type, count in reaction_type_count.items():
    print(f"{reaction_type:<25} | {count:>10} reactions")

print("=" * 50)

if reaction_type_count["LIPID"] > 0 or reaction_type_count["PLATELET"] > 0:
    print(f"Summary of Lipid/Platelet Reactions:")
    print(f"\tAll Lipid and Platelet Binding Sites: {unique_nbs}")
    print("=" * 50)

if reaction_type_count["FLOW"] > 0:
    print(f"Summary of Flow Reactions:")
    print(f"\tFlow Reactions occur at rate: {flowRate[0]}.")
    flow_species_list = sorted(flow_species_count.keys())
    print(f"\tSpecies that are involved in Flow: {', '.join(flow_species_list)}")
    print(f"\tAll Flow Related Modifiers: {unique_flow}")
    print("=" * 50)

print(f"End: Consistency Checking Biochemical Reactions")
print("-" * 50)

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

########
# DONE #
########
# (1) Funciton arguments can either be: previously defined parameter or NEW one
#   - Added an extra set of parameters for ones NOT defined elsewhere. (Default Value = 1)
# (2) For kflow reactions; fixed the reactions to have an in-out type.
# (3) Error Check: All FLOW rxns should have the same parameter (could easily allow for 2)
# (4) Error Check: All FLOW rxns should have either [1,0] or [0,1] species reactants


##########################################
# Step 3: Create Stoichiometric Matrices #
#   - Columns: Species                   #
#   - Rows:    Reactions                 #
##########################################

# Extract the PREFIX from the input file name
prefix, _  = os.path.splitext(sys.argv[1])

stoich = createBiochemicalMatrices(unique_species,parsed_reactions,flowRate,prefix,verbose)

print(f"Unique Species: {unique_species}")
print(f"Lipid Species: {parsedLipidSpecies}")
print(f"Platelet Species: {parsedPlateletSpecies}")
print(f"Platelet Sites: {parsedPlateletSites}")
print(f"Parameters: {parsedParameters}")  #Should also pull out the values for biochemical rates from rxns.

print(f"(Q): Do all Specified Species Occur in the Biochemical Reactions?")

# Check Lipid Species
for lipid in parsedLipidSpecies:
    if lipid not in unique_species:
        print(f"\tWarning: Lipid species '{lipid}' not found in biochemical reactions.")

# Check Platelet Species
for platelet in parsedPlateletSpecies:
    if platelet not in unique_species:
        print(f"\tWarning: Platelet species '{platelet}' not found in biochemical reactions.")

# Check Platelet Sites
for site in parsedPlateletSites:
    if site not in unique_species:
        print(f"\tWarning: Platelet site '{site}' not found in biochemical reactions.")

validate_parameters(parsedParameters,unique_species,verbose)

#####################################
# Step 4: Create Matlab File Output #
#####################################

print(f"Functions Defined: {functionsDefined}")

create_matlab_multipleFileOutput(sys.argv[1], prefix, stoich, parsed_reactions, unique_species, rates, unique_rates, unique_nbs, unique_flow, unique_ftnArgs, parsed_ICs, parsedParameters, specialSpecies, parsedFunctions, verbose)

print(f"FINISHED SUCCESSFULLY!")

#######################################################
## Step 5: Determine Conserved Quantities (If Needed) #
#######################################################
#
#if determineConserved:
#    s_count = 0
#    st_count = 0
#    sn_count = 0;
#
#    sLipidSites  = "Ls_Sites = ";
#    stLipidSites = "Lst_Sites = ";
#    snLipidSites = "Sn_Sites = ";
#
#    #Goal:
#    # For each species bound to lipid; determine the number of bindings sites.
#    # The binding site sizes we have are given in this list:
#    # unique_nbs
#
#    for spec in unique_species:
#        #Determine which strings in this list we contain a match to;
#        #matches = [nbs for nbs in unique_nbs if nbs in spec]
#        #print(f"Species {spec} and matches = {matches}")
#
#        # Count occurrences of "_st" first
#        local_st_count = spec.count("_st")
#    
#        # Replace "_st" with a placeholder to avoid counting it as "_s"
#        temp_spec = spec.replace("_st", "")
#    
#        # Count occurrences of "_st" first
#        local_sn_count = temp_spec.count("_sn")
#    
#        # Replace "_sn" with a placeholder to avoid counting it as "_s"
#        temp_spec = temp_spec.replace("_sn", "")
#    
#        # Count occurrences of "_s" in the modified string
#        local_s_count = temp_spec.count("_s")
#    
#        if local_st_count > 0:
#            if local_st_count > 1:
#                stLipidSites += f"+ {local_st_count} * {transform_string(spec)} "
#            else:
#                stLipidSites += f"+ {transform_string(spec)} "
#        if local_sn_count > 0:
#            if local_sn_count > 1:
#                snLipidSites += f"+ {local_sn_count} * {transform_string(spec)} "
#            else:
#                snLipidSites += f"+ {transform_string(spec)} "
#        if local_s_count > 0:
#            if local_s_count > 1:
#                sLipidSites += f"+ {local_s_count} * {transform_string(spec)} "
#            else:
#                sLipidSites += f"+ {transform_string(spec)} "
#    
#        st_count += local_st_count
#        sn_count += local_sn_count
#        s_count += local_s_count
#        
#    
#        if verbose: print(f"{spec}\t{local_st_count}\t{local_sn_count}\t{local_s_count}")
#
#    if verbose:
#        print(f'The string "_s" occurs {s_count} times (excluding "_st" and "_sn").')
#        print(f'The string "_st" occurs {st_count} times.')
#        print(f'The string "_sn" occurs {sn_count} times.')
#
#
#        print(f'The total bound lipid: {stLipidSites}')
#        print(f'The total bound ubBoundlipid: {sLipidSites}')
#        print(f'The total bound Sn Sites: {snLipidSites}')
#
#
#
