#!/usr/bin/env python3

import sys, re, io, os, textwrap, math, csv
import numpy as np
from collections import defaultdict


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
    -----------------------------------
    - Assumes  biochemical reactions (one per line) of the form:
            LHS <-> RHS, Rates, TYPE
    - Comments/Whitespace: Anything following "#" is ignored; Whitespace is ignored.
    - Supports forward and reversible reactions.
    - Reaction Operators allowed: '*', '+', '<->', '->'
    - Species names may contain: [0-9a-zA-Z_:]
    - Species names must start with a letter (no leading numbers).

    - 5 Reaction Types:
        * MASS_ACTION: Default (if no type given)
        * FLOW:        species entering or exiting reaction zone.
        * LIPID:       Binding on/off lipid (competition)
        * PLATELET:    Binding on/off platelet (no competition)
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

    Supported Features:
    -----------------------------------
      - Support for non-mass action lipid/platlet binding.
      - Outputs lipid/platlet binding sites as a separate parameter vector (nbs)
      - Can handle only mass action terms.
      - Outputs Species and rates output in input order.
      - Supports pure synthesis/degradation (e.g., "-> A", "B ->").
      - Consolidates duplicate kinetic rates.
      - Splits stoichiometric matrix for A + B -> A + C reactions.
      - Removes duplicate reactions (even one side of bidirectional ones).
      - Checks reaction rate dimensions for consistency.

    In-Progress Features (Should Check):
    -----------------------------------
      - Initial conditions can be set in the input file.
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

    Problems to Resolve:
    -----------------------------------
    - Check for valid reaction types (i.e., flow types must all have the same FLOW rate)
    - Currently doesn't handle the reaction: -> 3*A, kflow, Aup correctly

    Current Version:
        Suzanne Sindi, 03/25/2025

    """))
    sys.exit("Usage: python3 createCoagModel.py StaticCoag.txt")

########################
# Supporting Functions #
########################

#Matlab can not handle ":"'s in variable names!
def transform_string(s):
    #return s.replace(':', '_')
    return s.replace(':', 'b')

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

def parseReactions(reactions, verbose):

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
        result = parseEquation(reaction)
        
        # Ensure result is valid before unpacking
        if not result or all(val is None for val in result):
            print(f"\tError: Invalid equation format in reaction: {reaction}")
            continue  # Skip this reaction
            
        (
            reactionCount,
            reactants,
            reactantCoeffs,
            products,
            productCoeffs,
            rates,
            reaction_type
        ) = result
        if verbose:
            print(f"\tReaction Count = {reactionCount}")
            print(f"\tReactants = {reactants})")
            print(f"\tCoeffs = {reactantCoeffs})")
            print(f"\tProducts = {products}")
            print(f"\tProducts Coeffs = {productCoeffs}")
            print(f"\tRates = {rates}")
            print(f"\tReaction Type = {reaction_type}")
            print(f"**********")

        # Extract rate names and values
        # Note: There might be different rates & values;
        names = []
        values = []
        for rate in rates:
            if '=' in rate:
                name, value = rate.split('=')
                names.append(name.strip())
                values.append(float(value.strip()))
            else:
                names.append(rate.strip())
                values.append(-1)
                
        # Ensure that names and values have at least 3 elements, filling with defaults if needed
        while len(names) < (reactionCount+1):
            names.append("")  # Default empty string for missing names
        while len(values) < (reactionCount+1):
            values.append(-1)  # Default value of -1 for missing values
        
        if reactionCount == 1:  # We add only 1 case, easy
            # Create a Reaction object and add to the list
            reaction_obj = Reaction(
                equation=formatFwdReaction(reactants, reactantCoeffs, products, productCoeffs),
                rateName=names[0],
                rateValue=values[0],
                reactants=reactants,
                reactant_coeffs=reactantCoeffs,
                products=products,
                product_coeffs=productCoeffs,
                reactionType=reaction_type,
                rateModifier=names[1],
                rateModifierValue=values[1]
            )
            parsed_reactions.append(reaction_obj)
        
        elif reactionCount == 2:  # We add 2 objects for bi-directional reactions
            reaction_fwd = Reaction( equation=formatFwdReaction(reactants, reactantCoeffs, products, productCoeffs), rateName=names[0], rateValue=values[0], reactants=reactants, reactant_coeffs=reactantCoeffs, products=products, product_coeffs=productCoeffs, reactionType=reaction_type, rateModifier=names[2], rateModifierValue=values[2] )
            
            reaction_rev = Reaction( equation=formatFwdReaction(products, productCoeffs, reactants, reactantCoeffs), rateName=names[1], rateValue=values[1], reactants=products, reactant_coeffs=productCoeffs, products=reactants, product_coeffs=reactantCoeffs, reactionType=reaction_type, rateModifier=names[2], rateModifierValue=values[2] )
            parsed_reactions.append(reaction_fwd)
            parsed_reactions.append(reaction_rev)
        
        else:
            print(f"\tError: Invalid reaction count in: {reaction}")
            continue
    
    return parsed_reactions



def parseEquation(equation):
    # Remove comments if present
    equation = equation.split("#")[0].strip()

    # Split components based on commas
    parts = [p.strip() for p in equation.split(",")]

    # Ensure there are at least three components (reactants/products + at least one rate)
    if len(parts) < 2:
        print(f"Error: Malformed equation - {equation}")
        return None, None, None, None, None, None, None

    #Grab the Biochemical Equation
    equation_part = parts[0]

    # Extract reaction type (if specified) or default to "MASS_ACTION"
    # and reaction rates (either 1, 2 or 3);
    known_types = {"LIPID", "PLATELET", "FLOW", "MASS_ACTION","FUNCTION"}
    if parts[-1] in known_types:
        reaction_type = parts[-1]
        rates = parts[1:-1]  # capture rates only when reaction type is known
        
        # SPECIAL CASE FOR FUNCTION: Separate rate and rate modifier
        if reaction_type == "FUNCTION":
            if len(rates) < 2:
                print(f"Error: FUNCTION reaction must have a rate and a rateModifier - {equation}")
                return None, None, None, None, None, None, None
            
            ratePart = rates[0]  # First part is the actual rate
            rateModifier = ",".join(rates[1:])  # Join the rest into a single string
   
            rates = {ratePart, rateModifier}
        
            #if verbose:
            #    print(f"\tRate Part: {ratePart}")
            #    print(f"\tRate Modifier: {rateModifier}")
                    
    else:
        reaction_type = "MASS_ACTION"
        rates = parts[1:]  # treat all parts as rates when no reaction type is specified

    # Determine reaction direction and count
    if '<->' in equation_part:
        reactionCount = 2
        lhs, rhs = equation_part.split('<->')
        expected_rate_count = 2 if reaction_type == "MASS_ACTION" else 3
    elif '->' in equation_part:
        reactionCount = 1
        lhs, rhs = equation_part.split('->')
        expected_rate_count = 1 if reaction_type == "MASS_ACTION" else 2
    else:
        print(f"Error: Equation must contain '->' or '<->': {equation}")
        return None, None, None, None, None, None, None
        
    # Validate rate count
    if len(rates) != expected_rate_count:
        print(f"Error: Expected {expected_rate_count} rates for {reaction_type} reaction '{equation_part}', but got {len(rates)}.")
        return None, None, None, None, None, None, None

    # Split LHS and RHS on "+" and ignore white space
    reactants = [r.strip() for r in lhs.split('+')] if lhs.strip() else []
    products = [p.strip() for p in rhs.split('+')] if rhs.strip() else []

    # Call extract_coefficients to process reactants and products
    reactantCoeffs, reactants = extract_coefficients(reactants)
    productCoeffs, products = extract_coefficients(products)

    return reactionCount, reactants, reactantCoeffs, products, productCoeffs, rates, reaction_type

# Extract the Coefficients for a Set of Reactants
def extract_coefficients(terms):
    # Check if the input contains only an empty string
    if len(terms) == 1 and terms[0] == '':
        return [], []

    coeffs = []
    species = []
    for term in terms:
        # Match terms with coefficients followed by '*'
        match = re.match(r'(\d*)\s*\*\s*(.*)', term)
        if match:
            coeff_str, species_name = match.groups()
            coeff = int(coeff_str) if coeff_str else 1
            coeffs.append(coeff)
            species.append(species_name.strip())
        else:
            # Match terms with optional coefficients without '*'
            match = re.match(r'(\d*)\s*(.*)', term)
            if match:
                coeff_str, species_name = match.groups()
                coeff = int(coeff_str) if coeff_str else 1
                coeffs.append(coeff)
                species.append(species_name.strip())
            else:
                # Handle unexpected formats
                print(f"Error: Invalid term format: {term}")
                coeffs.append(1)
                species.append(term.strip())
    return coeffs, species

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
                    elif "= PLATELET_SITES" in line:
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
        "PLATELET_SITES": {
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
        if mode in ["LIPID", "PLATELET","PLATELET_SITES"]:
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
    duplicates = []  # List to store any duplicates found

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
    
        # Define headers for the initial conditions table
        p_headers = ["Parameter Name", "Value"]
    
        # Create a data list for parsed initial conditions
        p_data = [(p.name, p.value) for p in parsed_params]

        # Print table header for initial conditions
        print(f"\n{'=' * 50}")
        print(f"{'Step 1(f): Parse Specified Parameters':^30}")
        print(f"{'=' * 50}")
        print(f"{p_headers[0]:<30} | {p_headers[1]:>10}")
        print("-" * 50)

        # Print each parsed initial condition in the data list
        for name, value in p_data:
            if isinstance(value, float):  # Format floats to 2 decimal places
                print(f"{name:<30} | {value:>10.2f}")
            else:  # Print strings without decimal formatting
                print(f"{name:<30} | {value:>10}")
            
        # Print done message and finish with a line
        print(f"{'=' * 50}\n")

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
        print("Error: Flow list contains multiple unique values.")
        exit(-1)
    else:
        uniqueFlowRate = flowRate[0];

    if len(bad_rates) > 0:
        print(f"Error: We found some biochemical rates with different dimensions.")
        for rate, column_sum in bad_rates.items():
            if(not rate == uniqueFlowRate):
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
        self.name = name
        self.value = value

    def __repr__(self):
        return f"InitialCondition(name={self.name}, value={self.value})"

class Parameter:
    def __init__(self, name, value):
        self.name = name
        self.value = value
    
    def __repr__(self):
        return f"Parameter(name={self.name}, value={self.value})"

class Reaction:
    index_counter = 0  # Class variable to keep track of the index

    def __init__(self, equation, rateName, rateValue, reactants, reactant_coeffs, products, product_coeffs,reactionType,rateModifier,rateModifierValue):
        self.equation = equation
        self.rateName  = rateName
        self.rateValue = rateValue
        self.reactants = reactants
        self.reactant_coeffs = reactant_coeffs
        self.products = products
        self.product_coeffs = product_coeffs
        self.reactionType = reactionType
        self.rateModifier = rateModifier
        self.rateModifierValue = rateModifierValue
        self.index = Reaction.index_counter  # Assign the current index
        Reaction.index_counter += 1  # Increment the index for the next reaction

    def __repr__(self):
        return (f"Reaction(index={self.index}, equation={self.equation}, rateName={self.rateName}, "
                f"rateValue={self.rateValue}, reactants={self.reactants}, "
                f"reactant_coeffs={self.reactant_coeffs}, products={self.products}, "
                f"product_coeffs={self.product_coeffs}, rxn_type={self.reactionType}, "
                f"rate_modifier = {self.rateModifier}, rate_modifierValue = {self.rateModifierValue}")

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

def create_matlab_multipleFileOutput(input_file: str, outputPrefix: str,s: Stoich, parsed_reactions: list, species: list, rates: list, uniqueRates: list, uniqueNBS: list, uniqueFlow: list, parsed_ICs: list, parsed_params: list, verbose: bool = False):
    
    if verbose:
        print('-' * 50)
        print(f"Step 3: Creating Matlab Output Files")
    
    Ns = len(s.species)
    Nr = len(s.rates)
    
    NumReactions = len(parsed_reactions)
        
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
        fParam.write(" = 1; \n")                            #<- Should be able to set flow rates

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
            #fParam.write(str(uniqueNBS[i]))
            fParam.write(" = 1; \n")                          #<- Should be able to set flow rates
        
        fParam.write("\n")
        fParam.write("nbs = [ ")
        fParam.write(uniqueNBS[0])
        #fParam.write(str(uniqueNBS[0]))
        for i in range(1, len(uniqueNBS)):
            fParam.write(", ")
            fParam.write(uniqueNBS[i])
        fParam.write(" ];\n\n\n")

    f.write("% Set the Kinetic Parameters\n")
    f.write(f"{ParamPrefix}\n\n")

    ########################################################################################
    # Set Initial Conditions and Upstream Flow Values: Default value of 0 if not specified #
    ########################################################################################
    
    ic_dict = {ic.name: ic.value for ic in parsed_ICs}  # Convert list of objects to a dictionary

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
    
    if len(uniqueFlow) >=1:
        fIC.write("% Flow Rate Parameters \n")
        for i in range(len(uniqueFlow)):
            fIC.write(uniqueFlow[i])
            fIC.write(" = 1; \n")                    #<- Should be able to set flow rates
        
        fIC.write("\n")
        fIC.write("flow = [ ")
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
        args.append('flow')

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

    f.write("% ODEs from reaction equations \n\n")

    if verbose: print("\tWriting ODEs now....")

    for i in range(Ns): #For each Species
        speciesName = species[i]
        isLipidSpecies = "L_noTF" in speciesName or "L_TF" in speciesName
        isOnEmptyLipid = "_s" in speciesName and (not "_st" in speciesName) and (not "_sn" in speciesName)
        isOnTFLipid = "_st" in speciesName
        isOnSilica = "_sn" in speciesName
        f.write("% ")
        f.write(str(transform_string(species[i])))
        f.write("\n dy(")
        f.write(str(i+1))
        f.write(")  =")
        for j in range(Nr): #For Each Reaction Reaction
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
                    f.write(parsed_reactions[j].rateModifier)
                                        
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
                    f.write(parsed_reactions[j].rateModifier)
                if isLipidReaction and ( not isLipidSpecies ) and (isOnEmptyLipid or isOnTFLipid): #V_s
                    f.write("/")
                    f.write(parsed_reactions[j].rateModifier)
                    
                if isFlowReaction:
                    f.write(" * ")
                    f.write(parsed_reactions[j].rateModifier)

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
parsedPlateletSpecies = parseList(plateletSpecies,"PLATELET",verbose)
parsedPlateletSites   = parseList(plateletSites,"PLATELET_SITES",verbose)
parsedFunctions       = parseFunction(functionsDefined,verbose)
parsedParameters      = parseParameters(parameters, verbose)

###########################################
# Step 2: Parse the Biochemical Reactions #
###########################################

# (2a) Do initial parsing of reactions
parsed_reactions = parseReactions(biochemicalReactions, False)

# (2b) Check for & Remove Duplicates in reaction list
unique_reactions = set()
duplicates = [reaction for reaction in parsed_reactions if reaction in unique_reactions or unique_reactions.add(reaction)]

print("-" * 50)
print(f"Begin: Consistency Checking Biochemical Reactions:")

print("\t(1) Sanity Check: Are there duplicate Reactions?")
if duplicates:
    print("\t\t---> Suspicious Duplicate Reactions Found. Removing them:")
    print("\n".join(duplicates))
else:
    print("\t--->Sanity Check PASSED: No Duplicate Reactions Found.")

# (2c) Determine the Number of reaction types, List of Unique Species, Unique Reaction Rates
species  = []
rates    = []
nbs      = []
flow     = []
flowRate = []

# Dictionary to track the flow reactions for each species
flow_species_count  = {}
badFlowReactions    = []

# Dictionary counting occurrences of each reaction type
reaction_type_count = {"LIPID": 0, "PLATELET": 0, "MASS_ACTION": 0, "FLOW": 0, "FUNCTION": 0}

# Iterate through each reaction (now unidirectional after parsing)
for reaction in parsed_reactions:
    rates.append(reaction.rateName)
    
    # Increment the count for the reaction type
    if reaction.reactionType in reaction_type_count:
        reaction_type_count[reaction.reactionType] += 1
     
    if(reaction.reactionType=="LIPID" or reaction.reactionType=="PLATELET"):
        nbs.append(reaction.rateModifier)
    
    if(reaction.reactionType=="FLOW"):
        flow.append(reaction.rateModifier)
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

#Things to check:
# (1) If Cup = C_IN symbolically; we need to make sure that variable EXISTS!
# (2) Pull all the in-line parameters from the equations (if the values do NOT equal to the numerical value)
#(3) All specified parameters should occur somewhere! (Rxn rates, in a function).
#(4) All function arguments should be
#   - Specified Parameters    (go at the end of parameters list in p variable, real values) p[extra]
#   - Unspecified Parameters  (go at the end of parameters list in p variable, default) p[extra]
#   - Given kinetic rates     (already listed as a kinetic rate in another equation) use it's p[i]
#   - One of the biochemical species unique_species (use it's y[i] value)
#(5) Specified Parameters might be:
#    - Arguments to functions (end of parameters)  <- p
#    - NBS or Platelet binding rates               <- NBS
#    - Flow rate (this is potentially special      <- Flow goes separate
#(6) Make it possible to be able to in-line declare a biochemical rate:
#    - Check that any rate for a parameter is the same (error check)
#(7) DONE: For kflow reactions; do we always want an in-and-out; then maybe make it 1 equation:
#(8) DONE: All the FLOW reactions should have the same parameter (this is an error check)
#(9) DONE: All Flow reactions should have either [1,0] or [0,1] species reactants; should check.

#Allow this one too:
# Perhaps Remake THE FLOW Equation Form:    C <-> C, kflow, C_up, FLOW

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
print(f"Parameters: {parsedParameters}")

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

create_matlab_multipleFileOutput(sys.argv[1], prefix, stoich, parsed_reactions, unique_species, rates, unique_rates, unique_nbs, unique_flow, parsed_ICs, parsedParameters, verbose)

######################################################
# Step 5: Determine Conserved Quantities (If Needed) #
######################################################

if determineConserved:
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



