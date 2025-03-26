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

    - 4 Reaction Types:
        * MASS_ACTION: Default (if no type given)
        * LIPID: Binding on/off lipid (competition)
        * PLATELET: Binding on/off platelet (no competition)
        * FLOW: species entering or exiting reaction zone.
        * DILUTION: concentration changing due to reaction zone

    Example Reactions:
        A + 2 * B -> C , k_1 #Forward
        A + 2 * B <-> C , k_1 , k_2 #Reversible

    Lipid Binding Support:
        - Allows non-mass action terms for lipid binding:
            L_TF + II <-> II_st, kon_ii, koff_ii, nbs_ii, LIPID
        - kon_ii units: 1/(concentration * time * binding sites)

    Supported Features:
    -----------------------------------
      - Support for non-mass action lipid/platlet binding.
      - Outputs Species and rates output in input order.
      - Supports pure synthesis/degradation (e.g., "-> A", "B ->").
      - Consolidates duplicate kinetic rates.
      - Splits stoichiometric matrix for A + B -> A + C reactions.
      - Removes duplicate reactions (even one side of bidirectional ones).
      - Checks reaction rate dimensions for consistency.
      - Checks for valid reaciton types (i.e., flow types must all have the same FLOW rate)

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
    - Currently requires there to be some nbs <- Add condition to check if no nbs
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

def parseReactions(reactions):
    parsed_reactions = []

    for reaction in reactions:
        print(f"Processing Reaction: '{reaction}'")
        
        # Parse the equation
        result = parseEquation(reaction)
        if result is None:
            print(f"\tError: Invalid equation format in reaction: {reaction}")
            continue

        (
            reactionCount,
            reactants,
            reactantCoeffs,
            products,
            productCoeffs,
            rates,
            reaction_type
        ) = result
        
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
        return None, None, None, None, None, None

    #Grab the Biochemical Equation
    equation_part = parts[0]

    # Extract reaction type (if specified) or default to "MASS_ACTION"
    # and reaction rates (either 1, 2 or 3);
    known_types = {"LIPID", "PLATELET", "FLOW", "MASS_ACTION"}
    if parts[-1] in known_types:
        reaction_type = parts[-1]
        rates = parts[1:-1]  # capture rates only when reaction type is known
    else:
        reaction_type = "MASS_ACTION"
        rates = parts[1:]  # treat all parts as rates when no reaction type is specified

    # Determine if the equation has <-> or ->
    if '<->' in equation_part:
        reactionCount = 2
        lhs, rhs = equation_part.split('<->')
    elif '->' in equation_part:
        reactionCount = 1
        lhs, rhs = equation_part.split('->')
    else:
        print(f"Error: Equation must contain '->' or '<->': {equation}")
        return None, None, None, None, None, None

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

def parseInitialConditions(initialConditions):
    parsed_ICs = []
    
    for ic_str in initialConditions:
        # Split the string at the '=' character and strip any leading/trailing whitespace
        try:
            name, value_str = ic_str.split("=")
            name = name.strip()
            value = float(value_str.strip())  # Convert value to float (could be int or float)
                        
            # Create InitialCondition object and append to the result list
            parsed_ICs.append(InitialCondition(name, value))
        
        except ValueError:
            print(f"Error: The string '{ic_str}' could not be parsed.")
        except Exception as e:
            print(f"Unexpected error while parsing '{ic_str}': {e}")
    
    return parsed_ICs

#################
# Class Objects #
#################

class InitialCondition:
    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __repr__(self):
        return f"InitialCondition(name={self.name}, value={self.value})"

    def __eq__(self, other):
        if isinstance(other, InitialCondition):
            # Two InitialCondition objects are equal if they have the same name
            return self.name == other.name and self.value == other.value
        return False

    def __hash__(self):
        # We want the hash based on the 'name' and 'value' attributes.
        return hash((self.name, self.value))

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

def create_matlab_multipleFileOutput(input_file: str, outputPrefix: str,s: Stoich, parsed_reactions: list, species: list, rates: list, uniqueRates: list, uniqueNBS: list, v: bool = False):
    Ns = len(s.species)
    Nr = len(s.rates)
    
    NumReactions = len(parsed_reactions)
        
    #if Nr==NumReactions:
    #    print(f"Yay!")
    #    exit(-1)
    #else:
    #    print(f"Nooooo!")
    #    exit(-1)
    
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
        fParam.write(" = 1; \n")

    fParam.write("\np = [ ")
    fParam.write(uniqueRates[0])
    for i in range(1, len(uniqueRates)):
        fParam.write(", ")
        fParam.write(uniqueRates[i])
    fParam.write(" ];\n\n\n")
    
    fParam.write("% Binding Site Parameters \n")
    for i in range(len(uniqueNBS)):
        fParam.write(uniqueNBS[i])
        #fParam.write(str(uniqueNBS[i]))
        fParam.write(" = 1; \n")
        
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

    fIC.write("% Initial Conditions \n")
    for i in range(Ns):
        fIC.write(transform_string(species[i]))
        fIC.write("_IC")
        fIC.write(" = 0; \n")

    ##Line that's too long;
    fIC.write("\ninit_cond = [ ")
    fIC.write(transform_string(species[0]))
    fIC.write("_IC")
    for i in range(1,Ns):
        fIC.write(", ")
        fIC.write(transform_string(species[i]))
        fIC.write("_IC")
    fIC.write(" ];\n\n\n")
    
    f.write("% Set the Initial Conditions\n")
    f.write(f"{ICFilePrefix}\n\n")

    f.write("options = odeset('RelTol',1e-12,'AbsTol',1e-23);\n\n\n")

    f.write("%------------------------- Main Solve ----------------------%\n")
    f.write("[time,y] = ode15s(@(t,y)RHS(t,y,p,nbs), t_start:1:t_final, init_cond, options);\n")
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


    f.write("function dy = RHS(t,y,p,nbs)\n\n")
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

    f.write("% Rename Binding Site \n")
    for i in range(len(uniqueNBS)):
        f.write(str(uniqueNBS[i]))
        f.write(" = nbs(")
        f.write(str(i+1))
        f.write(");  \n")
    f.write("\n\n")

    f.write("% ODEs from reaction equations \n\n")

    if v: print("Writing ODEs now....\n")

    for i in range(Ns):
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
        for j in range(Nr):
            isLipidReaction = (parsed_reactions[j].reactionType=="LIPID")
            
            #If the reaction j reduces the amount of species i;
            if (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) < 0):
                f.write("  -  ")
                f.write(str(rates[j])) #Reaction Rate
                #Writing the Reaction Rate isLipidSpecies (true/false); isLipidReaction (true/false)
                #if isLipidReaction: #This must be a koff
                # V_s -> V + L; L is correct
                # (V_s)_dt = -rate
                # (V)_dt   = +rate
                # L_dt     = +rate*nbs
                if isLipidReaction and not (isLipidSpecies or isOnEmptyLipid or isOnTFLipid): #V_s
                    f.write("/")
                    f.write(parsed_reactions[j].rateModifier)
                for k in range(Ns):
                    if (not math.isnan(s.reactants[k][j])) and (int(s.reactants[k][j]) <= 0):
                        f.write(" * ")
                        f.write(transform_string(species[k]))
                        if (abs(int(s.reactants[k][j])) != 1) and (s.reactants[k][j] != 0):
                            f.write("^")
                            f.write(str(abs(int(s.reactants[k][j]))))
            #If the reaction j increases the amount of species i
            elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) > 0):
                f.write("  +  ")
                f.write(str(rates[j])) #Reaction Rate
                #Writing the Reaction Rate isLipidSpecies (true/false); isLipidReaction (true/false)
                if isLipidSpecies: #This must be a koff
                    f.write("*")
                    f.write(parsed_reactions[j].rateModifier)
                if isLipidReaction and ( not isLipidSpecies ) and (isOnEmptyLipid or isOnTFLipid): #V_s
                    f.write("/")
                    f.write(parsed_reactions[j].rateModifier)
                for k in range(Ns):
                    ##Double check that this makes sense; perhaps should be reactants.
                    if (not math.isnan(s.stoich[k][j])) and (int(s.stoich[k][j]) <= 0):
                        f.write(" * ")
                        if (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) > 1):
                            f.write(str(int(s.stoich[i][j])))
                            f.write(" * ")
                        f.write(transform_string(species[k]))
                        if (abs(int(s.reactants[k][j])) != 1) and (s.reactants[k][j] != 0):
                            f.write("^")
                            f.write(str(abs(int(s.reactants[k][j]))))
            
            #If the reaction j leaves species i untouched.
            elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) == 0):
                f.write("  +  ")
                f.write(" 0 ")
        if v: print(species[i]," complete")
        f.write(";\n\n")


    f.write("\n\n\n\n")
    f.write("end")
    
    
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


#############
# Main Code #
#############

verbose        = True
specialVerbose = True
previewVal     = 10

######################################
# Step 0: Preprocessing              #
#   - Read in Txt File               #
#   - Remove white space             #
#   - Store the biochemicalReactions #
######################################

if verbose:
    print('-' * 50)
    print(f"Step 0: Preprocessing")

# Initialize the initial biochemical arrays and counters
biochemicalReactions = []
initialConditions    = []

numReactionsReadIn         = 0
numInitialConditionsReadIn = 0

try:
    with open(sys.argv[1], 'r') as file:
        for line in file:
            line = line.rstrip()  # Remove trailing newline characters
            
            # Ignore blank lines or lines starting with '#'
            if not line or line.startswith('#'):
                continue
            
            # Remove everything after the first '#' (if any)
            line = line.split('#')[0].rstrip()  # Keep part before '#', remove trailing spaces

            # Split the line by commas
            fields = line.split(',')
            
            # Process as biochemical reaction if there are multiple comma-separated values
            if len(fields) > 1:
                biochemicalReactions.append(line)
                numReactionsReadIn += 1
            else:
                # Handle single-field lines (initial conditions)
                initialConditions.append(line)
                numInitialConditionsReadIn += 1
                print(f"\t\tInitial Condition:{line}")
                    
except IOError:
    sys.exit(f"Couldn't open {sys.argv[1]}")

# Get the input file name
input_file = sys.argv[1]
    
# Extract the PREFIX from the input file name
prefix, _ = os.path.splitext(input_file)

# Output the counts to the screen
if verbose:
    print(f"\t\tNumber of Biochemical Equations: {numReactionsReadIn}")
    print(f"\t\tNumber of Species with Initial Conditions: {numInitialConditionsReadIn}")
    print(f"DONE: Step 0 Preprocessing of {sys.argv[1]}")
    print('-' * 50)

#########################################
# Step 1: Parse the Initial Conditions  #
#########################################

if verbose:
    print(f"Step 1: Parse the Intial Conditions")

# Parse the reactions
parsed_ICs = parseInitialConditions(initialConditions)

if verbose:
    # Print parsed objects
    for ic in parsed_ICs:
        print(f"\t{ic}")
    print(f"DONE: Step 1 Parsed ")
    print('-' * 50)

# Define the output file names
stoich_output_file   = f"{prefix}Stoich.csv"
reactant_output_file = f"{prefix}Reactants.csv"
product_output_file  = f"{prefix}Products.csv"

#matlab_output_file = f"{prefix}Matlab.m"
#IC_output_file     = f"{prefix}ICs.m"
#param_output_file  = f"{prefix}Params.m"
#rename_output_file = f"{prefix}Rename.m"

if verbose:
    print(f"\tNumber of biochemical reactions: {len(biochemicalReactions)}")
    print(f"\tSome of the Biochemical Reactions:")
    num_reactions_to_print = min(previewVal, len(biochemicalReactions))
    for i in range(num_reactions_to_print):
        print(f"\t\t" + biochemicalReactions[i])
    print('-' * 50)

###########################################
# Step 2: Parse the Biochemical Reactions #
###########################################

if verbose:
    print(f"Step 2: Parsing the Biochemical Reacitons")

# Parse the reactions
parsed_reactions = parseReactions(biochemicalReactions)

#exit(-1);

########
# Error: Check for Duplicates in reaction list
#######

reaction_list = parsed_reactions

# Use a set to find duplicates
unique_reactions = set()
duplicates = []

for reaction in reaction_list:
    if reaction in unique_reactions:
        duplicates.append(reaction)
    else:
        unique_reactions.add(reaction)

# Output the results
if duplicates:
    if verbose: print("Duplicate Reactions Found. Removing them:")
    for duplicate in duplicates:
        if verbose: print(duplicate)
else:
    if verbose: print("No Duplicate Reactions Found.")
##End Search for Duplicates

# Determine the Number of Reactions, List of Unique Species, Unique Reaction Rates
species = []
rates = []
nbs = []

# Create a dictionary to count the occurrences of each reaction type
reaction_type_count = {"LIPID": 0, "PLATELET": 0, "MASSACTION": 0, "FLOW": 0}

# Iterate through each reaction
for reaction in parsed_reactions:
    #All reactions are now (after parsing) unidirectional
    rates.append(reaction.rateName)
    
    # Increment the count for the reaction type
    if reaction.reactionType in reaction_type_count:
        reaction_type_count[reaction.reactionType] += 1
     
    if(reaction.reactionType=="LIPID" or reaction.reactionType=="PLATELET"):
        nbs.append(reaction.rateModifier)
    # Append Reactants and Products if they are not empty
    if reaction.reactants:
        species.extend(reaction.reactants)
    if reaction.products:
        species.extend(reaction.products)

# Get unique entries
#unique_species = unique_entries_only(species)
#unique_rates   = unique_entries_only(rates)

# Output the counts of each reaction type in a table format
print(f"-----------")
print(f"{'Reaction Type':<15} {'Count'}")
print(f"-----------")

for reaction_type, count in reaction_type_count.items():
    print(f"{reaction_type:<15} {count:>5} reactions")

print(f"-----------")

unique_species = unique_entries_in_order(species)
unique_rates   = unique_entries_in_order(rates)
unique_nbs     = unique_entries_in_order(nbs)

print(f"All BindingSites: {unique_nbs}")
#exit(-1)
# Output Duplicate Rates

#Note: We will need BOTH rates and unique_rates.
#      The rates will be needed for the stoichiometric matrix.

#print(f"All Rates: {rates}")

# Output the partial results:
if specialVerbose:
    #print(f"\tDONE. Successfully Parsed Reactions")
    #print(f"Number of Total Reactions: {len(parsed_reactions)}")
    #print(f"\t{list_to_string(rates)}")
    #print(f"Number of Unique Species: {len(unique_species)}")
    #print(f"\t{list_to_string(unique_species)}")
    #print(f"Number of Unique Rates: {len(unique_rates)}")
    #print(f"\t{list_to_string(unique_rates)}")
    # Output the parsed reactions to the screen with reaction numbers
    #print(f"{'Num':<3} {'Equation':<20} {'KFwd':<10} {'Reactants':<10} {'RCoeffs':<10} {'Products':<10} {'PCoeffs':<10}")

    print(f"{'Num':<3}\t{'Equation':<20}\t{'Reactants':<10} {'Products':<10}")
    #print('-' * 50)
    num_reactions_to_print = len(parsed_reactions) #min(previewVal, len(parsed_reactions))
    for i, reaction in enumerate(parsed_reactions[:num_reactions_to_print], start=1):
        #print(f"{i:<3} {reaction.equation:<20} {reaction.kfwd:<10} {list_to_string(reaction.reactants):<10} {list_to_string(reaction.reactant_coeffs):<10} {list_to_string(reaction.products):<10} {list_to_string(reaction.product_coeffs):<10}")
        print(f"{i:<3}\t{reaction.equation:<20}\t {list_to_string(reaction.reactants):<10}\t {list_to_string(reaction.products):<10}")
    #print('-' * 50)

######################################
# Step 2: Create Stochiometric Matrix
#   - Columns: Species
#   - Rows:    Reactions
######################################

if verbose:
    print(f"Step 2: Creating Stoichiometry Matrix")

#unique_species.sort() #For ease in outputting; put variables in sorted order.

s = Stoich(unique_species,parsed_reactions)

if verbose:
    s.to_csv(stoich_output_file,"S")
    s.to_csv(reactant_output_file,"R")
    s.to_csv(product_output_file,"P")
    print(f"DONE! Successfully created Stoichiometry Matrix (output file {stoich_output_file})")
    print('-' * 50)

######################################
# Step 3: Create Matlab File Output
######################################

if verbose:
    print(f"Step 3: Creating Matlab Output Files")

create_matlab_multipleFileOutput(sys.argv[1], prefix, s, parsed_reactions, unique_species, rates, unique_rates, unique_nbs, verbose)

if verbose:
    print(f"DONE! Successfully created Matlab Files")
    print('-' * 50)

########################################
# Step 4: Determine Conserved Quantity #
########################################

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

##################################
# Step 4: Checking the Dimension #
##################################

bad_rates = s.check_rates()

if len(bad_rates) > 0:
    if verbose: print(f"We have some rates with different dimensions.")
    for rate, column_sum in bad_rates.items():
        if verbose: print(f"Rate {rate}: Sum of columns = {column_sum}")
else:
    if verbose: print("All rates are correct dimension")
