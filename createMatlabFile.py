#!/usr/bin/env python3

import sys, re, io, os, textwrap, math, csv

from collections import defaultdict


###############################
# Error and Usage Information #
###############################

# Check if the correct number of arguments is provided
if len(sys.argv) < 2:
    print(textwrap.dedent("""
    
    Python code that generates a Matlab Code for Law of Mass action Reactions
    
    Use: python3 createMatlabFile.py StaticCoag.txt

    StaticCoag.txt:
    - Can contain comments that are proceeded with "#" (will be ignored)
    - Can contain whitespace (will be ignored)
    - Two types of structure are allowed for biochemical reactions: Forward, Reversible
        
    SINGLE_REACTION , RATE_VARIABLE
    ex:
        A + 2 * B -> C , k_1

    In the case of reversible reactions, two rates must be specified, the forward reaction
        is always FIRST, eg:

        A + 2 * B <-> C , k_1 , k_2
    
    **Warning: Any invalid reactions will be skipped (and output to the screen).
  
    Output File: StaticCoagMatlab.m (assumes prefix).

    Requirements:
        1) Only the operators '*', '+', '<->', and '->' are allowed
        2) Reactions should be pre-simplified, this code will NOT reduce algebra
        3) Species names can only contain [0-9a-zA-Z_:]
        4) Species names MUST begin with a letter - no leading numbers allowed
        
    
    New Features (11/07/2024):
        1) Allows for specifying the value of the kinetic rates in the reacitons:
            Allowable:
                A + B -> C, k1=0.1
                A + B <-> C, kon, koff=100
       
    MIGHT NOT BE VALIE
    
        New Features (08/19/2024):
        1) Allows for pure synthesis or pure degradation. ( -> A, B ->)
        2) Handles comments/white space in the biochemical equation file.
        3) Handles duplicate kinetic rates (i.e., p vector is consolidated)
        4) Handles A + B -> A + C reactions (splits stochiometric matrix)
        5) Checks for (and removes) duplicate reactions. (i.e., identical reaction)
            -> Even if the reaction is 1 side of a bi-directional reaction.
        6) Checks that the dimension of reactions rates is the SAME (?).
    
    Coag Specific Changes:
        1) Creates variables to track active and in-active bound lipid sites (s and _st)
            -> Counts the number of "s" and "st" in each species name.
        2) Separately creates outputs for initial conditions, parameters and code.
            -> Could do a --coag flag = false to do the original mode.
        
    Other Features to Consider:
        1) Keep the species/rates in the same order the StaticCoag.txt file was in <-Helpful
            -> Currently sorting species alphabetically, rates are not in any particular order.
        2) Set the rates of reactions and initial conditions based on a file. <- Helpful
            -> Could have a user input values (separate file) and use 1 for missing values.
        3) Original version allowed "=" operator. Not sure what this was for.
        4) Text wrap for the long lines in Matlab. (Currently half-implemented)

    Based on previous Python Code by:
    Michael Stobb (Originally written on 8/4/2023)
      
    Current Version:
    Suzanne Sindi
    11/07/2024

    """))
    sys.exit("Usage: python3 createMatlabFile.py StaticCoag.txt")

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

##Parses the reactions from the input file.
##Checks for valid format, skips any invalid lines (outputs them as warnings).
#def parseReactions(reactions):
#    parsed_reactions = []
#
#    for reaction in reactions:
#        # Split the reaction line by commas
#        parts = [part.strip() for part in reaction.split(',')]
#        
#        # Check if we have at least the equation and one rate constant
#        if len(parts) < 2:
#            print(f"\tError: Invalid format {reaction} (not enough parts)")
#            continue
#
#        # Validate the direction and number of rate constants
#        if len(parts) == 2:
#            if '->' not in parts[0] or '<->' in parts[0]:
#                print(f"\tError: Expected '->' in reaction: {reaction}")
#                continue
#            direction = '->'
#            equation = parts[0]
#            kfwd = parts[1]
#            krev = 'N/A'
#        elif len(parts) == 3:
#            if '<->' not in parts[0]:
#                print(f"\tError: Expected '<->' in reaction: {reaction}")
#                continue
#            direction = '<->'
#            equation = parts[0]
#            kfwd = parts[1]
#            krev = parts[2]
#        else:
#            print(f"\tError: Invalid format {reaction}")
#            continue
#
#        # Parse the equation
#        result = parseEquation(equation)
#        if result is None:
#            print(f"\tError: Invalid equation format in reaction: {reaction}")
#            continue
#
#        reactionCount, reactants, reactantCoeffs, products, productCoeffs = result
#                
#        if reactionCount == 1: #We add only 1 case, easy
#            # Create a Reaction object and add to the list
#            reaction_obj = Reaction(
#                equation=equation,
#                kfwd=kfwd,
#                reactants=reactants,
#                reactant_coeffs=reactantCoeffs,
#                products=products,
#                product_coeffs=productCoeffs
#            )
#            parsed_reactions.append(reaction_obj)
#        
#        elif reactionCount == 2: #We add 2 objects;
#            reaction_fwd = Reaction(
#                equation=formatFwdReaction(reactants, reactantCoeffs, products, productCoeffs),
#                kfwd=kfwd,
#                reactants=reactants,
#                reactant_coeffs=reactantCoeffs,
#                products=products,
#                product_coeffs=productCoeffs
#            )
#            
#            reaction_rev = Reaction(
#                equation=formatFwdReaction(products, productCoeffs,reactants, reactantCoeffs),
#                kfwd=krev,
#                reactants=products,
#                reactant_coeffs=productCoeffs,
#                products=reactants,
#                product_coeffs=reactantCoeffs
#            )
#            parsed_reactions.append(reaction_fwd)
#            parsed_reactions.append(reaction_rev)
#        else:
#            print(f"\tError: Invalid equation format in reaction: {reaction}")
#            continue
#
#    return parsed_reactions

def parseReactions(reactions):
    parsed_reactions = []

    for reaction in reactions:
        # Split the reaction line by commas
        parts = [part.strip() for part in reaction.split(',')]
        
        # Check if we have at least the equation and one rate constant
        if len(parts) < 2:
            print(f"\tError: Invalid format {reaction} (not enough parts)")
            continue

        # Separate the equation and rate constants
        equation = parts[0]
        rate_constants = parts[1:]

        # Sanity check for number of rate constants based on the reaction direction
        if '<->' in equation:
            # For uni-directional reactions, we expect exactly 1 rate constant
            if len(rate_constants) != 2:
                print(f"\tError: Bi-directional reaction '{reaction}' must have exactly 2 rate constants.")
                continue
        elif '->' in equation:
            # For bi-directional reactions, we expect exactly 2 rate constants
            if len(rate_constants) != 1:
                print(f"\tError: Uni-directional reaction '{reaction}' must have exactly 1 rate constant.")
                continue
        else:
            print(f"\tError: Invalid reaction format (must include '->' or '<->'): {reaction}")
            continue
            
        # Parse the equation
        result = parseEquation(equation)
        if result is None:
            print(f"\tError: Invalid equation format in reaction: {reaction}")
            continue

        reactionCount, reactants, reactantCoeffs, products, productCoeffs = result
        
        # Parse the rate constants;
        names  = []
        values = []
        for rateString in rate_constants:
            # Check if there's an '=' symbol, which indicates a rate constant with a value
            if '=' in rateString:
                name, value = rateString.split('=')
                names.append(name.strip())
                values.append(float(value.strip()))  # Convert the value to a float
            else:
                # If no '=' symbol, it's just a symbolic rate constant (e.g., "kon")
                names.append(rateString.strip())
                values.append(-1)  # Unrealistic value so we know this wasn't set to a number

        if reactionCount == 1:  # We add only 1 case, easy
            # Create a Reaction object and add to the list
            reaction_obj = Reaction(
                equation=equation,
                rateName=names[0],
                rateValue=values[0],
                reactants=reactants,
                reactant_coeffs=reactantCoeffs,
                products=products,
                product_coeffs=productCoeffs
            )
            parsed_reactions.append(reaction_obj)
        
        elif reactionCount == 2:  # We add 2 objects for bi-directional reactions
            reaction_fwd = Reaction(
                equation=formatFwdReaction(reactants, reactantCoeffs, products, productCoeffs),
                rateName=names[0],
                rateValue=values[0],
                reactants=reactants,
                reactant_coeffs=reactantCoeffs,
                products=products,
                product_coeffs=productCoeffs
            )
            
            reaction_rev = Reaction(
                equation=formatFwdReaction(products, productCoeffs, reactants, reactantCoeffs),
                rateName=names[1],
                rateValue=values[1],
                reactants=products,
                reactant_coeffs=productCoeffs,
                products=reactants,
                product_coeffs=reactantCoeffs
            )
            parsed_reactions.append(reaction_fwd)
            parsed_reactions.append(reaction_rev)
        else:
            print(f"\tError: Invalid equation format in reaction: {reaction}")
            continue

    return parsed_reactions


# Adjusted parseEquation to handle None case
def parseEquation(equation):
    # Determine if the equation has <-> or ->
    if '<->' in equation:
        reactionCount = 2
        lhs, rhs = equation.split('<->')
    elif '->' in equation:
        reactionCount = 1
        lhs, rhs = equation.split('->')
    else:
        print(f"Error: Equation must contain '->' or '<->': {equation}")
        return None, None, None, None, None

    # Split LHS and RHS on "+" and ignore white space
    reactants = [r.strip() for r in lhs.split('+')]
    products = [p.strip() for p in rhs.split('+')]

    # Call extract_coefficients to process reactants and products
    reactantCoeffs, reactants = extract_coefficients(reactants)
    productCoeffs, products = extract_coefficients(products)

    return reactionCount, reactants, reactantCoeffs, products, productCoeffs

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

#################
# Class Objects #
#################

class Reaction:
    def __init__(self, equation, rateName, rateValue, reactants, reactant_coeffs, products, product_coeffs):
        self.equation = equation
        self.rateName  = rateName
        self.rateValue = rateValue
        self.reactants = reactants
        self.reactant_coeffs = reactant_coeffs
        self.products = products
        self.product_coeffs = product_coeffs

    def __repr__(self):
        return (f"Reaction(equation={self.equation}, rateName={self.rateName}, "
                f"rateValue={self.rateValue}, "
                f"reactants={self.reactants}, "
                f"reactant_coeffs={self.reactant_coeffs}, products={self.products}, "
                f"product_coeffs={self.product_coeffs})")
                
    def __eq__(self, other):
        if isinstance(other, Reaction):
            sorted_self_reactants = sorted(zip(self.reactants, self.reactant_coeffs))
            sorted_other_reactants = sorted(zip(other.reactants, other.reactant_coeffs))
            sorted_self_products = sorted(zip(self.products, self.product_coeffs))
            sorted_other_products = sorted(zip(other.products, other.product_coeffs))
            
            return (self.rateName == other.rateName and
                    self.rateValue == other.rateValue and
                    sorted_self_reactants == sorted_other_reactants and
                    sorted_self_products == sorted_other_products)
        return False

    def __hash__(self):
        sorted_reactants = tuple(sorted(zip(self.reactants, self.reactant_coeffs)))
        sorted_products = tuple(sorted(zip(self.products, self.product_coeffs)))
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
    
        self.N = len(self.species)
        self.M = len(self.rates)
        
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

def create_matlab_multipleFileOutput(input_file: str, outputPrefix: str, s: Stoich, species: list, rates: list, uniqueRates: list, v: bool = False):

    Ns = len(s.species)
    Nr = len(s.rates)
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
    f.write(outputPrefix)
    f.write("(t_final,t_start)\n")
    f.write("% Solves a system of ODEs from t=t_start to t=t_final \n")
    f.write("% If no start time is given, then t_start = 0 \n")
    f.write("% If no start or final time is given, then t_start = 0, t_final = 1 \n")
    f.write("%\n")
    f.write("%\n")
    f.write("% This file was created by issuing command: \n")
    f.write("%     python createMatlabFile.py ")
    f.write(input_file)
    f.write("\n")
    f.write("%\n")
    f.write("\n")
    f.write("\nif nargin == 1\n")
    f.write("     t_start = 0;  % Default start time is 0 \n")
    f.write("elseif nargin == 0\n")
    f.write("     t_start = 0; % Default start time is 0\n")
    f.write("     t_final = 1; % Default final time is 1\n")
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
    f.write("[time,y] = ode15s(@(t,y)RHS(t,y,p), [t_start t_final], init_cond, options);\n")
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


    f.write("function dy = RHS(t,y,p)\n\n")
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

    f.write("% ODEs from reaction equations \n\n")

    if v: print("Writing ODEs now....\n")

    for i in range(Ns):
        f.write("% ")
        f.write(str(transform_string(species[i])))
        f.write("\n dy(")
        f.write(str(i+1))
        f.write(")  =")
        for j in range(Nr):
            #If the reaction j reduces the amount of species i;
            if (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) < 0):
                f.write("  -  ")
                f.write(str(rates[j]))
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
                f.write(str(rates[j]))
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
    

def create_matlab_output(input_file: str, output_file: str, s: Stoich, species: list, rates: list, uniqueRates: list, v: bool = False):

    Ns = len(s.species)
    Nr = len(s.rates)
   #species = [i.name for i in s.species.keys()]
   #rates = [i for i in s.rates.keys()]
    if v: print('\nOutput File: \n', output_file)
    if v: print("\n")
    #f = open(output_file,'w')
    f = io.StringIO() #Let's us format the code to wrap after 80 characters;

    f.write("function [time,y] = ")
    f.write(output_file.strip(".m"))
    f.write("(t_start,t_final)\n")
    f.write("% Solves a system of ODEs from t=t_start to t=t_final \n")
    f.write("% If no start time is given, then t_start = 0 \n")
    f.write("% If no start or final time is given, then t_start = 0, t_final = 1 \n")
    f.write("%\n")
    f.write("%\n")
    f.write("% This file was created by issuing command: \n")
    f.write("%     python createMatlabFile.py ")
    f.write(input_file)
    f.write("\n")
    f.write("%\n")
    f.write("\n")
    f.write("\nif nargin == 1\n")
    f.write("     t_start = 0;  % Default start time is 0 \n")
    f.write("elseif nargin == 0\n")
    f.write("     t_start = 0; % Default start time is 0\n")
    f.write("     t_final = 1; % Default final time is 1\n")
    f.write("end\n\n\n")

    f.write("% Kinetic Parameters \n")
    for i in range(len(uniqueRates)):
        f.write(uniqueRates[i])
        f.write(" = 1; \n")

    f.write("\np = [ ")
    f.write(uniqueRates[0])
    for i in range(1, len(uniqueRates)):
        f.write(", ")
        f.write(uniqueRates[i])
    f.write(" ];\n\n\n")

    
    f.write("% Initial Conditions \n")
    for i in range(Ns):
        f.write(transform_string(species[i]))
        f.write("_IC")
        f.write(" = 0; \n")

    ##Line that's too long;
    f.write("\ninit_cond = [ ")
    f.write(transform_string(species[0]))
    f.write("_IC")
    for i in range(1,Ns):
        f.write(", ")
        f.write(transform_string(species[i]))
        f.write("_IC")
    f.write(" ];\n\n\n")

    f.write("options = odeset('RelTol',1e-12,'AbsTol',1e-23);\n\n\n")

    f.write("%-------------------------------- Main Solve -----------------------------%\n")
    f.write("[time,y] = ode15s(@(t,y)RHS(t,y,p), [t_start t_final], init_cond, options);\n")
    f.write("%-------------------------------------------------------------------------%")
    f.write("\n\n\n")

    f.write("% Rename solution components\n") ##Modify
    for i in range(Ns):
        f.write(transform_string(species[i]))
        f.write(" = y(:,")
        f.write(str(i+1))
        f.write("); \n")

    f.write("\n\n\n")

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


    f.write("function dy = RHS(t,y,p)\n\n")
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


    f.write("\n\n")

    f.write("% ODEs from reaction equations \n\n")

    if v: print("Writing ODEs now....\n")

    for i in range(Ns):
        f.write("% ")
        f.write(str(transform_string(species[i])))
        f.write("\n dy(")
        f.write(str(i+1))
        f.write(")  =")
        for j in range(Nr):
            if (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) < 0):
                f.write("  -  ")
                f.write(str(rates[j]))
                for k in range(Ns):
                    if (not math.isnan(s.stoich[k][j])) and (int(s.stoich[k][j]) <= 0):
                        f.write(" * ")
                        f.write(transform_string(species[k]))
                        if (abs(int(s.stoich[k][j])) != 1) and (s.stoich[k][j] != 0):
                            f.write("^")
                            f.write(str(abs(int(s.stoich[k][j]))))
            elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) > 0):
                f.write("  +  ")
                f.write(str(rates[j]))
                for k in range(Ns):
                    if (not math.isnan(s.stoich[k][j])) and (int(s.stoich[k][j]) <= 0):
                        f.write(" * ")
                        if (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) > 1):
                            f.write(str(int(s.stoich[i][j])))
                            f.write(" * ")
                        f.write(transform_string(species[k]))
                        if (abs(int(s.stoich[k][j])) != 1) and (s.stoich[k][j] != 0):
                            f.write("^")
                            f.write(str(abs(int(s.stoich[k][j]))))
            elif (not math.isnan(s.stoich[i][j])) and (int(s.stoich[i][j]) == 0):
                f.write("  +  ")
                f.write(" 0 ")
        if v: print(species[i]," complete")
        f.write(";\n\n")


    f.write("\n\n\n\n")
    f.write("end")
    
    formatted_code = f.getvalue()
    f.close()
    
    formatted_code_with_continuations = add_line_continuations(formatted_code)

    # Write the formatted code to the file
    with open(output_file, 'w') as file:
        file.write(formatted_code_with_continuations)

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

verbose        = False
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
            
            if len(line) > 1 and '#' not in line:  # Check if line is not empty and does not contain '#'
                # Split the line by commas
                fields = line.split(',')
                
                # Check if there are more than one field (i.e., multiple comma-separated values)
                if len(fields) > 1:
                    # Process as a biochemistry equation (append the whole line to the list)
                    biochemicalReactions.append(line)
                    numReactionsReadIn += 1
                else:
                    # Handle the case where there is only one field
                    # For now, you can just append it separately or process as needed
                    initialConditions.append(line)
                    numInitialConditionsReadIn += 1
                    print(f"Single field line (not parsed as equation): {line}")
                    
except IOError:
    sys.exit(f"Couldn't open {sys.argv[1]}")

# Output the counts to the screen
if verbose:
    print(f"\nNumber of Biochemical Equations: {numReactionsReadIn}")
    print(f"\nNumber of Species with Initial Conditions: {numInitialConditionsReadIn}")

# Get the input file name
input_file = sys.argv[1]
    
# Extract the PREFIX from the input file name
prefix, _ = os.path.splitext(input_file)
    
# Define the output file names
stoich_output_file   = f"{prefix}Stoich.csv"
reactant_output_file = f"{prefix}Reactants.csv"
product_output_file  = f"{prefix}Products.csv"

#matlab_output_file = f"{prefix}Matlab.m"
#IC_output_file     = f"{prefix}ICs.m"
#param_output_file  = f"{prefix}Params.m"
#rename_output_file = f"{prefix}Rename.m"

if verbose:
    print(f"\tDONE. Successfully processed {sys.argv[1]}")
    print(f"\tNumber of biochemical reactions: {len(biochemicalReactions)}")
    print(f"\tSome of the Biochemical Reactions:")
    num_reactions_to_print = min(previewVal, len(biochemicalReactions))
    for i in range(num_reactions_to_print):
        print(f"\t\t" + biochemicalReactions[i])
    print('-' * 50)


###########################################
# Step 1: Parse the Biochemical Reactions #
###########################################

if verbose:
    print(f"Step 1: Parsing the Biochemical Reacitons")

# Parse the reactions
parsed_reactions = parseReactions(biochemicalReactions)

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

# Iterate through each reaction
for reaction in parsed_reactions:
    #All reactions are now (after parsing) unidirectional
    rates.append(reaction.rateName)
        
    # Append Reactants and Products if they are not empty
    if reaction.reactants:
        species.extend(reaction.reactants)
    if reaction.products:
        species.extend(reaction.products)

# Get unique entries
unique_species = unique_entries_only(species)
unique_rates   = unique_entries_only(rates)

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

unique_species.sort() #For ease in outputting; put variables in sorted order.

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

#Deprecated - we should remove this.
#create_matlab_output(sys.argv[1], matlab_output_file, s, unique_species, rates, unique_rates, verbose)

create_matlab_multipleFileOutput(sys.argv[1], prefix, s, unique_species, rates, unique_rates, verbose)

if verbose:
    print(f"DONE! Successfully created Matlab Files")
    print('-' * 50)

########################################
# Step 4: Determine Conserved Quantity #
########################################

s_count = 0
st_count = 0

sLipidSites = "Ls_Sites = ";
stLipidSites = "Lst_Sites = ";

for spec in unique_species:
    # Count occurrences of "_st" first
    local_st_count = spec.count("_st")
    
    # Replace "_st" with a placeholder to avoid counting it as "_s"
    temp_spec = spec.replace("_st", "")
    
    # Count occurrences of "_s" in the modified string
    local_s_count = temp_spec.count("_s")
    
    if local_st_count > 0:
        if local_st_count > 1:
            stLipidSites += f"+ {local_st_count} * {transform_string(spec)} "
        else:
            stLipidSites += f"+ {transform_string(spec)} "
    if local_s_count > 0:
        if local_s_count > 1:
            sLipidSites += f"+ {local_s_count} * {transform_string(spec)} "
        else:
            sLipidSites += f"+ {transform_string(spec)} "
    
    st_count += local_st_count
    s_count += local_s_count
        
    
    if verbose: print(f"{spec}\t{local_st_count}\t{local_s_count}")

if verbose:
    print(f'The string "_s" occurs {s_count} times (excluding "_st").')
    print(f'The string "_st" occurs {st_count} times.')

    print(f'The total bound lipid: {stLipidSites}')
    print(f'The total bound ubBoundlipid: {sLipidSites}')

##################################
# Step 4: Checking the Dimension #
##################################

bad_rates = s.check_rates()

if len(bad_rates) > 0:
    for rate, column_sum in bad_rates.items():
        if verbose: print(f"Rate {rate}: Sum of columns = {column_sum}")
else:
    if verbose: print("All rates are correct dimension")
