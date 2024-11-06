# BioMassODE
Converts a biochemical system of equations with the law of mass action to a system of ODEs.

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

New Features (08/19/2024):
    1) Allows for pure synthesis or pure degradation. ( -> A, B ->)
    2) Handles comments/white space in the biochemical equation file.
    3) Handles duplicate kinetic rates (i.e., p vector is consolidated)
    4) Handles A + B -> A + C reactions (splits stochiometric matrix)
    5) Checks for (and removes) duplicate reactions. (i.e., identical reaction)
        -> Even if the reaction is 1 side of a bi-directional reaction.
    6) Checks that the dimension of reactions rates is the SAME.

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
08/19/2024


Usage: python3 createMatlabFile.py StaticCoag.txt
