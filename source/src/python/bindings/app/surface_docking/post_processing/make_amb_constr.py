#!/usr/bin/python2.6
# Script to help make ambiguous constraint files
# See https://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/file_constraints.html for guideline
# Emily Koo

def get_atom2_num(pdb_file, atom2):
    i = 0
    for line in pdb_file:
        if line.split()[2] == atom2:
            i += 1
    return i

def get_func_def(func_type):
    func_type_l = func_type.lower()
    if func_type_l == "bounded":
        lb = raw_input("Enter lower bound: ")
        ub = raw_input("Enter upper bound: ")
        sd = raw_input("Enter standard deviation: ")
        rswitch = raw_input("Enter rswitch: ")
        tag = raw_input("Enter tag: ")
        return (lb, ub, sd, rswitch, tag)    
    elif func_type_l == "harmonic" or func_type_l == "circular harmonic":
        x = raw_input("Enter x0: ")
        sd = raw_input("Enter standard deviation: ")
        return (x, sd)

def format_constr_type(constr_type):
    constr_type_l = constr_type.lower()
    if constr_type_l == "atompair":
        constr_type_format = "AtomPair"
    elif constr_type_l == "dihedral":
        constr_type_format = "Dihedral"

    return constr_type_format

def main():
    input_file = raw_input("Enter name of PDB file with ext: ")
    pdb_file = open(input_file, 'r')

    constr_type_raw = raw_input("Enter constraint type: ")
    constr_type = format_constr_type(constr_type_raw)

    atom1_name = raw_input("Enter name of atom 1 (e.g. C, CZ, CD): ")
    atom1_res = raw_input("Enter residue number of atom 1: ")

    atom2_name = raw_input("Enter name of atom 2 (e.g. C, CZ, CD): ")
    atom2_num = get_atom2_num(pdb_file, atom2_name.upper())
    func_type = raw_input("Enter function type: ")
    func_def = get_func_def(func_type)
    
    output_file = raw_input("Enter desired name of output file without ext: ")
    
    cst_file = open(output_file + ".cst", 'a')
    
    for n in range(1, atom2_num + 1):
        cst_file.write("AmbiguousConstraint " + constr_type + " " + atom1_name.upper() + " " + atom1_res + " " + atom2_name.upper() + " " + str(n) + " " + func_type.upper() + " ")
        for m in func_def:
            cst_file.write(m + " ")
        if n != atom2_num:
            cst_file.write("\n")

main()
