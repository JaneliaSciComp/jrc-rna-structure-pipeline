rna_atom_groups = {
    "A": {
        "all": ["P", "OP1", "OP2", "O5'","O3'"]+\
               ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"]+\
               ["N9", "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N6"]
    },
    "U": {
        "all": ["P", "OP1", "OP2", "O5'","O3'"]+\
               ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"]+\
               ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
    },
    "G": {
        "all": ["P", "OP1", "OP2", "O5'","O3'"]+\
                     ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"]+\
                     ["N9", "N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", ]
    },
    "C": {
        "all": ["P", "OP1", "OP2", "O5'","O3'"]+\
                     ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"]+\
                     ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
    },
    "N": {
        "all": ["P", "OP1", "OP2", "O5'","O3'"]+\
                     ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"]+\
                     ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
    }
}
atom_names=set()

for nt in rna_atom_groups.values():
    for atom_group in nt.values():
        for atom in atom_group:
            atom_names.add(atom)
atom_names=list(atom_names)
atom_names.sort()
print("atom_names=", atom_names)
print("number of atom names:", len(atom_names))


