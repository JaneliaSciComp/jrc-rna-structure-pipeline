from Bio.PDB import MMCIFParser, MMCIFIO, PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import os
from glob import glob
import numpy as np
from rna_reference import modified_to_unmodified

def split_cif_by_chains(input_cif, output_dir):
    """
    Splits a CIF file containing multiple chains into separate CIF files, saving only RNA chains,
    and writes the sequences of each chain to a FASTA file.

    Parameters:
        input_cif (str): Path to the input CIF file.
        output_dir (str): Directory to save the individual chain CIF files and FASTA files.

    Returns:
        list: Paths to the generated CIF and FASTA files.
    """


    for key in modified_to_unmodified:
        modified_to_unmodified[key] = modified_to_unmodified[key].strip()
    
    # Initialize the parser and parse the CIF file
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", input_cif)

    output_files = []
    ID = input_cif.name.split(".")[0]  # to take care of .cif.gz too

    for model in structure:
        for chain in model:
            # Replace modified nucleotides with their unmodified counterparts
            to_delete=[]
            for residue in chain.get_residues():
                if residue.resname in modified_to_unmodified:
                    residue.resname = modified_to_unmodified[residue.resname]

                #print(residue.type)

                # if residue.id[0] != " ":
                #     print(residue.resname)
                #exclude heteroatoms and non canonical residues
                # FIXME: That will leave fragments, make sure we want that
                if residue.resname not in {"A", "U", "G", "C"} or residue.id[0] != " ":
                    to_delete.append(residue.id)
            #exit()
            # print(chain.child_dict.keys())
            # exit()
            #print(to_delete)
            # Remove unwanted residues
            for res_id in to_delete:
                #print(chain.child_dict[res_id])
                chain.detach_child(res_id)
                #del chain.child_dict[res_id]

            #print(list(chain.child_dict.keys()))
            # FIXME: configurable minimum length 
            # continue if the chain has at least 10 residues left
            if len(chain.child_dict) < 10:
                continue

            # Check if 0.5 or more of the residues are RNA
            # FIXME: configurable threshold
            # TODO:idealy that should depend on the interactions, ie. it is fine to keep 
            #      RNA from protein complex if the RNA is interacting in limited fashion with the protein
            #      alternatively assign score for further filtering
            is_rna = np.array([residue.resname in ["A", "C", "G", "U"] for residue in chain.get_residues()]).mean() >= 0.5
            # print([residue.resname for residue in chain.get_residues()])  
            # print(chain)
            # print(is_rna)
            if is_rna:
                chain_id = chain.id

                # Create a new structure object with only this chain
                chain_structure = structure.__class__("{}_chain_{}".format(structure.id, chain_id))
                model_copy = model.__class__(model.id)
                model_copy.add(chain.copy())
                chain_structure.add(model_copy)

                # Write the chain to a new CIF file
                cif_output_path = os.path.join(output_dir, f"{ID}_{chain_id}.cif")
                io = MMCIFIO()
                io.set_structure(chain_structure)
                io.save(cif_output_path)
                output_files.append(cif_output_path)

                # Write the chain to a new PDB file
                try:
                    pdb_output_path = os.path.join(output_dir, f"{ID}_{chain_id}.pdb")
                    pdb_io = PDBIO()
                    pdb_io.set_structure(chain_structure)
                    pdb_io.save(pdb_output_path)
                    output_files.append(pdb_output_path)
                except:
                    pass

                # Extract the sequence and write to a FASTA file
                sequence = "".join(residue.resname for residue in chain.get_residues())
                fasta_output_path = os.path.join(output_dir, f"{ID}_{chain_id}.fasta")
                seq_record = SeqRecord(Seq(sequence), id=f"{ID}_{chain_id}", description=f"Chain {chain_id}")
                with open(fasta_output_path, "w") as fasta_file:
                    write(seq_record, fasta_file, "fasta")
                output_files.append(fasta_output_path)

    return output_files

if __name__ == "__main__":
    import argparse
    from functools import partial
    from pathlib import Path
    from utils import parallel_process

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir", type=str, help="Path to the directory with mmCIF files to check"
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Path to the directory where filtered files will be saved (symlinked)",
    )
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cif_files = list(input_dir.glob("*.cif"))
    # Use multiprocessing to process CIF files in parallel
    parallel_process(
        partial(split_cif_by_chains, output_dir=output_dir),
        cif_files,
        desc="Splitting CIF files by chains",
    )