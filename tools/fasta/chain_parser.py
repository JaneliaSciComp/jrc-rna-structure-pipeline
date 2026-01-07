import re
from typing import Dict, List, Tuple


def parse_fasta(fasta_content: str) -> Dict[str, Tuple[str, List[str]]]:
    """
    Parse FASTA content into dictionary.

    Args:
        fasta_content: Multi-line FASTA string with format:
        >1A1T_1|Chain A[auth B]|SL3 STEM-LOOP RNA|
        or
        >104D_1|Chains A[auth A], B[auth B]|DNA/RNA (...)|

    Returns:
        Dictionary mapping auth chain_id to (sequence, list_of_auth_chain_ids)
        Example: {"A": ("ACGT", ["A", "B"]), "C": ("UGCA", ["C"])}
        The key is the auth chain ID, and the list contains all auth chain IDs for this sequence
    """
    result = {}
    lines = fasta_content.strip().split("\n")

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith(">"):
            # Parse new format header: >104D_1|Chains A[auth A], B[auth B]|...| or >1A1T_1|Chain A[auth B]|...|
            # Extract the chains part (between first | and second |)
            parts = line.split("|")
            if len(parts) < 2:
                print("Warning: Malformed FASTA header:", line)
                auth_chain_ids = []
                chains_part = ""
            else:
                chains_part = parts[1].strip()

                # Extract auth chain IDs from patterns like "Chain A[auth B]" or "Chains A[auth A], B[auth B] or just "Chain A" or "Chains A, B"
                auth_chain_ids = []
                replaced_chains_part = re.sub(r"^Chains? ", "", chains_part)
                chains = replaced_chains_part.split(",")
                for chain in chains:
                    auth_match = re.search(r"\[auth ([^\]]+)\]", chain)
                    if auth_match:
                        auth_chain_ids.append(auth_match.group(1).strip())
                    else:
                        c = chain.strip()
                        if c:
                            auth_chain_ids.append(c)

            # # Find all Y[auth X] or just Y patterns
            # auth_matches = re.findall(r"[^ ,]+(?:\[auth ([^\]]+)\])?", chains_part)

            # auth_chain_ids = [match.strip() for match in auth_matches]

            if not auth_chain_ids:
                print("Warning: Empty chains part:", chains_part)
                primary_auth_chain = None
            else:
                # Use the first auth chain ID as the key
                primary_auth_chain = auth_chain_ids[0]

            # Read sequence (next lines until next header or end)
            sequence = ""
            while (i + 1) < len(lines) and lines[i + 1].startswith(">") is False:
                sequence += lines[i + 1].strip()
                i += 1
            for chain in auth_chain_ids:
                result[chain] = (sequence, auth_chain_ids)

        i += 1

    return result
