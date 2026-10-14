from pathlib import Path
from typing import Dict, List, Tuple, Union


def read_fasta(filepath: Union[str, Path]) -> List[Tuple[str, str]]:
    """Reads a FASTA file and returns a list of (header, sequence) tuples."""
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records = []
    current_header = None
    current_seq_chunks = []

    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_seq_chunks)))
            current_header = line[1:].strip()
            current_seq_chunks = []
        else:
            current_seq_chunks.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_seq_chunks)))

    return records


def write_fasta(records: List[Tuple[str, str]], filepath: Union[str, Path], line_width: int = 60) -> None:
    """Writes a list of (header, sequence) tuples to a FASTA file."""
    path = Path(filepath)
    lines = []
    for header, seq in records:
        lines.append(f">{header}")
        for i in range(0, len(seq), line_width):
            lines.append(seq[i : i + line_width])
    path.write_text("\n".join(lines) + "\n")


def read_msa(filepath: Union[str, Path]) -> Dict[str, str]:
    """Reads a multiple sequence alignment file, preserving aligned columns and gaps."""
    records = read_fasta(filepath)
    if not records:
        return {}

    lengths = [len(seq) for _, seq in records]
    if len(set(lengths)) > 1:
        raise ValueError(f"Inconsistent sequence lengths in MSA file {filepath}: {lengths}")

    return {header: seq for header, seq in records}
