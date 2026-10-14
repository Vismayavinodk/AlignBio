from pathlib import Path

folder = Path("data/raw/bench1.0/sabrem/ref")
files = sorted(path for path in folder.iterdir() if path.is_file())

total_sequences = 0
alignment_lengths = []
sequence_counts = []
invalid_files = []

for file in files:
    try:
        lines = file.read_text(errors="replace").splitlines()
    except OSError:
        invalid_files.append(file.name)
        continue

    sequences = []
    current_sequence = ""
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_sequence:
                sequences.append(current_sequence)
            current_sequence = ""
        else:
            current_sequence += line
    if current_sequence:
        sequences.append(current_sequence)

    if not sequences:
        invalid_files.append(file.name)
        continue

    lengths = [len(seq) for seq in sequences]
    # A valid multiple sequence alignment has equal row lengths.
    if len(set(lengths)) != 1:
        invalid_files.append(file.name)
        continue

    total_sequences += len(sequences)
    sequence_counts.append(len(sequences))
    alignment_lengths.append(lengths[0])

print(f"Total files: {len(files)}")
print(f"Valid alignment files: {len(sequence_counts)}")
print(f"Invalid files: {len(invalid_files)}")
if invalid_files:
    print(f"List of invalid files: {invalid_files}")
print(f"Total sequences: {total_sequences}")
if sequence_counts:
    print(f"Sequences per file: min={min(sequence_counts)}, max={max(sequence_counts)}, avg={sum(sequence_counts)/len(sequence_counts):.2f}")
if alignment_lengths:
    print(f"Alignment lengths: min={min(alignment_lengths)}, max={max(alignment_lengths)}, avg={sum(alignment_lengths)/len(alignment_lengths):.2f}")
