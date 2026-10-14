from pathlib import Path

data_folder = Path("data/raw")
print("Files found in the dataset folder:\n")
for path in sorted(data_folder.rglob("*")):
    if path.is_file():
        print(f"{path}  ({path.stat().st_size:,} bytes)")
