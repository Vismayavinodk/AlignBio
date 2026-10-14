from pathlib import Path

folder = Path("data/raw/bench1.0/sabrem/ref")

files = sorted(folder.iterdir())

print("Number of files:", len(files))

for file in files[:3]:
    if file.is_file():
        print(f"\n--- {file.name} ---")
        print(file.read_text(errors="replace")[:1500])
