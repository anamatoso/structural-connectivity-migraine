from pathlib import Path
import glob

# Select all text files
list_all = glob.glob('*.txt')

# Replace suffic .txt with .node
for filename in list_all:
    f = Path(filename)
    f.rename(f.with_suffix('.node'))
