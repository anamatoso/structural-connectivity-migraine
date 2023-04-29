from pathlib import Path
import glob
import sys

# Select all text files
list_all = glob.glob(sys.argv[1]+'/*.txt')

# Replace suffic .txt with .node
for filename in list_all:
    f = Path(filename)
    f.rename(f.with_suffix('.node'))
