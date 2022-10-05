from pathlib import Path
import glob

list_all = glob.glob('*.txt')

for filename in list_all:
    # creating a variable and storing the text that we want to search
    # search_text = "_"

    # # creating a variable and storing the text that we want to add
    # replace_text = "-"

    # # Opening our text file in read only mode using the open() function
    # with open(filename, 'r') as file:

    # 	# Reading the content of the file
    # 	# using the read() function and storing
    # 	# them in a new variable
    # 	data = file.read()

    # 	# Searching and replacing the text
    # 	# using the replace() function
    # 	data = data.replace(search_text, replace_text)

    # # Opening our text file in write only mode to write the replaced content
    # with open(filename, 'w') as file:

    # 	# Writing the replaced data in our
    # 	# text file
    # 	file.write(data)

    f = Path(filename)
    f.rename(f.with_suffix('.node'))
