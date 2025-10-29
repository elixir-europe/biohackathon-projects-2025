import sys
import os
import re
#tree_directory = sys.argv[1]
alignment_directory =  sys.argv[1]

#tree_list = os.listdir(tree_directory)
alignment_list = os.listdir(alignment_directory)

for file in alignment_list:
    file = os.path.join(alignment_directory,file)
    alignment = open(file,"r")
    alignment_out_file = file.replace(".faa.filtered.mafft.clipkit","_modified_names.faa.filtered.mafft.clipkit")
    alignment_out = open(alignment_out_file,"w")

    for line in alignment:
        if line.startswith(">"):
            final_line = line.rstrip()
            final_line = final_line.replace(":","_")
            if final_line.endswith("+"):
                final_line = final_line[:-1] + "_"

            #            first_half = line.split('|')[0]
#            splitted_half = first_half.split("_")
#            filtered_list = [x for x in splitted_half if not any(char.isdigit() for char in str(x))]
#            final_line = "_".join(filtered_list)
            alignment_out.write(final_line)
            alignment_out.write("\n")
        else:
            alignment_out.write(line)
    alignment.close()
    alignment_out.close()

'''
for file in tree_list:
    file = os.path.join(tree_directory, file)
    tree = open(file,"r")
    tree_out_file = file.replace(".faa.mafft.clipkit.treefile","_modified_names.faa.mafft.clipkit.treefile")
    tree_out = open(tree_out_file,"w")
    for line in tree:
        def replacer(match):
            full = match.group(0)
            species = re.match(r'([^\(\),:|]+_[^\(\),:|]+)', full).group(1)
            return species

        # Replace each taxon with the cleaned species name
        cleaned_tree = re.sub(r'([^\(\),:|]+(?:_[^\(\),:|]+)*?)_\d+at\d+[^\(\),:|]*\|[^\(\),:|]*', replacer,
                              line)
        tree_out.write(cleaned_tree)
    tree.close()
    tree_out.close()
'''


