from ete3 import Tree
import matplotlib.pyplot as plt

def count_unique_tips(tree_line):
    t = Tree(tree_line)
    unique_tips = set(tip.rsplit('_', 1)[0] for tip in t.get_leaf_names())
    return len(unique_tips)

def process_trees(input_file, output_file, save_path):
    values = []
    with open(input_file, 'r') as infile, open(output_file + "/num_species_gt.txt", 'w') as outfile:
        for line in infile:
            count = count_unique_tips(line)
            values.append(count)
            outfile.write(f"{count}\n")

    x_values = range(1, len(values) + 1)
    
    plt.figure(figsize=(10, 6))
    plt.bar(x_values, values)
    plt.title('Number of species in gene trees')
    plt.xlabel('Gene trees')
    plt.ylabel('Number of species')
    plt.xticks(x_values)  
    plt.tight_layout()
    
    # Save the plot to a file
    plt.savefig(save_path + "/num_species_gt.png")
    plt.close() 

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
save_path = sys.argv[3]

process_trees(input_file_path, output_file_path, save_path)
