# Function to extract gene IDs and their corresponding gene trees from the file
def extract_gene_info(file_name, interval_limit, total_genes):
    gene_info = {}
    with open(file_name, "r") as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) <= 2:
                continue  # Skip lines with only gene IDs
            gene_id = int(parts[0].strip())
            gene_tree = ",".join(parts[1:]).strip()
            gene_info[gene_id] = gene_tree

    x = gene_info.keys()
    y = gene_info.values()

    for i in range(interval_limit, total_genes + 1, interval_limit):
        gene_count = []
        for j in range(len(x)):
            if int(list(gene_info.keys())[j]) <= i:
                gene_count.append(list(gene_info.values())[j])
        with open(f"gene_trees_{i}.nwk", "w") as f:
            for item in gene_count:
                f.write(str(item) + "\n")


# Replace 'input_file.txt' with the name of your input file
# Replace the '50' in the function call with the user-provided interval limit
interval_limit = 50  # int(input("Enter the interval limit: "))
total_genes = 10000
extract_gene_info(
    "../converge_10kgenes_balancedmode/run_00/results/genetrees/original_list.txt", interval_limit, total_genes
)
