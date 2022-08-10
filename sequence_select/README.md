Use following code to generate seeding of the file
`python3 seeding.py -k [number of regions] --output [output].csv --path [input path directory]`
Then, Use following code to generate output selected regions fasta file. 
`./get_seq.sh [input seeding file].csv [input pathdirectory] [output file].fa [length of one region]`

If filename or length or number of regions is not specified, it is random choices or default sample file. 

Output is at out.fasta