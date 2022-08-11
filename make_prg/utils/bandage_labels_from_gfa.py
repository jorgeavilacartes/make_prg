import sys
from pathlib import Path
from collections import namedtuple

path_gfa = str(sys.argv[1])
assert Path(path_gfa).suffix == ".gfa", "Input must be a gfa file"
def bandage_labels(path_gfa):
    """
    Writes a CSV file with labels for Bandage
    """
    BandageLabels = namedtuple("BandageLables",["node_name","first_base","last_base", "seq"])
    info_labels = []
    
    with open(path_gfa) as fp:
        for line in fp:
            if line.startswith("S"):
                split_line  = line.split("\t")
                node_name = split_line[1]
                seq = split_line[2]
                first_base=seq[0]
                last_base=seq[-1]
                info_labels.append(BandageLabels(node_name, first_base, last_base, seq))
            
    # Write info in csv file
    colnames = ",".join(BandageLabels._fields)+"\n"

    filename_no_extension = Path(path_gfa).with_suffix("")

    with open(f"{filename_no_extension}_bandage_labels.csv","w") as fp:
        fp.write(colnames)
        for label in info_labels:
            fp.write(f"{label.node_name},{label.first_base},{label.last_base},{label.seq}\n")

bandage_labels(path_gfa)