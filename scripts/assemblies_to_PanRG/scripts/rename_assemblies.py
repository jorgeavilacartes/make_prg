input_dir="assemblies"
output_dir="assemblies_with_filenames_fixed"

import os
import shutil

os.makedirs(output_dir)
for index, file in enumerate(os.listdir(input_dir)):
    print(f"{input_dir}/{file} -> {output_dir}/{index}.fa")
    shutil.copy(f"{input_dir}/{file}", f"{output_dir}/{index}.fa")
