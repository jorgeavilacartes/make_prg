from pathlib import Path

def output_files_already_exist(output_prefix):
    return Path(output_prefix + "_tmp").exists() or \
           Path(output_prefix + ".prg.fa").exists() or \
           Path(output_prefix + ".log").exists() or \
           Path(output_prefix + ".pickle").exists()
