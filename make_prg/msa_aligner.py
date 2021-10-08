from abc import ABC, abstractmethod
import uuid
import shutil
from pathlib import Path
from typing import List
import time
import subprocess
from loguru import logger
import os
from Bio.AlignIO import MultipleSeqAlignment
from Bio import SeqIO
from make_prg.io_utils import load_alignment_file


class TempDirAlreadyExistsError(Exception):
    pass


class NotAValidExecutableError(Exception):
    pass


class MSAAligner(ABC):
    def _set_executable(self, executable: str) -> None:
        is_valid_executable = shutil.which(executable) is not None
        if not is_valid_executable:
            raise NotAValidExecutableError(f"Given MSA executable {executable} does not work or is invalid")
        self._executable = executable

    def _set_tmpdir(self, tmpdir: Path) -> None:
        tmpdir.mkdir(parents=True, exist_ok=True)
        self._tmpdir = tmpdir

    def __init__(self, executable: str, tmpdir: Path = Path(".")):
        self._set_executable(executable)
        self._set_tmpdir(tmpdir)

    @abstractmethod
    def get_updated_alignment(self, current_alignment: MultipleSeqAlignment, new_sequences: List[str]) -> MultipleSeqAlignment:
        pass

    @abstractmethod
    @classmethod
    def get_aligner_name(cls):
        pass

    def _run_aligner(self, args: str, env=None):
        start = time.time()
        process = subprocess.Popen(
            args,
            stderr=subprocess.PIPE,
            encoding="utf-8",
            shell=True,
            env=env,
        )
        exit_code = process.wait()

        if exit_code != 0:
            raise RuntimeError(
                f"Failed to execute {self.__class__.get_aligner_name()} for arguments {args} due to the following "
                f"error:\n{process.stderr.read()}"
            )
        stop = time.time()
        runtime = stop - start
        logger.debug(f"{self.__class__.get_aligner_name()} runtime for arguments {args} in seconds: {runtime:.3f}")

    @staticmethod
    def _create_new_sequences_file(directory: Path, new_sequences: List[str]) -> Path:
        new_sequences_filepath = directory / f"new_sequences.fa"
        with open(new_sequences_filepath, "w") as new_sequences_handler:
            for index_new_seq, new_seq in enumerate(new_sequences):
                print(
                    f">Denovo_path_{index_new_seq}_random_id_{uuid.uuid4()}",
                    file=new_sequences_handler,
                )
                print(new_seq, file=new_sequences_handler)
        return new_sequences_filepath


class MAFFT(MSAAligner):
    def __init__(self, executable: str, tmpdir: Path = Path(".")):
        super().__init__(executable, tmpdir)

    @classmethod
    def get_aligner_name(cls):
        return "MAFFT"

    def _prepare_run_tmpdir(self) -> Path:
        run_tmpdir = self._tmpdir / f"temp.MSA.{uuid.uuid4()}"
        if run_tmpdir.exists():
            raise TempDirAlreadyExistsError(f"{str(run_tmpdir)} already exists, cannot be created, "
                                            f"this is very unlikely.")
        run_tmpdir.mkdir(parents=True)
        return run_tmpdir

    def _cleanup_run(self, run_tmpdir) -> None:
        shutil.rmtree(run_tmpdir)

    def get_updated_alignment(self, current_alignment: MultipleSeqAlignment, new_sequences: List[str]) -> MultipleSeqAlignment:
        # setup
        run_tmpdir = self._prepare_run_tmpdir()

        current_msa_filepath = run_tmpdir / "previous_msa.fa"
        with open(current_msa_filepath, "w") as current_msa_handler:
            SeqIO.write(current_alignment, current_msa_handler, "fasta")

        new_sequences_filepath = self._create_new_sequences_file(run_tmpdir, new_sequences)
        new_msa = run_tmpdir / f"updated_msa.fa"

        # run
        args = " ".join(
            [
                self._executable,
                "--auto",
                "--quiet",
                "--thread",
                "1",
                "--add",
                str(new_sequences_filepath),
                str(current_msa_filepath),
                ">",
                str(new_msa),
            ]
        )
        env = os.environ
        env["TMPDIR"] = str(run_tmpdir)
        self._run_aligner(args, env)

        # load the updated alignment
        updated_alignment = load_alignment_file(str(new_msa), "fasta")

        # clean up
        self._cleanup_run(run_tmpdir)

        return updated_alignment

