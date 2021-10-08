from abc import ABC, abstractmethod
import uuid
import shutil
from pathlib import Path
from typing import List
import time
import subprocess
from loguru import logger
import os


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

    def _set_tmpdir(self, temp_prefix: Path) -> None:
        tmpdir = temp_prefix / f"temp.MSA.{uuid.uuid4()}"
        if tmpdir.exists():
            raise TempDirAlreadyExistsError(f"{str(tmpdir)} already exists, cannot be created.")
        tmpdir.mkdir(parents=True)
        self._tmpdir = tmpdir

    def __init__(self, executable: str, temp_prefix: Path = Path("."), env=None):
        self._set_executable(executable)
        self._set_tmpdir(temp_prefix)
        self._env = env

    @abstractmethod
    def get_updated_alignment(self, previous_alignment: Path, new_sequences: List[str]) -> Path:
        pass

    def _create_new_sequences_file(self, new_sequences: List[str]) -> Path:
        new_sequences_filepath = self._tmpdir / f"new_sequences_{uuid.uuid4()}.fa"
        with open(new_sequences_filepath, "w") as new_sequences_handler:
            for index_new_seq, new_seq in enumerate(new_sequences):
                print(
                    f">Denovo_path_{index_new_seq}_random_id_{uuid.uuid4()}",
                    file=new_sequences_handler,
                )
                print(new_seq, file=new_sequences_handler)
        return new_sequences_filepath

    @abstractmethod
    @classmethod
    def get_aligner_name(cls):
        pass

    def run_aligner(self, args: str, remove_temp_dir_after_running: bool = True):
        start = time.time()
        process = subprocess.Popen(
            args,
            stderr=subprocess.PIPE,
            encoding="utf-8",
            shell=True,
            env=self._env,
        )
        exit_code = process.wait()

        if remove_temp_dir_after_running:
            shutil.rmtree(self._tmpdir)

        if exit_code != 0:
            raise RuntimeError(
                f"Failed to execute {self.__class__.get_aligner_name()} for arguments {args} due to the following "
                f"error:\n{process.stderr.read()}"
            )
        stop = time.time()
        runtime = stop - start
        logger.debug(f"{self.__class__.get_aligner_name()} runtime for arguments {args} in seconds: {runtime:.3f}")


class MAFFT(MSAAligner):
    def __init__(self, executable: str, temp_prefix: Path = Path(".")):
        env = os.environ
        super().__init__(executable, temp_prefix, env)

    @classmethod
    def get_aligner_name(cls):
        return "MAFFT"

    def get_updated_alignment(self, previous_alignment: Path, new_sequences: List[str]) -> Path:
        new_sequences_filepath = self._create_new_sequences_file(new_sequences)
        new_msa = self._tmpdir / f"updated_msa_{uuid.uuid4()}.fa"

        args = " ".join(
            [
                self._executable,
                "--auto",
                "--quiet",
                "--thread",
                "1",
                "--add",
                str(new_sequences_filepath),
                str(previous_alignment),
                ">",
                str(new_msa),
            ]
        )
        self._env["TMPDIR"] = str(self._tmpdir)
        self.run_aligner(args, remove_temp_dir_after_running=True)

        return new_msa

