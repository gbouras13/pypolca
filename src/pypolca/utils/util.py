"""MISC FUNCTIONS
You shouldn't need to tweak these much if at all
"""

import os
import shutil
import subprocess as sp
import sys
import time
from pathlib import Path

import click
from loguru import logger


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


def pypolca_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(pypolca_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def print_splash():
    click.echo(
        """\b
                         _           
 _ __  _   _ _ __   ___ | | ___ __ _ 
| '_ \| | | | '_ \ / _ \| |/ __/ _` |
| |_) | |_| | |_) | (_) | | (_| (_| |
| .__/ \__, | .__/ \___/|_|\___\__,_|
|_|    |___/|_| 

"""
    )


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as lo:
            lo.write(msg)


def print_citation():
    with open(pypolca_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

"""
begin and end functions
"""


def begin_pypolca(params) -> None:
    """
    begins pypolca
    params: params a dictionary of params
    returns start time
    """
    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(params["--output"], f"pypolca_{start_time}.log")
    # adds log file
    logger.add(log_file)
    logger.add(lambda _: sys.exit(1), level="ERROR")
    print_splash()
    logger.info(
        f"pypolca: Standalone Python implementation of the POLCA polisher from MaSuRCA"
    )
    logger.info(f"You are using pypolca version {get_version()}")
    logger.info("Repository homepage is https://github.com/gbouras13/pypolca.")
    logger.info(
        "Written by George Bouras: george.bouras@adelaide.edu.au adapting the original POLCA code by Aleksey Zimin."
    )

    logger.info(f"Listing input parameters.")
    for key, value in params.items():
        logger.info(f"Parameter: {key} {value}.")

    return start_time


def end_pypolca(start_time) -> None:
    """
    finishes pypolca
    """

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("pypolca has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")


def copy_file(source_path: Path, destination_path: Path) -> None:
    """
    Copy a file from the source path to the destination path.

    Args:
        source_path (Path): The path to the source file to be copied.
        destination_path (Path): The path to the destination where the file will be copied to.

    Raises:
        FileNotFoundError: If the source file does not exist.
        IsADirectoryError: If the source path is a directory.
        PermissionError: If the user does not have permission to copy the file.
        shutil.Error: If an error occurs during the copy operation.

    Returns:
        None
    """
    # Check if the source file exists and is not a directory
    if not source_path.is_file():
        raise FileNotFoundError(f"Source file '{source_path}' not found.")

    # Use shutil.copy to perform the file copy
    try:
        shutil.copyfile(source_path, destination_path)
    except IsADirectoryError:
        raise IsADirectoryError(
            f"Source path '{source_path}' is a directory, not a file."
        )
    except PermissionError:
        raise PermissionError(
            f"You do not have permission to copy '{source_path}' to '{destination_path}'."
        )
    except shutil.Error as e:
        raise shutil.Error(f"An error occurred while copying the file: {e}")


def remove_directory(dir_path: Path) -> None:
    """
    Remove a directory and its contents if it exists.

    Args:
        dir_path (Path): The Path object representing the directory to be removed.

    Returns:
        None: This function does not return a value.

    Raises:
        FileNotFoundError: If the specified directory does not exist.
        NotADirectoryError: If the specified path is not a directory.

    Example:
        remove_directory(Path('/path/to/directory'))
    """
    if dir_path.exists() and dir_path.is_dir():
        shutil.rmtree(dir_path)
