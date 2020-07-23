#!/usr/bin/env python3
#
# Copyright 2018-2019 PSB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""helper_general.py

This module contains functions which are useful for a multitude of programs,
and which are hard to be categorized.
"""

# IMPORTS
# External modules
import json
import os
import pickle
import sys
from typing import Any, Dict, List


# CONSTANT SECTION
# As no useful UniProt API access libraries were found, the programs use the UniProt URL API directly.
# UNIPROT_UPLOADLISTS_URL represents the UniProt ID conversion search URL. The paremeters for the
# conversion are appendet after this string.
UNIPROT_UPLOADLISTS_URL = "https://www.uniprot.org/uploadlists"
# UNIPROT_URL represents the UniProtKB database search URL. After this string, the parameters are appended
# (e.g. a UniProt ID to refer to the protein's entry).
UNIPROT_URL = "https://www.uniprot.org/uniprot"


# PUBLIC FUNCTIONS SECTION
def check_argument(argument: str, command: str) -> str:
    """Checks wheter the given mandatory argument variable is None or not. If none, the program stops.

    It is checked wheter the argument is None as argparser sets an argument variable to None
    if the argument is not set.

    Arguments
    ----------
    * argument: str ~ The argument that shall be checked and which shall not be None.
    * command: str ~ The argument long name (i.e. the name after --) which is displayed
      in an error message if the argument is None.
    """
    if argument is None:
        print(f"Mandatory argument --{command} not given. Stopping program.")

        # Exit program with error code (not 0).
        sys.exit(-1)

    return argument


def ensure_folder_existence(folder: str) -> None:
    """Checks if the given folder exists. If not, the folder is created.

    Argument
    ----------
    * folder: str ~ The folder whose existence shall be enforced.
    """
    if os.path.isdir(folder):
        return
    os.makedirs(folder)


def get_files(path: str) -> List[str]:
    """Returns the names of the files in the given folder as a list of strings.

    Arguments
    ----------
    * path: str ~ The path to the folder of which the file names shall be returned
    """
    files: List[str] = []
    for (_, _, filenames) in os.walk(path):
        files.extend(filenames)
    return files


def get_float_cell_value(cell_value) -> float:
    """Returns the value of an openpyxl cell value.

    This function is used in oder to operate with spreadsheets
    which use the comma as the decimal mark instead of a period.

    Arguments
    ----------
    * cell_value ~ The openpyxl cell value
    """
    if type(cell_value) not in (int, float):
        cell_value = cell_value.replace(",", ".")
        cell_value = float(cell_value)
    return cell_value


def is_fitting_ec_numbers(ec_number_one: str, ec_number_two: str, wildcard_level: int) -> bool:
    """Check whether the EC numbers are the same under the used wildcard level.

    Arguments
    ----------
    * ec_number_one: str ~ The first given EC number.
    * ec_number_two: str ~ The second given EC number.
    * wildcard_level: int ~ The wildcard level.
    """
    if wildcard_level == 0:
        ec_number_one_full_numbers = ec_number_one.split(".")
        ec_number_two_full_numbers = ec_number_two.split(".")
    else:
        ec_number_one_full_numbers = ec_number_one.split(".")[:-wildcard_level]
        ec_number_two_full_numbers = ec_number_two.split(".")[:-wildcard_level]

    if ec_number_one_full_numbers == ec_number_two_full_numbers:
        return True
    else:
        return False


def json_load(path: str) -> Dict[Any, Any]:
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary


def json_write(path: str, dictionary: Dict[Any, Any]) -> None:
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written
    * dictionary: Dict[Any, Any] ~ The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)


def pickle_load(path: str) -> Any:
    """Returns the value of the given pickle file.

    Arguments
    ----------
    * path: str ~ The path to the pickle file.
    """
    pickle_file = open(path, 'rb')
    pickled_object = pickle.load(pickle_file)
    pickle_file.close()
    return pickled_object


def pickle_write(path: str, pickled_object: Any) -> None:
    """Writes the given object as pickled file with the given path

    Arguments
    ----------
    * path: str ~ The path of the pickled file that shall be created
    * pickled_object: Any ~ The object which shall be saved in the pickle file
    """
    pickle_file = open(path, 'wb')
    pickle.dump(pickled_object, pickle_file)
    pickle_file.close()


def mkdir(path: str) -> None:
    """Creates a directory if the given directory does not exist yet.

    Uses os.mkdir and os.path.exists internally.

    Argument
    ----------
    * path: str ~ The directory's path.
    """
    if not os.path.exists(path):
        os.mkdir(path)


def resolve_pathway_ids(pathway_ids: str, pathways: List[Any]) -> List[int]:
    """Returns a list index integer list of the given pathway_ids for the given pathway list.

    It gets 'all' and returns all indices of pathways or a comma-separated list of pathway numbers.

    Example
    ----------
    If an organism has 5 pathways in total, this function would return for 'all' as pathway_ids
    [0, 1, 2, 3, 4], and for '1,3,4' it would return [0, 2, 3]

    Arguments
    ----------
    * pathway_ids ~ The IDs of the pathways in a string. Can be 'all' to get all indices, or
      a comma-separated list of numbers.
    * pathways ~ The list of pathways of which an index integer list shall be created.
    """
    # Check if 'all' is given...
    if pathway_ids == "all":
        return [i for i in range(len(pathways))]
    # ...if not, it should be a comma-separated list of numbers.
    else:
        # Return a list with each number as integer and substracted by 1.
        # The substraction is done as e.g. the 1st pathway of an organism
        # has the index 0 in a Python pathway list.
        pathway_ids_list = pathway_ids.split(",")
        return [int(i) - 1 for i in pathway_ids_list]


def sanitize_path(text: str) -> str:
    """Replaces all invalid characters for a path with valid ones.

    E.g. useful for creating files with the name of KEGG pathways,
    as these names may contain invalid characters.

    Argument
    ----------
    * text: str ~ The string that may contain invalid characters.
    """
    return text.replace("\\", "_").replace("/", "_").\
        replace(":", "_").replace("*", "_").\
        replace("<", "_").replace(">", "_").\
        replace("|", "_")


def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder
