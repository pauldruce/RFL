import itertools
import json
import random
import sys


def create_combination_selection(combinations, random_seed, num_selected):
    random.seed(random_seed)
    return random.sample(combinations, num_selected)


def create_new_matrix(subset):
    json_versions = []
    for version in subset:
        json_version = {
            "BUILD_TYPE": version[0],
            "OS": version[1],
            "ARMA_VERSION": version[2],
        }
        json_versions.append(json_version)

    return json.dumps(json_versions, separators=(",", ":"))


def main(random_seed, num_selected, select_all):
    build_type = ["Release"]
    os_versions = [
        "ubuntu-latest",
        "ubuntu-20.04",
        "macos-13",
        "macos-latest"
    ]
    armadillo_version = ["10.8.2", "11.4.4"]

    list_of_lists = [build_type, os_versions, armadillo_version]

    combinations = [p for p in itertools.product(*list_of_lists)]

    if not select_all:
        subset = create_combination_selection(combinations, random_seed, num_selected)
    else:
        subset = combinations

    matrix = create_new_matrix(subset)
    print(matrix)
    return 0


if __name__ == "__main__":
    command_line_arguments = sys.argv
    if len(command_line_arguments) > 3:
        sys.exit(""""
        Too many arguments passed in. 
        Valid arguments are a positive integer representing the seed for random selection,
        followed by another positive integer representing the number of elements to select.
        I.e. "python -u combination_selection.py 1234 5" would seed the random selection
        using the value 1234 and would select 5 extra combinations.
                 """)

    select_all = False
    random_seed = None
    num_selected = 3

    if len(command_line_arguments) == 2:
        if command_line_arguments[1] == "all":
            select_all = True
        else:
            random_seed = int(command_line_arguments[1])

    if len(command_line_arguments) == 3:
        random_seed = int(command_line_arguments[1])
        num_selected = int(command_line_arguments[2])

    sys.exit(main(random_seed, num_selected, select_all))
