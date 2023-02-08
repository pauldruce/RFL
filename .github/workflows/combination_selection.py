import itertools
import random

default_config = {}

build_type = ["Release", "Debug"]
os_versions = [
    "ubuntu-22.04",
    "ubuntu-20.04",
    "ubuntu-18.04",
    "macos-12",
    "macos-11",
    "macos-10.15",
]
armadillo_version = ["10.8.2", "11.2.3"]

list_of_lists = [build_type, os_versions, armadillo_version]

combinations = [p for p in itertools.product(*list_of_lists)]

# print(combinations)
# print(len(combinations))


def create_combination_selection():
    random.seed(1)
    return random.sample(combinations, 3)


def create_new_matrix(subset):
    import json

    json_versions = []
    for version in subset:
        json_version = {
            "BUILD_TYPE": version[0],
            "OS": version[1],
            "ARMA_VERSION": version[2],
        }
        json_versions.append(json_version)

    return json.dumps(json_versions, separators=(",", ":"))


subset = create_combination_selection()
matrix = create_new_matrix(subset)
print(matrix)


# Proving that we can recreate the combinations by setting the random.seed
# for i in range(10):
#     random.seed(1)
#     print(random.sample(combinations, 3))
