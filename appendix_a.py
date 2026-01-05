__author__ = "Jason C. Klima"


import collections
import numpy as np
import pyrosetta

# Import locally defined helper methods
from appendix_module import (
    get_coord_digit_agreement,
    get_random_sequence,
    random_rotation,
    random_translation,
    roundtrip,
)


# Initialize PyRosetta
pyrosetta.init(options="-run:constant_seed 1 -mute all", silent=True)

# Set numpy RNG
rng = np.random.default_rng(111)

# Set number of round-trip repeats
num_repeats = 10

# Round-trip coordinates
file_format_digits = collections.defaultdict(list)
for _ in range(num_repeats):
    sequence = get_random_sequence(rng)
    pose_start = pyrosetta.pose_from_sequence(sequence)
    pose_start = random_translation(pose_start, rng)
    pose_start = random_rotation(pose_start, rng)
    for file_format in ("pdb", "mmcif", "mmtf", "silent", "pkl_pose", "b64_pose", "init"):
        pose_final = roundtrip(pose_start, file_format)
        for digits in get_coord_digit_agreement(pose_start, pose_final):
            file_format_digits[file_format].append(digits)

# Print results
for file_format, digits in file_format_digits.items():
    print(f"{file_format}: {np.mean(digits):.3g} ± {np.std(digits):.3g}")
