__author__ = "Jason C. Klima"


import collections
import json
import numpy as np
import os
import pyrosetta

# Import locally defined helper methods
from appendix_module import (
    apply_moves,
    calc_res_rmsd,
    random_rotation,
    random_translation,
    roundtrip,
    save_rg_state,
)


# Initialize PyRosetta
pyrosetta.init(options="-run:constant_seed 1 -mute all", silent=True)

# Save Rosetta RG internal state
rg_state = save_rg_state()

# Set numpy RNG
rng = np.random.default_rng(222)

# Generate starting `Pose` object with reference coordinates
pose_start = pyrosetta.pose_from_sequence("TEST-CHANGING-THE-SIGNIFICANT-DECIMAL-DIGITS")
pose_start = random_translation(pose_start, rng)
pose_start = random_rotation(pose_start, rng)

# Move `Pose` object from reference coordinates, restoring Rosetta's RG state
pose_final = apply_moves(pose_start.clone(), rg_state)
# Dump `Pose` object for visualization
pose_final.dump_pdb(f"./pose_moved.pdb")

# Run moves on round-trip coordinates
file_format_res_rmsds = collections.defaultdict(dict)
file_format_poses = {}
for file_format in ("pdb", "mmcif", "mmtf", "silent", "pkl_pose", "b64_pose", "init"):
    # Round-trip coordinates through save-and-reload operation
    pose_roundtrip = roundtrip(pose_start, file_format)
    # Move `Pose` object from round-trip coordinates, restoring Rosetta's RG state
    pose_roundtrip_final = apply_moves(pose_roundtrip, rg_state)
    # Compute results between reference and round-trip coordinates
    for res in range(1, pose_final.size() + 1):
        rmsd = calc_res_rmsd(pose_final, pose_roundtrip_final, res)
        file_format_res_rmsds[file_format][res] = rmsd
    # Save `Pose` object for visualization
    file_format_poses[file_format] = pose_roundtrip_final

# Save data
with open("./file_format_res_rmsds.json", "w") as f:
    json.dump(file_format_res_rmsds, f)

# Dump `Pose` objects for visualization
for file_format, pose in file_format_poses.items():
    pose.dump_pdb(f"./pose_{file_format}_roundtrip_moved.pdb")
