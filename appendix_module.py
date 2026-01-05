__author__ = "Jason C. Klima"


import numpy as np
import os
import pyrosetta
import pyrosetta.distributed.io as io
import sys

from pyrosetta import Pose
from pyrosetta.rosetta.core.io.mmcif import dump_cif
from pyrosetta.rosetta.core.io.mmtf import dump_mmtf
from pyrosetta.rosetta.numeric.random import rg
from pyrosetta.rosetta.std import istringstream, ostringstream
from scipy.spatial.transform import Rotation
from toolz import curry
from typing import Generator, Optional


def coords_are_identical(pose1: Pose, pose2: Pose) -> bool:
    """
    Test whether the atomic coordinates of the input `Pose` objects
    are identical.
    """
    assert pose1.size() == pose2.size(), "Input `Pose` object sizes are not equal."
    for i in range(1, pose1.size() + 1):
        residue1 = pose1.residue(i)
        residue2 = pose2.residue(i)
        for atom in range(1, residue1.natoms() + 1):
            if list(residue1.xyz(atom)) != list(residue2.xyz(atom)):
                return False
    return True

is_close_to_zero: curry = curry(np.isclose)(
    b=0.0, rtol=0.0, atol=pow(10.0, -sys.float_info.dig)
)

def get_significant_digits(a: float, b: float) -> Optional[float]:
    """
    Return the significant digit agreement between two `float` values.
    If either of the `float` values is close to zero, return `None`.
    """
    if any(map(is_close_to_zero, (a, b))):
        return None
    rel = abs(a - b) / abs(a)
    if rel == 0.0:
        return float(sys.float_info.dig)
    else:
        return max(-np.log10(rel), 0.0)

def get_coord_digit_agreement(
   pose1: Pose, pose2: Pose
) -> Generator[float, None, None]:
    """
    Yield the significant digit agreement of atomic coordinates
    between two input `Pose` objects.
    """
    assert pose1.size() == pose2.size(), "Input `Pose` object sizes are not equal."
    for res in range(1, pose1.size() + 1):
        residue1 = pose1.residue(res)
        residue2 = pose2.residue(res)
        for atom in range(1, residue1.natoms() + 1):
            for value1, value2 in zip(residue1.xyz(atom), residue2.xyz(atom)):
                digits = get_significant_digits(value1, value2)
                if digits is not None:
                    yield digits

def roundtrip(pose_start: Pose, file_format: str) -> Pose:
    """
    Perform a round-trip save-and-reload operation on
    an input `Pose` object through an intermediate file format.
    """
    assert file_format in (
        "pdb", "mmcif", "mmtf", "silent", "pkl_pose", "b64_pose", "init"
    )

    file = f"./roundtrip.{file_format}"
    if os.path.isfile(file):
        os.remove(file)

    if file_format == "pdb":
        pose_start.dump_pdb(file)
        pose_final = pyrosetta.pose_from_file(file)
        assert coords_are_identical(pose_start, pose_final) is False
    elif file_format == "mmcif":
        dump_cif(pose_start, file)
        pose_final = pyrosetta.pose_from_file(file)
        assert coords_are_identical(pose_start, pose_final) is False
    elif file_format == "mmtf":
        dump_mmtf(pose_start, file)
        pose_final = pyrosetta.pose_from_file(file)
        assert coords_are_identical(pose_start, pose_final) is False
    elif file_format == "silent":
        pyrosetta.io.poses_to_silent(pose_start, file)
        pose_final = next(iter(pyrosetta.io.poses_from_silent(file)))
        assert coords_are_identical(pose_start, pose_final) is False
    elif file_format == "pkl_pose":
        io.dump_pickle(pose_start, file)
        pose_final = io.pose_from_file(file).pose
        assert coords_are_identical(pose_start, pose_final) is True
    elif file_format == "b64_pose":
        io.dump_base64(pose_start, file)
        pose_final = io.pose_from_file(file).pose
        assert coords_are_identical(pose_start, pose_final) is True
    elif file_format == "init":
        pyrosetta.dump_init_file(file, poses=pose_start, verbose=False)
        pose_final = io.pose_from_init_file(file).pose
        assert coords_are_identical(pose_start, pose_final) is True

    return pose_final

def get_random_sequence(rng: np.random.Generator) -> str:
    """Return a pseudorandom amino acid sequence."""
    size = rng.integers(1, high=1_000)
    sequence = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=size))
  
    return sequence

def random_translation(pose: Pose, rng: np.random.Generator) -> Pose:
    """Return the input `Pose` object pseudorandomly translated."""
    translation_vector = rng.uniform(-1_000, 1_000, size=3)
    pose.translate(translation_vector)
  
    return pose

def random_rotation(pose: Pose, rng: np.random.Generator) -> Pose:
    """Return the input `Pose` object pseudorandomly rotated."""
    rotation_matrix = Rotation.random(rng=rng).as_matrix()
    pose.rotate(rotation_matrix)

    return pose

def save_rg_state() -> str:
    """Save Rosetta RandomGenerator internal state."""
    stream = ostringstream()
    rg().saveState(stream)

    return stream.str()

def restore_rg_state(rng_state: str) -> None:
    """Restore Rosetta RandomGenerator internal state."""
    stream = istringstream(rng_state)
    rg().restoreState(stream)

def apply_moves(pose: Pose, rng_state: str) -> Pose:
    """
    Restore Rosetta RandomGenerator state, then iteratively apply
    `SmallMover` followed by `MinMover` to amplify lever-arm effects.
    """
    restore_rg_state(rng_state)

    small_mover = pyrosetta.rosetta.protocols.simple_moves.SmallMover()
    small_mover.angle_max(35.0)
    small_mover.nmoves(100)

    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    min_mover.score_function(pyrosetta.get_fa_scorefxn())
    min_mover.min_type("dfpmin_armijo")
  
    for _ in range(1_000):
        small_mover.apply(pose)
        min_mover.apply(pose)

    return pose

def calc_res_rmsd(pose1: Pose, pose2: Pose, res: int) -> float:
    """
    Return the per-residue heavy atom RMSD for a residue number
    in two `Pose` objects.
    """
    assert pose1.size() == pose2.size(), "Input `Pose` object sizes are not equal."
    residues = pyrosetta.rosetta.utility.vector1_unsigned_long()
    residues.append(res)

    return pyrosetta.rosetta.core.scoring.all_atom_rmsd_nosuper(
        pose1=pose1,
        pose2=pose2,
        pose1_residues=residues,
        pose2_residues=residues,
    )
