# The PyRosetta initialization file format for precise reconstruction of Pose objects: Appendix material

Code and runtime environment for the appendix materials in the preprint publication:

<b>Klima JC. The PyRosetta initialization file format for precise reconstruction of Pose objects. Zenodo; 2026. doi:10.5281/zenodo.18113110</b>

Repository contents:

 - `appendix_a.py`: Python script used to generate data in Table 1.
 - `appendix_b.py`: Python script used to generate data in Figure 2.
 - `appendix_module.py`: Python module file containing helper functions imported by the Python scripts.
 - `requirements.txt`: uv project specifying the Python runtime environment used to run the Python scripts.
 - `LICENSE`: MIT license file associated with the code.
 - `.gitignore`, `README.md`, `CITATION.cff`: GitHub repository configuration files.

The PyRosetta build signature used to execute the Python scripts is the following:
`PyRosetta4.Release.python312.ubuntu.cxx11thread.serialization[extras:cxx11thread+serialization]2025.49+release.552382f5ef84f3c058ac922b717ceee014073be6`

The uv (Version 0.9.17) project specifying the Python runtime environment used to run the Python scripts was exported with the following command:
`uv export --no-hashes --no-annotate --no-header`

## Citation

If you use this code, please cite the Zenodo software archive:

```
@software{klima_2026_18157982,
  author       = {Klima, Jason C.},
  title        = {The PyRosetta initialization file format for
                   precise reconstruction of Pose objects: Appendix
                   material
                  },
  month        = jan,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18157982},
  url          = {https://doi.org/10.5281/zenodo.18157982},
}
```

If you reference the associated preprint manuscript, please cite the Zenodo preprint archive:

```
@misc{klima_2026_18113110,
  author       = {Klima, Jason C.},
  title        = {The PyRosetta initialization file format for
                   precise reconstruction of Pose objects
                  },
  month        = jan,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18113110},
  url          = {https://doi.org/10.5281/zenodo.18113110},
}
```
