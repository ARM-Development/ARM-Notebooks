# July 13, 2023 University of Puerto Rico Tutorials
These notebooks and presentations are from a tutorial session with students from the University of Puerto Rico.   These tutorials covered ARM's JupyterHub resources, basic python, and the [Atmospheric data Community Toolkit (ACT)](https://github.com/ARM-DOE/ACT).  The ACT tutorials were focused around exploring aerosol and lidar data from a dust event during TRacking Aerosol Convection interations ExpeRiment (TRACER) that occurred from July 16-19, 2022.

### JupyterHub and Jupyter Notebooks
[2-jupter_intro.ipynb](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/2-jupyter_intro.ipynb) is a notebook to get started with ARM's JupyterHub resources and navigate what JupyterHub notebooks are and some useful shortcuts.

### Python Tutorials
[3-Python-Basics.ipynb](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/3-Python-Basics.ipynb) is a what it seems, a notebook to go over beginner Python skills.  For those wanting to explore some more complex libraries in Python, checkout the other notebooks
- [3a-scientific_libraries_numpy.ipynb](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/3a-scientific_libraries_numpy.ipynb)
- [3b-scientific_libraries_pandas.ipynb](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/3b-scientific_libraries_pandas.ipynb)
- [3c-scientific_libraries_xarray.ipynb](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/3c-scientific_libraries_xarray.ipynb)

### Atmospheric data Community Toolkit (ACT) Tutorials
- [5-ACT-Basics-2023.ipynb](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/5-ACT-Basics-2023.ipynb) is a base tutorial using data from the TRACER Particle Soot Absorption Photometer (PSAP) to explore some of the base functionality of ACT including visualizations and working with quality control information.

- [6-ACT-TRACER-Dust.ipyn](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/6-ACT-TRACER-Dust.ipynb) is a notebook to bring a variety of datasets together for visualization and performing some of the more complex operations in ACT to produce a figure like the below.

![Output from the advanced ACT tutorial](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/images/output.png)
Figure 1. From the top, micropulse lidar linear depolarization ratio, aerodynamic particle sizer concentration, scanning mobility particle sizer concentration, aerosol chemical speciation monitor chemical compositions, single particle soot photometer black carbon concentration, PSAP aerosol absorption coefficient in the blue channel, and the surface meteorology station wind direction.

### Additional Tutorials
There are two additional tutorials for an [Introduction to GitHub and git](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/optional_github_intro.md) and how to perform [Branching and Pull Requests](https://github.com/ARM-Development/ARM-Notebooks/blob/main/Tutorials/UPR-ACT-2023/optional_github_branching.md) which are all very important when contributing work back into these open-source packages.
