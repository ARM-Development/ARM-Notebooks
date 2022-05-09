
# Earth System Model Aerosol-Cloud Diagnostics Package

This document describes the version 2.0.0 of Earth System Model (ESM) aerosol-cloud diagnostics package (ESMAC Diags) that facilitate routine evaluation of aerosols, clouds and aerosol-cloud interactions simulated by the Department of Energy’s (DOE) Energy Exascale Earth System Model (E3SM). This package focuses on comparing simulated aerosols, clouds and aerosol-cloud interactions with in-situ and remote-sensing measurements from surface, aircraft, ship and satellite platforms. Various types of diagnostics and metrics are performed to assess how well E3SM represents observed aerosol properties and aerosol-cloud interactions across spatial scales. 

More information can be found in README_ESMAC_Diags_v2.0.0.pdf

## To install
This code is best run using a conda virtual environment. To install the required environment one can do
```bash
conda env create -f environment.yml
```
to set up a esmac_diags environment. Note if running this on a HPC system, you may need to load the appropriate module for anaconda. 

Once the environment has been created you can activate it with 
```
conda activate esmac_diags
``` 
and then this code can be installed with
```bash
pip install -e .
```
Which will install the code as editable allowing you to make changes to the codebase and it be reflected in the installed package. 


# Test run
To verify the package, enter scripts/ directory and run 
```bash
python run_testcase.py
```
Then go to the directory in testcase/figures/. Compare the output figure with the plot testcase/figures_verify/timeseries_organic_HISCALE.png. If the two figures look the same, the testcase is successfully run.


