# Accessing the Compute Environment

## Register for Access to the ARM Data Workbench

Participants who do not already have one should first register for an ARM account. When filling out the form, indicate that you are interested in using ARM's notebook examples.

You can sign up for an ARM account using [this link](https://adc.arm.gov/armuserreg/#/new).

If you already have an ARM account, email clustersupport@arm.gov that you need to be granted access to the workshop or course materials.

## Running the Notebooks
You can either run the notebook using [Binder](https://mybinder.org/) or on your local machine.

### Running on Jupyter

The simplest way to interact with a Jupyter Notebook is through the
[ARM Jupyter](https://jupyterhub.arm.gov), which enables the execution of a
[Jupyter Book](https://jupyterbook.org) on ARM infrastructure. The details of how this works are not
important for now. Navigate your mouse to
the top right corner of the book chapter you are viewing and click
on the rocket ship icon, (see figure below), and be sure to select
“launch Jupyterhub”. After a moment you should be presented with a
notebook that you can interact with. I.e. you’ll be able to execute
and even change the example programs. You’ll see that the code cells
have no output at first, until you execute them by pressing
{kbd}`Shift`\+{kbd}`Enter`. Complete details on how to interact with
a live Jupyter notebook are described in [Getting Started with
Jupyter](https://foundations.projectpythia.org/foundations/getting-started-jupyter.html).

### Running on Your Own Machine
If you are interested in running this material locally on your computer, you will need to follow this workflow:

1. Clone the `https://github.com/ARM-Development/ARM-Notebooks` repository:

   ```bash
    git clone https://github.com/ARM-Development/ARM-Notebooks
    ```  
1. Move into the `ARM-Notebooks` directory
    ```bash
    cd ARM-Notebooks
    ```  
1. Create and activate your conda environment from the `environment.yml` file
    ```bash
    conda env create -f environment.yml
    conda activate arm-tutorial-dev
    ```  