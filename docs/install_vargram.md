# VARGRAM

## Quick install

Get the latest version with [pip](https://pip.pypa.io/en/stable/):
```
pip install vargram
```
It is recommended that you install VARGRAM on a virtual environment to avoid dependency conflicts. If you do not have Python installed yet or do not know how to create an environment, follow the instructions below.

## Installing Python 

Make sure that you have Python version ≥ 3.11 installed. Navigate to your terminal emulator (Terminal on Mac and Command Prompt on Windows) and determine the current version on the command line by entering `python --version`. 

If you do not have Python installed yet, you can download Python from the official [website](https://www.python.org/downloads/). You can also instead [download Python with Anaconda](https://anaconda.com/download), which comes with preinstalled libraries, environment management with Conda and other features. Installing Python this way is recommended if you are a beginner.

## Creating an environment

As an optional step, you may create a virtual environment, which is essentially an isolated folder structure where you can install and run VARGRAM. In this way, you can avoid running into package dependency conflicts. If you downloaded the Anaconda distribution, you may create an environment called, for e.g., `myvargram` with Python version e.g. 3.12 using `conda`:
```
conda create --name myvargram python=3.12
```
To access this environment, simply run `conda activate myvargram`. On the command line, you should see e.g. `(myvargram)` before your username, indicating that you are in the `myvargram` environment. To exit the environment, run `conda deactivate`.

In your created environment, install VARGRAM with pip:
```
pip install vargram
```
