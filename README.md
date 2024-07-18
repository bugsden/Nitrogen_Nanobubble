# diffusion-flux

## To Install MD-Analysis 

conda create -n myenv python=3.8  #Creating the environment

conda config --add channels conda-forge
conda install -c conda-forge mdanalysis

python -c "import MDAnalysis; print(MDAnalysis.__version__)"
conda install mdanalysistests 

conda activate myenv  # Activating Conda environment before running Python scripts

## Calling MD-Analysis

# Links to refer 
https://userguide.mdanalysis.org/stable/installation.html

