# Welcome!

This repository contains both the source code and a demo to show how to compute the **Fragmentation Environmental Index**, as devised in: [L. Gisolfi Master's Thesis](https://thesis.unipd.it/retrieve/b00fb71a-4118-444b-bb0e-3ab77846ce05/Gisolfi_Luigi.pdf.pdf). 

You're welcome to use it for your own analysis, suggest improvements and/or blame the author for his silly mistakes.

## Repo Content
- **README.md**: (this file)

- **Demo_FEI.ipynb**: Jupyter Notebook showing how to get to the FEI.

- **Demo_FEI.py**: Python file corresponding to the Demo_FEI_ipynb. 

- **fei_library.py**: Python file containing the source code


- **dens_mean_2023.dat**: ASCII file containing the density of each altitude shell. \
  To get information on how the spatial density of objects varies as a function of the altitude, a simulation of the evolution of the space debris environment spanning 200 years was carried out. This considered the population of objects larger than 10 cm from the MASTER 2009 population


- **background_pop.dat.5cm**: ASCII file containing the background population of objects. \
  The background environment at the initial fragmentation epoch is derived from the MASTER population, while the environment at different epochs in the future is obtained by evolving the MASTER population with the SDM model.


## Link to Download the Data I used for my Analysis

- **clouds**: [Google Drive Link to the "clouds" folder](https://drive.google.com/drive/folders/1wuXl7Wdsg7ubhrNmb5JuKKdAsFx2YuL_?usp=drive_link), to be put in ./FEI/clouds. This contains the original data for each considered cloud (as for the example, only the two fragmentation altitudes: 450 km, 1200 km are given). The cloud folder also contains plots of the cumulative indexes. Users are encouraged to check the validity of the results. \
  The cloud_XXX_.fla files constitute the output of Space Debris Mitigation Tool simulation.

P.S. In case you were wondering..."nube" = cloud in italian :) 

Cheers! \
Luigi