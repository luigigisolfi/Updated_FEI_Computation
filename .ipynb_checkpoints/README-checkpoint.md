# Welcome!

This repository contains both the source code and a demo to show how to compute the **Fragmentation Environmental Index**, as devised in: [L. Gisolfi Master's Thesis](https://thesis.unipd.it/retrieve/b00fb71a-4118-444b-bb0e-3ab77846ce05/Gisolfi_Luigi.pdf.pdf). 

You're welcome to use it for your own analysis, suggest improvements and/or blame the author for silly mistakes.

## Content

- Demo_FEI.ipynb: Jupyter Notebook showing how to get to the FEI.
- fei_library.py: Python file containing the source code
- dens_mean_2023.dat: ASCII file containing the density of each altitude shell.
  To get information on how the spatial density of objects varies as a function of the altitude, a simulation of the evolution of the space debris environment spanning 200 years was carried out. This considered the population of objects larger than 10 cm from the MASTER 2009 population

- background_pop.dat.5cm ASCII file containing the background population of objects.
  The background environment at the initial fragmentation epoch is derived from the MASTER population, while the environment at different epochs in the future is obtained by evolving the MASTER population with the SDM model.
  
- clouds: folder containing the data and outputs of each considered cloud (as for the example, only the two fragmentation altitudes: 450 km, 1200 km are given). The cloud folder also contains plots of the cumulative indexes. Users are encouraged to check the validity of the results.
  This is the output of the [Space Debris Mitigation Tool](https://conference.sdo.esoc.esa.int/proceedings/sdc5/paper/48/SDC5-paper48.pdf)

Cheers!
Luigi