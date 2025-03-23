#!/usr/bin/env python
# coding: utf-8

# # FRAGMENTATION INDEX COMPUTATION 
# 
# ## Objectives
# This code shows how to use functions in the `fei_library.py` file in order to compute the (Upgraded) Fragmentation Environmental Index, as devised in [L. Gisolfi Master's Thesis](https://thesis.unipd.it/retrieve/b00fb71a-4118-444b-bb0e-3ab77846ce05/Gisolfi_Luigi.pdf.pdf). 
# 
# **Author**: Luigi Gisolfi
# 
# **Year**: 2025
# 
# 
# Please, create the folder nube_XXX_km with the all the files nube_XXX_km/data/cloud_XXX.fla 
# Make sure you have access to the file dens_mean_no_weights_2023.dat
# Note: For the background, only MASTER population objects > 5 cm is considered (file: background_pop.dat.5cm) has to be used as input

# ## Import Statements
# Let's start by importing some relevant python modules that will be useful for the computation.

# In[1]:


from fei_library import ComputeFEI
import os
object = ComputeFEI()


# ## Set Parent Objects Properties
# 
# The `SetParent` class allows users to store the Parent's eccentricity and inclination (in degrees), as well as its mass (in kg).
# Of course, these must match the properties of the object that underwent fragmentation. 

# In[2]:


parent_mass = object.SetParent.mass(mass = 2000)
parent_e = object.SetParent.eccentricity(parent_e = 0.00003)
parent_inc = object.SetParent.inclination(parent_inc = 80.3)


# ## Retrieve h_frag
# 
# Each cloud folder is named after the given fragmentation altitude (nube_XXX_km, where XXX indicates the fragmentation altitude. Also, "nube" is the word for "cloud" in italian!), so we use this fact to retrieve the correspoding fragmentation altitude, $h_{frag}$.

# In[3]:


cloud_name = 'nube_1200_km'
h_frag = float(cloud_name.split('_')[1])
h_frag = object.SetFragmentation.fragmentation_altitude(h_frag = h_frag)


# ## Observing Network Settings
# 
# The `SetObservingNetwork` class allows to set the minimum observable size $s_{min}$ at a given reference altitude $h_{max}$. 
# Setting the network elevation is also possible and achieved via the `SetObservingNetwork.constant_elevation` metho (the elevation could in prinicple be time dependent, but no time dependent elevation method has been developed yet...)

# In[4]:


s_min = object.SetObservingNetwork.s_min(s_min = 5)
elevation = object.SetObservingNetwork.constant_elevation(elevation = 30)


# ## Background Population File
# 
# At this point, the background population file path has to be specified. 

# In[5]:


background_population_file = 'background_pop.dat.5cm'
clouds_folder_path = 'clouds'


# ## Network Type Automatic Selection
# 
# As explained in [L. Gisolfi Master's Thesis](https://thesis.unipd.it/retrieve/b00fb71a-4118-444b-bb0e-3ab77846ce05/Gisolfi_Luigi.pdf.pdf), we assume radar can probe LEO from the lowest shells up to (excluded) $h_{max} = 1200$ km, while optical networks are active between $h = 1200$ km and $h_{max} = 2000$ km.

# In[6]:


if h_frag < 1200:
    network_type = 'radar'
else:
    network_type = 'optical'


# ## Computing Relevant Quantities for Subsequent Analysis
# We are now ready to compute the global csi, the cloud csi and other quantities. These will allow us to compute the FEI. 
# We:
# - compute the weighted and unweighted CSI contributions to the background population (as it was before the fragmentation epioch, given the background population file defined above)
# - compute the weighted and unweighted CSI contributions of all objects after the fragmentation.
# 
# Relevant values are saved to files in the `nube_XXX_km/output` folders. 

# In[7]:


if network_type == 'radar':
    color_plot = 'red'
    h_max = object.SetObservingNetwork.h_max(h_max = 1200)
    pieces_of_strings = ['_5cm_1200km','_15cm_1200km','_20cm_1200km' ]
else:
    color_plot = 'black'
    h_max = object.SetObservingNetwork.h_max(h_max = 2000)
    pieces_of_strings = ['_5cm_2000km','_15cm_2000km','_20cm_2000km' ]

piece_of_string = '_' + str(round(s_min)) + 'cm_' + str(round(h_max)) + 'km'
cloud_folder_path = os.path.join(clouds_folder_path, cloud_name)
output_folder_path = os.path.join(cloud_folder_path, 'output')
figures_folder = os.path.join(output_folder_path, 'figures')
global_csi_folder = os.path.join(output_folder_path, 'global_csi')
csi_cloud_only_folder = os.path.join(output_folder_path, 'cloud_only_csi')

if network_type == 'radar':

    object.radar_background(cloud_name,
                            h_frag, 
                            piece_of_string, 
                            background_population_file
                           )
    
    day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list,global_csi_list_no_weights, ratios_list = object.radar(parent_mass, 
                                                                                                                                                        parent_e, 
                                                                                                                                                        parent_inc, 
                                                                                                                                                        cloud_name, 
                                                                                                                                                        h_frag, 
                                                                                                                                                        elevation, 
                                                                                                                                                        s_min, 
                                                                                                                                                        h_max, 
                                                                                                                                                        piece_of_string)
else:
    object.optical_background(cloud_name, h_frag, piece_of_string, background_population_file)
    day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list,global_csi_list_no_weights, ratios_list = object.optical(parent_mass, 
                                                                                                                                                          parent_e, 
                                                                                                                                                          parent_inc, 
                                                                                                                                                          cloud_name, 
                                                                                                                                                          h_frag, 
                                                                                                                                                          elevation, 
                                                                                                                                                          s_min, 
                                                                                                                                                          h_max, 
                                                                                                                                                          piece_of_string)

print('Done processing.')   


# ## Plotting the Results
# Finally, the CSI and FEI values are plotted (using the `plotter_csi()`, `plotter_FEI()` and `multi_potter_csi()` methods), and the **cumulative index** file corresponding to the given fragmentation and network is created via `get_cumulative_index_files()`.

# In[8]:
object.plotter_csi(network_type, color_plot, day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list, s_min, clouds_folder_path, cloud_name, piece_of_string)
object.plotter_FEI(network_type,color_plot, clouds_folder_path, cloud_name, piece_of_string, h_frag)
object.multi_plotter_csi(pieces_of_strings, cloud_name, network_type, figures_folder, global_csi_folder, csi_cloud_only_folder)
object.get_cumulative_index_files(output_folder_path, network_type, piece_of_string)


# ## Cumulative Index: How To Plot It 
# After running the above cells for all fragmentations altitudes (i.e. 450 km, 800 km, 1200 km, 1800 km) and for different networks (for instance, setting $s_{min} = 5,10,15,20,25$ cm) the user would be able to get the cumulative index plots via the following (commented for now) lines:

# In[ ]:


#cloud_names = ['nube_1200_km', 'nube_1800_km']
#network_type = 'optical'
#object.plot_cumulative_indexes(clouds_folder_path, network_type, cloud_names)

