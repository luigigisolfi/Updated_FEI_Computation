from fei_library import ComputeFEI
import os
object = ComputeFEI()

parent_mass = object.SetParent.mass(mass = 2000)
parent_e = object.SetParent.eccentricity(parent_e = 0.00003)
parent_inc = object.SetParent.inclination(parent_inc = 80.3)

h_frag = object.SetFragmentation.fragmentation_altitude(h_frag = 1200)

s_min = object.SetObservingNetwork.s_min(s_min = 20)
h_max = object.SetObservingNetwork.h_max(h_max = 2000)
elevation = object.SetObservingNetwork.constant_elevation(elevation = 30)
background_population_file = 'background_pop.dat.5cm'
cloud_name = 'nube_1200_km'
piece_of_string = '_' + str(round(s_min)) + 'cm_' + str(round(h_max)) + 'km'
network_type = 'optical'
clouds_folder_path = 'clouds'
color_plot = 'black'

pieces_of_strings = ['_5cm_2000km','_15cm_2000km','_20cm_2000km' ]
cloud_folder_path = os.path.join(clouds_folder_path, cloud_name)
output_folder_path = os.path.join(cloud_folder_path, 'output')
figures_folder = os.path.join(output_folder_path, 'figures' + f'/{network_type}' + piece_of_string)
global_csi_folder = os.path.join(output_folder_path, 'global_csi')
csi_cloud_only_folder = os.path.join(output_folder_path, 'cloud_only_csi')
object.multi_plotter_csi(pieces_of_strings, cloud_name, network_type, figures_folder, global_csi_folder, csi_cloud_only_folder)
exit()
object.optical_background(cloud_name, h_frag, piece_of_string, background_population_file)
day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list,global_csi_list_no_weights, ratios_list = object.optical(parent_mass, parent_e, parent_inc, cloud_name, h_frag, elevation, s_min, h_max, piece_of_string)

#object.plotter_csi(network_type, color_plot, day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list, s_min, clouds_folder_path, cloud_name, piece_of_string)
object.plotter_FEI(network_type,color_plot, clouds_folder_path, cloud_name, piece_of_string, h_frag)
object.multi_plotter_csi(pieces_of_strings, cloud_name, network_type, figures_folder, global_csi_folder, csi_cloud_only_folder)

