from fei_library import ComputeFEI
import os
object = ComputeFEI()

parent_mass = object.SetParent.mass(mass = 2000)
parent_e = object.SetParent.eccentricity(parent_e = 0.00003)
parent_inc = object.SetParent.inclination(parent_inc = 80.3)


cloud_name = 'nube_450_km'
h_frag = float(cloud_name.split('_')[1])
h_frag = object.SetFragmentation.fragmentation_altitude(h_frag = h_frag)
s_min = object.SetObservingNetwork.s_min(s_min = 25)
elevation = object.SetObservingNetwork.constant_elevation(elevation = 30)
background_population_file = 'background_pop.dat.5cm'
clouds_folder_path = 'clouds'

if h_frag < 1200:
    network_type = 'radar'
else:
    network_type = 'optical'

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

object.multi_plotter_csi(pieces_of_strings, cloud_name, network_type, figures_folder, global_csi_folder, csi_cloud_only_folder)
exit()
if network_type == 'radar':

    object.radar_background(cloud_name, h_frag, piece_of_string, background_population_file)
    day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list,global_csi_list_no_weights, ratios_list = object.radar(parent_mass, parent_e, parent_inc, cloud_name, h_frag, elevation, s_min, h_max, piece_of_string)

else:
    object.optical_background(cloud_name, h_frag, piece_of_string, background_population_file)
    day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list,global_csi_list_no_weights, ratios_list = object.optical(parent_mass, parent_e, parent_inc, cloud_name, h_frag, elevation, s_min, h_max, piece_of_string)

object.plotter_csi(network_type, color_plot, day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list, s_min, clouds_folder_path, cloud_name, piece_of_string)
object.plotter_FEI(network_type,color_plot, clouds_folder_path, cloud_name, piece_of_string, h_frag)

