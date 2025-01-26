from fei_library import ComputeFEI

object = ComputeFEI()

parent_mass = object.SetParent.mass(mass = 2000)
parent_e = object.SetParent.eccentricity(parent_e = 0.00003)
parent_inc = object.SetParent.inclination(parent_inc = 80.3)

h_frag = object.SetFragmentation.altitude(h_frag = 450)

s_min = object.SetObservingNetwork.s_min(s_min = 5)
h_max = object.SetObservingNetwork.h_max(h_max = 1200)
elevation = object.SetObservingNetwork.constant_elevation(elevation = 30)

cloud_name = 'nube_450_km'
piece_of_string = '_' + str(round(s_min)) + 'cm_' + str(round(h_max)) + 'km'

day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list,global_csi_list_no_weights, ratios_list = object.radar(parent_mass, parent_e, parent_inc, cloud_name, h_frag, elevation, s_min, h_max, piece_of_string)

print(day_list)
