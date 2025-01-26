# # FRAGMENTATION INDEX COMPUTATION
#
# **Author**: Luigi Gisolfi
#
# **Year**: 2023
#
# This code computes the (Upgraded) Fragmentation Environmental Index, as devised in [L. Gisolfi Master's Thesis](https://thesis.unipd.it/retrieve/b00fb71a-4118-444b-bb0e-3ab77846ce05/Gisolfi_Luigi.pdf.pdf).
#
# Please, create the folder nube_XXX_km with the all the files cloud_XXX.fla
# Also, create the empty folders: weights, figures, csi_out (containing csi0.out) INSIDE the folder named nube_XXX_km before running this code.
# Make sure you have access to the file dens_mean_no_weights_2023.dat
# Note: For the background, only MASTER population objects > 5 cm is considered (file: background_pop.dat.5cm) has to be used as input
#

# ## Import Statements
# Let's start by importing some relevant python modules that will be useful for the computation.

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
import os
import math


# # Auxiliary Functions - A library
# In what follows, we define some functions that we will need for
# * the computation of both radar and optical weights to be applied to the Criticality of Spacecraft Index (csi)
# * the computation of the csi
#
# ## Radar Range equation
# The observable given by a Resident Space Object is often given in terms of **range**, but the csi definition rather relies on its **altitude**. The functions _h_to_rho_ and *rho_to_h* allow the user to go back and forth between the two quantities.
#
# An object of the Earth's surface (with radius $r_{e} = 6378$ km) has an altitude of $0$ km.

# In[5]:

class ComputeFEI:

    # Assumes folder structure:
    # clouds/cloud_XXXkm/cloud_XXX.fla
    def __init__(self, elevation = None, r_e = 6378, mean_density_file_path = f'/Users/luigigisolfi/dens_mean_2023.dat', clouds_folder_path = 'clouds'):
        """
        Initialize the ComputeFEI class.

        Parameters:
        r_e (float): Radius of the Earth.
        elevation (float): Elevation angle of the telescope (in radians).
        """
        self.SetObservingNetwork = self.SetObservingNetwork(self)
        self.SetFragmentation = self.SetFragmentation(self)
        self.SetParent = self.SetParent(self)

        self.r_e = r_e
        self.elevation = elevation
        self.mean_density_file_path = mean_density_file_path
        self.clouds_folder_path = clouds_folder_path

    class SetParent:
        def __init__(self, outer_instance):
            """
            Initialize the SetParent class.

            Parameters:
            outer_instance (ComputeFEI): Reference to the outer class instance.
            """
            self.outer_instance = outer_instance
        def mass(self, mass):

            self.parent_mass = mass

            return(self.parent_mass)

        def eccentricity(self, parent_e):

            self.parent_e = parent_e

            return(self.parent_e)

        def semi_major_axis(self, h_frag):

            self.parent_a = h_frag + self.outer_instance.r_e

            return(self.parent_a)

        def inclination(self, parent_inc):

            self.parent_inc = parent_inc

            return(self.parent_inc)

    class SetFragmentation:
        def __init__(self, outer_instance):
            self.outer_instance = outer_instance

        def altitude(self, h_frag):
            self.h_frag = h_frag

            return(self.h_frag)

    class SetObservingNetwork:
        def __init__(self, outer_instance):
            self.outer_instance = outer_instance

        def constant_elevation(self, elevation):
            self.elevation = elevation
            self.outer_instance.elevation = elevation
            return(self.elevation)

        def h_max(self, h_max):
            self.h_max = h_max
            return(self.h_max)

        def s_min(self, s_min):
            self.s_min = s_min
            return(self.s_min)


    def rho_to_h(self, rho):

        """
        Compute the height (h) of the target given its range (rho).

        Parameters:
        rho (float): Target range.

        Returns:
        float: Height of the target above Earth's surface.

        """
        r_e = self.r_e
        elevation = self.elevation
        r = np.sqrt(r_e**2 + rho**2 - 2 * r_e * rho * np.cos(np.pi / 2 + elevation))
        self.h = r - r_e
        return(self.h)

    def h_to_rho(self, h):
        r_e = self.r_e
        elevation = self.elevation
        self.rho = r_e*(np.sqrt(((h+r_e)/r_e)**2 - np.cos(elevation)**2) - np.sin(elevation))
        return(self.rho)


    # ## Magnitude and Optical signature
    #
    # The magnitude of an RSO, together with its optical signature are computed as underlined in [Shell et al, 2010.](https://amostech.com/TechnicalPapers/2010/Systems/Shell.pdf)

    # In[7]:


    def m_obj(self,s,rho):
        self.m_obj = -26.732 - 2.5*np.log10(((s/100)**2/(rho*1000)**2)*0.175*(0.25 + 2/(3*np.pi))) #s in cm, rho in km #diffuse and reflected specular component considered
        return (self.m_obj)

    def e_rso(self,s,rho):

        self.e_rso = (5.6*10**10)*10**(-0.4*self.m_obj(s,rho)) #in photons/s/m^2
        return(self.e_rso)


    # ## Threshold Values at a Given Altitude Fragmentation $h_{frag}$
    #
    # In our main work, the performance of an optical sensor is set by defining the **faintest magnitude** of an object that is able to trigger a detection. Since, for both optical and radar sensors, we are talking about **reflective objects**, we assume that **the bigger the size (diameter), the brighter the object**. Therefore, the _capability_ of a sensor can be set by defining the minimum size of an object that can be detected at a maximum altitude. This is computed by the two threshold and threshold_radar functions.
    #
    # For each given size i (in $cm$, starting from $0.01$ $cm$), the ratio between the ratio of the magnitude of the object of size i at the collision altitude verus the least bright detectable object (which defines the capability) is computed. Due to the nature of magnitudes ( lower magnitude = brighter object), as soon as this ratio gets smaller than one, i defines the size of the smallest object that can be detected at $h_{frag}$.
    #
    # A similar threshold, involving the ratio between cross sections, allows to determine the minimum detectable size for radar sensors in _threshold_radar_

    # In[9]:


    def optical_threshold(self, h_frag,s_min,h_max):
        rho_frag = self.h_to_rho(h_frag)
        rho_max = self.h_to_rho(h_max)

        self.optical_threshold_value = 0
        ratio = 2
        while ratio >= 1:
            self.optical_threshold_value += 0.01
            ratio = self.m_obj(self.optical_threshold_value,rho_frag)/self.m_obj(s_min,rho_max)

        return(np.round(self.optical_threshold_value,1))

    def radar_threshold(self, h_frag,s_min,h_max):
        self.radar_threshold_value = 0.01
        ratio = 0.00001
        rho_frag = self.h_to_rho(h_frag)

        sigma_min = (np.pi/4)*s_min**2
        rho_max = self.h_to_rho(h_max)

        while ratio <= 1:
            self.radar_threshold_value += 0.01
            sigma_fragment = (np.pi/4)*self.radar_threshold_value**2
            ratio = (sigma_fragment/rho_frag**4)/(sigma_min/rho_max**4)

        return(np.round(self.radar_threshold_value, 2))


    # ## Optical Weight Computation
    # The optical weight is computed taking two factors into account:
    # 1) what SNR the RSO produces ($\omega_{t_{sig}}$)
    # 2) how fast the RSO is in the FOV (linked to the RSO altitude, $\omega_{E_{RSO}}$)
    #
    #
    # The function _get_omega_optical_sum_ allows to perform a weighted sum of the two optical weights

    # In[11]:


    def w_e_rso(self, s_fragment, a_fragment):
        h_fragment = a_fragment - r_e
        rho_fragment = self.h_to_rho(h_fragment)
        rho_max = self.h_to_rho(h_max)
        if self.m_obj(s_fragment,rho_fragment) <= self.m_obj(s_min,rho_max): #check visibility

            if s_fragment <= s_min:
                self.omega_e_rso = 1 - (self.e_rso(s_fragment,rho_fragment)/self.e_rso(s_min,rho_fragment))
            else:
                self.omega_e_rso = 0
        else:
            self.omega_e_rso = 1

        return (self.omega_e_rso)

    def w_t_sig(self,h_coll): #does not depend on fragment, only depends on h_coll of collision

        vel_HL = 500
        vel_ML = 1000
        vel_LL = 2000

        if (0 <= h_coll <= 500):
            self.w_t_sig = 1 - vel_HL/vel_LL
        elif(500 < h_coll <= 1000):
            self.w_t_sig = 1 - vel_HL/vel_ML
        else:
            self.w_t_sig = 0

        return(self.w_t_sig)

    def get_omega_optical_sum(self, omega_i_elem, omega_j_elem):

        if omega_i_elem == 1:
            self.sum_optical_weights = omega_i_elem
            return(self.sum_optical_weights)

        else:
            if (0<= h_coll <= 500):
                A = 0.1
                self.sum_optical_weights = (omega_i_elem*A + omega_j_elem)
            elif (500< h_coll <= 1000):
                A = 0.5
                self.sum_optical_weights = (omega_i_elem*A + omega_j_elem)
            else:
                self.sum_optical_weights = omega_i_elem

            return(self.sum_optical_weights)


    # ## Radar Weight
    #
    # The radar weight computation is more straightforward, as it ultimately only depends on the ratio between two radar cross sections.

    # In[18]:


    def w_radar(self, s_fragment,a_fragment):

        h_fragment = a_fragment - r_e
        sigma_fragment = (np.pi/4)*s_fragment**2
        rho_fragment_radar = self.h_to_rho(h_fragment)

        sigma_min_radar = (np.pi/4)*s_min**2
        rho_max_radar = self.h_to_rho(h_max)

        if (sigma_fragment/rho_fragment_radar**4) >= (sigma_min_radar/rho_max_radar**4):

            if s_fragment<= s_min:
                self.w_radar = 1 - sigma_fragment/sigma_min_radar #the two rho_frag cancel out

            else:
                self.w_radar = 0

        else:
            self.w_radar = 1

        return(self.w_radar)


    #

    # ## Computing the fractional csi
    #
    # The (modified) fractional csi, as defined in [Bombardelli et al](https://www.sciencedirect.com/science/article/abs/pii/S0273117717302491) and here incorporating the optical or radar weight, is computed as follows.
    #
    # The two functions are slightly different. The first one is optimized to compute the _fractional_csi_ for multiple objects at the same time. The second one, _parent_fractional_csi_ is used for single objects, and we will use it later to compute the contribution to the csi of the parent object.

    # In[22]:


    def fractional_csi(self, mass,a,e,inc,weight, r_in, r_out):

        h_in = r_in - r_e
        h_fragment = a - r_e

        phi = self.get_phi(a,e, r_in,r_out)
        mass_norm = mass/10000
        dens = self.h_to_dens(h_in)/(6.8*10**(-8))
        f = (13-3*np.cos(inc*np.pi/180))/16
        l = list(map(self.life, h_fragment))
        self.f_csi = phi*mass_norm*dens*l*f*weight

        return(self.f_csi)

    def parent_fractional_csi(self,mass,a,e,inc,weight,r_in, r_out):

        h_in = r_in - r_e
        h_fragment = a - r_e

        phi = self.get_phi_parent(a,e,r_in,r_out)

        mass_norm = mass/10000
        dens = self.h_to_dens(h_in)/(6.8*10**(-8))
        parent_f = (13-3*np.cos(inc*np.pi/180))/16
        parent_l = self.life(h_fragment)
        self.f_csi = phi*mass_norm*dens*parent_l*parent_f*weight

        return(self.f_csi)


    # ## csi Elements Computation
    #
    # The csi of an object depends on
    # * object lifetime
    # * object density of the crossed altitude shells
    # * time spent in each altitude shell
    # * object mass
    # * object orbital inclination
    #
    # The following auxiliary functions: _life_, _get_phi_ and _h_to_dens_ allow to compute all the necessary terms to compute the csi.
    # Please note that, _get_phi_parent_ is used later on for the parent object.

    # In[25]:


    def life(self,h_fragment):

        #coefficients from Bombardelli

        h_ref = 1031.5 #in km
        a_coeff= 6.5215
        b_coeff=0.2583
        c_coeff=-33.8481

        life_fragment = math.exp(a_coeff*(h_fragment**b_coeff) + c_coeff)
        life_ref = math.exp(a_coeff*(h_ref**b_coeff) + c_coeff)
        self.l = life_fragment/life_ref

        if self.l >= 1:
            self.l = 1
            return(self.l)
        else:
            return(self.l)

    def get_phi(self, a,e, r_in, r_out):

        peri = a*(1-e)
        apo = a*(1+e)
        E_out = np.arccos((a - r_out)/(a*e)) #E_out of selected objects only
        E_in = np.arccos((a - r_in)/(a*e))#E_in of selected objects only

        condition1_peri = peri > r_in
        condition1_apo = apo < r_out
        condition2_peri = peri < r_in
        condition2_apo = apo > r_out
        condition3_peri = peri < r_in
        condition31_apo =  apo < r_out
        condition32_apo = apo > r_in
        condition41_peri = peri < r_out
        condition42_peri = peri > r_in
        condition4_apo = apo > r_out
        condition0_peri = peri > r_out
        condition0_apo = apo < r_in

        pos0 = np.where(condition0_peri|condition0_apo)
        pos1 = np.where(condition1_peri & condition1_apo)
        pos2 = np.where(condition2_peri & condition2_apo)
        pos3 = np.where(condition3_peri & (condition31_apo & condition32_apo))
        pos4 = np.where((condition41_peri & condition42_peri) & condition4_apo)

        self.phi_array = np.zeros(len(peri))
        if len(pos0[0]) != 0:
            self.phi_array[pos0] = 0
        if len(pos1[0]) != 0:
            self.phi_array[pos1] = 1
        if len(pos2[0]) != 0:
            self.phi_array[pos2] = (E_out[pos2] - E_in[pos2] - e[pos2]*(np.sin(E_out[pos2]) - np.sin(E_in[pos2])))/np.pi
        if len(pos3[0]) != 0:
            self.phi_array[pos3] = 1 - (E_in[pos3] - e[pos3]*(np.sin(E_in[pos3])))/np.pi
        if len(pos4[0]) != 0:
            self.phi_array[pos4] = (E_out[pos4] - e[pos4]*(np.sin(E_out[pos4])))/np.pi
        else:
            self.phi_array = self.phi_array

        return(self.phi_array)

    def h_to_dens(self, h):

        altitude, dens = np.loadtxt(self.mean_density_file_path, unpack = True, usecols = (0,1))
        altitude = np.round(altitude)
        pos = np.where(altitude == h)[0]
        pos_1 = pos +1
        self.density = (dens[pos] + dens[pos_1])/2
        return(self.density)

    def get_phi_parent(self, a,e, r_in, r_out):

        E_out_par = np.arccos((a - r_out)/(a*e))
        E_in_par = np.arccos((a - r_in)/(a*e))

        peri_par = a*(1-e)
        apo_par = a*(1+e)

        cond0_peri_par = peri_par > r_out
        cond0_apo_par = apo_par < r_in
        condition1_peri = peri_par > r_in
        condition1_apo = apo_par < r_out
        condition2_peri = peri_par < r_in
        condition2_apo = apo_par > r_out
        condition3_peri = peri_par < r_in
        condition31_apo =  apo_par < r_out
        condition32_apo = apo_par > r_in
        condition41_peri = peri_par < r_out
        condition42_peri = peri_par > r_in
        condition4_apo = apo_par > r_out

        if (cond0_peri_par|cond0_apo_par):
            self.phi = 0
            return(self.phi)

        if (condition1_peri and condition1_apo):
            self.phi = 1
            return(self.phi)
        elif (condition2_peri and condition2_apo):
            self.phi = (E_out_par - E_in_par -e *(np.sin(E_out_par) - np.sin(E_in_par)))/np.pi
            return(self.phi)
        elif (condition3_peri and (condition31_apo and condition32_apo)):
            return(self.phi)
        elif ((condition41_peri and condition42_peri) and condition4_apo):
            self.phi = (E_out_par - e*(np.sin(E_out_par)))/np.pi
            return(self.phi)
        else:
            print('I dont really know what to print...')
            return(self.phi)


    # ## Radar Main
    #
    # This function is called as a main when a pure radar network is assumed to be in place.
    # In it, all the above defined functions are used, so as to retrieve the needed information and compute the FEI.
    #
    # It takes as inputs:
    # * the cloud file (from MASTER)
    # * the fragmentation altitude $h_{frag}$
    # * the minimum detectable size at $h_{max}$
    # * the maximum altitude at which the use of radar sensors is considered to be effective
    # * a piece of string, created with s_min and h_max, so it is in principle not needed...
    #
    # It gives as outputs:
    # * day_list (list of epoch after fragmentation, in Days)
    # * global_csi_cloud_only_list (weights are considered)
    # * global_csi_cloud_only_list_no_weights (weights are all set to one)
    # * global_csi_list (cloud + background csi, weights are considered)
    # * global_csi_list_no_weights (cloud + background csi, weights are all set to one)
    # * ratios_list (percentage FEI values)

    # In[21]:


    def radar(self, parent_mass, parent_e, parent_inc, cloud_name, h_frag, elevation, s_min, h_max, piece_of_string):
        r_e = self.r_e
        elevation = self.SetObservingNetwork.constant_elevation(elevation)
        clouds_folder_path = self.clouds_folder_path

        radar_threshold_value = self.radar_threshold(h_frag,s_min,h_max) #minimum detectable size for a given rho_frag

        print(f'The minimum detectable size for a fragment at {h_frag} km is: {radar_threshold_value} cm\n')

        self.global_csi_list = []
        self.global_csi_list_no_weights = []
        self.global_csi_cloud_only_list = []
        self.global_csi_cloud_only_list_no_weights = []

        self.day_list = []
        self.ratios_list = []
        self.ratios_no_weights_list = []

        parent_mass = self.SetParent.mass(2000) # parent_mass = 2000 kg
        parent_a = self.SetParent.semi_major_axis(h_frag)
        parent_e =  self.SetParent.eccentricity(parent_e) #parent_e = 0.00003
        parent_inc = self.SetParent.inclination(parent_inc) # parent_inc = 80.3

        cloud_folder_path = os.path.join(clouds_folder_path, cloud_name)
        csi_background_folder = os.path.join(cloud_folder_path, 'csi_out')
        background_filename = 'csi0_radar' + piece_of_string + '.out'
        figures_folder_path = os.path.join(clouds_folder_path, 'figures')
        weights_folder_path = os.path.join(clouds_folder_path, 'weights')

        csi_background, csi_background_no_weights = np.loadtxt(os.path.join(csi_background_folder, background_filename), unpack = True, usecols = (0,1))
        csi_background_array = np.array(csi_background)
        csi_background_no_weights_array = np.array(csi_background_no_weights)

        if not os.path.isdir(figures_folder_path + '/radar' + piece_of_string):
            os.mkdir(figures_folder_path + '/radar' + piece_of_string)

        if not os.path.isdir(weights_folder_path + '/radar' + piece_of_string):
            os.mkdir(weights_folder_path + '/radar' + piece_of_string)

        for filename in os.listdir(cloud_folder_path):
            print(f'Processing cloud file: {filename}')
            f = os.path.join('/Users/luigigisolfi/' + str(cloud_name), filename)
            if (os.path.isdir(f) == True):
                print(f'File: {f} is a directory. Skipping...\n')
                continue
            else:
                if (filename[6:9].isnumeric() and float(filename[6:9]) > 100): # process only the first 100 days.
                    continue
                elif (f[-4:] == '.fla' and filename != 'cloud_init.fla' and filename != 'fragment.fla' and filename != 'delta_v_1200_km.fla'):
                    np.set_printoptions(threshold=np.inf)
                    epoch, n, area, mass, a, e, inc, Omega, omega, M = np.loadtxt(f, unpack = True, usecols = (0,1,2,3,4,5,6,7,8,9))
                    f_weights = weights_folder_path + '/radar' + piece_of_string + '/' + filename[:-4] + '_weights.fla'
                else:
                    continue

            sizes = np.sqrt(4*area/np.pi) #area in m^2, sizes in m
            sizes = sizes*100 #in cm
            sizes_cond = sizes >=1
            sensor_cond = a-r_e < 1200
            combined = sizes_cond & sensor_cond
            e_cond = e <= 0.5
            combined = combined & e_cond
            a = a[combined]
            sizes = sizes[combined]
            e = e[combined]
            inc = inc[combined]
            mass = mass[combined]

            weights = list(map(self.w_radar, sizes,a))

            csi_shell_list = []
            csi_shell_list_no_weights = []
            csi_parent_list = []
            csi_parent_list_no_weights = []
            csi_post_list= []
            csi_post_list_no_weights = []
            csi_pre_list = []

            no_weights = np.ones(len(weights))

            np.set_printoptions(threshold=np.inf)

            array = np.transpose(np.array([sizes, weights])) #array with sizes and associated weights

            with open(weights_folder_path + '/radar' + piece_of_string + '/' + str(filename[:-4]) + '_weights.fla', 'w') as fw:

                for line in array:
                    fw.writelines(str(line)[1:-1] + '\n')

                for r_in in range(200 + r_e, 1200 + r_e, 50): # define the shells in the range 200 - 1200 km (radar shell range)

                    r_out = r_in + 50

                    f_csi_shell = self.fractional_csi(mass,a,e,inc,weights, r_in, r_out) #all objects weighted fractional csi on a single shell
                    f_csi_shell_array = np.array(f_csi_shell) #make it into array so as to perform operations on it later on

                    f_csi_shell_no_weights =self.fractional_csi(mass,a,e,inc,no_weights, r_in, r_out)#all objects fractional csi (w=1) on a single shell
                    f_csi_shell_no_weights_array = np.array(f_csi_shell_no_weights)

                    f_csi_shell_parent = 0
                    f_csi_shell_parent_no_weights = self.parent_fractional_csi(parent_mass,parent_a,parent_e,parent_inc,1,r_in,r_out) #parent object csi (w =1)

                    total_csi_shell = np.sum(f_csi_shell_array)  #computing sum on j of (Post_j-Pre_j) = (Frag_j - Parent_j) on every shell j
                    total_csi_shell_no_weights = np.sum(f_csi_shell_no_weights_array) #same but w = 1

                    csi_shell_list.append(total_csi_shell) #append the result into a list made of each day's results
                    csi_shell_list_no_weights.append(total_csi_shell_no_weights) #same as above with w = 1

                    csi_parent_list.append(f_csi_shell_parent)
                    csi_parent_list_no_weights.append(f_csi_shell_parent_no_weights)

                csi_shell_list_array = np.array(csi_shell_list) #make list into array so as to perform operations later on
                csi_shell_list_no_weights_array = np.array(csi_shell_list_no_weights)  #make list into array so as to perform operations later on

                csi_parent_list_array = np.array(csi_parent_list)
                csi_parent_list_no_weights_array = np.array(csi_parent_list_no_weights)

                csi_background_array = csi_background_array + csi_parent_list_array #made of master + parent
                csi_background_no_weights_array = csi_background_no_weights_array + csi_parent_list_no_weights_array #made of master + parent

                csi_post = csi_background_array + csi_shell_list_array - csi_parent_list_array #add cloud, subtract parent
                csi_pre = csi_background_array

                csi_post_no_weights = csi_background_no_weights_array + csi_shell_list_no_weights_array - csi_parent_list_no_weights_array #add cloud, subtract parent
                csi_pre_no_weights = csi_background_no_weights_array

                csi_post_list.append(csi_post)
                csi_post_list_no_weights.append(csi_post_no_weights)
                csi_pre_list.append(csi_pre)

                ratios = (csi_post - csi_pre)/csi_pre
                ratios_no_weights = (csi_post_no_weights - csi_pre_no_weights)/csi_pre_no_weights

                self.ratios_list.append(ratios) #put them in a list (every element contains the ratio for each altitude at each day)
                self.ratios_no_weights_list.append(ratios_no_weights)

            array_csi_post_pre = np.transpose(np.array([csi_post_list, csi_pre_list]))

            if float(filename[6:9]) <= 100 and float(filename[6:9]) > 1:
                with open(filename, 'w') as fw:
                    for line in array_csi_post_pre:
                        fw.writelines(str(line)[2:-2] + '\n')

            global_csi_cloud_only = np.sum(csi_shell_list_array) #sum of the cumulative csi on all shells (only Frag - Parent)
            global_csi_cloud_only_no_weights = np.sum(csi_shell_list_no_weights_array) #same

            global_csi = np.sum(np.array(csi_post_list)) #sum of the cumulative csi on all shells
            global_csi_no_weights = np.sum(np.array(csi_post_list_no_weights)) #same

            #outputs to be used in plotter function
            self.global_csi_cloud_only_list.append(global_csi_cloud_only)
            self.global_csi_cloud_only_list_no_weights.append(global_csi_cloud_only_no_weights)
            self.global_csi_list.append(global_csi) #append it for each day
            self.global_csi_list_no_weights.append(global_csi_no_weights)
            #day_list.append(float(filename[6:9])) #list of days
            self.day_list.append(float(filename[6:9])) #list of days

        array_cumulative_cloud_csi = np.transpose(np.array([self.day_list,self.global_csi_cloud_only_list]))
        array_cumulative_cloud_csi_no_weights = np.transpose(np.array([self.day_list,self.global_csi_cloud_only_list_no_weights]))
        array_global_csi = np.transpose(np.array([self.day_list, self.global_csi_list]))

        with open(cloud_folder_path + '/radar_array_global_csi_cloud_only' + piece_of_string, 'w') as fw:
            for line in array_cumulative_cloud_csi:
                fw.writelines(str(line)[1:-1] + '\n')

        with open(cloud_folder_path + '/radar_array_global_csi_cloud_only_list_no_weights' + piece_of_string, 'w') as fw:
            for line in array_cumulative_cloud_csi_no_weights:
                fw.writelines(str(line)[1:-1] + '\n')

        with open(cloud_folder_path +  '/radar_array_global_csi' + piece_of_string, 'w') as fw:
            for line in array_global_csi:
                fw.writelines(str(line)[1:-1] + '\n')

        pos_1 = self.day_list.index(1) #at day 1
        pos_100 = self.day_list.index(100) #at day 100

        ratios_100 = self.ratios_list[pos_100]
        ratios_1 = self.ratios_list[pos_1]
        shells = np.arange(250,1250,50)

        array_shells_ratios_100 = np.transpose(np.array([shells,ratios_100]))
        array_shells_ratios_1 = np.transpose(np.array([shells,ratios_1]))

        with open(cloud_folder_path +  '/array_shells_ratios_100_radar' + piece_of_string, 'w') as fw:
            for line in array_shells_ratios_100:
                fw.writelines(str(line)[1:-1] + '\n')
        with open(cloud_folder_path + '/array_shells_ratios_1_radar' + piece_of_string, 'w') as fw:
            for line in array_shells_ratios_1:
                fw.writelines(str(line)[1:-1] + '\n')

        return(self.day_list, self.global_csi_cloud_only_list, self.global_csi_cloud_only_list_no_weights, self.global_csi_list, self.global_csi_list_no_weights, self.ratios_list)


    # ## Optical Main
    #
    # This function is called as a main when a pure optical network is assumed to be in place.
    # In it, all the above defined functions are used, so as to retrieve the needed information and compute the FEI.
    #
    # It takes as inputs:
    # * the cloud file (from MASTER)
    # * the fragmentation altitude $h_{frag}$
    # * the minimum detectable size at $h_{max}$
    # * the maximum altitude at which the use of radar sensors is considered to be effective
    # * a piece of string, created with s_min and h_max, so it is in principle not needed...
    #
    # It gives as outputs:
    # * day_list (list of epoch after fragmentation, in Days)
    # * global_csi_cloud_only_list (weights are considered)
    # * global_csi_cloud_only_list_no_weights (weights are all set to one)
    # * global_csi_list (cloud + background csi, weights are considered)
    # * global_csi_list_no_weights (cloud + background csi, weights are all set to one)
    # * ratios_list (percentage FEI values)

    # In[23]:


    def optical_main(self, nube, h_coll, s_min,h_max, piece_of_string):

        optical_threshold_value = self.optical_threshold(h_coll,s_min,h_max) #minimum detectable size for a given rho_frag
        global_csi_list = []
        global_csi_list_no_weights = []
        global_csi_cloud_only_list = []
        global_csi_cloud_only_list_no_weights = []

        day_list = []
        ratios_list = []
        ratios_no_weights_list = []

        parent_mass = 2000
        parent_a = h_coll + r_e
        parent_e = 0.00003
        parent_inc = 80.3

        csi_background, csi_background_no_weights = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name) + '/csi_out/csi0_optical' + piece_of_string + '.out', unpack = True, usecols = (0,1))
        csi_background_array = np.array(csi_background)
        csi_background_no_weights_array = np.array(csi_background_no_weights)

        omega_j_elem = w_t_sig(h_coll) #associated w_t_sig weight (depending on wether we are in LOW LEO, MED LEO or HIGH LEO)

        if os.path.isdir('/Users/luigigisolfi/' + str(cloud_name) + '/figures/optical' + piece_of_string):
            print('path for figures already exists!')
        else:
            os.mkdir('/Users/luigigisolfi/' + str(cloud_name) + '/figures/optical' + piece_of_string)

        if os.path.isdir('/Users/luigigisolfi/' + str(cloud_name)+ '/weights/optical' + piece_of_string):
            print('path for weights already exists!')
        else:
            os.mkdir('/Users/luigigisolfi/' + str(cloud_name)+ '/weights/optical' + piece_of_string)

        for filename in os.listdir('/Users/luigigisolfi/' + str(cloud_name)):
            f = os.path.join('/Users/luigigisolfi/' + str(cloud_name), filename)
            if (os.path.isdir(f) == True):
                print(str(f) + ' is a directory!')
                continue
            else:
                if (filename[6:9].isnumeric() and float(filename[6:9]) > 100):
                    continue
                elif (f[-4:] == '.fla' and filename != 'cloud_init.fla' and filename != 'fragment.fla'):
                    print(filename)
                    np.set_printoptions(threshold=np.inf)
                    epoch, n, area, mass, a, e, inc, Omega, omega, M = np.loadtxt(f, unpack = True, usecols = (0,1,2,3,4,5,6,7,8,9))
                    f_weights = os.path.join('/Users/luigigisolfi/' + str(cloud_name) + '/weights/optical' + piece_of_string + '/' + filename[:-4] + '_weights.fla')

                else:
                    continue

                sizes = np.sqrt(4*area/np.pi) #area in m^2, sizes in m
                sizes = sizes*100 #in cm
                sizes_cond = sizes >=1
                sensor_cond = a-r_e >=1200
                combined = sizes_cond & sensor_cond
                e_cond = e <= 0.5
                combined = combined & e_cond
                a = a[combined]

                sizes = sizes[combined]
                e = e[combined]
                inc = inc[combined]
                mass = mass[combined]

                omega_i_map = list(map(w_e_rso, sizes,a))
                omega_i_map = np.array(omega_i_map)

                omega_j_map = np.full(shape= len(omega_i_map), fill_value=omega_j_elem,dtype=float)

                weights = list(map(get_omega_optical_sum, omega_i_map, omega_j_map))

                np.set_printoptions(threshold=np.inf)

                array = np.transpose(np.array([sizes, weights])) #array with sizes and associated weights

                with open('/Users/luigigisolfi/' + str(cloud_name)+ '/weights/optical' + piece_of_string + '/' + str(filename[:-4]) + '_weights.fla', 'w') as fw:
                    for line in array:
                        fw.writelines(str(line)[1:-1] + '\n')

                csi_shell_list = []
                csi_shell_list_no_weights = []
                csi_parent_list = []
                csi_parent_list_no_weights = []
                csi_post_list= []
                csi_post_list_no_weights = []
                csi_pre_list = []
                no_weights = np.ones(len(weights))


                for r_in in range(1200 + 6378,2000 + 6378,50):

                    r_out = r_in + 50

                    f_csi_shell = fractional_csi(mass,a,e,inc, weights, r_in, r_out) #all objects weighted fractional csi on a single shell
                    f_csi_shell_array = np.array(f_csi_shell) #make it into array so as to perform operations on it later on

                    f_csi_shell_no_weights =fractional_csi(mass,a,e,inc, no_weights, r_in, r_out)#all objects fractional csi (w=1) on a single shell
                    f_csi_shell_no_weights_array = np.array(f_csi_shell_no_weights)

                    f_csi_shell_parent = 0
                    f_csi_shell_parent_no_weights = parent_fractional_csi(parent_mass,parent_a,parent_e,parent_inc,1, r_in,r_out) #parent object csi (w =1)

                    total_csi_shell = np.sum(f_csi_shell_array)  #computing sum on j of (Post_j-Pre_j) = (Frag_j - Parent_j) on every shell j
                    total_csi_shell_no_weights = np.sum(f_csi_shell_no_weights_array) #same but w = 1

                    csi_shell_list.append(total_csi_shell) #append the result into a list made of each day's results
                    csi_shell_list_no_weights.append(total_csi_shell_no_weights) #same as above with w = 1

                    csi_parent_list.append(f_csi_shell_parent)
                    csi_parent_list_no_weights.append(f_csi_shell_parent_no_weights)

                csi_shell_list_array = np.array(csi_shell_list) #make list into array so as to perform operations later on
                csi_shell_list_no_weights_array = np.array(csi_shell_list_no_weights)  #make list into array so as to perform operations later on

                csi_parent_list_array = np.array(csi_parent_list)
                csi_parent_list_no_weights_array = np.array(csi_parent_list_no_weights)

                csi_background_array = csi_background_array + csi_parent_list_array #made of master + parent
                csi_background_no_weights_array = csi_background_no_weights_array + csi_parent_list_no_weights_array #made of master + parent

                csi_post = csi_background_array + csi_shell_list_array - csi_parent_list_array #add cloud, subtract parent
                csi_pre = csi_background_array

                csi_post_no_weights = csi_background_no_weights_array + csi_shell_list_no_weights_array - csi_parent_list_no_weights_array #add cloud, subtract parent
                csi_pre_no_weights = csi_background_no_weights_array

                csi_post_list.append(csi_post)
                csi_post_list_no_weights.append(csi_post_no_weights)
                csi_pre_list.append(csi_pre)

                ratios = (csi_post - csi_pre)/(csi_pre)

                ratios_no_weights = (csi_post_no_weights - csi_pre_no_weights)/(csi_pre_no_weights)

                ratios_list.append(ratios) #put them in a list (every element contains the ratio for each altitude at each day)
                ratios_no_weights_list.append(ratios_no_weights)

            array_csi_post_pre = np.transpose(np.array([csi_post_list, csi_pre_list]))

            if filename[6:9] == '100' or filename[6:9] == '001':
                with open('/Users/luigigisolfi/' + str(cloud_name)+ '/optical_csi_post_pre_' + filename[6:9] +piece_of_string, 'w') as fw:
                    for line in array_csi_post_pre:
                        fw.writelines(str(line)[2:-2] + '\n')

            global_csi_cloud_only = np.sum(csi_shell_list_array) #sum of the cumulative csi on all shells (only Frag - Parent)
            global_csi_cloud_only_no_weights = np.sum(csi_shell_list_no_weights_array) #same

            global_csi = np.sum(np.array(csi_post_list)) #sum of the cumulative csi on all shells
            global_csi_no_weights = np.sum(np.array(csi_post_list_no_weights)) #same

            #outputs to be used in plotter function
            global_csi_cloud_only_list.append(global_csi_cloud_only)
            global_csi_cloud_only_list_no_weights.append(global_csi_cloud_only_no_weights)
            global_csi_list.append(global_csi) #append it for each day
            global_csi_list_no_weights.append(global_csi_no_weights)
            day_list.append(float(filename[6:9])) #list of days

        array_cumulative_cloud_csi = np.transpose(np.array([day_list,global_csi_cloud_only_list]))
        array_cumulative_cloud_csi_no_weights = np.transpose(np.array([day_list,global_csi_cloud_only_list_no_weights]))
        array_global_csi = np.transpose(np.array([day_list, global_csi_list]))

        with open('/Users/luigigisolfi/' + str(cloud_name)+ '/optical_array_global_csi_cloud_only' + piece_of_string, 'w') as fw:
            for line in array_cumulative_cloud_csi:
                fw.writelines(str(line)[1:-1] + '\n')

        with open('/Users/luigigisolfi/' + str(cloud_name)+ '/optical_array_global_csi_cloud_only_list_no_weights' + piece_of_string, 'w') as fw:
            for line in array_cumulative_cloud_csi_no_weights:
                fw.writelines(str(line)[1:-1] + '\n')

        with open('/Users/luigigisolfi/' + str(cloud_name)+  '/optical_array_global_csi' + piece_of_string, 'w') as fw:
            for line in array_global_csi:
                fw.writelines(str(line)[1:-1] + '\n')


        pos_1 = day_list.index(1) #at day 1
        pos_100 = day_list.index(100) #at day 100

        ratios_100 = ratios_list[pos_100]
        ratios_1 = ratios_list[pos_1]
        shells = np.arange(1250,2050,50)

        array_shells_ratios_100 = np.transpose(np.array([shells,ratios_100]))
        array_shells_ratios_1 = np.transpose(np.array([shells,ratios_1]))

        with open('/Users/luigigisolfi/' + str(cloud_name)+ '/array_shells_ratios_100_optical' + piece_of_string, 'w') as fw:
            for line in array_shells_ratios_100:
                fw.writelines(str(line)[1:-1] + '\n')
        with open('/Users/luigigisolfi/' + str(cloud_name)+ '/array_shells_ratios_1_optical' + piece_of_string, 'w') as fw:
            for line in array_shells_ratios_1:
                fw.writelines(str(line)[1:-1] + '\n')

        return(day_list, global_csi_cloud_only_list, global_csi_cloud_only_list_no_weights, global_csi_list, global_csi_list_no_weights, ratios_list)


# ## Background Population FEI Computation
# As done for the fragmentaiton cloud, we compute the FEI for each of the objects that were present in the atmosphere before the fragmentation event. These constitute the background population.
# As done above, this is computed for optical and/or radar.

# In[25]:


def radar_background(nube, h_coll, s_min,h_max, piece_of_string, background_pop):

    if not os.path.isfile(background_pop):
        print('Could not find background population file. Aborting...')
        exit()
    n, mass, sizes, area, a, e, inc, Omega, omega, M = np.loadtxt(background_pop, unpack = True, usecols = (0,1,2,3,4,5,6,7,8,9))
    sensor_action_range = np.where(a - r_e >= 1200)

    sizes = np.sqrt(4*area/np.pi) #area in m^2, sizes in m
    sizes = sizes*100 #in cm
    sizes_cond = sizes >=1
    sensor_cond = a-r_e < 1200
    #index_sizes = np.where(sizes >= 1)
    combined = sizes_cond & sensor_cond
    e_cond = e <= 0.5
    combined = combined & e_cond
    a = a[combined]
    sizes = sizes[combined]

    weights = list(map(w_radar, sizes,a))

    np.set_printoptions(threshold=np.inf)

    array = np.transpose(np.array([sizes, weights])) #array with sizes and associated weights

    with open('/Users/luigigisolfi/' f'weights_background_radar{piece_of_string}.fla', 'w') as fw:
        for line in array:
            fw.writelines(str(line)[1:-1] + '\n')

    e = e[combined]
    inc = inc[combined]
    mass = mass[combined]

    no_weights = np.ones(len(weights))

    with open('/Users/luigigisolfi/' + str(cloud_name) + '/csi_out' + '/csi0_radar' + piece_of_string + '.out', 'w') as fw_csi0:

        for r_in in range(200 + 6378,1200 + 6378,50):

            r_out = r_in + 50

            f_csi_shell = fractional_csi(mass,a,e,inc,weights, r_in, r_out) #all objects weighted fractional csi on a single shell
            f_csi_shell_array = np.array(f_csi_shell) #make it into array so as to perform operations on it later on

            f_csi_shell_no_weights =fractional_csi(mass,a,e,inc,no_weights, r_in, r_out)#all objects fractional csi (w=1) on a single shell
            f_csi_shell_no_weights_array = np.array(f_csi_shell_no_weights)

            total_csi_shell = np.sum(f_csi_shell_array)  #computing sum on j of (Post_j-Pre_j) = (Frag_j - Parent_j) on every shell j
            total_csi_shell_no_weights = np.sum(f_csi_shell_no_weights_array) #same but w = 1


            fw_csi0.writelines(str(total_csi_shell) + ' ' + str(total_csi_shell_no_weights) + '\n')

def optical_main_background(nube, h_coll, s_min,h_max, piece_of_string, background_pop):

    omega_j_elem = w_t_sig(h_coll) #associated w_t_sig weight (depending on wether we are in LOW LEO, MED LEO or HIGH LEO)

    if not os.path.isfile(background_pop):
        print('Could not find background population file. Aborting...')
        exit()
    n, mass, sizes, area, a, e, inc, Omega, omega, M = np.loadtxt(background_pop, unpack = True, usecols = (0,1,2,3,4,5,6,7,8,9))
    sensor_action_range = np.where(a - r_e >= 1200)

    sizes = np.sqrt(4*area/np.pi) #area in m^2, sizes in m
    sizes = sizes*100 #in cm
    sizes_cond = sizes >=1
    sensor_cond = a-r_e >=1200
    #index_sizes = np.where(sizes >= 1)
    combined = sizes_cond & sensor_cond
    e_cond = e <= 0.5
    combined = combined & e_cond
    a = a[combined]
    sizes = sizes[combined]

    omega_i_list =[] #initialize w_e_rso list
    omega_i_map = list(map(w_e_rso, sizes,a))
    omega_i_map = np.array(omega_i_map)

    omega_j_map = np.full(shape= len(omega_i_map), fill_value=omega_j_elem,dtype=float)

    weights = list(map(get_omega_optical_sum, omega_i_map, omega_j_map))

    np.set_printoptions(threshold=np.inf)

    array = np.transpose(np.array([sizes, weights])) #array with sizes and associated weights

    with open('/Users/luigigisolfi/' f'weights_background_optical{piece_of_string}.fla', 'w') as fw:
        for line in array:
            fw.writelines(str(line)[1:-1] + '\n')

    e = e[combined]
    inc = inc[combined]
    mass = mass[combined]

    no_weights = np.ones(len(weights))

    with open('/Users/luigigisolfi/' + str(cloud_name) + '/csi_out' + f'/csi0_optical{piece_of_string}' + '.out', 'w') as fw_csi0:

        for r_in in range(1200 + 6378,2000 + 6378,50):

            r_out = r_in + 50

            f_csi_shell = fractional_csi(mass,a,e,inc,weights, r_in, r_out) #all objects weighted fractional csi on a single shell
            f_csi_shell_array = np.array(f_csi_shell) #make it into array so as to perform operations on it later on

            f_csi_shell_no_weights =fractional_csi(mass,a,e,inc,no_weights, r_in, r_out)#all objects fractional csi (w=1) on a single shell
            f_csi_shell_no_weights_array = np.array(f_csi_shell_no_weights)

            total_csi_shell = np.sum(f_csi_shell_array)  #computing sum on j of (Post_j-Pre_j) = (Frag_j - Parent_j) on every shell j
            total_csi_shell_no_weights = np.sum(f_csi_shell_no_weights_array) #same but w = 1

            fw_csi0.writelines(str(total_csi_shell) + ' ' + str(total_csi_shell_no_weights) + '\n')


        # ## Visualizing the Results
#
# Some functions are written to properly visualize and make sense of the various results we obtained from the simulations:
#
# 1) plotter_csi
# 2) plotter_FEI
# 3) multi_plotter_csi
# 4) modulated_FEI

# In[27]:


def plotter_csi(network_type, c):

    plt.plot(day_list, global_csi_cloud_only_list, 'o', ms = 3, color = c)
    plt.title(f'Cumulative Cloud csi ({network_type}, {size})')
    plt.xlabel('Days From Collision')
    plt.ylabel('Cumulative Cloud csi')
    plt.tight_layout()
    plt.savefig('/Users/luigigisolfi/' + str(cloud_name) + f'/figures/{network_type}{piece_of_string}' + '/Cumulative_Cloud_csi')
    plt.show()


    plt.plot(day_list, global_csi_cloud_only_list_no_weights, 'o', ms = 3, color = 'blue')
    plt.title('Cumulative Cloud csi (w_tr = 1)')
    plt.xlabel('Days From Collision')
    plt.ylabel('Cumulative Cloud csi')
    plt.tight_layout()
    plt.savefig('/Users/luigigisolfi/' + str(cloud_name) + f'/figures/{network_type}{piece_of_string}' + '/Cumulative_Cloud_csi_no_track')
    plt.show()


    plt.plot(day_list, global_csi_list, 'o', ms = 3, color = c)
    plt.title(f'Cumulative csi  ({network_type}, {size})')
    plt.xlabel('Days From Collision')
    plt.ylabel('Cumulative csi')
    plt.tight_layout()
    plt.savefig('/Users/luigigisolfi/' + str(cloud_name) + f'/figures/{network_type}{piece_of_string}' + '/Cumulative_csi')
    plt.show()


# In[28]:


def plotter_FEI(network_type,c):
    size = piece_of_string.split('_')[1]
    shells, ratios_1 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/array_shells_ratios_1_{network_type}' + piece_of_string, unpack= True, usecols = (0,1))
    shells, ratios_100 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/array_shells_ratios_100_{network_type}' + piece_of_string, unpack= True, usecols = (0,1))

    post_1, pre_1 =  np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_csi_post_pre_001' + piece_of_string, unpack= True, usecols = (0,1))
    post_100, pre_100 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_csi_post_pre_100' + piece_of_string, unpack= True, usecols = (0,1))

    diff_1 = post_1 - pre_1
    diff_100 = post_100 - pre_100

    plt.plot(shells, ratios_1, label = 'Day 1', color = c)
    plt.plot(shells, ratios_100, label = 'Day 100',linestyle = '--', color = c)
    plt.axvline(h_coll,0,linestyle = '-.',color = 'silver',label = 'Collision Altitude')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Percentage FEI')
    plt.legend(loc = 'lower right', prop={'size': 6})
    plt.title(f'Percentage FEI 1 and 100 Days After Collision, ({network_type}, {size})')
    plt.yscale('log')
    plt.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}{piece_of_string}' + '/Percentage_FEI_T0_T100')
    plt.show()

    plt.plot(shells, diff_1, label = 'Day 1', color = c)
    plt.plot(shells, diff_100, label = 'Day 100', linestyle = '--', color = c)
    plt.axvline(h_coll,0,linestyle = '-.',color = 'silver',label = 'Collision Altitude')
    plt.xlabel('Altitude (km)')
    plt.ylabel('csi_post - csi_pre')
    plt.legend(loc = 'lower right', prop={'size': 6})
    plt.title(f'FEI 1 and 100 Days After Collision ({network_type}, {size})')
    plt.yscale('log')
    plt.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}{piece_of_string}' + '/Diff_FEI_T0_T100')
    plt.show()

    plt.plot(shells, diff_1*ratios_1, label = 'Day 1', color = c)
    plt.plot(shells, diff_100*ratios_100, label = 'Day 100',linestyle = '--', color = c)
    plt.axvline(h_coll,0,linestyle = '-.',color = 'silver',label = 'Collision Altitude')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Perc_FEI * Diff')
    plt.legend(loc = 'lower right', prop={'size': 6})
    plt.title(f'Modulated Perc FEI ({network_type}, {size})')
    plt.yscale('log')
    plt.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}{piece_of_string}' + '/Modulated_FEI_T0_T100')
    plt.show()


# In[29]:


def multi_plotter_csi(pieces_of_strings, nube, network_type):

    colors = ['black', 'grey', 'red']
    if len(nube) == 11:
        h_coll = round(float(nube[5:8]))
    elif len(nube) == 12:
        h_coll = round(float(nube[5:9]))

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()

    for piece_of_string, c in zip(pieces_of_strings, colors):
        print(piece_of_string, c)
        day_list, global_csi_cloud_only_list= np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi_cloud_only' + f'{piece_of_string}', unpack= True, usecols = (0,1))
        global_csi_list = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi' + f'{piece_of_string}', unpack= True, usecols = 1)
        cloud_percentage_contribution = 100*global_csi_cloud_only_list/global_csi_list

        if piece_of_string[0:4] == '_5cm':
            ax1.plot(day_list, global_csi_cloud_only_list, 'o', markersize = 3, label = f'{piece_of_string[1:4]} ' + f'{piece_of_string[5:]}' , color = c)
            ax2.plot(day_list, cloud_percentage_contribution, 'o', markersize = 3, label = f'{piece_of_string[1:4]} ' + f'{piece_of_string[5:]}' , color = c)
            ax3.plot(day_list, global_csi_list, 'o', markersize = 3, label = f'{piece_of_string[1:4]} ' + f'{piece_of_string[5:]}' , color = c)

        else:
            ax1.plot(day_list, global_csi_cloud_only_list, 'o', markersize = 3, label = f'{piece_of_string[1:5]} ' + f'{piece_of_string[6:]}' , color = c)
            ax2.plot(day_list, cloud_percentage_contribution, 'o', markersize = 3, label = f'{piece_of_string[1:5]} ' + f'{piece_of_string[6:]}' , color = c)
            ax3.plot(day_list, global_csi_list, 'o', markersize = 3, label = f'{piece_of_string[1:5]} ' + f'{piece_of_string[6:]}' , color = c)

    day_list, global_csi_cloud_only_list_no_weights = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi_cloud_only_list_no_weights' + f'{piece_of_string}', unpack= True, usecols = (0,1))
    ax1.plot(day_list, global_csi_cloud_only_list_no_weights, 'o', markersize = 3, label = 'no weights' , color = 'blue')
    ax1.set(xlabel = 'Time From Collision (Days)', ylabel = 'Cloud csi')
    ax1.legend(loc = 'upper right', prop={'size': 6})
    ax1.set_title("Cloud csi (on different " + f"{network_type}" + " networks)", fontsize = 'small')
    plt.tight_layout()
    fig1.show()
    fig1.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}' + '_Cumulative_Cloud_csi_comparison',  bbox_inches="tight")

    ax2.plot(day_list, global_csi_cloud_only_list_no_weights, 'o', markersize = 3, label = 'no weights' , color = 'blue')
    ax2.set(xlabel = 'Time From Collision (Days)', ylabel = 'Cloud csi Contribution to Global csi (%)')
    ax2.legend(loc = 'upper right', prop={'size': 6})
    ax2.set_title("Cloud's Contribution to Global csi (on different " + f"{network_type}" + " networks)", fontsize = 'small')
    plt.tight_layout()
    fig2.show()
    fig2.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}' + '_Cumulative_Cloud_csi_comparison_perc',  bbox_inches="tight")

    ax3.set(xlabel = 'Time From Collision (Days)', ylabel = 'Global csi')
    ax3.legend(loc = 'upper right', prop={'size': 6})
    ax3.set_title("Global csi (on different " + f"{network_type}" + " networks)", fontsize = 'small')
    plt.tight_layout()
    fig3.show()
    fig3.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}' + '_Global_csi_comparison',  bbox_inches="tight")


    if network_type == 'radar':
        piece_of_string_0 = '_5cm_1200km'
        piece_of_string_1 = '_15cm_1200km'
        day_list_0, global_csi_list_0 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi' + f'{piece_of_string_0}', unpack= True, usecols = (0,1))
        day_list_1, global_csi_list_1 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi' + f'{piece_of_string_1}', unpack= True, usecols = (0,1))
        ratio_performance = global_csi_list_0/global_csi_list_1
        percentage_ratio_performance = (1 - global_csi_list_0/global_csi_list_1)*100

        print(global_csi_list_0)
        print(global_csi_list_1)
        print(ratio_performance)

        fig, axs = plt.subplots(3, 1)

        axs[0].plot(day_list_0, global_csi_list_1,'o', markersize = 3)
        axs[0].set_title(f'Radar 15 cm, Coll. Altitude = {h_coll} km', fontsize = 'small')
        axs[1].plot(day_list_0, global_csi_list_0, 'o', markersize = 3)
        axs[1].set_title(f'Radar 5 cm, Coll. Altitude = {h_coll} km',  fontsize = 'small')
        axs[2].plot(day_list_0, percentage_ratio_performance, 'o', markersize = 3)
        axs[2].set_title('Performance Comparison', fontsize = 'small')


        axs[0].set(xlabel='Time From Collision (Days)', ylabel='Global csi')
        axs[1].set(xlabel='Time From Collision (Days)', ylabel='Global csi')
        axs[2].set(xlabel='Time From Collision (Days)', ylabel='Risk Reduction (%)')

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()

        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)
        fig.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}' + '_Performance_Ratios',  bbox_inches="tight")

    elif network_type == 'optical':
        piece_of_string_0 = '_5cm_2000km'
        piece_of_string_1 = '_20cm_2000km'
        day_list_0, global_csi_list_0 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi' + f'{piece_of_string_0}', unpack= True, usecols = (0,1))
        day_list_1, global_csi_list_1 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi' + f'{piece_of_string_1}', unpack= True, usecols = (0,1))
        ratio_performance = global_csi_list_0/global_csi_list_1
        percentage_ratio_performance = (1 - global_csi_list_0/global_csi_list_1)*100

        print(global_csi_list_0)
        print(global_csi_list_1)
        print(ratio_performance)

        fig, axs = plt.subplots(3, 1)
        axs[0].plot(day_list_0, global_csi_list_1,'o', markersize = 3)
        axs[0].set_title(f'Optical 20 cm, Coll. Altitude = {h_coll} km', fontsize = 'small')
        axs[1].plot(day_list_0, global_csi_list_0, 'o', markersize = 3)
        axs[1].set_title(f'Optical 5 cm, Coll. Altitude = {h_coll} km',  fontsize = 'small')
        axs[2].plot(day_list_0, percentage_ratio_performance, 'o', markersize = 3)
        axs[2].set_title('Performance Comparison', fontsize = 'small')


        axs[0].set(xlabel='Time From Collision (Days)', ylabel='Global csi')
        axs[1].set(xlabel='Time From Collision (Days)', ylabel='Global csi')
        axs[2].set(xlabel='Time From Collision (Days)', ylabel='Risk Reduction (%)')

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()

        plt.ticklabel_format(useOffset=False)
        plt.tight_layout()
        fig.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}' + '_Performance_Ratios',  bbox_inches="tight")


# In[30]:


def modulated_FEI(pieces_of_strings, nube, h_coll, network_type):

    for piece_of_string in pieces_of_strings:
        size = piece_of_string.split('_')[1]
        shells, ratios_1 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/array_shells_ratios_1_{network_type}' + piece_of_string, unpack= True, usecols = (0,1))
        shells, ratios_100 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/array_shells_ratios_100_{network_type}' + piece_of_string, unpack= True, usecols = (0,1))
        csi_post_100_base, csi_pre_100_base = np.loadtxt('/Users/luigigisolfi/' + nube + f'/{network_type}_csi_post_pre_' + '100' + f'{piece_of_string}', unpack = True, usecols = (0,1))
        csi_post_1_base, csi_pre_1_base = np.loadtxt('/Users/luigigisolfi/' + nube + f'/{network_type}_csi_post_pre_' + '001' + f'{piece_of_string}', unpack = True, usecols = (0,1))

        plt.plot(shells, ratios_1*(csi_post_1_base-csi_pre_1_base), label = 'Day 1', color = 'grey')
        plt.plot(shells, ratios_100*(csi_post_1_base-csi_pre_1_base), linestyle = '--', label = 'Day 100', color = 'grey')
        plt.axvline(h_coll,0,linestyle = '-.',color = 'silver',label = 'Collision Altitude')
        plt.title(f'Modulated Perc FEI ({network_type}, {size})')
        plt.xlabel('Altitude (km)')
        plt.ylabel('Modulated Perc FEI')
        plt.legend()
        plt.tight_layout
        plt.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}{piece_of_string}' + '/Modulated_FEI_T0_T100', bbox_inches = 'tight')
        plt.show()

# def multi_plotter_FEI(pieces_of_strings, nube, h_coll, network_type):

#     colors = ['black', 'grey', 'red']
#     for piece_of_string, c in zip(pieces_of_strings, colors):
#         print(piece_of_string, c)
#         shells, ratios_1 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/array_shells_ratios_1_{network_type}' + piece_of_string, unpack= True, usecols = (0,1))
#         shells, ratios_100 = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/array_shells_ratios_100_{network_type}' + piece_of_string, unpack= True, usecols = (0,1))
#         plt.plot(shells, ratios_1, label = f'{piece_of_string[1:5]} ' + f'{piece_of_string[6:]}' , color = c)
#         plt.plot(shells, ratios_100, linestyle = '--', color = c)

#     plt.axvline(h_coll,0,linestyle = '-.',color = 'silver',label = 'Collision Altitude')
#     plt.xlabel('Altitude (km)')
#     plt.ylabel('Percentage FEI')
#     plt.legend(loc = 'lower right', prop={'size': 6})
#     plt.title('Percentage FEI 1 and 100 Days After Collision (different ' + f'{network_type}' + ' networks)')
#     plt.yscale('log')
#     plt.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}' + '_Percentage_FEI_T0_T100_multi_plot')
#     plt.show()

# if network_type == 'optical':
#     csi_post_100_base, csi_pre_100_base = np.loadtxt('/Users/luigigisolfi/' + nube + f'/{network_type}_csi_post_pre_' + '100' + '_20cm_2000km', unpack = True, usecols = (0,1))
#     csi_post_1_base, csi_pre_1_base = np.loadtxt('/Users/luigigisolfi/' + 'nube_1800_km' + f'/{network_type}_csi_post_pre_' + '001' + '_20cm_2000km', unpack = True, usecols = (0,1))
#     days, csi_cloud_base = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi_cloud_only' + '_20cm_2000km', unpack= True, usecols = (0,1))

# elif network_type == 'radar':
#     csi_post_100_base, csi_pre_100_base = np.loadtxt('/Users/luigigisolfi/' + nube + f'/{network_type}_csi_post_pre_' + '100' + '_15cm_1200km', unpack = True, usecols = (0,1))
#     csi_post_1_base, csi_pre_1_base = np.loadtxt('/Users/luigigisolfi/' + nube + f'/{network_type}_csi_post_pre_' + '001' + '_15cm_1200km', unpack = True, usecols = (0,1))
#     days, csi_cloud_base = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi_cloud_only' + '_15cm_1200km', unpack= True, usecols = (0,1))

# pos_100_base = np.where(days == 100)
# pos_1_base = np.where(days == 1)
# csi_cloud_100_base = csi_cloud_base[pos_100_base]
# csi_cloud_1_base = csi_cloud_base[pos_1_base]
# csi_parent_pre_100 = csi_pre_100_base + csi_cloud_100_base - csi_post_100_base
# print(csi_parent_pre_100, ' at 100')
# csi_parent_pre_1 = csi_pre_1_base + csi_cloud_1_base - csi_post_1_base
# print(csi_parent_pre_1, ' at 1')

# for piece_of_string, c in zip(pieces_of_strings, colors):
#     days, csi_cloud = np.loadtxt('/Users/luigigisolfi/' + str(cloud_name)+ f'/{network_type}_array_global_csi_cloud_only' + piece_of_string, unpack= True, usecols = (0,1))
#     pos_100 = np.where(days == 100)
#     pos_1 = np.where(days == 1)
#     csi_cloud_100 = csi_cloud[pos_100]
#     csi_cloud_1 = csi_cloud[pos_1]
#     csi_post_100 = csi_pre_100_base + csi_cloud_100 - csi_parent_pre_100
#     csi_post_1 = csi_pre_1_base + csi_cloud_1 - csi_parent_pre_1

#     print('csi parent 100 is ', csi_parent_pre_100)
#     print('csi parent 1 is ', csi_parent_pre_1)
#     print(f'{piece_of_string} has csi_post_100 as ', csi_post_100)
#     print(f'{piece_of_string} has csi_cloud_100 as ', csi_cloud_100)
#     print(f'{piece_of_string} has csi_post_1 as ', csi_post_1)
#     print(f'{piece_of_string} has csi_cloud_1 as ', csi_cloud_1)

#     ratios_100_base = (csi_post_100 - csi_pre_100_base)/csi_pre_100_base
#     ratios_1_base = (csi_post_1 - csi_pre_1_base)/csi_pre_1_base
#     print(piece_of_string, c)
#     plt.plot(shells, ratios_1_base, label = f'{piece_of_string[1:5]} ' + f'{piece_of_string[6:]}' , color = c)
#     plt.plot(shells, ratios_100_base, linestyle = '--', color = c)

# plt.axvline(h_coll,0,linestyle = '-.',color = 'silver',label = 'Collision Altitude')
# plt.xlabel('Altitude (km)')
# plt.ylabel('Percentage FEI')
# plt.legend(loc = 'lower right', prop={'size': 6})
# plt.title('Percentage FEI 1 and 100 Days After Collision (different ' + f'{network_type}' + ' networks)')
# plt.yscale('log')
# plt.savefig('/Users/luigigisolfi/' + str(cloud_name)+ f'/figures/{network_type}' + '_Percentage_FEI_T0_T100_multi_plot_base')
# plt.show()