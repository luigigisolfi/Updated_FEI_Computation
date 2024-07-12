import numpy as np
import matplotlib.pyplot as plt  
import os
import math

################################################################################################
# THIS CODE COMPUTES THE J2 INDICATOR FOR DECISION MAKERS AS DESCRIBED IN:
#  "Quantifying the Effect of an In-Orbit Fragmentation of Space Objects: 
#  An Environmental Index for the Space Surveillance and Tracking (SST) Services"
# Author: Luigi Gisolfi, Year: 2024
################################################################################################

################################################################################################
# PLEASE NOTE:
# for the Inc Rel formula, we assume:
# Omega_debris = mean_Omega and i_debris = mean_I
################################################################################################

################################################################################################
# INPUT PARAMETERS: 
# - Cloud's Fragments Orbital Elements file, of the type "nube_xxx_km", xxx = h_collision
# - Asset Altitude 
# - Asset i, Ω, ω, M (Degrees)
#
# OUTPUTS:
# - J2 Indicator Value
# - Histogram Plot of Potentially Involved Fragments 
# - Fragments Distribution and Satellite Orbit 3D Plot
# PLEASE NOTE:
# Every output is given at:
# 3 different epochs (1, 10 and 100 Days from the fragmentation epoch) and
# for 3 different asset's eccentricity values (0.001,0.01,0.1);
# --> Maximum Number of Output Images: 9 (can be less than 9 if some computations are skipped)
################################################################################################

def J2_histogram(nube, a_satellite, e_satellite, f):
    np.set_printoptions(threshold=np.inf)
    epoch, n, area, mass, a, e, inc, Omega, omega, M = np.loadtxt(f, unpack = True, usecols = (0,1,2,3,4,5,6,7,8,9))

    #epoch, n, area, mass, a, e, inc, Omega, omega, M = np.loadtxt(f, unpack = True, usecols = (0,1,2,3,5,6,7,8,9,10))
    mean_inc = np.radians(np.mean(inc))
    mean_Omega = np.radians(np.mean(Omega))

    h = a - r_e

    combined = np.where((h_perigee_satellite < h) & (h < h_apogee_satellite))
    #print(list(combined))

    sizes = np.sqrt(4*area/np.pi) #area in m^2, sizes in m
    sizes = sizes*100 #in cm

    h_involved = h[combined]
    sizes = sizes[combined]
    orbital_elements = np.array([a[combined],e[combined],np.radians(inc[combined]),np.radians(Omega[combined]),np.radians(omega[combined]),np.radians(M[combined])])
    #print(orbital_elements)
    sizes_10cm = sizes>=10   
    sizes_5cm = sizes >=5 
    print('h_involved is ', len(h_involved))
    h_involved_10cm = h_involved[sizes_10cm]
    h_involved_5cm = h_involved[sizes_5cm]
    print('h_involved_5cm is ', len(h_involved_5cm))
    print('h_involved_10cm is ', len(h_involved_10cm))
    count, edges, plot = plt.hist(h, bins=np.linspace(200,2000,num = 36), histtype='step' )
    count_involved, edges_involved, plot_involved = plt.hist(h_involved, bins=np.linspace(200,2000,num = 36), color= 'cornflowerblue', label = 'All Involved')
    count_involved_10cm, edges_involved_10cm, plot_involved_10cm = plt.hist(h_involved_10cm, bins=np.linspace(200,2000,num = 36), color= 'red', alpha = 0.7, label = '> 10 cm')
    count_involved_5cm, edges_involved_10cm, plot_involved_5cm = plt.hist(h_involved_5cm, bins=np.linspace(200,2000,num = 36), color= 'black', alpha = 0.5, label = '> 5 cm')
    percentage_of_fragments = np.sum(count_involved)/np.sum(count)

    print('count involved', count_involved)
    print('count', count)
    if a_satellite >= 10000:
        mock_ones = np.where(h<10000)
        h= h[mock_ones]
    print(f'Percentage of fragments: {percentage_of_fragments} %')
    print('=============================')
    h_coll = nube.split('_')[1]
    plt.title('Involved Fragments, $h_{coll}$ =' f'{h_coll} km,' + f' Day {day}')
    plt.xlabel(f'Altitude Shell')
    plt.ylabel(f'Number of Objects')
    plt.legend()
    plt.savefig(f'/Users/luigigisolfi/{nube}/percentage_{nube}_a_{a_satellite}_e_{e_sat}_day_{day}.png', dpi=300)
    #plt.show()
    plt.close()
    return(h_involved, percentage_of_fragments, mean_inc, mean_Omega, orbital_elements,h)


def get_phi_new(a,e, r_in, r_out):
    
    E_out_par = np.arccos((a - r_out)/(a*e))
    E_in_par = np.arccos((a - r_in)/(a*e))
    
    peri_par = a*(1-e)
    apo_par = a*(1+e)
    
    cond0_peri_par = peri_par > r_out
    cond0_apo_par = apo_par < r_in
    condizione1_peri = peri_par > r_in
    condizione1_apo = apo_par < r_out
    condizione2_peri = peri_par < r_in
    condizione2_apo = apo_par > r_out
    condizione3_peri = peri_par < r_in
    condizione31_apo =  apo_par < r_out
    condizione32_apo = apo_par > r_in
    condizione41_peri = peri_par < r_out
    condizione42_peri = peri_par > r_in
    condizione4_apo = apo_par > r_out
    
    if (cond0_peri_par|cond0_apo_par): 
        phi = 0
        return(phi) 
    
    if (condizione1_peri and condizione1_apo):
        phi = 1
        return(phi)
    elif (condizione2_peri and condizione2_apo):
        phi = (E_out_par - E_in_par -e *(np.sin(E_out_par) - np.sin(E_in_par)))/np.pi
        return(phi)
    elif (condizione3_peri and (condizione31_apo and condizione32_apo)):
        phi = 1 - (E_in_par - e*(np.sin(E_in_par)))/np.pi
        return(phi)
    elif ((condizione41_peri and condizione42_peri) and condizione4_apo):
        phi = (E_out_par - e*(np.sin(E_out_par)))/np.pi
        return(phi)
    else: 
        print('i dont really know what to print')
        phi = phi
        return(phi)

# Convert elliptical orbital elements with mean anomaly to Cartesian coordinates
def elliptical_elements_with_mean_anomaly_to_cartesian(a, e, i, Ω, ω, M):
    # Convert mean anomaly to eccentric anomaly using Kepler's equation
    E = M
    while True:
        E_next = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        if np.abs(E_next - E) < 1e-8:
            break
        E = E_next

    # Calculate true anomaly from eccentric anomaly
    ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    # Elliptical orbital elements to Cartesian coordinates conversion
    r = a * (1 - e**2) / (1 + e * np.cos(ν))
    x = r * (np.cos(Ω) * np.cos(ω + ν) - np.sin(Ω) * np.sin(ω + ν) * np.cos(i))
    y = r * (np.sin(Ω) * np.cos(ω + ν) + np.cos(Ω) * np.sin(ω + ν) * np.cos(i))
    z = r * (np.sin(i)) * np.sin(ω + ν)
    return x, y, z

def get_Tb(mean_delta_v, mean_inc, h_coll):
    #parent_omega = float(input('Please, input the parent object omega in degrees.'))
    J_two = 0.001082
    beta_omega = np.arctan(5*np.sin(2*mean_inc)*np.cos(np.radians(345))/(14*(2 - 2.5*(np.sin(mean_inc))**2)))
    beta_Omega = np.arctan(1/7*(np.tan(mean_inc)*np.cos(np.radians(345))))

    delta_omega_dot = -1.5*J_two*(r_e**2/(h_coll + r_e)**3)*(7*(2 - 2.5*(np.sin(mean_inc))**2)*np.cos(beta_omega)+ 2.5*np.sin(2*mean_inc)*np.sin(beta_omega)*np.cos(np.radians(345)))*mean_delta_v
    delta_Omega_dot = -1.5*J_two*(r_e**2/(h_coll + r_e)**3)*(7*np.cos(mean_inc)*np.cos(beta_Omega) + np.sin(mean_inc)*np.sin(beta_Omega)*np.cos(np.radians(345)))*mean_delta_v
    
    T_omega = np.pi/np.abs(delta_omega_dot)
    T_Omega = np.pi/np.abs(delta_Omega_dot)
    
    print('beta omega', beta_omega, 'beta Omega', beta_Omega)
    print(f'T omega Ash at {np.degrees(mean_inc)} is {T_omega/86400}, T Omega Ash is,  {T_Omega/86400}')

    return(max(T_omega, T_Omega))

def eccentric_anomaly_from_mean_anomaly(M, e):
    # Solving Kepler's equation for eccentric anomaly (E) using Newton's method
    E_guess = M  # Initial guess
    tol = 1e-8  # Tolerance for convergence
    max_iter = 100  # Maximum number of iterations

    for _ in range(max_iter):
        E_next = E_guess - (E_guess - e * np.sin(E_guess) - M) / (1 - e * np.cos(E_guess))
        if abs(E_next - E_guess) < tol:
            break
        E_guess = E_next

    return E_next

def orbit(a, e, i, Omega, omega, M, mu=398600.4418):
    n = np.sqrt(mu / a**3)
    t = np.linspace(0, 2 * np.pi / n, 1000)

    M = n * t
    E_list = []
    for m in M: 
        E = eccentric_anomaly_from_mean_anomaly(m, e)
    
        E_list.append(E)

    x = a * (np.cos(E_list[:]) - e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(np.array(E_list)[:] / 2), np.sqrt(1 - e) * np.cos(np.array(E_list)[:] / 2))

    r = a * (1 - e**2) / (1 + e * np.cos(nu))
    x = r * (np.cos(Omega) * np.cos(omega + nu) - np.sin(Omega) * np.sin(omega + nu) * np.cos(i))
    y = r * (np.sin(Omega) * np.cos(omega + nu) + np.cos(Omega) * np.sin(omega + nu) * np.cos(i))
    z = r * (np.sin(i)) * np.sin(omega + nu)


    return x, y, z

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else np.zeros(len(d))
    return data[s<m]

# Function to create a 3D sphere with inclination
def plot_inclined_sphere(ax, radius=6378, inclination=0, color='b', alpha=0.8):
    phi, theta = np.mgrid[0.0:2.0*np.pi:1000j, 0.0:np.pi:500j]
    
    # Tilt the sphere by the specified inclination angle
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    
    # Rotate the sphere around the x-axis to simulate inclination
    x, y, z = rotate_x(x, y, z, inclination)
    
    ax.plot_surface(x, y, z, color=color, alpha=alpha, antialiased=False)

# Function to rotate coordinates around the x-axis
def rotate_x(x, y, z, angle):
    angle_rad = np.radians(angle)
    new_y = y * np.cos(angle_rad) - z * np.sin(angle_rad)
    new_z = y * np.sin(angle_rad) + z * np.cos(angle_rad)
    return x, new_y, new_z

def write_file(mass_proj, coll_vel, mass_parent, a_par, e_par, i_par, omega_par, Omega_par, M_par, h_coll, mean_inc, mean_delta_v, h_sat, e_sat, i_sat, omega_sat, Omega_sat, M_sat, percentage_involved, percentage_period, I_weight, J2_indicator, T_b):
    if os.path.isfile(f'/Users/luigigisolfi/{nube}/J2_indicator_report_temp.txt'):                
        with open(f'/Users/luigigisolfi/{nube}/J2_indicator_report_temp.txt', 'a') as J2_report:
            a_sat = h_sat + r_e
            J2_report.write(f'Asset Object Orbital Elements at Fragmentation Epoch: (a,e,i,ω,Ω,M) = {a_sat,e_sat, np.degrees(i_sat), np.degrees(omega_sat), np.degrees(Omega_sat), np.degrees(M_sat)}\n')
            J2_report.write(f'\n')
            J2_report.write(f'Percentage of Produced Fragments Involved: {percentage_involved*100} %\n')
            J2_report.write(f'Percentage of Asset Orbital Period Involved: {percentage_period*100} %\n')
            J2_report.write(f'Contribution Due to PFO-Asset Relative Inclination: {I_weight}\n')    
            J2_report.write(f'Estimated Spreading Timescale: {T_b/(86400)}\n')           
            J2_report.write(f'J2 Indicator Value: {J2_indicator}\n\n')
            J2_report.write(f'============================================================================================\n\n')       

    else: 
        with open(f'/Users/luigigisolfi/{nube}/J2_indicator_report_temp.txt', 'w') as J2_report:
            a_sat = h_sat + r_e
            J2_report.write(f'================================= Details and Outputs =====================================\n\n')
            J2_report.write(f'                  Results Refer to a Time t = 1 Day After the Fragmentation Event           \n')
            J2_report.write(f'              Angles: Degrees, Distances: Km, Masses: Kg, Velocities: km/s, Times: Days     \n\n')    
            J2_report.write(f'============================================================================================\n\n')                
            J2_report.write(f'Projectile Mass: {mass_proj}\n')           
            J2_report.write(f'Fragmented Object Mass: {mass_parent}\n')
            J2_report.write(f'Collisional Velocity (km/s): {coll_vel}\n')
            J2_report.write(f'Fragmented Object Orbital Elements (a,e,i,ω,Ω,M): {a_par,e_par,  np.degrees(i_par), int(h_coll),  np.degrees(omega_par),  np.degrees(Omega_par),  np.degrees(M_par)}\n')
            J2_report.write(f'\n')
            J2_report.write(f'Collison Altitude (km): {h_coll}\n')
            J2_report.write(f'Mean Fragments Inclination : {np.degrees(mean_inc)}\n')    
            J2_report.write(f'Mean Fragments Delta V: {mean_delta_v}\n')      
            J2_report.write(f'\n\n')
            J2_report.write(f'============================================================================================\n\n')         
            J2_report.write(f'Asset Object Orbital Elements at Fragmentation Epoch: (a,e,i,ω,Ω,M) = {a_sat,e_sat, np.degrees(i_sat), np.degrees(omega_sat), np.degrees(Omega_sat), np.degrees(M_sat)}\n')
            J2_report.write(f'\n')
            J2_report.write(f'Percentage of Produced Fragments Involved: {percentage_involved*100} %\n')
            J2_report.write(f'Percentage of Asset Orbital Period Involved: {percentage_period*100} %\n')
            J2_report.write(f'Contribution Due to PFO-Asset Relative Inclination: {I_weight}\n')
            J2_report.write(f'Estimated Spreading Timescale: {T_b/(86400)}\n')           
            J2_report.write(f'J2 Indicator Value: {J2_indicator}\n\n')
            J2_report.write(f'============================================================================================\n\n')                


if __name__ == '__main__':

    r_e = 6378
    nube = input('Nube? ')
    h_coll = nube.split('_')[1]
    h_satellite = float(input('Satellite Altitude?'))
    i_satellite = float(input('Satellite Inclination (Degrees)?'))
    i_satellite = np.radians(i_satellite)
    Omega_satellite = float(input('Satellite Omega (Degrees)?'))
    Omega_satellite = np.radians(Omega_satellite)
    omega_satellite = float(input('Satellite omega (Degrees)?'))
    omega_satellite = np.radians(omega_satellite)
    M_satellite = float(input('Satellite M (Degrees)?'))
    M_satellite = np.radians(M_satellite)

    f1 = '/Users/luigigisolfi/' + str(nube) + '/cloud_001.fla'
    f2 = '/Users/luigigisolfi/' + str(nube) + '/cloud_050.fla'
    f3 =  '/Users/luigigisolfi/' + str(nube) + '/cloud_100.fla'
    f4 = '/Users/luigigisolfi/' + str(nube) + '/cloud_415.fla'
    a_satellite = r_e + h_satellite
    
    e_satellite_list = [0.001,0.01,0.1]

    for e_sat in e_satellite_list:

        perigee_satellite = a_satellite*(1- e_sat)
        print('satellite perigee is ', perigee_satellite)
        apogee_satellite = a_satellite*(1 + e_sat)

        h_perigee_satellite = perigee_satellite - r_e
        h_apogee_satellite = apogee_satellite - r_e

        print('Satellite Altitude is ', h_satellite)
        print(f'Assuming e_sat = {e_sat}')
        if (h_perigee_satellite < 100) or (h_apogee_satellite <= 100):
            print(f'For e_sat = {e_sat}, we get Perigee = {h_perigee_satellite}, Apogee = {h_apogee_satellite}. One of them is negative. Skipping... \n')
            continue

        for f in [f1,f2,f3]:
            if f==f1:
                day = 1
            elif f==f2:
                day = 50
            elif f==f3:
                day = 100
            elif f==f4:
                day = 488
            print('==========================================================')
            print('                 RUNNING FILE: ', f, '...                \n')
            h_involved, percentage_of_fragments, mean_inc, mean_Omega, orbital_elements, h = J2_histogram(nube, a_satellite, e_sat, f)
            delta_v = np.loadtxt(f'/Users/luigigisolfi/{nube}/delta_v_{h_coll}_km.txt', usecols = (-1), skiprows= 1,  unpack = True)
            delta_v = reject_outliers(delta_v)
            mean_delta_v = np.mean(delta_v)
            print('average fragments delta v is ', mean_delta_v)
            print('average fragments inclination is', np.degrees(mean_inc))
            T_b_mean_inc = get_Tb(mean_delta_v, mean_inc, int(h_coll))
            T_b_90 = get_Tb(mean_delta_v, np.pi/2, int(h_coll))
            T_b_0 = get_Tb(mean_delta_v, 0, int(h_coll))

            #print(f'T_b is {T_b_mean_inc}, {T_b_90}, {T_b_0}')
            print(f'T_b in days: {T_b_mean_inc/(86400)}, {T_b_90/(86400)}, {T_b_0/(86400)}')

            if len(h_involved) < 1:
                h_involved = h
                percentage_orbital_period = 0

            else: 
                r_in = np.min(h_involved) + r_e
                r_out =np.max(h_involved) + r_e

                percentage_op = get_phi_new(a_satellite, e_sat, r_in, r_out)
                percentage_orbital_period = percentage_op


            print(f'Percentage of satellite Orbital Period: {percentage_orbital_period*100} %')
            square_root = np.sqrt(np.sin(0.5*(mean_inc - i_satellite))**2 + np.sin(i_satellite)*np.sin(mean_inc)*(np.sin(0.5*(mean_Omega - Omega_satellite)))**2) 
            #print('Square Root Inclination part', np.sin(0.5*(mean_inc - i_satellite))**2)
            #print('Square Root Omega part',np.sin(0.5*(mean_Omega - Omega_satellite))**2)

            I = 2*math.asin(square_root)
            
            if abs(I) > np.pi:
                I = I - np.pi

            I_weight = I/np.pi
            print('I weight: ', I_weight)
            indicator = (percentage_orbital_period) * (percentage_of_fragments) * I_weight * T_b_mean_inc/T_b_90

            print('J2 Indicator:', indicator )

            mass_par = 2000 #kg
            a_par = int(h_coll) + r_e
            e_par = 0.00003
            inc_par = np.radians(80.3)  
            Omega_par = np.radians(24)
            omega_par = np.radians(345)
            M_par = np.radians(32)

            mass_proj = 15 #kg 
            coll_vel = 10 #km/s

            if f == f1:
                write_file(mass_proj, coll_vel, mass_par, a_par, e_par, inc_par, omega_par, Omega_par, M_par, h_coll, mean_inc, mean_delta_v, h_satellite, e_sat, i_satellite, omega_satellite, Omega_satellite, M_satellite, percentage_of_fragments, percentage_orbital_period, I_weight, indicator, T_b_mean_inc)
            # Create a 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plot the Earth's sphere with inclination
            inclination_angle = 0
            plot_inclined_sphere(ax, radius=6378, inclination=inclination_angle, color='skyblue', alpha=0.2)

            # Plot the equator (circle in the xy-plane)
            theta_eq = np.linspace(0, 2 * np.pi, 100)
            x_eq = 6378*np.cos(theta_eq)
            y_eq = 6378*np.sin(theta_eq)
            z_eq = np.zeros_like(theta_eq)

            # Rotate the equator to match the inclination
            x_eq, y_eq, z_eq = rotate_x(x_eq, y_eq, z_eq, inclination_angle)

                # Calculate coordinates for the meridian
            phi_green = np.linspace(0, 2 * np.pi, 100)
            theta_green = np.linspace(-np.pi / 2, np.pi / 2, 100)
            
            # Rotate the meridian to match the inclination
            y_green = 6378* np.cos(phi_green)
            x_green = np.zeros_like(theta_green)
            z_green = 6378*np.sin(phi_green)

            x_green, y_green, z_green = rotate_x(x_green, y_green, z_green, inclination_angle)

            ax.plot(x_eq, y_eq, z_eq, color='r', linestyle = '--', label= "Equator",alpha = 0.7)
            ax.plot(x_green,y_green,z_green, color = 'black', linestyle = '--', label = "Greenwich Meridian", alpha = 0.7)

            orbital_elements = [tuple(row) for row in np.transpose(orbital_elements)]
            # Plot each object's elliptical orbit with mean anomaly
            altitude_ticks = []
            for elements in orbital_elements:
                x, y, z = elliptical_elements_with_mean_anomaly_to_cartesian(*elements)
                ax.plot([x], [y], [z], marker='o', markersize=0.5, color = 'grey')

            x_sat, y_sat, z_sat = elliptical_elements_with_mean_anomaly_to_cartesian(a_satellite, e_sat, i_satellite, Omega_satellite, omega_satellite, M_satellite)
            ax.plot([x_sat], [y_sat], [z_sat], marker='o', markersize=6, color = 'dodgerblue', alpha = 0.5, label = 'Asset (x,y,z) at Fragmentation Epoch')

            x_eci, y_eci, z_eci = orbit(a_satellite,e_sat,i_satellite,Omega_satellite,omega_satellite,M_satellite)
            ax.plot(x_eci, y_eci, z_eci, label = 'Asset Orbit', color = 'dodgerblue')
            # Set labels
            # Set equal scaling for all 3 axes
            ax.set_box_aspect([1, 1, 1])
            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.set_zlabel('Z (km)')
            ax.set_aspect('equal')
            ax.set_title("Cloud's 3D Plot $h_{coll}$" + f'= {h_coll} km' + f', Day {day}')
            
            # Add legend
            # Show the plot
            ax.legend(loc='center left', bbox_to_anchor=(0.8, 0.9), fontsize=5)
            plt.savefig(f'/Users/luigigisolfi/{nube}/3D_{nube}_h_{a_satellite - r_e}_e_{e_sat}_day_{day}.png', dpi=300, bbox_inches='tight')
            #plt.show()
            plt.close()

    os.system('mv ' + f'/Users/luigigisolfi/{nube}/J2_indicator_report_temp.txt ' + f'/Users/luigigisolfi/{nube}/J2_indicator_report_{h_satellite}_km_{i_satellite}.txt')
    print('Done.')
    print('==========================================================')

    exit()



