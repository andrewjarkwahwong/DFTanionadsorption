#%% READMe
"""
This is a script that calculates equilibrium potentials of specific anion adsorption using the aGC-DFT framework.

The following are obtained in this script:
- Calculates the equilibrium adsorption potential of anions at a given set of EDL properties
- Plots the U0 at a given value of the metal descriptors, which are Work function, O* binding, and OH* binding
- Given the plots, the R^2 is also calculated and the "error bar" based on different EDL parameters
- Plots Electrosorption Valency across metal descriptors and has the same features as above.

The main inputs are: 
1. Excel Sheet (provided a template)
2. Set of EDL properties
3. Valency of the anion 
4. pKa and pH

The following blog post tutorials are provided to inform about the fundamental theory:
Basics of pH, pKa, and equilibrium potential of specific anion adsorption: https://andrewjarkwahwong.github.io/blog/2024/specificanionads/

@author: Andrew Jark-Wah Wong
Contact: ajwongphd@gmail.com
Date: July 27th 2024

"""
#%% Import Packages + Functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate as tb
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator
from sklearn.linear_model import LinearRegression


# Function to calculate capacitance terms
def calculate_capacitance(er_val, d_val, a, e_vac, dm_A, dm_in, diff_dm, u_prime, i):
    C = er_val * a[i] * e_vac / d_val
    diff_dm_sq = (dm_A[i] ** 2 - dm_in[i] ** 2)
    C_0 = -0.5 * (diff_dm_sq) / (C * d_val ** 2)
    C_const_1 = diff_dm[i] / d_val
    C_1 = C_const_1 * u_prime[i]
    c_total = C_0 + C_1
    return c_total, C_0, C_const_1

# Function to calculate dipole-field terms
def calculate_dipole_field(C_0, u_prime, C_const_1, i):
    dm_0 = 2 * C_0
    dm_1 = u_prime[i] * C_const_1
    dm_total = dm_0 + dm_1
    return dm_total

# Function to calculate polarizability terms
def calculate_polarizability(er_val, e_vac, a, d_val, diff_dm_polar_sq, u_prime, diff_a_dm, diff_polar, i):
    p_0_d = 2 * ((er_val * e_vac) ** 2) * (a[i] ** 2) * (d_val ** 2)
    p_0 = diff_dm_polar_sq[i] / p_0_d
    p_1_d = er_val * e_vac * a[i] * d_val ** 2
    p_1 = -1 * u_prime[i] * diff_a_dm[i] / p_1_d
    p_2 = 0.5 * (u_prime[i] ** 2) * (diff_polar[i]) * (1 / d_val) * (1 / d_val)
    p_total = p_0 + p_1 + p_2
    return p_total

# Function to calculate U0 values
def calculate_u0(g_values, u_prime, order):
    TL = np.polyfit(u_prime, g_values, order)
    TL_poly = np.poly1d(TL)
    if order == 1:
        u0 = -TL[1] / TL[0]
        if np.isreal(u0):
            return u0
        else:
            return np.nan
    elif order == 2:
        roots = np.roots(TL_poly)
        valid_roots = [root.real for root in roots if np.isreal(root) and -30 <= root.real <= 10]
        return valid_roots[0] if valid_roots else np.nan

def calc_ev(diff_dm, diff_polar, diff_a_dm, u_eV, er_val, d_val, a, i):
    # Initialize the list to store ev_2c values
 
    
    # Iterate over the list of u values

    ev_2c = ( -1+2 * diff_dm[i] / d_val - diff_a_dm[i] / (er_val * e_vac * a[i] * d_val**2) + (diff_polar[i] * u_eV[i] / d_val**2))

    # Calculate the average of ev_2c values
    ev_2c_avg = np.mean(ev_2c)
    
    return ev_2c_avg
    

#%% Inputs by extracting Excel data

# Input Sheet Name and Path below
sheet = '111' # 110 or 111
path = "C:/Users/dreww/OneDrive/Documents/ajwongphd/OneDrive/Ph.D/Papers in Progress/Acetate_Anion/Acetate_O_OH_Metals.xlsx"

#Code: Imports Data
df = pd.read_excel(path, sheet_name=sheet)

# Initial State:
m = df['metals'].tolist()
G_Bare = df['bare_energy'].tolist()
a = df['area'].tolist()
upzc = df['upzc'].tolist()
wf = df['wf'].tolist()
dm_in = df['dm_bare'].tolist()
polar_bare=df['polar_bare']
#Energy:
G_O = df['G_O'].tolist()
G_OH = df['G_OH'].tolist()
G_A = df['G_Acetate'].tolist()

#Dipole Moments
dm_O = df['dm_O'].tolist()
dm_OH = df['dm_OH'].tolist()
dm_A = df['dm_Acetate'].tolist()


#Polarizability
polar_O = df['polar_O'].tolist()
polar_OH = df['polar_OH'].tolist()
polar_A = df['polar_Acetate'].tolist()


#Constants
G_Acetate = df['Acetate'][0].tolist()
G_AcetaticAcid = df['Hacetate'][0].tolist()
G_H2O = df['H2O'][0].tolist()
G_H2 = df['H2'][0].tolist()

c

# Simple math and rearrangement
G_in =[(G_Bare[i] + G_AcetaticAcid) for i in range(len(G_Bare))]

#Binding Energy

obind = [(G_O[i] - G_Bare[i] - G_H2O  + G_H2) for i in range(len(G_Bare))]
ohbind = [(G_OH[i] - G_Bare[i] - G_H2O  + 0.5*G_H2) for i in range(len(G_Bare))]

#Dipole moment and Polarizability Changes

polar_fin =  [(polar_A[i]-polar_bare[i]) for i in range(len(G_Bare))]
polar_in = [(polar_bare[i]-polar_bare[i]) for i in range(len(G_Bare))]

diff_dm = [i - y for i,y in zip(dm_A,dm_in)]
diff_polar = [i - y for i,y in zip(polar_fin,polar_in)]
diff_dm_polar_sq = [y * x * x - z * a * a for x,y,z,a in zip(dm_A,polar_fin,dm_in,polar_in)]
diff_a_dm = [x * y - z* a for x,y,z,a in zip(dm_A,polar_fin,dm_in,polar_in)]




#%% Math: aGC-DFT to predict U0 values

# Values of pH
pH=0 #pH to calculate shift 
pKa = 3.76 #Reference pKa 
n = 1 #Valency

# Example Input for pH level
if pH  > pKa:
    pH_level = "high"  # "high" or "low" relative to the pKa of interest
else:
    pH_level = "low"  # "high" or "low" relative to the pKa of interest  


#Potential Range of Interest
u_low = -1 #Lower Potential V-SHE that you input
u_high = 2 #Upper Potential V-SHE 
u = np.linspace(u_low,u_high,25)

u_prime = [u - upzc[i] for i in range(len(m))]
u_eV = [u - 0*i for i in range(len(m))]
#EDL Properties

#Inputs 
er = [1,2,8,13,78.4] #Relative permittivity (Dielectric Constant)
d = [3,4.5,6,10] #Helmholtz EDL Width in Angstrom

#Constants
e_vac = 0.00553 #vacuum permittivity
vac_nhe = 4.6 #abs to SHE

# Define lists of values for er and d
er_values = er  # Relative permittivity (Dielectric Constant)
d_values = d  # Helmholtz EDL Width in Angstrom

# Initialize results dictionary
resultsEDL = {}

# Nested loops to iterate over er and d values
for er_val in er_values:
    for d_val in d_values:
        for i in range(len(m)):
            # Calculate g_1a for the current iteration\
            if pH_level== "high":
                G_in =[(G_Bare[i] + G_Acetate) for i in range(len(G_Bare))]
                g_1a = G_A[i] - G_in[i] - n*(upzc[i]+vac_nhe) + g_solv 
                
                # Calculate g_1b for the current iteration
                g_1b = g_1a - np.array(u_prime[i]) 
            elif pH_level== "low":
                
                G_in =[(G_Bare[i] + G_AcetaticAcid) for i in range(len(G_Bare))]
                # Calculate g_1a and g_1b
                g_1a = G_A[i]- G_in[i]   + n*(0.5 * G_H2  -  (wf[i] - vac_nhe) + 0.059 * pH) + g_solv
                g_1b = g_1a - np.array(u_prime[i])
            else:
                raise ValueError("Invalid input. Please use 'high' or 'low' relative to the pKa.")
                
            # Calculate capacitance terms
            c_total, C_0, C_const_1 = calculate_capacitance(er_val, d_val, a, e_vac, dm_A, dm_in, diff_dm, u_prime, i)
            
            # Calculate g_2a
            g_2a = g_1b + c_total

            # Calculate dipole-field terms
            dm_total = calculate_dipole_field(C_0, u_prime, C_const_1, i)
            g_2b = g_2a + dm_total

            # Calculate polarizability termser
            p_total = calculate_polarizability(er_val, e_vac, a, d_val, diff_dm_polar_sq, u_prime, diff_a_dm, diff_polar, i)
            g_2c = g_2b + p_total

            # Calculate EDL total
            EDL_total = p_total + dm_total + c_total

            # Calculate U0 values
            u0_1b = calculate_u0(g_1b, u_prime[i], 1)
            u0_2c = calculate_u0(g_2c, u_prime[i], 2)
            
            # calculate eV
            eV_avg = calc_ev(diff_dm, diff_polar, diff_a_dm, u_eV, er_val, d_val, a, i)

            # Store the result in the dictionary
            resultsEDL[(er_val, d_val), m[i]] = {'wf': wf[i],  'u0_2c': u0_2c, 'u0_1b': u0_1b,'eV':eV_avg,'diffdm':diff_dm[i],'diffpolar':diff_polar[i]}


#%% Plots: U0 vs metal descriptors (WF, O* and OH* binding)
# Create an empty dictionary to store min and max values for each key[1]
minmax_dict = {}

# Iterate over unique m[i] values
for unique_m_value in set(key[1] for key in resultsEDL):
    # Filter values for the current m[i]
    u0_2c_values = [resultsEDL[(er_val, d_val), m_value]['u0_2c'] for (er_val, d_val), m_value in resultsEDL.keys() if m_value == unique_m_value]

    # Convert the list to a numpy array for min and max calculation
    u0_2c_values = np.array(u0_2c_values)

    # Calculate and store min and max values for the current m[i]
    min_u0_2c = np.min(u0_2c_values)
    max_u0_2c = np.max(u0_2c_values)
    minmax_dict[unique_m_value] = {'min_u0_2c': min_u0_2c, 'max_u0_2c': max_u0_2c}

# Sort the minmax_dict based on the order
minmax_dict = {key: minmax_dict[key] for key in m}

# Extract data for plotting
elements = list(minmax_dict.keys())
plot_values = [resultsEDL[(2, 3), element]['u0_2c'] for element in elements] #Plot the value of U0 at er = 2 and d=3, which you can change to whatever value
min_values = [abs(minmax_dict[element]['min_u0_2c'] - plot_values[i]) for i, element in enumerate(elements)]
max_values = [abs(minmax_dict[element]['max_u0_2c'] - plot_values[i]) for i, element in enumerate(elements)]

# Specify Inputs
datasets = [
    ('Work Function (eV)', wf, (round((min(wf)-1),0), round((max(wf)+1),0)), (round((min(plot_values)-1),0), round((max(plot_values)+1),0)), 'R$^2$ = ', 'Max_u0_2c'),
    ('O* Binding Energy (eV)',obind, (round((min(obind)-1),0), round((max(obind)+1),0)), (round((min(plot_values)-1),0), round((max(plot_values)+1),0)), 'R$^2$ = ', 'Max_u0_2c'),
    ('OH* Binding Energy (eV)',ohbind, (round((min(ohbind)-1),0), round((max(ohbind)+1),0)),  (round((min(plot_values)-1),0), round((max(plot_values)+1),0)), 'R$^2$ = ', 'Max_u0_2c')
]

# Make the plot 
colors = sns.color_palette('deep', n_colors=8)
fig, axes = plt.subplots(1, 3, figsize=(28, 8),sharey=True)

#Iterate through data set to calculate 
for ax, (xlabel, x, xlim, ylim, rsquared_text, legend_label) in zip(axes, datasets):
    ax.scatter(x, y=plot_values, label=legend_label, marker='o', color='white', s=1000, zorder=300, edgecolors=colors[0], linewidth=5, alpha=1)
    ax.errorbar(x, y=plot_values, yerr=[min_values, max_values], fmt='none', capsize=10, capthick=5, color=colors[0], elinewidth=5)
    
    #Regression line and R^2 Calculation
    r_squared = np.corrcoef(x, plot_values)[0, 1] ** 2
    print(r_squared)
    slope, intercept = np.polyfit(x, plot_values, 1)
    linear_fit = np.poly1d((slope, intercept))
    x_fit = np.linspace(min(x)-1, max(x)+1, 100)
    y_fit = linear_fit(x_fit)
    
    #Plot R^2 and linear Fit
    ax.plot(x_fit, y_fit, linestyle='-.', alpha=0.15, color='black', linewidth=6, zorder=0)
    ax.text(x=min(x), y=max(plot_values)+0.25, s=rsquared_text + str(round(r_squared, 2)), fontsize=36, color=colors[0])
    
    #Annotate each point to be a metal
    for i, txt in enumerate(elements):
        ax.annotate(txt, (x[i], plot_values[i]), fontsize=16, fontweight='bold', ha='center', va='center', xytext=(0, 0),
                    textcoords='offset points', zorder=1000, color=colors[0])

    ax.set_xlabel(xlabel, fontsize=40)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='minor', length=10, width=4, direction='in')
    ax.tick_params(axis='both', which='major', length=18, width=4, direction='in', labelsize=30, pad=30)
    ax.yaxis.grid(False, which='minor')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.grid(False, which='minor')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_facecolor('white')
    ax.patch.set_edgecolor('black')
    ax.patch.set_linewidth(5)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        
axes[0].set_ylabel(r'U$^0$ (V-SHE)', fontsize=40)

plt.tight_layout()
plt.show()

#%% Plots: eV vs metal descriptors (WF, O* and OH* binding)
# Create an empty dictionary to store min and max values for each key[1]
minmax_dict = {}

# Iterate over unique m[i] values
for unique_m_value in set(key[1] for key in resultsEDL):
    # Filter values for the current m[i]
    u0_2c_values = [resultsEDL[(er_val, d_val), m_value]['eV'] for (er_val, d_val), m_value in resultsEDL.keys() if m_value == unique_m_value]

    # Convert the list to a numpy array for min and max calculation
    u0_2c_values = np.array(u0_2c_values)

    # Calculate and store min and max values for the current m[i]
    min_u0_2c = np.min(u0_2c_values)
    max_u0_2c = np.max(u0_2c_values)
    minmax_dict[unique_m_value] = {'min_eV': min_u0_2c, 'max_eV': max_u0_2c}

# Sort the minmax_dict based on the order
minmax_dict = {key: minmax_dict[key] for key in m}

# Extract data for plotting
elements = list(minmax_dict.keys())
plot_values = [resultsEDL[(2, 3), element]['eV'] for element in elements] #Plot the value of U0 at er = 2 and d=3, which you can change to whatever value
min_values = [abs(minmax_dict[element]['min_eV'] - plot_values[i]) for i, element in enumerate(elements)]
max_values = [abs(minmax_dict[element]['max_eV'] - plot_values[i]) for i, element in enumerate(elements)]

# Specify Inputs
datasets = [
    ('Work Function (eV)', wf, (round((min(wf)-1),0), round((max(wf)+1),0)), (round((min(plot_values)-1),0), round((max(plot_values)+1),0)), 'R$^2$ = ', 'Max_u0_2c'),
    ('O* Binding Energy (eV)',obind, (round((min(obind)-1),0), round((max(obind)+1),0)), (round((min(plot_values)-1),0), round((max(plot_values)+1),0)), 'R$^2$ = ', 'Max_u0_2c'),
    ('OH* Binding Energy (eV)',ohbind, (round((min(ohbind)-1),0), round((max(ohbind)+1),0)),  (round((min(plot_values)-1),0), round((max(plot_values)+1),0)), 'R$^2$ = ', 'Max_u0_2c')
]

# Make the plot 
colors = sns.color_palette('deep', n_colors=8)
fig, axes = plt.subplots(1, 3, figsize=(28, 8),sharey=True)

#Iterate through data set to calculate 
for ax, (xlabel, x, xlim, ylim, rsquared_text, legend_label) in zip(axes, datasets):
    ax.scatter(x, y=plot_values, label=legend_label, marker='o', color='white', s=1000, zorder=300, edgecolors=colors[0], linewidth=5, alpha=1)
    ax.errorbar(x, y=plot_values, yerr=[min_values, max_values], fmt='none', capsize=10, capthick=5, color=colors[0], elinewidth=5)
    ax.hlines(xmin=xlim[0],xmax=xlim[1],y=-1,linewidth=5,alpha=0.5,linestyle='--',color='black')
    #Regression line and R^2 Calculation
    r_squared = np.corrcoef(x, plot_values)[0, 1] ** 2
    print(r_squared)
    slope, intercept = np.polyfit(x, plot_values, 1)
    linear_fit = np.poly1d((slope, intercept))
    x_fit = np.linspace(min(x)-1, max(x)+1, 100)
    y_fit = linear_fit(x_fit)
    
    #Plot R^2 and linear Fit
    ax.plot(x_fit, y_fit, linestyle='-.', alpha=0.15, color='black', linewidth=6, zorder=0)
    ax.text(x=min(x), y=max(plot_values)+0.1, s=rsquared_text + str(round(r_squared, 2)), fontsize=36, color=colors[0])
    
    #Annotate each point to be a metal
    for i, txt in enumerate(elements):
        ax.annotate(txt, (x[i], plot_values[i]), fontsize=16, fontweight='bold', ha='center', va='center', xytext=(0, 0),
                    textcoords='offset points', zorder=1000, color=colors[0])

    ax.set_xlabel(xlabel, fontsize=40)
    ax.set_xlim(xlim)
    ax.set_ylim(-1.4,-0.6)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='minor', length=10, width=4, direction='in')
    ax.tick_params(axis='both', which='major', length=18, width=4, direction='in', labelsize=30, pad=30)
    ax.yaxis.grid(False, which='minor')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.grid(False, which='minor')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_facecolor('white')
    ax.patch.set_edgecolor('black')
    ax.patch.set_linewidth(5)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        
axes[0].set_ylabel(r'Electrosorption Valency', fontsize=36)

plt.tight_layout()
plt.show()

