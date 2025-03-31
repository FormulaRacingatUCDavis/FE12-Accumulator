# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 23:32:43 2025

@author: laptop
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import pandas as pd
from scipy.interpolate import interp1d

data_raw = pd.read_csv('FE11 Endurance Full Data V2.csv',header = 1)
data = np.array(data_raw,dtype=float)
time = data[:,0]
voltage = data[:,1]
current = data[:,2]
temp = data[:,3]
SOC_data = data[:,6]

for i in range(len(current)):
    if current[i] < 0:
        current[i] = 0

plt.plot(time/60,current)
plt.plot(time/60,temp)
# plt.show()

D_cell = 21.55e-3; # [m] Cell diameter
R_cell = D_cell/2; # [m] Cell radius
L_cell = 70.15e-3; # [m] Cell length/height
SA_cell = 0.5 * (L_cell-12.7e-3)*2*np.pi*R_cell; # [m^2] Surface Area of cell, length minus total cell holders 

# R = 0.011 #[Ohms]
m = 0.07 # [kg]
# htc = 7 # [W/m2 K]
Cp_cell = 1360 # [J/kg K]


#htc natural convection
beta = 3.245e-3
g = 9.8
Tinf = temp[0]
alpha = 2.34e-5
kin_vis = 1.655e-5
k_air = 0.027

temp_range = np.linspace(35,55,100)
def htc_natural(x):
    htc_array = np.zeros(len(temp_range))
    for i in range(len(x)):
        htc_array[i] = 0.5 * ( (g * beta * (x[i] - Tinf) * D_cell**3) / (kin_vis * alpha) )**(1/4) * k_air / D_cell
    return htc_array

def interpolate_htc(x, target):
    linear_interp = interp1d(x,htc_natural(x),kind='linear', fill_value='extrapolate')
    return linear_interp(target)


# # Air Properties at 35 C
rho = 1.1459 # [kg/m^3]
# Cp = 1006.7 # [J/kg K]
# mu = 18.915e-6 # [N s/m^2]
# k = 26.71e-3 # [W/m K] 
# Pr = mu*Cp/k

# External htc
#Single Cylinder in Crossflow configs
# C= 0.683
# m_constant = 0.466
# V_max = 3.5    #[m/s]
# s_t = 24e-3    #[m]
# k_cell = 2.21  #[W /m K] radial thermal conductivity

# D = s_t-D_cell
# Re = rho*V_max*D_cell/mu
# Nu = C*(Re**m_constant) * (Pr**(1/3))
# # htc_external = Nu*k_cell/D_cell


# Internal Resistance
SOC = np.linspace(10,100,10)
IRdata_raw = pd.read_csv('IR_data.csv')
IR_data = np.array(IRdata_raw,dtype=float)
IR_P42A_298 = IR_data[:,0]
IR_P42A_303 = IR_data[:,1]
IR_P42A_313 = IR_data[:,2]
IR_P45B_303_2C = IR_data[:,3]
IR_P45B_308_2C = IR_data[:,4]
IR_P45B_313_2C = IR_data[:,5]
IR_P45B_318_2C = IR_data[:,6]
IR_P45B_323_2C = IR_data[:,7]
IR_P45B_303_4C = IR_data[:,3]
IR_P45B_308_4C = IR_data[:,4]
IR_P45B_313_4C = IR_data[:,5]
IR_P45B_318_4C = IR_data[:,6]
IR_P45B_323_4C = IR_data[:,7]

interp_IR_P42A_298 = interp1d(SOC,IR_P42A_298,kind='cubic', fill_value='extrapolate')
interp_IR_P42A_303 = interp1d(SOC,IR_P42A_303,kind='cubic', fill_value='extrapolate')
interp_IR_P42A_313 = interp1d(SOC,IR_P42A_313,kind='cubic', fill_value='extrapolate')

interp_IR_P45B_303_2C = interp1d(SOC,IR_P45B_303_2C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_308_2C = interp1d(SOC,IR_P45B_308_2C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_313_2C = interp1d(SOC,IR_P45B_313_2C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_318_2C = interp1d(SOC,IR_P45B_318_2C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_323_2C = interp1d(SOC,IR_P45B_323_2C,kind='cubic', fill_value='extrapolate')

interp_IR_P45B_303_4C = interp1d(SOC,IR_P45B_303_4C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_308_4C = interp1d(SOC,IR_P45B_308_4C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_313_4C = interp1d(SOC,IR_P45B_313_4C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_318_4C = interp1d(SOC,IR_P45B_318_4C,kind='cubic', fill_value='extrapolate')
interp_IR_P45B_323_4C = interp1d(SOC,IR_P45B_323_4C,kind='cubic', fill_value='extrapolate')

SOC_new = np.linspace(SOC.min(), SOC.max(), 100)

# plt.plot(SOC_new,interp_IR_P42A_298(SOC_new),label = '25 C')
# plt.plot(SOC_new,interp_IR_P42A_303(SOC_new),label = '30 C')
# plt.plot(SOC_new,interp_IR_P42A_313(SOC_new),label = '35 C')
# plt.legend()
# plt.show()

diff_313_303 = interp_IR_P42A_313(SOC_new)-interp_IR_P42A_303(SOC_new)       # interpolate IR for missing temperatures based on difference 
interp_IR_P42A_323 = interp_IR_P42A_313(SOC_new) - diff_313_303
interp_IR_P42A_333 = interp_IR_P42A_323 - diff_313_303    


def interpolated_resistance(target_temp,cell,current_index,num_parallel):
    for i in range(100):  
        if cell == 'P42A':
            temperatures = [25, 30, 40, 50, 60]
            resistances = [interp_IR_P42A_298(SOC_new)[i], interp_IR_P42A_303(SOC_new)[i], interp_IR_P42A_313(SOC_new)[i], interp_IR_P42A_323[i], interp_IR_P42A_333[i]]
        else:
            temperatures = [30, 35, 40, 45, 50]
            
            # check criteria for 2C (0 - 13.5 A) - if not use 4C parameters
            if 0 < current[current_index]/num_parallel <= 13.5: 
                resistances = [interp_IR_P45B_303_2C(SOC_new)[i], interp_IR_P45B_308_2C(SOC_new)[i], interp_IR_P45B_313_2C(SOC_new)[i], interp_IR_P45B_318_2C(SOC_new)[i], interp_IR_P45B_323_2C(SOC_new)[i]]
            else:
                resistances = [interp_IR_P45B_303_4C(SOC_new)[i], interp_IR_P45B_308_4C(SOC_new)[i], interp_IR_P45B_313_4C(SOC_new)[i], interp_IR_P45B_318_4C(SOC_new)[i], interp_IR_P45B_323_4C(SOC_new)[i]]
        f = interp1d(temperatures, resistances, kind='quadratic') 
    return f(target_temp)

# plotting IR trend at each temperature in target_temepratures
def plot_resistance(cell):
    target_temperatures = np.linspace(25,50,51)
    def interpolated_resistance_plotting(target_temp):
        new_resistances = []       
        for i in range(100):  
            if cell == 'P42A':
                temperatures = [25, 30, 40, 50, 60]
                resistances = [interp_IR_P42A_298(SOC_new)[i], interp_IR_P42A_303(SOC_new)[i], interp_IR_P42A_313(SOC_new)[i], interp_IR_P42A_323[i], interp_IR_P42A_333[i]]
            else:
                temperatures = [30, 35, 40, 45, 50]
                if cell == 'P45B_2C':
                    resistances = [interp_IR_P45B_303_2C(SOC_new)[i], interp_IR_P45B_308_2C(SOC_new)[i], interp_IR_P45B_313_2C(SOC_new)[i], interp_IR_P45B_318_2C(SOC_new)[i], interp_IR_P45B_323_2C(SOC_new)[i]]
                else:
                    resistances = [interp_IR_P45B_303_4C(SOC_new)[i], interp_IR_P45B_308_4C(SOC_new)[i], interp_IR_P45B_313_4C(SOC_new)[i], interp_IR_P45B_318_4C(SOC_new)[i], interp_IR_P45B_323_4C(SOC_new)[i]]
            f = interp1d(temperatures, resistances, kind='quadratic') 
            new_resistances.append(f(target_temp))
        return np.array(new_resistances)
            
    for temp in target_temperatures:
        plt.plot(SOC_new, interpolated_resistance_plotting(target_temperatures), label=f'{temp:.2f} C') 
    
    plt.xlabel('SOC (%)')
    plt.ylabel('Interpolated Internal Resistance (Ohms)')
    plt.title(f'Interpolated Internal Resistance vs SOC, {cell}')
    plt.grid(True)
    plt.show()
    
    
# Entropic Heat

EH_raw = pd.read_csv('Entropic Heat Data.csv')
EH_data = np.array(EH_raw,dtype = float)
EH_P42A = EH_data[:,0]
EH_P45B_2C = EH_data[:,1]
EH_P45B_4C = EH_data[:,2]

interp_EH_P42A = interp1d(SOC,EH_P42A,kind='cubic',fill_value = 'extrapolate')
interp_EH_P45B_2C = interp1d(SOC,EH_P45B_2C,kind='cubic',fill_value = 'extrapolate')
interp_EH_P45B_4C = interp1d(SOC,EH_P45B_4C,kind='cubic',fill_value = 'extrapolate')
# plt.plot(SOC_new,interp_EH_P42A(SOC_new),label = 'P42A')
# plt.plot(SOC_new,interp_EH_P45B_2C(SOC_new),label = 'P45B 2C')
# plt.plot(SOC_new,interp_EH_P45B_4C(SOC_new),label = 'P45B 4C')
# plt.xlabel('SOC (%)')
# plt.ylabel('Entropic Heat Coefficient (mV/K)')
# plt.legend()
# plt.show()

def entropic_heat(targetSOC,cell,current_index,num_parallel):
    for i in range(100):
        if cell == 'P42A':
            coeff = interp_EH_P42A(SOC_new)
        else: # check criteria for 2C (0 - 13.5 A) - if not use 4C parameters
            if 0 <= current[current_index]/num_parallel <= 13.5:        # current_index implemented in cell_temp_ode
                coeff = interp_EH_P45B_2C(SOC_new)
            else:
                coeff = interp_EH_P45B_2C(SOC_new)
        dOCVdt = interp1d(SOC_new, coeff,kind = 'cubic')
    return dOCVdt(targetSOC)
    
# P_avg = 18.6; # [kW] Average power 
# V_nom = 432 # [V]
# I = (P_avg*1000/V_nom)/4; # [A]


# Interpolating current data

# sorted_indices = np.argsort(time)
# xx_sorted = time[sorted_indices]
# yy_sorted = current[sorted_indices]

# xx_unique, indices = np.unique(xx_sorted, return_index=True)
# yy_unique = yy_sorted[indices]

# current_interp = interp1d(xx_unique, yy_unique, kind='cubic', fill_value='extrapolate')
# calculated_interp = current_interp(current)
# plt.plot(time/60,current_interp(time))


def cell_temp_ode(Y, t, cell, num_parallel):
    current_index = int(t) 

    if Y[0] > 55:  # P45B internal resistance data stops at 50 C, could manually extend it to 60 C
           return None
       
    if SOC_data[current_index] > 100 or SOC_data[current_index] < 0:
        SOC_data[current_index] = 11
    
    # To quickly test constant values, comment the corresponding lines out
    # R= 0.011 #Ohms      
    # htc = htc_external
    
    R = interpolated_resistance(Y[0],cell,current_index,num_parallel)                   # find interpolated R based on cell temp
    qgen = (current[current_index]/num_parallel)**2 * R #+ (current[current_index]/num_parallel * (Y[0]+273.15) * entropic_heat(SOC_data[current_index],cell,current_index,num_parallel)/1000)     # using raw current data
    # qgen = I**2 * R                                   # using average power
    # qgen = calculated_interp[current_index]**2 * R    # using interpolated current data
    
    target = Y[0]
    htc = interpolate_htc(temp_range, target)

    # htc = 7
    if 0 <= current_index < len(current):
        C1 = qgen / (m * Cp_cell)
        C2 = htc * SA_cell / (m * Cp_cell)
        
        return [C1 - C2 * (Y[0]-temp[0])]
    
    return [0]


# initialize time array, heavily discretized 
h = 0.5
h = np.arange(0,2400+h,h)

Y0 = [temp[0]]
# Y0 = [35]
solution = si.odeint(cell_temp_ode,Y0, h, args = ('P45B',3))

plt.plot(h/60,solution[:],'r-')

print('Final cell temperature:', solution[-1])
plt.xlabel('time (min)')
plt.ylabel('temperature (degC)')
plt.show()

#%% 

# General ODE solver

D_cell = 21.55e-3; # [m] Cell diameter
R_cell = D_cell/2; # [m] Cell radius
L_cell = 70.15e-3; # [m] Cell length/height
SA_cell = 0.5*(L_cell-12.7e-3)*2*np.pi*R_cell; # [m^2] Surface Area of cell, length minus total cell holders 
Cp_cell = 1360     # [J/kg K]

#Air Properties at 35 C
rho = 1.1459; # [kg/m^3]

# main inputs

R = 0.011     #[Ohms]
m = 0.07      # [kg]
htc = 7       # [W/m2 K]
P_avg = 18.6; # [kW] Average power 

def transient_solver(power, parallel):
    V_nom = 432 # [V]
    I = (P_avg*1000/V_nom)/parallel; # [A]
    
    qgen = I**2*R
    
    C1 = qgen/(m*Cp_cell) 
    C2 = htc*SA_cell/(m*Cp_cell)
    
    F = lambda Y,t: [C1-C2*(Y[0]-35)]
    
    Y0 = [35]
    L = 1800 # seconds
    dx = 0.05
    x = np.arange(0,L+dx,dx)
    
    solution = si.odeint(F,Y0, x)
    
    plt.plot(x/60,solution[:,0],'r-')
    plt.xlabel('time (min)')
    plt.ylabel('temperature (degC)')
    plt.title('Transient, Constant Power = 18.6 kW')
    plt.show()

transient_solver(18.6,4)

#%%
#coupled system of Tcell and Tamb

# def final_temp_ode(Y,t):  
#     current_index = int(t)
#     Vair = (23e-3*23e-3*57.45e-3) - (np.pi*R_cell**2*L_cell) # Volume of air around cell minus volume of cell

#     if Y[1] > 50:      # P45B internal resistance data stops at 50 C, could manually extend it to 60 C
#            return [0,0]
       
#     if SOC_data[current_index] > 100 or SOC_data[current_index] < 0:
#         SOC_data[current_index] = 11
        
#     R = interpolated_resistance(Y[1])   
#     qgen = current[current_index]**2 * R + (current[current_index]/4 * (Y[1]+273.15) * entropic_heat(SOC_data[current_index])/1000)
    
#     # qgen = I**2 * R 

#     target = Y[1]
#     htc = interpolate_htc(temp_range, target)
    
    # if 0 <= current_index < len(current):
    #     ambT = (qgen - htc*SA_cell * (Y[1]-Y[0])) / (rho * Vair * Cp)
    #     cellT = (qgen + htc*SA_cell * (Y[0]-Y[1])) / (m * Cp_cell)
        
    #     # Y[0] = Tamb, Y[1] = Tcell
    #     return [ ambT , cellT ]

# Y0 = [temp[0],temp[0]]

# solution = si.odeint(final_temp_ode,Y0,time,args = ('P42A',))

# h = 0.5
# h = np.arange(0,2400+h,h)

# plt.plot(time/60,solution[:,0],'r-',label = 'ambient')
# plt.plot(time/60,solution[:,1],'b-',label = 'cell')
# # solution_array = solution[0] #extract the array from the tuple
# # plt.plot(h/60, solution_array[:, 0], 'r-', label='ambient')
# # plt.plot(h/60, solution_array[:, 1], 'b-', label='cell')

# plt.xlabel('time (min)')
# plt.ylabel('temperature (degC)')
# plt.legend()
# plt.show()

# Y0 = [temp[0]]
#Y0 = [35]

