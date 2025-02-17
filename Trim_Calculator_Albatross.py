"""
Aircraft Trim Model
Case Study 1 - Albatross over the ocean
AERO32202 Flight Dynamics
Nicolas Arroyo 11091029
"""
"""
0. Library Imports and Functions
"""
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def output(results):
    print('|      Variable |      Value |      Units |')
    for var, val, unit in results:
        print(f'| {var:>13} | {round(val, 3):>10} | {unit:>10} |')

"""
1. Aircraft Flight Condition
"""
ht_ft = 6000 # Altitude (ft) from case study
ht = ht_ft * 0.3048 # Altitude (m)
m = 6.126 # Aircraft Dry Mass + Payload Mass (kg) exp
mtow = 10 # Maximum Take-Off Weight (kg) from datasheet
h = 0.085 # CG Position from root leading edge (m) from datasheet page
gamma_e_deg = 0 # Flight Path Angle (deg) MIGHT CHANGE
gamma_e = gamma_e_deg / 57.3 # Flight Path Angle (rad)
g = 9.81 # Gravity Constant (ms^-2)

"""
2. Air Density Calculation
"""
R = 287.05 # Gas Constant (Nm kg^-1 K^-1)
Ir = -0.0065 # Lapse Rate (Km^-1)
temp = 288.16 + (Ir * ht) # Temperature (K)
rho = 1.225 * (temp / 288.16)**(-((g / (Ir * R)) + 1)) # Air Density (kgm^-3)
sigma = rho / 1.225 # Density ratio

"""
Check Results 1
"""
results1 = [
    ('temp', temp, 'K'),
    ('rho', rho, 'kg m^-3'),
    ('sigma', sigma, '-')
]
print('==Output=1=================================')
output(results1)

"""
3. Set up Velocity Range for Computations
(note that true airspeed (TAS) is assumed unless otherwise stated)
"""
# See Section 10. Basic Performance Parameters

"""
4. Aircraft Geometry
"""
# Wing Geometry
b = (3.01 + 2.96) / 2 # Wing Span (m) exp averaged
c_tip = 0.13 # Tip Chord Length (m) exp
c_root = 0.29 # Root Chord Length (m) exp
h = h / c_root # CG Position (% of root chord) from datasheet
c_w = (((c_tip + c_root) / 2) + 0.223) / 2 # Wing Mean Aerodynamic Chord (m) exp derived averaged
S = b * c_w # Wing Area (m^2) exp derived
Ar = b**2 / S # exp derived
lambda_ = 0 # Wing Quarter Chord Sweep (deg) exp
z_w = 0 # Z-Coordinate of Quarter Chord (m) exp
alpha_w_r_deg = (11.539 + 9.46) / 2 # Wing Rigging Angle (deg) exp averaged
alpha_w_r = alpha_w_r_deg / 57.3 # Wing Rigging Angle (rad)

# Tailplane Geometry
s_T = (0.38 + (0.63 / 2)) / 2 # Tailplane Semi-Span (m) exp averaged
tau_T_deg = 36.1 # Tailplane Dihedral (deg) exp
c_MAC_T = 0.14 # Tailplane Mean Aerodynamic Chord (m) exp
b_T = 2 * s_T * np.cos(np.radians(tau_T_deg)) # Tailplane Span (m) exp derived
S_T = (0.11424 + (c_MAC_T * b_T)) / 2 # Tailplane Area (m^2) exp derived averaged
Ar_T = (3.47 + (b_T**2 / S_T)) / 2 # Tailplane Aspect Ratio exp derived averaged
l_t = 1.04 # Tail Arm, Quarter Chord Wing to Quarter Chord Tail (m) exp
lambda_T = 14.04 # Tailplane Sweep Angle (deg)
z_T = -0.35 # Quarter Chord Z-Coordinate (m) exp
eta_T_deg = 0 # Tailplane Setting Angle (deg) exp
eta_T = eta_T_deg / 57.3 # Tailplane Setting Angle (rad)

# General Geometry
F_d = 0.25 # Fuselage Diameter or Width (m) exp
z_tau = -0.07 # Thrust Line Z-Coordinate (m) exp
kappa_deg = 0 # Engine Thrust Line Angle (deg) exp
kappa = kappa_deg / 57.3 # Engine Thrust Line Angle (rad)

"""
5. Wing-Body Aerodynamics
"""
a = 2 * np.pi * Ar / (2 + Ar) # Wing-body CL-alpha (rad^-1) aero notes
C_L_max = 1.27 # Maximum Lift Coefficient 3rd year project +- 10%
C_m_0 = -0.5 # Zero Lift Pitching Moment Coefficient 3rd year project
C_D0 = 0.032 # Zero Lift Drag Coefficient 3rd year project
alpha_w0_deg = -2 # Zero Lift Angle of Attack (deg) 3rd year project
alpha_w0 = alpha_w0_deg / 57.3 # Zero Lift Angle of Attack (rad)
h_0 = 0.25 # Wing-Body Aero Centre 3rd year project

"""
6. Tailplane Aerodynamics
"""
a1_numerator = 2 * np.pi * Ar_T
beta = (1 - 0.2**2)**0.5 # Mach Number Parameter
kappa_ratio = 1 # Ratio of 2D Lift Curve Slope to 2pi (assumed to be perfect, i.e. 1)
term1 = (Ar_T * beta / kappa_ratio)**2
term2 = np.tan(np.radians(lambda_T))**2 / beta**2
a1_denominator = 2 + np.sqrt(term1 * (1 + term2) + 4)
a1 = a1_numerator / a1_denominator
a2 = 0.26 * a1 # Elevator CL-eta (rad^-1) aero notes
epsilon_0_deg = 2 # Zero Lift Downwash Angle (deg) aero notes
epsilon_0 = epsilon_0_deg / 57.3 # Zero Lift Downwash Angle (rad)

"""
7. Wing and Tailplane Calculations
"""
s = b / 2 # Wing Semi-Span (m)
l_T = l_t - c_w * (h - 0.25) # Tail Arm cg to Tail Quarter Chord (m)
V_T = (S_T * l_T) / (S * c_w) # Tail Volume Coefficient

"""
Check Results 2
"""
results2 = [
    ('Ar', Ar, '-'),
    ('s', s, 'm'),
    ('l_T', l_T, 'm'), 
     ('V_T', V_T, '-')
]
print('==Output=2=================================')
output(results2)

"""
8. Downwash at Tail
"""
x = l_t / b
z = (z_w - z_T) / b
A = a / ((np.pi**2) * Ar)
final = np.pi / 180

num3 = x
den3 = (x**2) + (z**2)

total = 0
for fi in range(5, 175):
    num1 = (0.5 * np.cos(fi * np.pi / 180))**2
    den1 = np.sqrt((x**2) + ((0.5 * np.cos((fi * np.pi / 180)))**2) + (z**2))
    num2 = x + np.sqrt((x**2) + ((0.5 * np.cos((fi * np.pi / 180)))**2) + (z**2))
    den2 = ((0.5 * np.cos((fi * np.pi / 180)))**2) + (z**2)
    total = total + (num1 / den1) * ((num2 / den2) + (num3 / den3))

d_epsilon_alpha = A * total * final # Tail Position Relative to Wing (% of span)

"""
Check Results 3
"""
results3 = [
    ('d_ε_α', d_epsilon_alpha, '% of span')
]
print('==Output=3=================================')
output(results3)

"""
9. Induced Drag Factor
"""
C = F_d / b
S_d = (0.9998 + (0.0421 * C)) - (2.6286 * C**2) + (2 * C**3) # Fuselage drag factor
k_D = (-3.333 * 10**(-4) * lambda_**2) + (6.667 * 10**(-5) * lambda_) + 0.38 # Empirical Constant
e = 1 / (np.pi * Ar * k_D * C_D0 + (1 / (0.99 * S_d))) # Oswald Efficiency Factor
# verified reasonable 'e' value with Fundamentals of Flight 2nd edition, page 187
# tailplane contribution is not considered?

K = 1 / (np.pi * Ar * (e))

"""
Check Results 4
"""
results4 = [
    ('S_d', S_d, '-'),
    ('k_D', k_D, '-'),
    ('e', e, '-'), 
     ('K', K, '-')
]
print('==Output=4=================================')
output(results4)

"""
10. Basic Performance Parameters
"""
VCL = np.sqrt(2 * m * g / (rho * S)) # Useful expression of V * CL from lift equation
V_md = (VCL * (K / C_D0)**0.25) / 0.515 # Minimim Drag Speed (knots)
V_md_eas = V_md * np.sqrt(sigma) # Equivalent Minimum Drag Speed (knots)
V_stall = (VCL / np.sqrt(C_L_max)) / 0.515 # Stall Speed (knots)
V_stall_eas = V_stall * np.sqrt(sigma) # Equivalent Stall Speed (knots)
h_n = h_0 + V_T * (a1 / a) * (1 - d_epsilon_alpha) # Neutral Points - controls fixed
K_n = h_n - h # Static Margin - controls fixed

V_max_km_hr = 129 # Maximum Airspeed (kmhr^-1) from datasheet
V_max_knots = V_max_km_hr / 3.6 / 0.515 # Maximum Airspeed (knots)
V_max_i = V_max_knots * 0.515 # Maximum Airspeed (ms^-1)

V_knots = [] # True Airspeed (knots)
V_i = [] # True Airspeed (ms^-1)
V_eas = [] # Equivalent Airspeed (knots)
array_min = 0.95 * V_stall
array_max = 1.05 * V_max_knots
intervals = 50

V_knots = np.linspace(array_min, array_max, intervals)
V_i = V_knots * 0.515
V_eas = V_knots * np.sqrt(sigma)

"""
11. Trim Calculation
"""
def equations(vars): # Defining the system of equations
    C_L, C_LW, C_D, C_tau, alpha_e, C_LT = vars  # Unpack variables
    eq1 = 2 * m * g / (rho * vel**2 * S) * np.cos(alpha_e + gamma_e) - \
        (C_L * np.cos(alpha_e) + C_D * np.sin(alpha_e) + C_tau * np.sin(kappa))
    eq2 = 2 * m * g / (rho * vel**2 * S) * np.sin(alpha_e + gamma_e) - \
        (C_tau * np.cos(kappa) - C_D * np.cos(alpha_e) + C_L * np.sin(alpha_e))
    eq3 = - C_D + C_D0 + K * C_L**2
    eq4 = - C_LW + a * (alpha_e + alpha_w_r - alpha_w0)
    eq5 = C_m_0 + (h - h_0) * C_LW - V_T * C_LT + C_tau * z_tau / c_w
    eq6 = - C_LT + (C_L - C_LW) * S / S_T
    return [eq1, eq2, eq3, eq4, eq5, eq6]

initial_guesses = [0.7, 0.5, 0.02, 0.4, 0.1, 0.1]

"""
12. Trim Variables Calculation
"""
C_L_i = np.zeros(shape=len(V_i))
C_LW_i = np.zeros(shape=len(V_i))
C_D_i = np.zeros(shape=len(V_i))
C_tau_i = np.zeros(shape=len(V_i))
alpha_e_i = np.zeros(shape=len(V_i))
C_LT_i = np.zeros(shape=len(V_i))

for i in range(0, len(V_i)):
    vel = V_i[i]
    solution, info, ier, msg = fsolve(equations, initial_guesses, full_output=True)
    C_L_i[i], C_LW_i[i], C_D_i[i], C_tau_i[i], alpha_e_i[i], C_LT_i[i] = solution

alpha_w_i = alpha_e_i + alpha_w_r # Wing Incidence (rad)
eta_e_i = C_LT_i / a2 - a1 / a2 * (alpha_w_i * (1 - d_epsilon_alpha) + eta_T - alpha_w_r - epsilon_0) # Trim Elevator Angle (rad)

index = 4
results_check = [
    (f'C_LT_i_{index}/a2', C_LT_i[index] / a2, ''), 
    ('a1/a2', a1 / a2, ''), 
    (f'α_w_i{index}', alpha_w_i[index], '-'), 
    ('ηT-αwr-ε0', eta_T - alpha_w_r - epsilon_0, ''),
    ('(1 - d_ε_α)', (1 - d_epsilon_alpha), '')
]
output(results_check) # C_LT_i is too high, making the elevator very sensitive

theta_e_i = gamma_e + alpha_w_i - alpha_w_r # Pitch Attitude (rad)
alpha_T_i = alpha_w_i * (1 - d_epsilon_alpha) + eta_T - epsilon_0 - alpha_w_r # Tail Angle of Attach (rad)
LD_i = C_LW_i / C_D_i # Lift to Drag Ratio

"""
13. Conversions of Angles to Degrees
"""
alpha_w_i = alpha_w_i * 57.3
alpha_e_i = alpha_e_i * 57.3
theta_e_i = theta_e_i * 57.3
alpha_T_i = alpha_T_i * 57.3
eta_e_i = eta_e_i * 57.3
gamma_e = gamma_e * 57.3

"""
14. Total Trim Forces Acting on Aircraft
"""
L_i = []
D_i = []
T_i = []

for i in range(0, len(V_i)):
    vel = V_i[i]
    L_i.append(0.5 * rho * vel**2 * S * C_L_i[i]) # Total Lift Force (N)
    D_i.append(0.5 * rho * vel**2 * S * C_D_i[i]) # Total Lift Force (N)
    T_i.append(0.5 * rho * vel**2 * S * C_tau_i[i]) # Total Lift Force (N)

"""
15. Definition of Flight Condition
"""
results_flight_condition = [
    ('mg', m*g, 'N'), 
    ('ht_ft', ht_ft, 'ft'), 
    ('gamma_e', gamma_e, 'deg'),
    ('h', h, '% of chord'), 
    ('h_n', h_n, '-'), 
    ('V_md', V_md, 'knots'), 
    ('V_md_eas', V_md_eas, 'knots'), 
    ('V_stall', V_stall, 'knots'), 
    ('V_stall_eas', V_stall_eas, 'knots'), 
    ('K_n', K_n, '-')
]
print('==Definition=of=Flight=Condition===========')
output(results_flight_condition)

"""
16. Trim Conditions as a Function of Aircraft Velocity
"""
results_trim_conditions = [
    ('V_i', V_i[14], 'knots'), 
    ('C_L_i', C_L_i[14], '-'), 
    ('C_D_i', C_D_i[14], '-'),
    ('C_LW_i', C_LW_i[14], '-'), 
    ('C_LT_i', C_LT_i[14], '-'), 
    ('LD_i', LD_i[14], '-'), 
    ('C_tau_i', C_tau_i[14], '-'), 
    ('alpha_w_i', alpha_w_i[14], 'deg'), 
    ('alpha_e_i', alpha_e_i[14], 'deg'), 
    ('theta_e_i', theta_e_i[14], 'deg'), 
    ('alpha_T_i', alpha_T_i[14], 'deg'), 
    ('eta_e_i', eta_e_i[14], 'deg'), 
    ('L_i', L_i[14], 'N'), 
    ('D_i', D_i[14], 'N'), 
    ('T_i', T_i[14], 'N')
]

print('==Trim=Conditions==Demo=Values=============')
output(results_trim_conditions)

"""
17. Some Useful Trim Plots
"""
plt.rcParams.update({'font.size': 8})
plt.figure(figsize=(4,3))
plt.title('Velocity vs Lift to Drag Ratio')
plt.xlabel('Velocity (knots)')
plt.ylabel('Lift to Drag Ratio (-)')
plt.plot(V_knots, LD_i, linewidth=1.5)
plt.axvline(x=V_stall, color='r', label='V_Stall', linestyle='--', linewidth=0.8)
plt.axvline(x=V_max_knots, color='g', label='V_Max', linestyle='--', linewidth=0.8)
plt.axvline(x=V_md, color='b', label='Min Drag', linestyle='--', linewidth=0.8)
plt.grid()
plt.legend(fontsize='small')

plt.figure(figsize=(4,3))
plt.title('Velocity vs Elevator Angle')
plt.xlabel('Velocity (knots)')
plt.ylabel('Elevator Angle (deg)')
plt.plot(V_knots, eta_e_i, linewidth=1.5)
plt.axvline(x=V_stall, color='r', label='V_Stall', linestyle='--', linewidth=0.8)
plt.axvline(x=V_max_knots, color='g', label='V_Max', linestyle='--', linewidth=0.8)
plt.axvline(x=V_md, color='b', label='Min Drag', linestyle='--', linewidth=0.8)
plt.grid()
plt.legend(fontsize='small')

plt.figure(figsize=(4,3))
plt.title('Velocity vs Total Drag')
plt.xlabel('Velocity (knots)')
plt.ylabel('Total Drag (N)')
plt.plot(V_knots, D_i, linewidth=1.5)
plt.axvline(x=V_stall, color='r', label='V_Stall', linestyle='--', linewidth=0.8)
plt.axvline(x=V_max_knots, color='g', label='V_Max', linestyle='--', linewidth=0.8)
plt.axvline(x=V_md, color='b', label='Min Drag', linestyle='--', linewidth=0.8)
plt.grid()
plt.legend(fontsize='small')

plt.figure(figsize=(4,3))
plt.title('Drag Polar')
plt.xlabel('Drag Coefficient')
plt.ylabel('Lift Coefficient')
plt.plot(C_D_i, C_L_i, linewidth=1.5)
plt.errorbar(x=[min(C_D_i), max(C_D_i)], y=[C_L_max, C_L_max], yerr=0.1*C_L_max, 
             color='r', label='Max CL ± 10%', linestyle='--', linewidth=0.8, capsize=3)
plt.grid()
plt.legend(fontsize='small')

"""
References
==========
datasheet page - https://www.appliedaeronautics.com/albatross-uav#:~:text=The%20Albatross%20UAV%20offers%20robust,needing%20to%20land%20to%20recharge.
accounting for tailplane sweep https://www.sciencedirect.com/science/article/pii/B978012397308500009X#fd108
"""
