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

def make_plot(x, y, xlabel, ylabel, title, i):
    plt.rcParams.update({'font.size': 8})
    plt.figure(figsize=(4,3))
    # plt.subplot(1, 2, i) # plot two plots on the same line
    plt.plot(x, y, linewidth=0.4)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    if i == 2:
        plt.show()

"""
1. Aircraft Flight Condition
"""
ht_ft = 6000 # Altitude (ft) from case study
ht = ht_ft * 0.3048 # Altitude (m)
m = 6.126 # Aircraft Dry Mass + Payload Mass (kg) exp
m_max = 10 # Maximum Take-Off Weight (kg) from datasheet
h = (0.085 + 0.095) / 2 # CG Position (m from root leading edge) from datasheet page
gamma_e_deg = 0 # Flight Path Angle (deg)
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
b = (3.01 + 2.96) / 2 # Wing Span (m) exp
c_tip = 0.13 # Tip Chord Length (m) exp
c_root = 0.29 # Rood Chord Length (m) exp
h = h / c_root # CG Position (% of root chord) from datasheet page ASK
c_w = (((c_tip + c_root) / 2) + 0.223) / 2 # Wing Mean Aerodynamic Chord (m) exp derived averaged
S = b * c_w # Wing Area (m^2) exp derived
Ar = b ** 2 / S # exp derived
lambda_ = 0 # Quarter Chord Sweep (deg) exp
z_w = 0 # Z-Coordinate of Quarter Chord (m) exp
alpha_w_r_deg = (11.539 + 9.46) / 2 # Wing Rigging Angle (deg) exp averaged
alpha_w_r = alpha_w_r_deg / 57.3 # Wing Rigging Angle (rad)

# Tailplane Geometry
s_T = 0.38 # Tailplane Semi-Span (m) exp
tau_T_deg = 36.1 # Tailplane Dihedral (deg) exp
c_MAC_T = 0.14 # Tailplane Mean Aerodynamic Chord (m) exp
b_T = (0.63 + 2 * s_T * np.cos(np.radians(tau_T_deg))) / 2 # Tailplane Span (m) exp averaged
S_T = c_MAC_T * b_T # Tailplane Area (m^2)
l_t = 1.04 # Tail Arm Quarter Chord Wing to Quarter Chord Tail (m) exp
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
a = 5.19 # Wing-body CL-alpha (rad^-1)
C_L_max = 1.52 # Maximum Lift Coefficient
C_m_0 = -0.05 # Zero Lift Pitching Moment Coefficient
C_D0 = 0.041 # Zero Lift Drag Coefficient
alpha_w0_deg = -2.5 # Zero Lift Angle of Attack (deg)
alpha_w0 = alpha_w0_deg / 57.3 # Zero Lift Angle of Attack (rad)
h_0 = 0.25 # Wing-Body Aero Centre

"""
6. Tailplane Aerodynamics
"""
a1 = 3.2 # Tail plane CL-alpha (rad^-1)
a2 = 2.414 # Elevator CL-eta (rad^-1)
epsilon_0_deg = 2 # Zero Lift Downwash Angle (deg)
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
K = 1 / (np.pi * Ar * e)

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

V_max_km_hr = 129 # Maximum Airspeed (kmhr^-1)
V_max_knots = V_max_km_hr / 3.6 / 0.515 # Maximum Airspeed (knots)
V_max_i = V_max_knots * 0.515 # Maximum Airspeed (ms^-1)
interval = (V_max_knots - V_stall_eas) / 11

V_knots = [] # True Airspeed (knots)
V_i = [] # True Airspeed (ms^-1)
V_eas = [] # Equivalent Airspeed (knots)
for i in range(0, 15):
    V_knots.append(V_stall * 0.8 + (interval * i))
    V_i.append(V_knots[i] * 0.515)
    V_eas.append(V_knots[i] * np.sqrt(sigma))

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
    return [eq1,eq2,eq3,eq4,eq5,eq6]

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
    ('V_i', V_i[4], 'knots'), 
    ('C_L_i', C_L_i[4], '-'), 
    ('C_D_i', C_D_i[4], '-'),
    ('C_LW_i', C_LW_i[4], '-'), 
    ('C_LT_i', C_LT_i[4], '-'), 
    ('LD_i', LD_i[4], '-'), 
    ('C_tau_i', C_tau_i[4], '-'), 
    ('alpha_w_i', alpha_w_i[4], 'deg'), 
    ('alpha_e_i', alpha_e_i[4], 'deg'), 
    ('theta_e_i', theta_e_i[4], 'deg'), 
    ('alpha_T_i', alpha_T_i[4], 'deg'), 
    ('eta_e_i', eta_e_i[4], 'deg'), 
    ('L_i', L_i[4], 'N'), 
    ('D_i', D_i[4], 'N'), 
    ('T_i', T_i[4], 'N')
]

print('==Trim=Conditions==Demo=Values=============')
output(results_trim_conditions)

"""
17. Some Useful Trim Plots
"""
plt.rcParams.update({'font.size': 8})
plt.figure(figsize=(4,3))
plt.plot(V_knots, LD_i, linewidth=0.8)
plt.axvline(x=V_stall, color='r', 
            label='V_Stall', linestyle='--', linewidth=0.5)
plt.text(V_stall, max(LD_i) * 0.9, 'Stall Velocity', color='r', 
         rotation=0, ha='center', va='bottom', fontsize=8)
plt.xlabel('Velocity (knots)')
plt.ylabel('Lift to Drag Ratio (-)')
plt.title('Velocity vs Lift to Drag Ratio')
plt.grid()

plt.figure(figsize=(4,3))
plt.plot(V_knots, eta_e_i, linewidth=0.8)
plt.axvline(x=V_stall, color='r', 
            label='V_Stall', linestyle='--', linewidth=0.5)
plt.text(V_stall, max(eta_e_i) * 0.9, 'Stall Velocity', color='r', 
         rotation=0, ha='center', va='bottom', fontsize=8)
plt.xlabel('Velocity (knots)')
plt.ylabel('Elevator Angle (deg)')
plt.title('Velocity vs Elevator Angle')
plt.grid()

plt.figure(figsize=(4,3))
plt.plot(V_knots, D_i, linewidth=0.8)
plt.axvline(x=V_stall, color='r', 
            label='V_Stall', linestyle='--', linewidth=0.5)
plt.text(V_stall, max(D_i) * 0.9, 'Stall Velocity', color='r', 
         rotation=0, ha='center', va='bottom', fontsize=8)
plt.xlabel('Velocity (knots)')
plt.ylabel('Total Drag (N)')
plt.title('Velocity vs Total Drag')
plt.grid()

plt.figure(figsize=(4,3))
plt.plot(C_D_i, C_L_i, linewidth=0.8)
plt.xlabel('Drag Coefficient')
plt.ylabel('Lift Coefficient')
plt.title('Drag Polar')
plt.grid()

# =============================================================================
# make_plot(V_knots, LD_i, 'Velocity (knots)', 'Lift to Drag Ratio (-)', 'Velocity vs Lift to Drag Ratio', 1)
# make_plot(V_knots, eta_e_i, 'Velocity (knots)', 'Elevator Angle (deg)', 'Velocity vs Elevator Angle', 2)
# make_plot(V_knots, D_i, 'Velocity (knots)', 'Total Drag (N)', 'Velocity vs Total Drag', 1)
# make_plot(C_D_i, C_L_i, 'Drag Coefficient', 'Lift Coefficient', 'Drag Polar', 2)
# =============================================================================

"""
References
==========
datasheet page - https://www.appliedaeronautics.com/albatross-uav#:~:text=The%20Albatross%20UAV%20offers%20robust,needing%20to%20land%20to%20recharge.

"""
