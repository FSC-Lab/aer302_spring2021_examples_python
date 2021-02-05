from math import pi

W = 13122 # weight (N)
S = 16.17 # wing plan area (m^2)
b = 10.91 # wing span (m)
AR = b**2 / S # aspect ratio

# engine data
engine = 'prop'
Ps = 171580 # propeller max power (W) at sea level
eta = 0.8 # propulsion efficiency
# Fuel = 1119 # fuel capacity (gal)
# Wf = 1119 * 6.67 * 4.448 # gal * 6.67 (lb/gal of kerosene) * 4.48 (N/lb) = N
# TSFC = 0.6 # (lb of fuel) / (lb of thrust) (hr)

# aerodynamics property
CD0 = 0.025 # parasite drag coefficient
ee = 0.8 # Oswald efficiency factor
K = 1 / ( pi * ee * AR)
# c_t = TSFC / 3600 # consistent unit of second

# constants
rho_s = 1.2250 # kg/m^3
p_s = 1.01325e05 # N /m^2
T_s = 288.16 # K

# derivations
PAs = eta * Ps # sea-level power available
WS = W / S # wing load

# C_L_stall = 1.0 # maximum C_L on the ground
# V_stall = sqrt(2*W / (rho_s * S * C_L_stall)) # stall velocity during takeoff (airborne)
# C_L_LOF = C_L_stall / 1.21
# V_LOF = 1.1 * V_stall
# C_L_2 = C_L_stall / 1.44
# V_2 = 1.2 * V_stall
# wingh = 1.83 # (m) wing height on the ground
# phi = (16 * wingh / b)^2 / ( 1 +  (16 * wingh / b)^2)
# mu_t = 0.02 # friction between wheels and ground when takeoff