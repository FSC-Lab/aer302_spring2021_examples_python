g = 9.81 				# m/s^2
rho_s = 1.2250 			# kg/m^3
p_s = 1.01325e05 		# N /m^2
T_s = 288.16 			# K
F_s = 2 * Tengine 		# two turbojet engines

# geometry and weight
W = 88141 				# weight (N)
S = 29.54 				# wing plan area (m^2)
b = 16.25 				# wing span (m)
AR = b^2 / S 			# aspect ratio

# aerodynamics property
C_D_0 = 0.02 			# parasite drag coefficient
ee = 0.81 				# Oswald efficiency factor
K = 1 / ( pi * ee * AR)
C_L_max = 2.5
C_L_ground = 1.0 		# maximum C_L on the ground

# engine data
engine = 'prop'
Tengine = 16236 			# one of the two turbojet engines thrust (N)
Fuel = 1119 				# fuel capacity (gal)
Wf = 1119 * 6.67 * 4.448 	# gal * 6.67 (lb/gal of kerosene) * 4.48 (N/lb) = N
TSFC = 0.6 				# (lb of fuel) / (lb of thrust) (hr)
c_t = TSFC / 3600 		# consistent unit of second

# V_ground = sqrt(2*W / (rho_s * S * C_L_ground)) # stall velocity during takeoff (airborne)
# C_L_LOF = C_L_ground / 1.21
# V_LOF = 1.1 * V_ground
# C_L_2 = C_L_ground / 1.44
# V_2 = 1.2 * V_ground
# wingh = 1.83 			# (m) wing height on the ground
# phi = (16 * wingh / b)^2 / ( 1 +  (16 * wingh / b)^2)
# mu_t = 0.02 			# friction between wheels and ground when takeoff