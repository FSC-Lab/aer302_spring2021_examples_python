import argparse
import numpy as np
import matplotlib.pyplot as plt

import SAM
from aircraft import CP1, CJ1
from funcs import *

################## arguments from command line #################
ap = argparse.ArgumentParser("AER302 Performance example")
ap.add_argument('-plane', default='CP1')
ap.add_argument('-V', nargs='?', type=float)
ap.add_argument('-H', nargs='?', type=float, default=0)
ap.add_argument(
	'--endurance', 
	action='store_true',
	default=False,
	help='compute the endurance of the plane'
)
ap.add_argument(
	'--range', 
	action='store_true',
	default=False,
	help='compute the range of the plane'
)
ap.add_argument(
	'--power', 
	action='store_true',
	default=False,
	help='compute power required/available of the plane'
)
ap.add_argument(
	'--climb_rate', 
	action='store_true',
	default=False,
	help='compute the climb rate of the plane'
)
# the following analysis requires to provide velocity range
ap.add_argument(
	'--power_plot', 
	nargs='?',
	default='',
	help='plot power required/available vs. velocity \'lb ub\''
)
ap.add_argument(
	'--hodograph', 
	nargs='?',
	default='',
	help='plot the hodograph, please specify velocity range as \'lb ub\''
)
ap.add_argument(
	'--ceiling', 
	nargs='?',
	default='0 12000',
	help='plot (R/C)max vs. altitude and compute the ceiling of the plane\
	      please specify altitude range as \'lb ub\''
)

args = ap.parse_args()

if args.plane == 'CP1':
	plane = CP1
if args.plane == 'CJ1':
	plane = CJ1

h = args.H
if h == 0:
	air = SAM.air_sea
if 0 <= h < 11000:
	air = SAM.tropo(h)
if h > 11000:
	air = SAM.strato(h)

V = args.V
if V is None:
	V = 60 # TODO, compute Vminpr as default

RAD2DEG = 180/3.14

############## power available/required ###############
if V is not None and h is not None:
	print('%s performance analysis at V=%.2fm/s, H=%.0fm:'%(args.plane, V, h))
elif h is not None:
	print('%s performance analysis at H=%.0fm:'%(args.plane, h))
else:
	print('%s performance analysis at V=%.2fm:'%(args.plane, V))

if args.power:
	PR = power_required(plane, V, air)
	PA = power_available(plane, air)
	print('power required: %.2fW'%PR)
	print('power available: %.2fW'%PA)

if args.climb_rate:
	RC, gmm = climb_rate(plane, V, air)
	print('climb rate: %.2fm/s'%RC)
	print('flight path angle: %.2f degree'%(gmm*RAD2DEG))

if args.endurance:
	E = endurance(plane, h)
	print('endurance: %.2fh'%E)

if args.range:
	R = range(plane, h)
	print('range: %.2fkm'%(R/1000))

if args.ceiling:
	temp = args.ceiling.split(' ')
	lb = float(temp[0])
	ub = float(temp[1])
	hs = np.linspace(lb, ub, 100)
	RCs = climb_rate_max(plane, hs)
	# TODO: solve h from pyinverse

	# TODO: polish plot
	plt.plot(RCs, hs)
	plt.grid()
	plt.show()	

if args.power_plot:
	temp = args.power_plot.split(' ')
	lb = float(temp[0])
	ub = float(temp[1])
	Vs = np.linspace(lb, ub, 100)
	PR, PA, RC, Gmm = power_vs_v(plane, air, Vs)

	# TODO: polish plot
	plt.plot(Vs, PA)
	plt.plot(Vs, PR)
	plt.grid()
	plt.show()

if args.hodograph:
	temp = args.hodograph.split(' ')
	lb = float(temp[0])
	ub = float(temp[1])
	Vs = np.linspace(lb, ub, 100)
	PR, PA, RC, Gmm = power_vs_v(plane, air, Vs)
	Vx, Vy = hodograph(Vs, Gmm)

	plt.plot(Vx, Vy)
	plt.grid()
	plt.show()
