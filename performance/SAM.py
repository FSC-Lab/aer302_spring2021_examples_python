from math import exp
from ideal_gas import IdealGas
		
# parameters
a = -6.5e-03 # gradiant slope
g = 9.8
R = 287

# sea level constants
air_sea = IdealGas(rho=1.225, P=1.01325e5, T = 288.16)

def tropo(h): 
	# troposphere
	# input: h in meter
	# output: IdeaGas
	assert 0 <= h <= 11000
	T = air_sea.T + a*h

	r = T/air_sea.T
	k = g/(a*R)

	P = air_sea.P*r**(-k)
	rho = air_sea.rho*r**(-k-1)

	return IdealGas(rho=rho, P=P, T=T)

air_11k = tropo(11000)

def strato(h):
	# stratosphere
	# input: h in meter
	# output: IdeaGas
	# TODO: other layers of the atomosphere
	assert h >= 11000

	T = air_11k.T
	k = exp(-g/(R*T)*(h - 11000))

	P = k*air_11k.P
	rho = k*air_11k.rho

	return IdealGas(rho=rho, P=P, T=T)

# test
if __name__ == '__main__':
	air = tropo(100)
	print(air.P)