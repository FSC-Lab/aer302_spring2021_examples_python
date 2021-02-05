class IdealGas(object):
	""" air """
	def __init__(self, rho=None, P=None, T=None):

		self.R = 287

		N = sum([p is not None for p in [rho, P, T]])
		assert N >= 2

		if N == 3:
			assert abs(P - rho*self.R*T) < 20
		else:
			if rho is None:
				rho = P/(self.R*T)
			if P is None:
				P = rho*self.R*T
			if T is None:
				T = P/(rho*self.R)

		self.rho = rho
		self.P = P
		self.T = T