import numpy as np

# This is the standard material library

def formdsig(E,nu,nip):

	if nip == 4:

		c = E/(1.0-nu**2);

		dee = np.zeros(shape=(3,3))

		dee[0,0] = 1.0
		dee[0,1] = nu
		dee[0,2] = 0.0

		dee[1,0] = nu
		dee[1,1] = 1.0
		dee[1,2] = 0.0

		dee[2,0] = 0.0
		dee[2,1] = 0.0
		dee[2,2] = (1.0-nu)/2.0

		dee = c*dee

	if nip == 8:

		dee = np.zeros(shape=(6,6))

		dee[0,0] = 1.0 - nu
		dee[1,1] = 1.0 - nu
		dee[2,2] = 1.0 - nu
		dee[3,3] = (1.0 - 2.0*nu)/2.0
		dee[4,4] = (1.0 - 2.0*nu)/2.0
		dee[5,5] = (1.0 - 2.0*nu)/2.0

		dee[0,1] = nu
		dee[0,2] = nu
		dee[1,0] = nu
		dee[1,2] = nu
		dee[2,0] = nu
		dee[2,1] = nu


		dee = (E/((1.0 + nu)*(1.0-2*nu)))*dee


	return dee