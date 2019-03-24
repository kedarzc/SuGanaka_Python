import numpy as np

# This is the standard material library

def formdsig(E,nu):

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

	return dee