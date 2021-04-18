import numpy as np

def P2R(radii, angles):
    return radii * exp(1j*angles)

def R2P(x):

    return np.abs(x)*(np.cos(np.angle(x, deg=False)) + 1j*np.sin(np.angle(x, deg=False)))

class Knot:
	

	def __init__(self, ind, ref, amp, angle, genP, genQ, conP, conQ): 
		self.index = ind

		self.amp = amp
		self.angle = angle
		self.genP = genP
		self.genQ = genQ
		self.conP = conP
		self.conQ = conQ
		self.type = ref
		self.converged = False
		print("Knoten ", self.index, self.type, self.amp, self.angle, " erstellt")

	def __str__(self):
		return str(self.index) + " U: " + str(self.amp) + " Ang: " + str(self.angle) + " Gen P:" + str(self.genP) + " GenQ " + str(self.genQ) + " ConP " + str(self.conP) + " ConQ " + str(self.conQ) + " Type " + str(self.type) 
class Cable:

	def __init__(self, fr, to, admit, imp):
		self.from_knot = fr
		self.to_knot = to

		self.imp = imp
		self.admit = admit
		print(self.from_knot.index, "-->", self.to_knot.index, "imp: ", self.imp,"admit: ", self.admit)

# PAUSE BIS 10:32 !! Ich sags

def main():
	knots = []
	cables = []

	knots.append(Knot(0, 0, 1.05, 0, 0, 0, 0, 0))
	knots.append(Knot(1, 1, 1, 0, 0, 0, 96, 62))
	knots.append(Knot(2, 1, 1, 0, 0, 0, 35, 14))
	knots.append(Knot(3, 1, 1, 0, 0, 0, 16, 8))
	knots.append(Knot(4, 2, 1.02, 0, 48, 0, 24, 11))
	print("----------------------------------")
	cables.append(Cable(knots[0],knots[1],0.03j,0.02+0.1j))
	cables.append(Cable(knots[0],knots[4],0.02j,0.05+0.25j))
	cables.append(Cable(knots[1],knots[2],0.025j,0.04+0.2j))
	cables.append(Cable(knots[1],knots[4],0.02j,0.05+0.25j))
	cables.append(Cable(knots[2],knots[3],0.02j,0.05+0.25j))
	cables.append(Cable(knots[2],knots[4],0.01j,0.08+0.4j))
	cables.append(Cable(knots[3],knots[4],0.075j,0.1+0.5j))


	print("----------------------------------")
	knotenadmittanzmatrix = np.zeros((len(knots), len(knots)), dtype=complex)

	for i in range(len(knots)):
		y_ii = 0
		for cable in cables:
			if cable.from_knot.index == i or cable.to_knot.index == i:
				y_ii += cable.admit + 1/cable.imp
		knotenadmittanzmatrix[i,i] = y_ii

	for cable in cables:
		knotenadmittanzmatrix[cable.from_knot.index, cable.to_knot.index] = -1/cable.imp
		knotenadmittanzmatrix[cable.to_knot.index, cable.from_knot.index] = -1/cable.imp
	
	np.set_printoptions(precision=2)
	print(knotenadmittanzmatrix)
	print("----------------------------------")
	convergence_border = 0.005

	n_con = 1
	for i in range(30):
		if n_con == len(knots):
			print("Konvergenzkriterium erreicht")
			break
		for knot in knots:

			print("Knoten:",knot.index)
			U = 0
			if knot.type == 0 or knot.converged:

				continue

			if knot.type == 2:
				Q = np.conjugate(R2P(knot.amp * np.exp(1j *knot.angle)))
				A = 0
				for aknot in knots:

					A += knotenadmittanzmatrix[knot.index, aknot.index] * aknot.amp * np.exp(1j* aknot.angle)

				knot.genQ = - 100 * np.imag(A*Q) 
				print("Blindleistung von Knoten:", knot.index, knot.genQ)
			temp = knot.amp
			n = np.conjugate(R2P(knot.amp * np.exp(1j * knot.angle)))

			#U = (knot.conP - 1j * knot.conQ)/n
			U = (-knot.conP/100 + 1j * knot.conQ/100)/n
			for aknot in knots:
				if aknot.index == knot.index:
					continue

				U -= knotenadmittanzmatrix[knot.index, aknot.index] * aknot.amp * np.exp(1j* aknot.angle)

			U = U*1/knotenadmittanzmatrix[knot.index, knot.index]
			knot.amp = np.abs(U)
			knot.angle = np.angle(U, deg=False)
			if np.abs(knot.amp-temp) <= convergence_border:
				print("Knoten", knot.index, "konvergiert")
				knot.converged = True
				n_con += 1
			print("Iterationsergebnis:",knot.amp, np.angle(U, deg=True), "Konvergenz",knot.amp-temp)
			print("--------------------------")


	# rounding

	for knot in knots:
		knot.amp = np.around(knot.amp, decimals=3) 
		knot.angle = np.around(knot.angle, decimals=4) 
		knot.genQ = np.around(knot.genQ, decimals=2) 
		print(knot)




main()