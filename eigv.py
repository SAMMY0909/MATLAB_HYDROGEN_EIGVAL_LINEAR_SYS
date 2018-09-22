import numpy
import scipy.io

data=scipy.io.loadmat("Hmat.mat")

H=data['H']

print(sorted(numpy.linalg.eigvals(H[:-1,:]))[:10])
