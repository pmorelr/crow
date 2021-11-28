from cylinderModule import *

# shift value
sigma = 0.2-0.7j
# number of eigenmodes to compute
k = 10
kk = 1

print(cpp.la.linear_algebra_backends())
cylinder = LIAproblem2D(100,"./Results/DNS/DNS/")
# cylinder.steady_state()
print('teste')
cylinder.eigenvalues(sigma=sigma, k=k, kk=kk)
