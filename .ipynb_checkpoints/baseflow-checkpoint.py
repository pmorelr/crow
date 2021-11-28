from cylinderModule import *

problem = LIAproblem2D(100,"./Results/DNS/")
problem.nonlinearTimestepper_Euler(0.05,1)
