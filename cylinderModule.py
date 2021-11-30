from dolfin import *
from spaces import *
import sys
import os

import math
import numpy as np
import scipy.sparse as sps
import scipy.io as sio
import scipy.sparse.linalg as la
import time
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

######################################### Cartesian Homogeneous 2 Components ####################################
class LIAproblem2D():
   def __init__(self, Re, respath):
      
      self.Re = Re
      self.respath = respath
      os.makedirs(respath,exist_ok=True) # create dir if necessary
      
      #teste de git

      self.X = X   #from spaces.py
      self.X3 = X3 
      self.Uv = Uv
      self.Up = Up
      self.Uw = Uw
      
      print('Number of degrees of freedom = '+str(self.X.dim()))
      # functions that decide whether any mesh point x lies on a boundary, and on which:
      #def cylinder1(x, on_boundary):
      #   return ((x[0]+0.25)**2+x[1]**2)**0.5 < 0.01001 and on_boundary
      #def cylinder2(x, on_boundary):
      #   return ((x[0]-0.25)**2+x[1]**2)**0.5 < 0.01001 and on_boundary
      def left(x, on_boundary):       
         return x[0] < (-1.9999) and on_boundary
      #def outlet(x, on_boundary):         
         #return x[0] > (0.9999) and on_boundary
      def upper(x, on_boundary):
         return x[1] > 1.9999 and on_boundary 
      def lower(x, on_boundary):
         return x[1] < -1.9999 and on_boundary
      def symmetry(x, on_boundary):
         return x[0] > -0.0001 and on_boundary
      # define boundary conditions for baseflow & DNS
    
      
      #cyl1   = Expression(("1*x[1]","-1*x[0]"), degree=2)
      #cyl2   = Expression(("-1*x[1]","1*x[0]"), degree=2)
      
      #bcs_cyl1=DirichletBC(X.sub(0), cyl1, cylinder1)
      #bcs_cyl2=DirichletBC(X.sub(0), cyl2, cylinder2)
      #coreu  = Expression(("1/x[1]","-1/x[0]"), degree=2)
      #bcu    = DirichletBC(X.sub(0), coreu)
      
      #bcs_left = DirichletBC(X.sub(1), 0, left)
      #bcs_upper = DirichletBC(X.sub(1), 0, upper)
      #bcs_lower = DirichletBC(X.sub(1), 0, lower)
      bcs_symmetry = DirichletBC(X.sub(0).sub(0), Constant(0.0), symmetry)
      self.bcs = [bcs_symmetry]#bcs_left,bcs_upper,bcs_lower]
      
      # define boundary conditions for perturbations
     # comentados:
      #bcp_cyl=DirichletBC(X.sub(0), (0,0), cylinder)
      #bcp_left= DirichletBC(X.sub(1), 0, left)
     # bcp_upper = DirichletBC(X.sub(1), 0, upper)
     # bcp_lower = DirichletBC(X.sub(1), 0, lower)
      bcp_symmetry = DirichletBC(X.sub(0).sub(1), Constant(0.0), symmetry)
      self.bcp = [bcp_symmetry]#bcp_upper,bcp_lower,bcp_symmetry]



# solves steady state equations
   def steady_state(self):
      up = Function(self.X) # function object representing the solution
      #Uinit = Expression(("1.", "0.", "0."),degree=1)
      #up = interpolate(Uinit, self.X)
      
      dup = TrialFunction(self.X)
      vp = TestFunction(self.X) 
      
      Re = Constant(self.Re)
      # define the nonlinear problem
      u  = as_vector((up[0],up[1]))   # velocity
      p  = up[2]                      # pressure
      v  = as_vector((vp[0],vp[1]))   # velocity
      q  = vp[2]                      # pressure
      F   =  inner(grad(u)*u, v)*dx      \
          + 1/Re*inner(grad(u), grad(v))*dx \
          - div(v)*p*dx - q*div(u)*dx
      # define its Jacobian
      dF  = derivative(F, up, dup) 
      problem = NonlinearVariationalProblem(F, up, self.bcs, dF)
      # solve the problem with Newton
      solver  = NonlinearVariationalSolver(problem)
      solver.parameters['newton_solver']['linear_solver']  = 'petsc'
      solver.parameters['newton_solver']['relative_tolerance']  = 1e-14
      solver.parameters['newton_solver']['absolute_tolerance']  = 1e-14
      solver.parameters['newton_solver']['maximum_iterations']  = 10
      solver.solve()
      # write to file
      File(self.respath+"baseflow.xml") << up.vector()
      u,p = up.split()
      u.rename("v","velocity");  p.rename("p","pressure");
      print("Saving vtk files baseflow_u.pvd, baseflow_p.pvd ")
      File(self.respath+"baseflow_u.pvd") << u
      File(self.respath+"baseflow_p.pvd") << p
      
      cs=plot(u.sub(0))
      plt.show()
      
           
#---------------------------------------------------------------------------------------#      
   # Returns dof indices which are free
   # freeinds = free indices of velocity, temperature, pressure
   # pinds    = free indices of pressure
   def get_indices(self):
      # Collect all dirichlet boundary dof indices
      bcinds = []
      for b in self.bcp:
         bcdict = b.get_boundary_values()
         bcinds.extend(bcdict.keys())

      # total number of dofs
      N = self.X.dim()
      #bcp_left
      # indices of free nodes
      freeinds = np.setdiff1d(range(N),bcinds,assume_unique=True).astype(np.int32)

      # pressure indices
      pinds = self.X.sub(1).dofmap().dofs()

      return freeinds, pinds
#----------------------------------------------------------------------------------------
   # Compute k eigenvalues/vectors   
   def eigenvalues(self, sigma, k, kk):
    parameters['linear_algebra_backend'] = 'Eigen'
    
    Re = Constant(self.Re)
    self.kk=kk
    kk = Constant(self.kk)

    # load baseflow
    ups = Function(self.X)
    File(self.respath+"solution.xml") >> ups.vector()
    us= as_vector((ups[0],ups[1]));

    u,w,p = TrialFunctions(self.X3)
    v,j,q = TestFunctions(self.X3)

    # define RHS matrix
    Ma = assemble(inner(u,v)*dx + inner(w,j)*dx) ###mod
    # Convert to sparse format
    rows, cols, values = as_backend_type(Ma).data()
    Ma = sps.csr_matrix((values, cols, rows))
    
    # define LHS matrix
    Aform = - inner(grad(u)*us, v)*dx - inner(grad(us)*u, v)*dx \
          + div(v)*p*dx - 1/Re*(inner(grad(u), grad(v))-kk**2*dot(u,v))*dx \
          - inner(dot(grad(w),us),j)*dx - 1/Re*(inner(grad(w), grad(j))-kk**2*w*j)*dx - kk*p*j*dx \
          + (div(u)-kk*w)*q*dx
    Aa = assemble(Aform)
    # Convert to sparse format
    rows, cols, values = as_backend_type(Aa).data()
    Aa = sps.csr_matrix((values, cols, rows))

    # remove Dirichlet points from the system
    freeinds,pinds = self.get_indices() 
    M = Ma[freeinds,:][:,freeinds]
    A = Aa[freeinds,:][:,freeinds]
    
    
    if k>0:
      # Compute eigenvalues/vectors of (A,M)
      print("Computing eigenvalues/vectors ...")
      ncv = np.max([10,2*k]) # number of Krylov vectors
      vals, vecs = la.eigs(A, k=k, M=M, sigma=sigma,  ncv=ncv, maxiter=40, tol=10e-10) 
      
      file = open(self.respath+"evals.dat","w")
        
      #lista associando autovalores e autovetores  
      valvec = []  
      for jj in range(len(vals)):
            im_val = np.imag(vals[jj])
            valvec.append([jj, vals[jj], im_val])
            
      valvec.sort(key = lambda x: x[2])
      
            
      for val in valvec:
          print(val[2], val[1])
          file.write("%s\n" % val)
      file.close()   
      
      # only writing real parts of eigenvectors to file
      ua = Function(self.X)
      for i in range(0,k):
          ua.vector()[freeinds] = vecs[:,i]
          File(self.respath+"evec"+str(i+1)+".xml") << ua.vector()
          u,p  = ua.split()
          File(self.respath+"evec_u_"+str(i+1)+".pvd") << u
          File(self.respath+"evec_p_"+str(i+1)+".pvd") << p
      
    
    
       # Compute k eigenvalues/vectors   
   def eigenvalues_k(self, sigma, k, kk_list):
    im_val = []
    for kk in kk_list:
        parameters['linear_algebra_backend'] = 'Eigen'

        Re = Constant(self.Re)
        self.kk=kk
        kk = Constant(self.kk)

        # load baseflow
        ups = Function(self.X)
        File(self.respath+"solution.xml") >> ups.vector()
        us= as_vector((ups[0],ups[1]));

        u,w,p = TrialFunctions(self.X3)
        v,j,q = TestFunctions(self.X3)

        # define RHS matrix
        Ma = assemble(inner(u,v)*dx + inner(w,j)*dx) ###mod
        # Convert to sparse format
        rows, cols, values = as_backend_type(Ma).data()
        Ma = sps.csr_matrix((values, cols, rows))

        # define LHS matrix
        Aform = - inner(grad(u)*us, v)*dx - inner(grad(us)*u, v)*dx \
              + div(v)*p*dx - 1/Re*(inner(grad(u), grad(v))-kk**2*dot(u,v))*dx \
              - inner(dot(grad(w),us),j)*dx - 1/Re*(inner(grad(w), grad(j))-kk**2*w*j)*dx - kk*p*j*dx \
              + (div(u)-kk*w)*q*dx
        Aa = assemble(Aform)
        # Convert to sparse format
        rows, cols, values = as_backend_type(Aa).data()
        Aa = sps.csr_matrix((values, cols, rows))

        # remove Dirichlet points from the system
        freeinds,pinds = self.get_indices() 
        M = Ma[freeinds,:][:,freeinds]
        A = Aa[freeinds,:][:,freeinds]


        if k>0:
          # Compute eigenvalues/vectors of (A,M)
          ncv = np.max([10,2*k]) # number of Krylov vectors
          vals, vecs = la.eigs(A, k=k, M=M, sigma=sigma,  ncv=ncv, maxiter=40, tol=10e-10) 

          im_vals1 = [] 
          for jj in range(len(vals)):
                im_vals1.append(np.imag(vals[jj]))

          max_im_val = max(im_vals1)
        
        im_val.append(max_im_val)
        
    plt.figure(figsize=(7, 9))
    
    plt.title('Growth rate vs Wave number', fontdict=font)
    plt.xlabel('Wave number value', fontdict=font)
    plt.ylabel('Imaginary component of the eigenvalue', fontdict=font)
    plt.plot(kk_list, im_val, 'k')
 
    plt.savefig('k_im.png', dpi= 300)
        
    

   
    
          
#--------------------------------------------------------------------------

   def nonlinearTimestepper_Euler(self, dt, tend):
      dnspath = self.respath+"DNS/"
      Re = Constant(self.Re)
      u,p = TrialFunctions(self.X)
      v,q = TestFunctions(self.X)
      up = Function(self.X)

      # initialize
      t = 0.
      it = 0
      ifile = 0
            #U_init = Expression(("1.","0.05*exp(-(x[0]-0.5)*(x[0]-0.5)/0.2 - x[1]*x[1]/0.2)","0"),  degree=2)
      #U_init = Expression(("x[1]*(sqrt(x[0]*x[0]+x[1]*x[1])<0.5)+(0.5*0.5*x[1])/(x[0]*x[0]+x[1]*x[1])*(sqrt(x[0]*x[0]+x[1]*x[1])>0.5)","-x[0]*(sqrt(x[0]*x[0]+x[1]*x[1])<0.5)-(0.5*0.5*x[0])/(x[0]*x[0]+x[1]*x[1])*(sqrt(x[0]*x[0]+x[1]*x[1])>0.5)","0."),  degree=2)
      #U_init = Expression(("-om*x[1]*(sqrt((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])<a) -(om*a*a*x[1])/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])*(sqrt((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])>a) +om*x[1]*(sqrt((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])<a) +(om*a*a*x[1])/((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])*(sqrt((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])>a)","om*x[0]*(sqrt((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])<a) +(om*a*a*x[0])/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])*(sqrt((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])>a) -om*x[0]*(sqrt((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])<a) -(om*a*a*x[0])/((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])*(sqrt((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])>a)","0."),om=1., a=0.2, b=1., degree=2)
      #U_init = Expression(("-om*a*a*x[1]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])*(1-exp(-1*((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])/(a*a)) + om*a*a*x[1]/((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])*(1-exp(-1*((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])/(a*a)))","om*a*a*x[0]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])*(1-exp(-1*((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])/(a*a))) - om*a*a*x[0]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])*(1-exp(-1*((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])/(a*a)))","0."), om=1., a=0.1, b=1., degree=2)
      U_init = Expression(("-om*a*a*x[1]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1] + 1e-3)*(1-exp(-1*((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])/(a*a)))","om*a*a*x[0]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1]+ 1e-3)*(1-exp(-1*((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])/(a*a)))","0."), om=1., a=0.1, b=-1., degree=2)
         
      #U_init = Expression(("(-om*a*a*x[1]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1] + 1e-3)*(1-exp(-1*((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])/(a*a))) * ((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1] > 0.01) + om*a*a*x[1]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])*(1-exp(-1*((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])/(a*a))))*((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1] > 0.01)","om*a*a*x[0]/((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])*(1-exp(-1*((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1])/(a*a)))* ((x[0]-b/2)*(x[0]-b/2)+x[1]*x[1] > 0.01) -om*a*a*x[0]/((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])*(1-exp(-1*((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1])/(a*a)))*((x[0]+b/2)*(x[0]+b/2)+x[1]*x[1] > 0.01)","0."), om=1., a=0.1, b=1., degree=2)

      up = interpolate(U_init, X)
#     File("solution.xml") >> up.vector()
      upast = as_vector((up[0], up[1]))
      #print("Velocidades:",upast,up)
      plt.figure(figsize=(7, 9))
      plot((up[1]**2+up[0]**2)**(1/2))
      #plt.show()
      plt.savefig('u_mag_ini.png')
  
      F = 1./dt*inner(u - upast,v)*dx \
            + inner( grad(upast)*upast, v )*dx\
            - div(v)*p*dx \
            + 1./Re*inner(grad(u), grad(v))*dx \
            + div(u)*q*dx
      Aform = lhs(F)
      bform = rhs(F)
      A = assemble(Aform)
      [bc.apply(A) for bc in self.bcs]
      solver = LUSolver(A, 'petsc')

      while (t<tend + DOLFIN_EPS):
          print("time t = "+str(np.round(t+dt, decimals=4)))
          b = assemble(bform)
          [bc.apply(b) for bc in self.bcs]
          solver.solve(up.vector(), b)
          t = t+dt
          it = it+1
          if (it%50==0):
              ifile=ifile+1
              u,p = up.split()
              File(dnspath+"u"+f"{ifile:03d}"+".pvd") << u
          print("Saving xml file")
          File(dnspath+"solution{0}.xml".format(t)) << up.vector()
          u,p = up.split()
          print("Saving vtk files")
          File(dnspath+"solution_u{0}.pvd".format(t)) << u
          plt.figure(figsize=(7, 9))
          cs=plot((u.sub(0)**2+u.sub(1)**2)**(1/2))
          #plt.show()
          plt.savefig('u_mag{0}.png'.format(t))
          
          #fig.savefig('fig{0}.jpg'.format(t))
          
      File(dnspath+"solution.xml") << up.vector()
      u,p = up.split()
      print("Saving vtk files")
      File(dnspath+"solution_u.pvd") << u    
               
#----------------------------------------------------------------------------------------

   def nonlinearTimestepper_CN(self, dt, tend):
      Re = Constant(self.Re)
      u,p = TrialFunctions(self.X)
      v,q = TestFunctions(self.X)
      up = Function(self.X)
      up_pred  = Function(self.X)
      up_past2 = Function(self.X)

      # initialize
      t = 0.
      it = 0
      ifile = 0
#      U_init = Expression(("1","0","0"), degree=2)
      U_init = Expression(("1.","0.05*exp(-(x[0]-0.5)*(x[0]-0.5)/0.2 - x[1]*x[1]/0.2)","0"), degree=2)
      up = interpolate(U_init, X)
      
      u_past  = as_vector((up[0], up[1]))
      u_pred  = as_vector((up_pred[0], up_pred[1]))
      
      NSOp = 1./dt*inner(u - u_past,v)*dx \
            + inner( grad(u_pred)*u_pred, v )*dx\
            + 0.5/Re*inner(grad(u), grad(v))*dx \
            + 0.5/Re*inner(grad(u_past), grad(v))*dx \
            - div(v)*p*dx \
            + div(u)*q*dx
      Aform = lhs(NSOp)
      bform = rhs(NSOp)
      A = assemble(Aform)
      [bc.apply(A) for bc in self.bcs]
      CNsolver = LUSolver(A, 'petsc')

      while (t<tend + DOLFIN_EPS):
          print("time t = "+str(np.round(t+dt, decimals=4)))
          if (t<dt):
             up_pred.vector()[:] = up.vector()
#             up_pred = up
          else:
             up_pred.vector()[:] = 1.5*up.vector() - 0.5*up_past2.vector()
#             up_pred = 1.5*up - 0.5*up_past2
          up_past2.vector()[:] = up.vector()
#          up_past2 = up
          b = assemble(bform)
          [bc.apply(b) for bc in self.bcs]
          CNsolver.solve(up.vector(), b)
#          u,p = up.split(True)
          t = t+dt
          it = it+1
          if (it%50==0):
              ifile=ifile+1
              u,p = up.split()
              File("u"+f"{ifile:03d}"+".pvd") << u
          
      print("Saving xml file")
      File("solution.xml") << up.vector()
      u,p = up.split()
      print("Saving vtk files")
      File("solution_u.pvd") << u
#      File("solution_p.pvd") << p
      plot(u.sub(0))
      
           
