from dolfin import *

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

mesh = Mesh()
with XDMFFile("Mesh/Mesh.xdmf") as infile:
    infile.read(mesh)
    
#mesh = BoxMesh(Point (-1,-1,0), Point (1,1,-0.1), 10, 10, 1)

#mesh = UnitSquareMesh(100, 100)

# boundaries = MeshFunction("size_t", mesh, "./Mesh/Mesh_boundaries.xdmf")

###############################################################################

V=VectorElement("Lagrange",mesh.ufl_cell(),2)
W=FiniteElement("Lagrange",mesh.ufl_cell(),1)
Q=FiniteElement("Lagrange",mesh.ufl_cell(),1)

TH3=MixedElement([V, W, Q])
TH=MixedElement([V, Q])
X=FunctionSpace(mesh,TH)
X3=FunctionSpace(mesh,TH3)

Uv=FunctionSpace(mesh,V)
Uw=FunctionSpace(mesh,W)
Up=FunctionSpace(mesh,Q)
