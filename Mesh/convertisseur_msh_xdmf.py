# by Antoine Tatin 2020
import numpy
import meshio
import h5py

mesh=meshio.read("Mesh.msh")

##
def create_mesh(mesh, cell_type):
    cells = numpy.vstack([cell.data for cell in mesh.cells if cell.type==cell_type])
    data = numpy.hstack([mesh.cell_data_dict["gmsh:geometrical"][key]
                           for key in mesh.cell_data_dict["gmsh:geometrical"].keys() if key==cell_type])
    mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells},
                               cell_data={"name_to_read":[data]})
    return mesh



def prune_z_0(self, tol=1.0e-13):
        """Remove third (z) component of points if it is 0 everywhere (up to a
        tolerance).
        """
        if self.points.shape[1] == 3 and numpy.all(numpy.abs(self.points[:, 2]) < tol):
            self.points = self.points[:, :2]


triangle = create_mesh(mesh, "triangle")
facet_mesh = create_mesh(mesh, "line")

prune_z_0(triangle)
prune_z_0(facet_mesh)

meshio.write("Mesh.xdmf", triangle)
meshio.write("Mesh_boundaries.xdmf", facet_mesh)

