import sys
import shutil
import os
import netCDF4 as nc

def add_bounds(file, mesh_varname='Mesh2d', verbose=False):

    ncfile = nc.Dataset(file, 'a', format='NETCDF4')
    if verbose: print (ncfile)
    
    mesh_var = ncfile.variables[mesh_varname]
    if verbose: print(mesh_var)
    
    face_node_connectivity = ncfile.variables[mesh_var.face_node_connectivity]
    if verbose: print (face_node_connectivity)

    for face_coord,node_coord in zip(mesh_var.face_coordinates.split(" "),mesh_var.node_coordinates.split(" ")):

        face_coordvar = ncfile.variables[face_coord]
        if verbose: print (face_coordvar)

        node_coordvar = ncfile.variables[node_coord]
        if verbose: print (node_coordvar)

        face_coord_bnds = f'{face_coord}_bounds'
        face_coordvar.bounds = face_coord_bnds
        bnds = ncfile.createVariable(face_coord_bnds, face_coordvar.dtype, face_node_connectivity.dimensions)
        bnds[:] = node_coordvar[face_node_connectivity[:].flatten()-face_node_connectivity.start_index].reshape(bnds.shape)
        if verbose: print (bnds)

    ncfile.close()

if __name__ == "__main__":

    if len(sys.argv) == 2:
        infile = sys.argv[1]
        name, ext = os.path.splitext(infile)
        outfile = f'{name}_bnds{ext}'
    elif len(sys.argv) > 2:
        infile = sys.argv[1]
        outfile = sys.argv[2]

    shutil.copy2(infile,outfile)
    add_bounds(outfile, verbose=False)
