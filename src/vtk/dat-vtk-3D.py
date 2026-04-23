'''
Script to read input and output data and convert it to a vtk rectilinear data file .vtr
for post processing.
'''
import numpy as np
import vtk
from vtk.util import numpy_support as vtk_np

# list data outputs
sclrs = [r'Enthalpy (J/kg)']
time_series = True

# read data input and output
wsl_prfx = r'\\wsl.localhost\Ubuntu\home\zach_t\me6313\ARYA\_case'
file_in = r'\energy\unst_energy3D.in' # in (Nx, Ny, Nz, Lx, Ly, Lz, Area, k, BCN, BCE, BCS, BCW, max_it, res_tol)
files_out = r'\energy\ss_r_0.1' # output data - directory if its a series

post_prfx = r'\\wsl.localhost\Ubuntu\home\zach_t\me6313\post'
vtr_fname = r'\energy\ss_bad' # a directory if outputting a series

in_param_flags = ['Nx', 'Ny', 'Nz', 'Lx', 'Ly', 'Lz']

in_param = []

#sclr_data = []

n = 1

def read_file_in(filename):
    with open(wsl_prfx+filename, 'r') as file:
        for line in file:
            for flag in in_param_flags:
                if flag in line:
                    delim = '='
                    idx = line.index(delim)
                    in_param.append(float(line[idx+1:]))
                else:
                    pass

def read_file_out(sclr_data, out_file):
    with open(out_file, 'r') as file:
        for line in file:
            sclr_data.append(float(line))
    return sclr_data
    

def data_vtr_cnvrtr(in_file, out_file, main_it):
    sclr_data = []

    read_file_in(in_file)
    read_file_out(sclr_data, out_file)

    Nx = int(in_param[0])
    Ny = int(in_param[1])
    Nz = int(in_param[2])
    Lx = in_param[3]
    Ly = in_param[4]
    Lz = in_param[5]
    

    grid_data = np.empty((Nx+2)*(Ny+2)*(Nz+2))
    grid_data = np.array(sclr_data, dtype=np.float64)
    grid_data = np.round(grid_data, decimals=1)
    print(grid_data.shape)

    grid = vtk.vtkRectilinearGrid()

    grid.SetDimensions(int(Nx+2), int(Ny+2), int(Nz+2)) # number of points in each direction
  
    # build grid for scalar array in physical units (m) as numpy arrays
    sclrXgrid = np.zeros(Nx+2, dtype = float)
    sclrYgrid = np.zeros(Ny+2, dtype = float)
    sclrZgrid = np.zeros(Nz+2, dtype = float)

    expected = (Nx+2) * (Ny+2) * (Nz+2)
    actual = len(sclr_data)

    if actual != expected:
        raise ValueError(f"Data size mismatch: expected {expected}, got {actual}")


    sclrXgrid[0], sclrXgrid[-1]  = 0, Lx
    sclrYgrid[0], sclrYgrid[-1]  = 0, Ly
    sclrZgrid[0], sclrZgrid[-1]  = 0, Lz

    X0 = (Lx/Nx)/2
    Y0 = (Ly/Ny)/2
    Z0 = (Lz/Nz)/2
    for i in range(1, Nx+1):
    
        sclrXgrid[i] = X0
        X0 = X0 + Lx/Nx

    for i in range(1, Ny+1):
      
        sclrYgrid[i] = Y0
        Y0 = Y0 + Ly/Ny

    for i in range(1, Nz+1):

        sclrZgrid[i] = Z0
        Z0 = Z0 + Lz/Nz
    

    sclrXgrid = vtk_np.numpy_to_vtk(sclrXgrid)
    sclrYgrid = vtk_np.numpy_to_vtk(sclrYgrid)
    sclrZgrid = vtk_np.numpy_to_vtk(sclrZgrid)

    grid.SetXCoordinates(sclrXgrid)
    grid.SetYCoordinates(sclrYgrid)
    grid.SetZCoordinates(sclrZgrid)
    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(1) # this is 3 for a vector
    array.SetNumberOfTuples(grid.GetNumberOfPoints())
   
    vtk_array = vtk_np.numpy_to_vtk(np.array(grid_data), deep=True)
    vtk_array.SetName(sclrs[0])
    grid.GetPointData().AddArray(vtk_array)
    

    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetInputData(grid)
    writer.SetDataModeToBinary()  # write data in binary - issues with ASCII format
    writer.EncodeAppendedDataOff()
    
    if time_series:
        tmp = r'\data_{}.vtr'.format(main_it)
        print(post_prfx+vtr_fname+tmp)
        writer.SetFileName(post_prfx+vtr_fname+tmp)
        writer.Write()
        #n += 1
    else:
        writer.SetFileName(post_prfx+vtr_fname)
        writer.Write()

if __name__ == '__main__':
    if time_series == False:
        data_vtr_cnvrtr(file_in, wsl_prfx+files_out, 0) # pass output data file name
    else:
        import os
        os_list = sorted(os.listdir(wsl_prfx+files_out))
        print(os_list)
        for i in range(len(os_list)):
            out_file_tmp = os.path.join(wsl_prfx + files_out, os_list[i])
            data_vtr_cnvrtr(file_in, out_file_tmp, i+1)
