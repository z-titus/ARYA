'''
Script to read input and output data and convert it to a vtk rectilinear data file .vtr
for post processing.
'''
import numpy as np
import vtk
from vtk.util import numpy_support as vtk_np

# list data outputs
sclrs = [r'u (m/s)', r'v (m/s)']
dataflags = ['u', 'v']
time_series = False

# read data input and output
wsl_prfx = r'\\wsl.localhost\Ubuntu\home\zach_t\me6313\ARYA\_case'
file_in = r'\sliding_lid\lid.in' # in (Nx, Ny, Nz, Lx, Ly, Lz, Area, k, BCN, BCE, BCS, BCW, max_it, res_tol)
files_out = r'\sliding_lid\data\Re_400' # output data - directory if its a series

post_prfx = r'\\wsl.localhost\Ubuntu\home\zach_t\me6313\post'
vtr_fname = r'\sliding_lid\tests' # a directory if outputting a series

in_param_flags = ['Nx', 'Ny', 'Lx', 'Ly']

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
    

def data_vtr_cnvrtr(in_file, out_file, main_it, data):
    sclr_data = []

    read_file_in(in_file)
    read_file_out(sclr_data, out_file)

    Nx = int(in_param[0])
    Ny = int(in_param[1])
   
    Lx = in_param[2]
    Ly = in_param[3]

    
    if data == 'u':
        xdim = Nx+1
        ydim = Ny+2

        X0 = Lx/Nx
        Y0 = (Ly/Ny)/2
    if data == 'v':
        xdim = Nx+2
        ydim = Ny+1

        X0 = (Lx/Nx)/2
        Y0 = Ly/Ny
    if data == 'sclr':
        xdim = Nx+2
        ydim = Ny+2

        X0 = (Lx/Nx)/2
        Y0 = (Ly/Ny)/2

    grid_data = np.empty(xdim*ydim*1)

    grid_data = np.array(sclr_data, dtype=np.float64)
    grid_data = np.round(grid_data, decimals=1)
    print(grid_data.shape)

    grid = vtk.vtkRectilinearGrid()

    grid.SetDimensions(int(xdim), int(ydim), int(1)) # number of points in each direction
  
    # build grid for scalar array in physical units (m) as numpy arrays
    sclrXgrid = np.zeros(xdim, dtype = float)
    sclrYgrid = np.zeros(ydim, dtype = float)

    expected = xdim * ydim * 1
    actual = len(sclr_data)

    if actual != expected:
        raise ValueError(f"Data size mismatch: expected {expected}, got {actual}")


    sclrXgrid[0], sclrXgrid[-1]  = 0, Lx
    sclrYgrid[0], sclrYgrid[-1]  = 0, Ly

    for i in range(1, Nx+1):
    
        sclrXgrid[i] = X0
        X0 = X0 + Lx/Nx

    for i in range(1, Ny+1):
      
        sclrYgrid[i] = Y0
        Y0 = Y0 + Ly/Ny
    

    sclrXgrid = vtk_np.numpy_to_vtk(sclrXgrid)
    sclrYgrid = vtk_np.numpy_to_vtk(sclrYgrid)

    grid.SetXCoordinates(sclrXgrid)
    grid.SetYCoordinates(sclrYgrid)

    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(1) # this is 3 for a vector
    array.SetNumberOfTuples(grid.GetNumberOfPoints())
   
    vtk_array = vtk_np.numpy_to_vtk(np.array(grid_data), deep=True)
    vtk_array.SetName(sclrs[main_it])
    grid.GetPointData().AddArray(vtk_array)
    

    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetInputData(grid)
    writer.SetDataModeToBinary()  # write data in binary - issues with ASCII format
    writer.EncodeAppendedDataOff()
    
    tmp = r'\data_{}.vtr'.format(data)
    print(post_prfx+vtr_fname+tmp)
    writer.SetFileName(post_prfx+vtr_fname+tmp)
    writer.Write()

def vector_vtr_cnvrtr(in_file, out_file, main_it, data):
    u = []
    v = []

    read_file_in(in_file)
    read_file_out(u, out_file[0])
    read_file_out(v, out_file[1])

    Nx = int(in_param[0])
    Ny = int(in_param[1])
   
    Lx = in_param[2]
    Ly = in_param[3]

    # center both sets of data so they are colocated
    u = np.array(u)
    u = np.reshape(u, (Nx+1, Ny+2))
    u_c = 0.5 * (u[0:Nx, :] + u[1:Nx+1, :])

    v = np.array(v)
    v = np.reshape(u, (Nx+2, Ny+1))
    v_c = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny+1])

    N = Nx * Ny

    vector_data = np.zeros((N, 3))
    vector_data[:, 0] = u_c.flatten()
    vector_data[:, 1] = v_c.flatten()
    vector_data[:, 2] = 0.0


    grid_data = np.empty((Nx+2))*(Ny+2))*1)

    grid_data = np.array(sclr_data, dtype=np.float64)
    grid_data = np.round(grid_data, decimals=1)
    print(grid_data.shape)

    grid = vtk.vtkRectilinearGrid()

    vtk_vector = vtk_np.numpy_to_vtk(vector_data, deep=True)

    vtk_vector.SetNumberOfComponents(3)
    vtk_vector.SetName("U")  

    grid.SetDimensions(int(Nx+2), int(Ny+2), int(1)) # number of points in each direction
  
    # build grid for scalar array in physical units (m) as numpy arrays
    sclrXgrid = np.zeros(Nx+2, dtype = float)
    sclrYgrid = np.zeros(Ny+2, dtype = float)

 
    X0 = (Lx/Nx)/2
    Y0 = (Ly/Ny)/2

    sclrXgrid[0], sclrXgrid[-1]  = 0, Lx
    sclrYgrid[0], sclrYgrid[-1]  = 0, Ly

    for i in range(1, Nx+1):
    
        sclrXgrid[i] = X0
        X0 = X0 + Lx/Nx

    for i in range(1, Ny+1):
      
        sclrYgrid[i] = Y0
        Y0 = Y0 + Ly/Ny
    

    sclrXgrid = vtk_np.numpy_to_vtk(sclrXgrid)
    sclrYgrid = vtk_np.numpy_to_vtk(sclrYgrid)

    grid.SetXCoordinates(sclrXgrid)
    grid.SetYCoordinates(sclrYgrid)

    #array = vtk.vtkDoubleArray()
    vtk_vector = vtk_np.numpy_to_vtk(vector_data, deep=True)

    vtk_vector.SetNumberOfComponents(3)
    vtk_vector.SetName("U")  
    #array.SetNumberOfComponents(3) # this is 3 for a vector
    #array.SetNumberOfTuples(grid.GetNumberOfPoints())
   
    vtk_array = vtk_np.numpy_to_vtk(np.array(grid_data), deep=True)
    vtk_array.SetName('U')
    grid.GetPointData().SetVectors(vtk_vector)
    

    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetInputData(grid)
    writer.SetDataModeToBinary()  # write data in binary - issues with ASCII format
    writer.EncodeAppendedDataOff()
    
    tmp = r'\data_{}.vtr'.format(data)
    print(post_prfx+vtr_fname+tmp)
    writer.SetFileName(post_prfx+vtr_fname+tmp)
    writer.Write()
       

if __name__ == '__main__':
    import os
    os_list = sorted(os.listdir(wsl_prfx+files_out))
    print(os_list)
    for i in range(len(os_list)):
        out_file_tmp = os.path.join(wsl_prfx + files_out, os_list[i])
        #data_vtr_cnvrtr(file_in, out_file_tmp, main_it = i, data = dataflags[i])
