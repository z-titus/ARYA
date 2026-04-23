'''
Script to read input and output data and convert it to a vtk rectilinear data file .vtr
for post processing.
'''
import numpy as np
import vtk
from vtk.util import numpy_support as vtk_np

# list data outputs
sclrs = [r'Energy (J/kg)']
BC_add = True # add boundary conditions to the grid

# read data input and output
wsl_prfx = r'\\wsl.localhost\Ubuntu\home\zach_t\me6313\ARYA\_case'
file_in = r'\energy\energy.in' # in (Nx, Ny, Lx, Ly, Area, k, BCN, BCE, BCS, BCW, max_it, res_tol)
files_out = [r'\energy\energy.out'] # output data

post_prfx = r'\\wsl.localhost\Ubuntu\home\zach_t\me6313\post'
vtr_fname = r'\energy\energy_adv_diff_PeX_10_Y_1.vtr'

in_param_flags = ['Nx', 'Ny', 'Lx', 'Ly']
if 'Energy (J/kg)' in sclrs:
    in_param_flags.append('BCNT')
    in_param_flags.append('BCET')
    in_param_flags.append('BCST')
    in_param_flags.append('BCWT')

in_param = []

sclr_data = []

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

def read_file_out(out_files):
    for filename in out_files:
        with open(wsl_prfx+filename, 'r') as file:
            for line in file:
                sclr_data.append(float(line))

def bld_grid_data(Nx, Ny, BCN, BCE, BCS, BCW, sclr_data, BC_add, grid_data):
    if BC_add:
        # first row is south BC
        grid_data[0:Nx+2] = BCS
        n = 0
        
        for i in range(1,Ny+1):
            x_strt = Nx+n+2
            grid_data[x_strt] = BCW
            grid_data[x_strt+1:x_strt+Nx+1]= sclr_data[(i-1)*Nx:i*Nx]
            grid_data[x_strt+Nx+1] = BCE
            
            n += Nx + 2 
        grid_data[(Nx+2)*(Ny+2)-(Nx+2):(Nx+2)*(Ny+2)] = BCN
        
    else:
        grid_data = sclr_data
        pass

    return grid_data
    

def main():
    read_file_in(file_in)
    read_file_out(files_out)

    Nx = int(in_param[0])
    Ny = int(in_param[1])
    Lx = in_param[2]
    Ly = in_param[3]

    if 'Energy (J/kg)' in sclrs:
        BCN = (in_param[4]+273)*1005
        BCE = (in_param[5]+273)*1005
        BCS = (in_param[6]+273)*1005
        BCW = (in_param[7]+273)*1005

    grid_data = np.empty((Nx+2)*(Ny+2))
    grid_data = bld_grid_data(Nx, Ny, BCN, BCE, BCS, BCW, sclr_data, BC_add, grid_data)


    grid = vtk.vtkRectilinearGrid()

    grid.SetDimensions(int(Nx+2), int(Ny+2), int(1)) # number of points in each direction
  
    # build grid for scalar array in physical units (m) as numpy arrays
    sclrXgrid = np.zeros(Nx+2, dtype = float)
    sclrYgrid = np.zeros(Ny+2, dtype = float)

    sclrXgrid[0], sclrXgrid[-1]  = 0, Lx
    sclrYgrid[0], sclrYgrid[-1]  = 0, Ly

    X0 = (Lx/Nx)/2
    Y0 = (Ly/Ny)/2
    for i in range(1, Nx+1):
        # if i == 0:
        #     sclrXgrid[i] = 0
        # elif i == Nx+2:
        #     sclrXgrid = X0 + Lx/2
        # else:
        sclrXgrid[i] = X0

        X0 = X0 + Lx/Nx

    for i in range(1, Ny+1):
        # if i == 0:
        #     sclrYgrid[i] = 0
        # elif i == Nx+2:
        #     sclrYgrid = Y0 + Ly/2
        # else:
        sclrYgrid[i] = Y0

        Y0 = Y0 + Ly/Ny
    

    sclrXgrid = vtk_np.numpy_to_vtk(sclrXgrid)
    sclrYgrid = vtk_np.numpy_to_vtk(sclrYgrid)

    grid.SetXCoordinates(sclrXgrid)
    grid.SetYCoordinates(sclrYgrid)
    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(1) # this is 3 for a vector
    array.SetNumberOfTuples(grid.GetNumberOfPoints())
    for i in range(grid.GetNumberOfPoints()):
        array.SetValue(i, grid_data[i])
    
    grid.GetPointData().AddArray(array)
    print(grid)
    array.SetName(sclrs[0])

    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetFileName(post_prfx+vtr_fname)
    writer.SetInputData(grid)
    writer.Write()

if __name__ == '__main__':
    main()
