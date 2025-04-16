import h5py 
import imageio.v3 as iio
import matplotlib.pyplot as plt
import numpy as np
import os

path: str        = "../3_convert_to_hdf5/hdf5_files"
fileName: str    =  lambda index : path + "/" + f"naca4412_5deg_200k.{index:05d}.h5.uvw"
path_frames: str = "./frames"

# Load grid information
with h5py.File(path + "/naca4412_5deg_200k.grid.h5", 'r') as f:
    Rec    = np.array(f['Rec'])     # Reynolds number (nu = 1/Rec)
    utz    = np.array(f['utz'])     # Friction velocity
    xc_vec = np.array(f['xc_vec'])  # Chord-wise distribution
    xs_vec = np.array(f['xs_vec'])  # Wall-tangent distribution
    yn_vec = np.array(f['yn_vec'])  # Wall-normal distribution
    z_vec  = np.array(f['z_vec'])   # Spanwise distribution
    alpha  = np.array(f['alpha'])   # Angle of the profile for projection
    x_pts  = np.array(f['x_pts'])   # X points (absolute frame)
    y_pts  = np.array(f['y_pts'])   # Y points (absolute frame)
    z_pts  = np.array(f['z_pts'])   # Z points (absolute frame)
Nx = len(xs_vec)
Ny = len(yn_vec)
Nz = len(z_vec)
print(f"Ny:{Ny}, Nx:{Nx}, Nz:{Nz}")
XC, YN = np.meshgrid(xc_vec, yn_vec)

# Load velocity field
start_index = 1
num_files = 200
for i in range(num_files):
    index = start_index + i
    with h5py.File(fileName(index), 'r') as f:
        time = np.array(f['time'])
        vx   = np.array(f['vx'])
        vy   = np.array(f['vy'])
        vz   = np.array(f['vz'])
    
    # Check size of matrices
    _ny, _nx, _nz = vx.shape
    if (_nx != Nx or _ny != Ny or _nz != Nz):
        raise Exception("Error: dimensions of vx for index {index} do not match the loaded grid.")

    print(f"File {index}, {time*1000:.0f} ms")

    # Plot a frame
    plt.rcParams.update({'font.size': 28})
    
    fig = plt.figure(figsize=(12,3))
    h   = plt.contourf(XC, YN, vx[:,:,0], levels=36, vmin=0, vmax=1.4, cmap="turbo")
    plt.xlabel("x/c")
    plt.ylabel("n/c")
    plt.title("$u_t=U_t+u_t'$" + f", time: {time:.3f}s")
    plt.axis('scaled')
    plt.tight_layout()
    fig.savefig(path_frames+"/"+f"frame_{time*1000:.0f}ms.png",dpi=400)


# Create a movie
frame_files = [f for f in os.listdir(path_frames) if os.path.isfile(f"{path_frames}/{f}")]
frame_files.sort()
frames = []
for i_frame_file in frame_files:
    print(i_frame_file)
    frames.append(iio.imread(path_frames+"/"+i_frame_file))
iio.imwrite('movie.gif', frames, duration=100) #duration in ms
        
