import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.io as sio

from py_func.mat2dict import mat2dict

#-----------------------------------------------------------------------
# Define the case type
#-----------------------------------------------------------------------
case_type: str     = "suction"  # Options: suction, pressure, wake

path_output: str   = "./hdf5_files"
flag_save_h5: bool = True

#------------------------------------------------------------------------
# Import the interpolation grid
#------------------------------------------------------------------------
print("Imported interpolation grid:", flush=True)
load_mat: dict = sio.loadmat(f"../00_Generate_grid/int_stat_mesh_{case_type}.mat")
int_mesh = load_mat["int_mesh"]

Rec = float(np.squeeze(int_mesh["Rec"][0,0]))
utz = float(np.squeeze(int_mesh["utz"][0,0]))
nu  = 1 / Rec
print(f"{Rec = }")
print(f"{nu = }")
print(f"{utz = }")

x_pts   = np.squeeze(int_mesh["x_pts"][0,0])
y_pts   = np.squeeze(int_mesh["y_pts"][0,0])
alpha   = np.squeeze(int_mesh["alpha"][0,0])
xc      = np.squeeze(int_mesh["xc"][0,0])     # Chord-wise distribution
xs      = np.squeeze(int_mesh["xs"][0,0])     # Wall-tangent distribution
yn      = np.squeeze(int_mesh["yn"][0,0])     # Wall-normal distribution

Nx      = len(xc)
Ny      = len(yn)
Npoints = Nx * Ny
print(f"{Ny = }, {Nx = }", flush=True)

#-----------------------------------------------------------------------
# Load the data
#-----------------------------------------------------------------------
print("Loading statistics file...")
load_mat = sio.loadmat(f"../02_Full_statistics/output/stat_naca4412_200k_{case_type}.mat",squeeze_me=True)
_out = load_mat['out']
print(f"Extracting {case_type} case...")
sts, _ = mat2dict(_out)

# Print the keys of the dictionary
print(f"----------------")
print(f"Loaded statistics:\n{sts.keys()}")
print(f"----------------")

#------------------------------------------------------------------------
# Save in .h5 format
#------------------------------------------------------------------------
if flag_save_h5:
    print(f"----------------")
    print(f"Saving .h5 file...")
    print(f"----------------")
    _fileName: str       = f"naca4412_5deg_200k_{case_type}.sts.h5"
    fileName_interp: str = path_output + "/" + _fileName
    with h5py.File(fileName_interp,'w') as fInterp:
        fInterp.create_dataset("Umean", data=sts['U'])
        fInterp.create_dataset("Vmean", data=sts['V'])
        fInterp.create_dataset("Wmean", data=sts['W'])
        fInterp.create_dataset("uurms", data=np.sqrt(sts['uu']))
        fInterp.create_dataset("vvrms", data=np.sqrt(sts['vv']))
        fInterp.create_dataset("wwrms", data=np.sqrt(sts['ww']))
        for key in sts.keys():
            fInterp.create_dataset(key, data=sts[key])
    print(f"  * Saved!", flush=True)

#-----------------------------------------------------------------------
# Visualize the data
#-----------------------------------------------------------------------
fig_path = f"figures_{case_type}"
try:
    os.makedirs(fig_path)
except FileExistsError:
    pass

XX, YY = np.meshgrid(sts['xa'], sts['yn'][:,0])

# Friction coefficient
fig = plt.figure(figsize=(10, 5))
plt.plot(sts['xa'], sts['Cf'], 'k-')
plt.xlabel('x/c')
plt.ylabel('Cf')
plt.title(f'Friction coefficient for {case_type}')
plt.grid()
plt.tight_layout()
plt.savefig(f'{fig_path}/Cf.png', dpi=300)

# Friction velocity
fig = plt.figure(figsize=(10, 5))
plt.plot(sts['xa'], sts['ut'], 'k-')
plt.xlabel('x/c')
plt.ylabel('u_tau')
plt.title(f'Friction velocity for {case_type}')
plt.grid()
plt.tight_layout()
plt.savefig(f'{fig_path}/utau.png', dpi=300)

# Dissipation rate
fig = plt.figure(figsize=(10, 5))
h = plt.contourf(XX, YY, sts['Dk'], 100, cmap='jet')
plt.colorbar(h)
plt.contour(XX, YY, sts['Dk'], 10, colors='k', linewidths=0.5)
plt.xlabel('x/c')
plt.ylabel('Dk')
plt.title(f'Dissipation rate for {case_type}, avg: {np.mean(sts['Dk'])}')
plt.grid()
plt.axis('scaled')
plt.tight_layout()
plt.savefig(f'{fig_path}/Dk.png', dpi=300)

# Velocities
fig = plt.figure(figsize=(10, 5))
h = plt.contourf(XX, YY, sts['U'], 100, cmap='jet')
plt.colorbar(h)
plt.contour(XX, YY, sts['U'], 10, colors='k', linewidths=0.5)
plt.xlabel('x/c')
plt.ylabel('Umean')
plt.title(f'Umean for {case_type}')
plt.axis('scaled')
plt.tight_layout()
plt.savefig(f'{fig_path}/Umean.png', dpi=300)

fig = plt.figure(figsize=(10, 5))
h = plt.contourf(XX, YY, sts['V'], 100, cmap='jet')
plt.colorbar(h)
plt.contour(XX, YY, sts['V'], 10, colors='k', linewidths=0.5)
plt.xlabel('x/c')
plt.ylabel('Vmean')
plt.title(f'Vmean for {case_type}')
plt.axis('scaled')
plt.tight_layout()
plt.savefig(f'{fig_path}/Vmean.png', dpi=300)

fig = plt.figure(figsize=(10, 5))
h = plt.contourf(XX, YY, sts['W'], 100, cmap='jet')
plt.colorbar(h)
plt.contour(XX, YY, sts['W'], 10, colors='k', linewidths=0.5)
plt.xlabel('x/c')
plt.ylabel('Wmean')
plt.title(f'Wmean for {case_type}')
plt.axis('scaled')
plt.tight_layout()
plt.savefig(f'{fig_path}/Wmean.png', dpi=300)




