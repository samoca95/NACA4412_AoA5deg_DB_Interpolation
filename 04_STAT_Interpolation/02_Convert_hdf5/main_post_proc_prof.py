import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.io as sio

from load_int_file import load_int_file

#------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------
case_type: str = "suction"  # suction, pressure or wake, use none to avoid projection

path_int_fields: str = "../01_NEK/INTERP"
path_output: str     = "./hdf5_files"

flag_save_h5: bool   = True

#------------------------------------------------------------------------
# Import the interpolation grid
#------------------------------------------------------------------------
print("Imported interpolation grid:", flush=True)
load_mat: dict = sio.loadmat("../00_Generate_grid/int_stat_mesh_suction.mat")
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

#------------------------------------------------------------------------
# Load the fields
#------------------------------------------------------------------------
print("")
print("Loading input...", flush=True)
out = load_int_file(case_type, path_int_fields, "int_fld", x_pts, y_pts, alpha, yn, xc)

#------------------------------------------------------------------------
# Save in .h5 format
#------------------------------------------------------------------------
print("")
print("Saving input...", flush=True)
if flag_save_h5:
    _fileName: str       = f"naca4412_5deg_200k.sts.h5"
    fileName_interp: str = path_output + "/" + _fileName
    with h5py.File(fileName_interp,'w') as fInterp:
        fInterp.create_dataset("Rec", data=Rec)  # Reynolds number
        fInterp.create_dataset("utz", data=utz)  # Friction velocity at x/c=0.4

        fInterp.create_dataset("xc_vec", data=xc)  # Points in x in the chord-wise reference frame
        fInterp.create_dataset("xs_vec", data=xs)  # Points in x in the wall reference frame
        fInterp.create_dataset("yn_vec", data=yn)   # Points in y in the wall reference frame

        fInterp.create_dataset("Umean", data=out['U'])
        fInterp.create_dataset("Vmean", data=out['V'])
        fInterp.create_dataset("Wmean", data=out['W'])

        fInterp.create_dataset("Pmean", data=out['P'])

    print(f"  * Saved!", flush=True)

#------------------------------------------------------------------------
# Plot interpolation
#------------------------------------------------------------------------
XC, YN = np.meshgrid(xc, yn)

plt.rcParams.update({'font.size': 28})

# Velocity field
fig = plt.figure(figsize=(12,5))
h = plt.pcolormesh(XC, YN, out['U'], edgecolors=("k", 0), cmap="turbo")
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("mean U")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("Umean.png",dpi=400)

fig = plt.figure(figsize=(12,5))
h = plt.pcolormesh(XC, YN, out['V'], edgecolors=("k", 0), cmap="turbo")
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("mean V")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("Vmean.png",dpi=400)

fig = plt.figure(figsize=(12,5))
h = plt.pcolormesh(XC, YN, out['W'], edgecolors=("k", 0), cmap="turbo")
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("mean W")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("Wmean.png",dpi=400)

# Pressure
fig = plt.figure(figsize=(12,5))
h = plt.pcolormesh(XC, YN, out['P'], edgecolors=("k", 0), cmap="turbo")
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("mean P")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("Pmean.png",dpi=400)