import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.io as sio

from load_int_file import load_int_file

#------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------
case_type:str = "suction" # pressure, suction, wake
out_start_index: int = 1  # Index of the first file in the global database
in_start_index: int  = 1  # File to start processing (in case you do it by batches), normally 1
n_files: int = 10000      # Number of files to transform to hdf5 (total in the folder)

# path_int_fields: str = "../2_run_interpolation/INTERP"
path_int_fields: str = "../../../NACA4412_5deg_DB/2_INTERPOLATED_RAW_SS/00001_10000"
path_db_fields: str  = "./hdf5_files"

flag_suction_side: bool  = True
flag_pressure_side: bool = False
flag_velocity: bool      = True
flag_vorticity: bool     = True
flag_la2: bool           = False
flag_project: bool       = False
flag_save_h5: bool       = True

#------------------------------------------------------------------------
# Import the interpolation grid
#------------------------------------------------------------------------
print("Imported interpolation grid:", flush=True)
load_mat: dict = sio.loadmat(f"../0_generate_grid/int_mesh_{case_type}.mat")
int_mesh = load_mat["int_mesh"]

Rec = float(np.squeeze(int_mesh["Rec"][0,0]))
utz = float(np.squeeze(int_mesh["utz"][0,0]))
nu  = 1 / Rec
print(f"{Rec = }")
print(f"{nu = }")
print(f"{utz = }")

x_pts   = np.squeeze(int_mesh["x_pts"][0,0])
y_pts   = np.squeeze(int_mesh["y_pts"][0,0])
z_pts   = np.squeeze(int_mesh["z_pts"][0,0])
alpha   = np.squeeze(int_mesh["alpha"][0,0])
xc      = np.squeeze(int_mesh["xc"][0,0])     # Chord-wise distribution
xs      = np.squeeze(int_mesh["xs"][0,0])     # Wall-tangent distribution
yn      = np.squeeze(int_mesh["yn"][0,0])     # Wall-normal distribution
z       = np.squeeze(int_mesh["z"][0,0])      # Spanwise distribution

Nx      = len(xc)
Ny      = len(yn)
Nz      = len(z)
Npoints = Nx * Ny * Nz
print(f"{Ny = }, {Nx = }, {Nz = }", flush=True)

sinalpha = np.sin(alpha)
cosalpha = np.cos(alpha)

#------------------------------------------------------------------------
# Load the fields
#------------------------------------------------------------------------
# Save the grid information
if flag_save_h5:
    if not os.path.exists(path_db_fields):
        os.makedirs(path_db_fields)
    _fileName: str       = f"naca4412_5deg_200k_{case_type}.grid.h5"
    fileName_grid: str = path_db_fields + "/" + _fileName
    with h5py.File(fileName_grid,'w') as fGrid:
        fGrid.create_dataset("Rec", data=Rec)
        fGrid.create_dataset("utz", data=utz)
        fGrid.create_dataset("xc_vec", data=xc)     # Points in x in the chord-wise reference frame
        fGrid.create_dataset("xs_vec", data=xs)     # Points in x in the wall reference frame
        fGrid.create_dataset("yn_vec", data=yn)     # Points in y in the wall reference frame
        fGrid.create_dataset("z_vec",  data=z)      # Points in z in the wall reference frame       
        fGrid.create_dataset("alpha",  data=alpha)
        fGrid.create_dataset("x_pts",  data=x_pts)
        fGrid.create_dataset("y_pts",  data=y_pts)
        fGrid.create_dataset("z_pts",  data=z_pts)
    print("* Saved grid information", flush=True)

print("")
print("Loading input...", flush=True)
for i_file in range(in_start_index,n_files+1):
    in_index = i_file
    fHistoryname: str        = path_int_fields + "/" + "intHistory.f" + f"{in_index:05d}"
    history: np.ndarray[int] = np.loadtxt(fHistoryname, dtype=int)
    ncores: int              = history[-1,0] + 1
    print(f"- File {i_file:05d}, {ncores = }")

    if flag_velocity:
        variable_name: str = "intVX.f"
        bb_buff_vx, time = load_int_file(path_int_fields, variable_name, i_file, history, ncores)

        variable_name: str = "intVY.f"
        bb_buff_vy, _ = load_int_file(path_int_fields, variable_name, i_file, history, ncores)

        variable_name: str = "intVZ.f"
        bb_buff_vz, _ = load_int_file(path_int_fields, variable_name, i_file, history, ncores)

        int_vx = np.reshape(bb_buff_vx,(Ny,Nx,Nz), order='F')
        int_vy = np.reshape(bb_buff_vy,(Ny,Nx,Nz), order='F')
        int_vz = np.reshape(bb_buff_vz,(Ny,Nx,Nz), order='F')

    if flag_vorticity:
        variable_name: str = "intOX.f"
        bb_buff_ox, time = load_int_file(path_int_fields, variable_name, i_file, history, ncores)

        variable_name: str = "intOY.f"
        bb_buff_oy, _ = load_int_file(path_int_fields, variable_name, i_file, history, ncores)

        variable_name: str = "intOZ.f"
        bb_buff_oz, _ = load_int_file(path_int_fields, variable_name, i_file, history, ncores)

        int_ox = np.reshape(bb_buff_ox,(Ny,Nx,Nz), order='F')
        int_oy = np.reshape(bb_buff_oy,(Ny,Nx,Nz), order='F')
        int_oz = np.reshape(bb_buff_oz,(Ny,Nx,Nz), order='F')

    if flag_la2:
        variable_name: str = "intLA.f"
        bb_buff_la, time = load_int_file(path_int_fields, variable_name, i_file, history, ncores)

        int_la = np.reshape(bb_buff_la,(Ny,Nx,Nz), order='F')       

    #------------------------------------------------------------------------
    # Project the xy vectors to the wall reference frame
    #------------------------------------------------------------------------
    if flag_project:
        cos = cosalpha.reshape(1,-1,1)
        sin = sinalpha.reshape(1,-1,1)
        
        if flag_suction_side:
            if flag_velocity:
                _proj_vx = int_vx * cos + int_vy * sin
                _proj_vy = int_vy * cos - int_vx * sin
                int_vx   = _proj_vx
                int_vy   = _proj_vy
            if flag_vorticity:
                _proj_ox = int_ox * cos + int_oy * sin
                _proj_oy = int_oy * cos - int_ox * sin
                int_ox   = _proj_ox
                int_oy   = _proj_oy
        
        if flag_pressure_side:
            if flag_velocity:
                _proj_vx = int_vx * cos  + int_vy * sin
                _proj_vy = int_vx * sin  - int_vy * cos
                int_vx   = _proj_vx
                int_vy   = _proj_vy
            if flag_vorticity:
                _proj_ox = int_ox * cos  + int_oy * sin
                _proj_oy = int_ox * sin  - int_oy * cos
                int_ox   = _proj_ox
                int_oy   = _proj_oy

    #------------------------------------------------------------------------
    # Save in .h5 format
    #------------------------------------------------------------------------
    if flag_save_h5:
        out_index = out_start_index + i_file - 1
        if flag_velocity:  # Velocity vector
            _fileName: str       = f"naca4412_5deg_200k_{case_type}.{out_index:05d}.h5.uvw"
            fileName_interp: str = path_db_fields + "/" + _fileName
            with h5py.File(fileName_interp,'w') as fInterp:
                fInterp.create_dataset("time", data=time)
                fInterp.create_dataset("vx", data=int_vx)
                fInterp.create_dataset("vy", data=int_vy)
                fInterp.create_dataset("vz", data=int_vz)

        if flag_vorticity:  # Vorticity vector
            _fileName: str       = f"naca4412_5deg_200k_{case_type}.{out_index:05d}.h5.omega"
            fileName_interp: str = path_db_fields + "/" + _fileName
            with h5py.File(fileName_interp,'w') as fInterp:
                fInterp.create_dataset("time", data=time)
                fInterp.create_dataset("ox", data=int_ox)
                fInterp.create_dataset("oy", data=int_oy)
                fInterp.create_dataset("oz", data=int_oz)

        if flag_la2:  # Lambda2 scalar field
            _fileName: str       = f"naca4412_5deg_200k_{case_type}.{out_index:05d}.h5.la2"
            fileName_interp: str = path_db_fields + "/" + _fileName
            with h5py.File(fileName_interp,'w') as fInterp:
                fInterp.create_dataset("time", data=time)
                fInterp.create_dataset("la2", data=int_la)

    print(f"  * Saved {_fileName}, time: {time:.3f}s", flush=True)

#------------------------------------------------------------------------
# Plot interpolation
#------------------------------------------------------------------------
XC, YN = np.meshgrid(xc, yn)

plt.rcParams.update({'font.size': 28})

# # Interpolation grid
# fig = plt.figure(figsize=(12,5))
# plt.scatter(x_pts, y_pts, c=int_vx.T, marker=".", cmap="turbo")
# plt.xlabel("x/c")
# plt.ylabel("y/c")
# plt.axis('scaled')
# plt.tight_layout()
# fig.savefig("grid_u.png",dpi=400)

# Velocity field
fig = plt.figure(figsize=(12,5))
h = plt.pcolormesh(XC, YN, int_vx[:,:,0], edgecolors=("k", 0.5), cmap="turbo")#, vmin=0, vmax=1.4)
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("$u_t=U_t+u_t'$")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("u_grid.png",dpi=400)

fig = plt.figure(figsize=(12,5))
h   = plt.pcolormesh(XC, YN, int_vx[:,:,0], edgecolors=("k", 0), cmap="turbo")#, vmin=0, vmax=1.4)
cb  = plt.colorbar(h)
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("$u_t=U_t+u_t'$")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("u.png",dpi=400)

fig = plt.figure(figsize=(12,5))
h   = plt.pcolormesh(XC, YN, int_vy[:,:,0], edgecolors=("k", 0), cmap="turbo")#, vmin=-0.36, vmax=0.2)
cb  = plt.colorbar(h)
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("$v_t=V_t+v_t'$")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("v.png",dpi=400)

fig = plt.figure(figsize=(12,5))
h   = plt.pcolormesh(XC, YN, int_vz[:,:,0], edgecolors=("k", 0), cmap="turbo")#, vmin=-0.31, vmax=0.29)
cb  = plt.colorbar(h)
plt.xlabel("x/c")
plt.ylabel("n/c")
plt.title("$w_t=W_t+w_t'$")
plt.axis('scaled')
plt.tight_layout()
fig.savefig("w.png",dpi=400)