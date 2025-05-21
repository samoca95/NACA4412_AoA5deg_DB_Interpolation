import numpy as np

def mat2dict(input: np.ndarray) -> dict:
    """
    Convert the stat structure to a dictionary where vectors have shape (nx,) and arrays (ny,nx).
    """
    out     = {}
    out_dim = {}
    # Extract variable names
    vars = input.dtype.names
    for var in vars:
        print(f"Extracting {var}...")
        if isinstance(input[var][0], (int, float, complex)): # Vector of scalars
            ndim = 1
            nx = input[var].shape[0]
            mat = np.zeros((nx,))
            for ix in range(nx):  # This is inefficient, but its done to avoid type errors later
                mat[ix] = input[var][ix]
        else: # Array, vector of vectors  
            # Variables of shape nx x nx not implemented
            if np.ndim(input[var][0]) == 2:
                print(f"Warning: {var} not read.")
                continue
            else:
                ndim = 2
                nx = input[var].shape[0]
                ny = len(input[var][0])
                if ny == 0:
                    print(f"Warning: {var} has no data.")
                    continue
                mat = np.zeros((ny, nx))
                for ix in range(nx):
                    for iy in range(ny):
                        mat[iy, ix] = input[var][ix][iy]
        out[var]     = mat
        out_dim[var] = ndim
    return out, out_dim