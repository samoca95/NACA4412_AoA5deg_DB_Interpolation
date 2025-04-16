import numpy as np
import scipy.io as sio

def load_int_file(path_in: str, variable_name: str, 
                  i_file: int, history: np.ndarray, ncores: int) -> np.ndarray[float]:
    """
    Reads a binary data file and extracts relevant information.
    
    Parameters:
    - path_in (str): Input file path.
    - variable_name (str): Variable name prefix.
    - i_file (int): File index.
    - history (numpy array): Array containing history of data distribution.
    - ncores (int): Number of processing cores.
    
    Returns:
    - bb_buff (numpy array): Extracted data.
    """
    fname = path_in + "/" + variable_name + f"{i_file:05d}"
    with open(fname, "rb") as fid:
        # print(f"  - Reading {fname}...", flush=True)

        # Read header
        hdr     = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        F       = np.fromfile(fid, dtype=np.byte,    count=hdr)
        dum5    = np.fromfile(fid, dtype=np.float64, count=1)[0]
        Rer     = np.fromfile(fid, dtype=np.float64, count=1)[0]
        Domain  = np.fromfile(fid, dtype=np.float64, count=3)
        nel     = np.fromfile(fid, dtype=np.int32,   count=3)
        Poly    = np.fromfile(fid, dtype=np.int32,   count=3)
        time    = np.fromfile(fid, dtype=np.float64, count=1)[0]
        npoints = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        nfields = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        dum6    = np.fromfile(fid, dtype=np.float64, count=1)[0]

        # Read remaining data
        aa = np.fromfile(fid, dtype=np.float64)

        bb_buff = np.zeros(np.sum(history[:, 1]), dtype=np.float64)
        for i in range(ncores):
            start_ibuff = np.sum(history[:i, 1])
            end_ibuff   = np.sum(history[:i+1, 1])
            # print(f"{start_ibuff = }")
            # print(f"{end_ibuff = }")
            start_iaa   = start_ibuff + i
            end_iaa     = end_ibuff + i
            # print(f"{start_iaa = }")
            # print(f"{end_iaa = }")

            # Adjusting for processor gaps within aa
            bb_buff[start_ibuff:end_ibuff] = aa[start_iaa:end_iaa]

            # # Debugging
            # test_read = bb_buff[start_ibuff:end_ibuff]
            # minrankval = np.min(test_read)
            # maxrankval = np.max(test_read)
            # print(f"irank={i}, min={minrankval}, max={maxrankval}")
            # print(test_read[:10])
            # print(test_read[-10:])

    return bb_buff, time
