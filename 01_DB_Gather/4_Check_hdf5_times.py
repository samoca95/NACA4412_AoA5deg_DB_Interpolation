import os
import h5py
import numpy as np
import re
import sys

# Configuration
BASE_DIR: str      = "./3_INTERPOLATED_HDF5"   # root directory
EXT: str           = ".h5.uvw"                 # file extension
EXPECTED_DT: float = 0.006                     # expected timestep
EPSILON: float     = 1e-6                      # tolerance for dt check

def extract_time(filepath):
    with h5py.File(filepath, "r") as f:
        return float(f["time"][()])

def sorted_folder_list(base_dir):
    folder_pattern = re.compile(r"^\d{5}_\d{5}$")
    return sorted([
        os.path.join(base_dir, d)
        for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d)) and folder_pattern.match(d)
    ])

def check_folder_times(folder_path):
    files = sorted([
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.endswith(EXT)
    ])
    times = [extract_time(f) for f in files]
    filenames = [os.path.basename(f) for f in files]
    
    # Check intra-folder dt
    for i in range(1, len(times)):
        dt = times[i] - times[i-1]
        if abs(dt - EXPECTED_DT) > EPSILON:
            raise ValueError(
                f"[ERROR] Δt mismatch in {folder_path} between:\n"
                f"        {filenames[i-1]} and {filenames[i]} → Δt = {dt:.10f}"
            )
        else:
            print(f"time:{times[i]:.6f}s OK for {filenames[i]}")
    return times[0], times[-1], files[0], files[-1]

def main():
    folders = sorted_folder_list(BASE_DIR)
    prev_end_time = None
    prev_end_file = None

    for folder in folders:
        print(f"Checking folder: {folder}")
        try:
            t_start, t_end, f_start, f_end = check_folder_times(folder)
        except Exception as e:
            print(e)
            sys.exit(1)

        # Check inter-folder continuity
        if prev_end_time is not None:
            dt_between = t_start - prev_end_time
            if abs(dt_between - EXPECTED_DT) > EPSILON:
                raise ValueError(
                    f"[ERROR] Inter-folder discontinuity:\n"
                    f"        {prev_end_file} → {f_start} → Δt = {dt_between:.10f}"
                )
            else:
                print(f"[OK] Inter-folder continuity for {prev_end_file} → {f_start}")

        prev_end_time = t_end
        prev_end_file = f_end

    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(f"[OK] All intra- and inter-folder Δt ≈ {EXPECTED_DT} validated.")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

if __name__ == "__main__":
    main()

