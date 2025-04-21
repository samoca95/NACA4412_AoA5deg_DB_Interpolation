import numpy as np
import scipy.io as sio

def load_int_file(case_type: str,
                  path_in: str, binary_file: str, 
                  x_pts: np.ndarray[float], y_pts: np.ndarray[float],
                  alpha: np.ndarray[float], 
                  yn:np.ndarray[float], xc: np.ndarray[float]) -> np.ndarray[float]:
    """
    Reads a binary data file and extracts relevant information.
    """

    npoints = len(x_pts)

    # Lenght of wall-normal vector
    ny = len(yn)

    # Number of profiles
    nx = len(xc)

    # Read the file
    fname = path_in + "/" + binary_file
    with open(fname, "rb") as fid:
        # print(f"  - Reading {fname}...", flush=True)

        # Read header
        hdr     = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        F       = np.fromfile(fid, dtype='S1', count=hdr).tobytes().decode('utf-8')
        dum5    = np.fromfile(fid, dtype=np.float64, count=1)[0]
        Rer     = np.fromfile(fid, dtype=np.float64, count=1)[0]
        Domain  = np.fromfile(fid, dtype=np.float64, count=3)
        nel     = np.fromfile(fid, dtype=np.int32,   count=3)
        Poly    = np.fromfile(fid, dtype=np.int32,   count=3)
        nstat   = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        nderiv  = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        times   = np.fromfile(fid, dtype=np.float64, count=1)[0]
        timee   = np.fromfile(fid, dtype=np.float64, count=1)[0]
        atime   = np.fromfile(fid, dtype=np.float64, count=1)[0]
        DT      = np.fromfile(fid, dtype=np.float64, count=1)[0]
        nrec    = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        Tint    = np.fromfile(fid, dtype=np.float64, count=1)[0]
        npoints = np.fromfile(fid, dtype=np.int32,   count=1)[0]
        dum6    = np.fromfile(fid, dtype=np.float64, count=1)[0]

        # Define the structures and coordinates
        Mat = np.zeros((ny, nx))
        fields_F = {}
        fields_D = {}
        # D = F.copy()
        
        # Arrange stat fields in array form
        for ii in range(1, nstat+nderiv+1):
            
            # Skip 8 bytes after the first field
            if ii != 1:
                fid.seek(8, 1)  # '1' means relative to current position
            
            # Read one field of npoints floats
            MatR = np.fromfile(fid, dtype=np.float64, count=npoints)

            # Arrange current field in array form
            for i in range(nx):
                Mat[:, i] = MatR[ny * i : ny * (i + 1)]

            # Generate variable name
            if ii <= nstat:
                fields_F[f'F{ii}']         = Mat.copy()
            else:
                fields_D[f'D{ii - nstat}'] = Mat.copy()
    
    # Simulation parameters
    # Bulk Reynolds number Reb=Ub*c/nu
    Reb=Rer

    # Domain dimensions
    Lx = Domain[0]
    Ly = Domain[1]

    # Fluid density
    rho = 1

    # Kinematic viscosity. Both Ub and c are unit normalizing parameters
    nu = 1 / Reb

    # Molecular visosity
    mu = rho * nu

    # Fields in the binary record from stat files

    # Statistic fields
    # u,v,w,p are instantaneous quantities
    # averaged in time and homogeneous direction z

    # 1.  <u>           % F1
    # 2.  <v>           % F2
    # 3.  <w>           % F3
    # 4.  <p>           % F4

    # 5.  <uu>          % F5
    # 6.  <vv>          % F6
    # 7.  <ww>          % F7
    # 8.  <pp>          % F8

    # 9.  <uv>          % F9
    # 10. <vw>          % F10
    # 11. <uw>          % F11

    # 12. <pu>          % F12
    # 13. <pv>          % F13
    # 14. <pw>          % F14

    # 15. <pdudx>       % F15
    # 16. <pdudy>       % F16
    # 17. <pdudz>       % F17

    # 18. <pdvdx>       % F18
    # 19. <pdvdy>       % F19
    # 20. <pdvdz>       % F20

    # 21. <pdwdx>       % F21
    # 22. <pdwdy>       % F22
    # 23. <pdwdz>       % F23

    # 24. <uuu>         % F24
    # 25. <vvv>         % F25
    # 26. <www>         % F26
    # 27. <ppp>         % F27

    # 28. <uuv>         % F28
    # 29. <uuw>         % F29
    # 30. <vvu>         % F30
    # 31. <vvw>  	    % F31
    # 32. <wwu>         % F32
    # 33. <wwv>         % F33
    # 34. <uvw>         % F34

    # 35. <uuuu>        % F35
    # 36. <vvvv>        % F36
    # 37. <wwww>        % F37
    # 38. <pppp>        % F38

    # 39. <uuuv>        % F39
    # 40. <uuvv>        % F40
    # 41. <uvvv> 	    % F41

    # 42. e11: <((du/dx)^2+(du/dy)^2+(du/dz)^2)> % F42
    # 43. e22: <((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)> % F43
    # 44. e33: <((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)> % F44
    # 45. e12: <(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>  % F45
    # 46. e13: <(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>  % F46
    # 47. e23: <(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>  % F47

    # 48. <dw/dx*dw/dx> % F48
    # 49. <dw/dy*dw/dy> % F49
    # 50. <dw/dx*dw/dy> % F50

    # 51. <du/dx*du/dx> % F51
    # 52. <du/dy*du/dy> % F52
    # 53. <du/dx*du/dy> % F53

    # 54. <dv/dx*dv/dx> % F54
    # 55. <dv/dy*dv/dy> % F55
    # 56. <dv/dx*dv/dy> % F56

    # 57. <du/dx*dv/dx> % F57
    # 68. <du/dy*dv/dy> % F58
    # 59. <du/dx*dv/dy> % F59
    # 60. <du/dy*dv/dx> % F60

    # MA: from the forcing terms:
    # 61. <fx> % F61
    # 62. <fy> % F62
    # 63. <fz> % F63
    # 64. <fx*ux> % F64
    # 65. <fy*uy> % F65
    # 66. <fz*uz> % F66

    #Derivative fields
    # 1. dU/dx           % D1
    # 2. dU/dy           % D2
    # 3. dV/dx           % D3
    # 4. dV/dy           % D4

    # 5. dW/dx           % D5
    # 6. dW/dy           % D6
    # 7. dP/dx           % D7
    # 8. dP/dy           % D8

    # 9.  d<uu>/dx       % D9
    # 10. d<uu>/dy       % D10
    # 11. d<vv>/dx       % D11
    # 12. d<vv>/dy       % D12

    # 13. d<ww>/dx       % D13
    # 14. d<ww>/dy       % D14
    # 15. d<pp>/dx       % D15
    # 16. d<pp>/dy       % D16

    # 17. d<uv>/dx       % D17
    # 18. d<uv>/dy       % D18
    # 19. d<vw>/dx       % D19
    # 20. d<vw>/dy       % D20

    # 21. d<uw>/dx       % D21
    # 22. d<uw>/dy       % D22
    # 23. d<uuu>/dx      % D23
    # 24. d<uuu>/dy      % D24

    # 25. d<vvv>/dx      % D25
    # 26. d<vvv>/dy      % D26
    # 27. d<www>/dx      % D27
    # 28. d<www>/dy      % D28

    # 29. d<ppp>/dx      % D29
    # 30. d<ppp>/dy      % D30
    # 31. d<uuv>/dx      % D31
    # 32. d<uuv>/dy      % D32

    # 33. d<uuw>/dx      % D33
    # 34. d<uuw>/dy      % D34
    # 35. d<vvu>/dx      % D35
    # 36. d<vvu>/dy      % D36

    # 37. d<vvw>/dx      % D37
    # 38. d<vvw>/dy      % D38
    # 39. d<wwu>/dx      % D39
    # 40. d<wwu>/dy      % D40

    # 41. d<wwv>/dx      % D41
    # 42. d<wwv>/dy      % D42
    # 43. d<uvw>/dx      % D43
    # 44. d<uvw>/dy      % D44

    # 45. d2U/dx2        % D45
    # 46. d2U/dy2        % D46
    # 47. d2V/dx2        % D47
    # 48. d2V/dy2        % D48

    # 49. d2W/dx2        % D49
    # 50. d2W/dy2        % D50
    # 51. d2<uu>/dx2     % D51
    # 52. d2<uu>/dy2     % D52

    # 53. d2<vv>/dx2     % D53
    # 54. d2<vv>/dy2     % D54
    # 55. d2<ww>/dx2     % D55
    # 56. d2<ww>/dy2     % D56

    # 57. d2<uv>/dx2     % D57
    # 58. d2<uv>/dy2     % D58
    # 59. d2<uw>/dx2     % D59
    # 60. d2<uw>/dy2     % D60

    # 61. d2<vw>/dx2     % D61
    # 62. d2<vw>/dy2     % D62

    # 63. d<pu>/dx       % D63
    # 64. d<pu>/dy       % D64
    # 65. d<pv>/dx       % D65
    # 66. d<pv>/dy       % D66

    # 67. d<pw>/dx       % D67
    # 68. d<pw>/dy       % D68

    # Dictionary to store the fields
    out = {}
    out['mu'] = mu
    out['nu'] = nu
    out['rho'] = rho

    # Prepare for projection if needed
    sin = np.sin(alpha).reshape(1,-1)
    cos = np.cos(alpha).reshape(1,-1)
    
    # Mean velocities, tensors of Rank 1
    print("U V W")
    if case_type in ['none', 'wake']:
            out['U'] = fields_F['F1'].copy()
            out['V'] = fields_F['F2'].copy()
            out['W'] = fields_F['F3'].copy()
    else:
        out['U'] = np.zeros((ny, nx))
        out['V'] = np.zeros((ny, nx))
        _x = fields_F['F1'].copy()
        _y = fields_F['F2'].copy()
        if case_type == 'suction':
            _proj_x = _x * cos + _y * sin
            _proj_y = _y * cos - _x * sin
        elif case_type == 'pressure':
            _proj_x = _x * cos  + _y * sin
            _proj_y = _x * sin  - _y * cos
        else:
            raise ValueError("case_type must be 'suction', 'pressure', 'wake' or 'none'")
        out['U'] = _proj_x
        out['V'] = _proj_y
        out['W'] = fields_F['F3'].copy()
    
    # Mean pressure, scalar
    print("P")
    out['P'] = fields_F['F4'].copy()

    return out
    
    



    

    