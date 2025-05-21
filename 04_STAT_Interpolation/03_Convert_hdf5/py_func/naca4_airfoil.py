import numpy as np

def naca4_airfoil(naca_code, num_points=100, chord=1.0):
    """
    Generates coordinates for a NACA 4-digit airfoil.

    Args:
        naca_code (str): The 4-digit NACA code (e.g., '4412').
        num_points (int): Number of points to generate for each surface (upper/lower).
        chord (float): The chord length of the airfoil.

    Returns:
        tuple: (x_coords, y_coords) where x_coords and y_coords are numpy arrays
               representing the airfoil shape, starting from the trailing edge,
               going along the upper surface to the leading edge, and back along
               the lower surface to the trailing edge.
    """
    if len(naca_code) != 4:
        raise ValueError("NACA code must be 4 digits long.")

    # Extract parameters from NACA code
    m_percent = int(naca_code[0])
    p_tenth = int(naca_code[1])
    t_percent = int(naca_code[2:])

    # Convert parameters to fractions of chord
    m = m_percent / 100.0
    p = p_tenth / 10.0
    t = t_percent / 100.0

    # Generate x-coordinates with cosine spacing for better resolution
    # near leading and trailing edges
    # Ensure num_points is at least 2 for linspace
    if num_points < 2:
        num_points = 2
    beta = np.linspace(0, np.pi, num_points)
    # Use cosine spacing: denser points near LE and TE
    x = chord * (0.5 * (1.0 - np.cos(beta)))

    # Calculate thickness distribution (yt) - using the common approximation
    # Note: The last term -0.1015*x^4 is often used, -0.1036*x^4 provides a closed TE
    # Avoid potential division by zero or sqrt of zero at x=0 if chord=0
    x_c = x / chord
    # Handle potential sqrt of negative number for tiny floating point inaccuracies near 0
    x_c_safe = np.maximum(x_c, 0)
    yt = 5 * t * chord * (
        0.2969 * np.sqrt(x_c_safe)
        - 0.1260 * x_c
        - 0.3516 * x_c**2
        + 0.2843 * x_c**3
        - 0.1015 * x_c**4 # Use -0.1036 for a theoretically closed TE
        # - 0.1036 * x_c**4
    )

    # Handle the case where p = 0 (symmetric airfoil)
    if p == 0:
        # For symmetric airfoil, camber is zero
        yc = np.zeros_like(x)
        theta = np.zeros_like(x)
        # Upper surface (x, yt)
        xu = x
        yu = yt
        # Lower surface (x, -yt)
        xl = x
        yl = -yt
    else:
        # Calculate camber line (yc) and its gradient (dyc_dx)
        yc = np.zeros_like(x)
        dyc_dx = np.zeros_like(x)

        # Identify forward and aft sections based on x coordinates
        mask_forward = x <= p * chord
        mask_aft = x > p * chord

        # Calculate normalized x coordinates for each section
        x_f = x_c[mask_forward]
        x_a = x_c[mask_aft]

        # Calculate camber for forward section
        if x_f.size > 0: # Ensure there are points in the forward section
            yc[mask_forward] = m * chord * (x_f / p**2) * (2 * p - x_f)
            # Calculate gradient for forward section - *** FIX HERE ***
            dyc_dx[mask_forward] = (2 * m / p**2) * (p - x_f) # Use x_f directly

        # Calculate camber for aft section
        if x_a.size > 0: # Ensure there are points in the aft section
            yc[mask_aft] = m * chord * ((1 - x_a) / (1 - p)**2) * (1 + x_a - 2 * p)
            # Calculate gradient for aft section - *** FIX HERE ***
            dyc_dx[mask_aft] = (2 * m / (1 - p)**2) * (p - x_a) # Use x_a directly

        # Calculate angle theta (slope of the camber line)
        theta = np.arctan(dyc_dx)

        # Calculate upper and lower surface coordinates using camber and thickness
        xu = x - yt * np.sin(theta)
        yu = yc + yt * np.cos(theta)
        xl = x + yt * np.sin(theta)
        yl = yc - yt * np.cos(theta)

    # Combine coordinates for plotting: TE -> Upper -> LE -> Lower -> TE
    # Reverse lower surface points to go from LE to TE
    # Ensure xu/xl and yu/yl are numpy arrays before flipping/concatenating
    xu = np.asarray(xu)
    yu = np.asarray(yu)
    xl = np.asarray(xl)
    yl = np.asarray(yl)

    x_coords = np.concatenate((np.flip(xu), xl[1:])) # Avoid duplicating LE point
    y_coords = np.concatenate((np.flip(yu), yl[1:]))

    # Ensure closed trailing edge if needed (especially if using -0.1015 term)
    # Force the start and end points to be the same (trailing edge at (chord, 0))
    if x_coords.size > 0: # Check if array is not empty
        x_coords[0] = chord
        y_coords[0] = 0.0
        x_coords[-1] = chord
        y_coords[-1] = 0.0

    return x_coords, y_coords