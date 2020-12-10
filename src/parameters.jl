# These are parameters that affect the performance of the algorithms. They may need to be changed for different data sets.

# Where to expect the purkinje reflection within the cornea. The value is normalized by the radius of the cornea, zero is at the corneal center, and the value should be in the interval [-1,1]. This only affects one initialization of the detection optimization, and if it is wrong, no harm is done. 
PURKINJE_LOCATION = -1/3

# Which color channel is best for detection of the purkinje image? Use red, green, blue, or Gray.
PURKINJE_CHANNEL = green

# These are the angles that are included in the cornea detection optimization. They are relative to the corneal center, with zero angle at six o'clock. Make a vector of intervals to include. The angles are selected from those spaced by 1 degree around the full circle.
DETECTION_ANGLES = [ π*[-3/4,-1/4], π*[1/4,3/4] ]

# These determine bounds on the allowable center location and cornea radius. They are expressed as fractions of the image height.
BOUNDS_COORD = [0.05,0.95]    # both row and column location
BOUNDS_RADIUS = [0.25,1.4]    # radius limits

# Initial values to try for the optimizer. Returns a vector of (i,j,r) tuples.
INITIALIZATION(m,n) =  [(m/2.,n/2.,0.37m),
						(m/2.,n/2.,0.6m),
						(0.5m,0.65n,0.4m),
						(0.5m,0.35n,0.4m),
]

# The image is blurred with a gaussian having standard deviation equal to m/BLUR_WIDTH, where m is the row size.
BLUR_WIDTH = 80