# These are parameters that affect the performance of the algorithms. They may need to be changed for different data sets.

initialize(m,n) = [(0.5m,0.5n,0.37m),
					(0.5m,0.5n,0.5m),
					(0.35m,0.5n,0.5m),
					(0.65m,0.5n,0.5m),
					(0.5m,0.5n,0.7m),
					(0.5m,0.65n,0.5m),
					(0.5m,0.35n,0.5m),
				]

get_defaults() = Options(-1/3,
blue,
[ π*[-3/4,-1/4], π*[1/4,3/4] ],
80,
1/8000,
[0.05,0.95],[0.05,0.95],[0.25,0.8],
initialize, 
( NewtonTrustRegion(), NelderMead(), ),
)

mutable struct Options
	# Where to expect the purkinje reflection within the cornea. The value is normalized by the radius of the cornea, zero is at the corneal center, and the value should be in the interval [-1,1]. This only affects one initialization of the detection optimization, and if it is wrong, no harm is done. 
	purkinje_location 
	# Which color channel is best for detection of the purkinje image? Use red, green, blue, or Gray.
	purkinje_channel 
	# Minimum allowable area to designate the Purkinje. Expressed as a fraction of the total image's area. If zero, any size is allowed.
	purkinje_minarea
	# These are the angles that are included in the cornea detection optimization. They are relative to the corneal center, with zero angle at six o'clock. Make a vector of intervals to include. The angles are selected from those spaced by 1 degree around the full circle.
	angle_ranges
	# The image is blurred with a gaussian having standard deviation equal to m/BLUR_WIDTH, where m is the row size.
	blur_width 
	# These determine bounds on the allowable center location and cornea radius. 
	bounds_row   # fraction of the row size 
	bounds_column  # fraction of column size 
	bounds_radius  # fraction of the row size 
	# Initial values to try for the optimizer. Returns a vector of (i,j,r) tuples.
	initializer
	# Methods to try in the minimization.
	optimizers
end

function Options(;kwargs...)
	opt = get_defaults()
	for (name,val) in kwargs
		setfield!(opt,Symbol(name),val)
	end
	return opt 
end

function show(io::IO,opt::Options)
	map(fieldnames(Options)) do fn
		println(io,fn,":  ",getfield(opt,fn))
	end
end

function saveopt(options::Options,file="options.txt")
	open(file,"w") do io 
		show(io,options)
	end
	return file
end
