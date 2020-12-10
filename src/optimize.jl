"""
	initvals(img,thresh=0.5)
Try to find good initial guesses for the cornea position in image `img`, which should be full color. The `thresh` parameter is the fraction of the max in the green channel intensity that a pixel should have to be considered part of the sclera. Returns a vector of zero to two tuples of (center_i,center_j,radius)
""" 
function initvals(X::AbstractMatrix{T} where T <: AbstractRGB,thresh=0.5)
	m,n = size(X)
	purk = findpurkinje(X,rectangle=true)
	
	# vertically, use the purkinje if possible
	i_c = length(purk) > 50 ? median(i[1] for i in purk) : m/2

	# horizontal guess 1: use the sclera's brightness
	G = green.(X) 
	G[purk] .= 0  # ignore the purkinje 
	idx = findall(G.>thresh*maximum(G))
	
	## get distribution peaks
	j_c = r_c = jleft = jright = 0
	try
		jleft,jright = findpeaks(idx,(m,n),2)
		@debug "jleft,jright = $jleft,$jright"
		j_c = (jleft+jright)/2
		r_c = (jright-jleft)/2
	catch 
	end

	## reject if the radius is unrealistic
	@debug "init1 = ($i_c,$j_c,$r_c)"
	guess = BOUNDS_RADIUS[1] < r_c/m < BOUNDS_RADIUS[2] ? [(i_c,j_c,r_c)] : []

	# horizontal guess 2: use the purkinje location heuristically
	if length(purk) > 50
		j_p = median(i[2] for i in purk)
		r_c = (jright - j_p)/(1-PURKINJE_LOCATION)
		j_c = jright - r_c
		@debug "init2 = ($i_c,$j_c,$r_c)"
		push!(guess,(i_c,j_c,r_c))
	end

	return guess
end

"""
	totalgrad(icen,jcen,r,Z[,trange])
Computes the sum of the radial derivative in an image along the specified arc of a circle. `Z` must be a callable function that returns numerical values for real arguments (i.e., allow interpolation). If given, `trange` is a 2-vector defining the range of angles to be used, measured ccw from "straight down" in the usual image visualization (i.e., vertically flipped).
"""
totalgrad(ic,jc,r,Z,θ) = sum(circlegrad(ic,jc,r,Z,θ))
function circlegrad(ic,jc,r,Z,θ)
	g = zeros(size(θ))
	p2c = (r,t) -> (ic+r*cos(t),jc+r*sin(t))
    for (k,t) in enumerate(θ)
		Zplus = Z(p2c(r+2,t)...) 
		Zminus = Z(p2c(r-2,t)...) 
		g[k] = Zplus - Zminus
	end
    return g
end

"""
	fitcircle(imgfun,i0,j0,r0[,t_range][;options=missing])
Find the circle that optimizes the total r-difference criterion, given initial values for the center and radius. Optionally specify the range of angles to be used (measured from straight down in the image) and options for the `optimize` solver.
"""
function fitcircle(imgfun,m,n,i0,j0,r0,θ=[-π,π];options=missing)
	# don't let the center get too close to the edges, or the radius get too small or too close to the nearest edge
	function penalty(u)
		f = (x,lo,hi) -> -log( max(1e-12,(x-lo)*(hi-x)) )
		i,j,r = u
		bc = BOUNDS_COORD
		br = BOUNDS_RADIUS
		return f(i/m,bc...) + f(j/n,bc...) + f(r/m,br...)
	end
	
	if length(θ)==2  # given a range
		# space angles at 1 degree, regardless of the size of the range
		θ = filter( t->θ[1] ≤ t ≤ θ[2], 2π*(-180:180)/360)
	end

	objective = u -> -totalgrad(u...,imgfun,θ) + 0.4penalty(u) 
	if ismissing(options)
		options = (
			method = NewtonTrustRegion(),
			x_tol = 5e-2,
			f_tol = 1e-8
		)
	end
	res = optimize(objective,[i0,j0,r0];options...)
	if !Optim.converged(res)
		@warn "Optimization did not converge"
		#show(res)
	end
	return res.minimizer,res.minimum
end
