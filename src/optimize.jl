function isvalidcircle(i,j,r,m,n,options)
	return (options.bounds_row[1] < i/m < options.bounds_row[2]) && 
		(options.bounds_column[1] < j/n < options.bounds_column[2]) && 
		(options.bounds_radius[1] < r/m < options.bounds_radius[2])
end

"""
	initvals(img,thresh=0.5)
Try to find good initial guesses for the cornea position in image `img`, which should be full color. The `thresh` parameter is the fraction of the max in the green channel intensity that a pixel should have to be considered part of the sclera. Returns a vector of zero to two tuples of (center_i,center_j,radius)
""" 
function initvals(X::AbstractMatrix{T} where T <: AbstractRGB,purk=[];options=get_defaults())
	m,n = size(X)
	if isempty(purk)
		purk = findpurkinje(X;options)
	end
	
	# vertically, use the purkinje if possible
	i_c = length(purk) > 50 ? median(i[1] for i in purk) : m/2

	# horizontal guess 1: use the sclera's brightness
	G = green.(X) 
	G[purk] .= 0  # ignore the purkinje 
	
	idx = findall(G.>0.5*maximum(G))
	## get distribution peaks
	j_c = r_c = jleft = jright = wleft = wright = 0
	try
		(jleft,wleft),(jright,wright) = findpeaks(idx,(m,n),2)
		@debug "jleft,jright = $jleft,$jright"
		j_c = (jleft+jright)/2
		r_c = (jright-jleft)/2
	catch 
	end

	## reject if the radius is unrealistic
	@debug "init1 = ($i_c,$j_c,$r_c)"
	guess = isvalidcircle(i_c,j_c,r_c,m,n,options) ? [(i_c,j_c,r_c)] : []

	# horizontal guess 2: use the purkinje location heuristically
	if length(purk) > 50
		j_p = median(i[2] for i in purk)
		γ = options.purkinje_location
		if wleft > wright
			r_c = (jright - j_p)/(1-γ)
			j_c = jright - r_c
		else
			r_c = (j_p-jleft)/(1+γ)
			j_c = jleft + r_c 
		end
		@debug "init2 = ($i_c,$j_c,$r_c)"
		if isvalidcircle(i_c,j_c,r_c,m,n,options)
			push!(guess,(i_c,j_c,r_c))
		end
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
function fitcircle(imgfun,m,n,i0,j0,r0,method,θ=[-π,π];options=get_defaults())
	# don't let the center get too close to the edges, or the radius get too small or too close to the nearest edge
	function penalty(u)
		f = (x,lo,hi) -> log( max(1e-12,4(x-lo)*(hi-x)/(hi-lo)^2) )^4
		i,j,r = u
		b1 = options.bounds_row
		b2 = options.bounds_column
		b3 = options.bounds_radius
		return f(i/m,b1...) + f(j/n,b2...) + f(r/m,b3...)
	end
	
	if length(θ)==2  # given a range
		# space angles at 1 degree, regardless of the size of the range
		θ = filter( t->θ[1] ≤ t ≤ θ[2], 2π*(-180:180)/360)
	end

	objective = u -> -totalgrad(u...,imgfun,θ) + 0.1penalty(u) 
	res = optimize(objective,[i0,j0,r0],method,Optim.Options(x_tol = 1e-2,f_tol = 1e-8))
	if !Optim.converged(res)
		@warn "Optimization did not converge"
		#show(res)
	end
	return res.minimizer,res.minimum
end
