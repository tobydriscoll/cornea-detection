"""
	findpurkinje(img,thresh=0.5)
Return indices of pixels believed to be in the purkinje of the image `img`. Decrease `thresh` toward zero to include more pixels, or increase to one to include fewer. Returns a vector of `CartesianIndex`.
"""
function findpurkinje(img,thresh=0.5)
	B = blue.(img)
	idx = findall(B .> thresh*maximum(B))
	m,n = size(img)
	keep = i -> (5 < i[1] < 0.66m) && (5 < i[2] < 0.66n)
	return filter(keep,idx)
end

"""
	findpeaks(index,sz,dim)
Look for peaks in the distribution of a vector `index` of indexes. Tuple `sz` is the size of the image that the indexes are taken from, and `dim` is the dimension (1=rows,2=columns) to operate in. Returns a pair of real values at the peaks of the distributions in the two halves of the range of the index.
"""
function findpeaks(idx,sz,dim)
	# look for the peaks in a distribution of indices 
	vals = getindex.(idx,dim)

	edges = LinRange(1,sz[dim]+1,101)
	hist = fit(Histogram,vals,edges)
	w = hist.weights

	## search left/right sides separately 
	k = argmax(w[1:50])
	jleft = mean(edges[k:k+1])
	k = argmax(w[51:end])
	jright = mean(edges[50+k:50+k+1])

	return jleft,jright
end

"""
	initvals(img,thresh=0.5)
Try to find good initial guesses for the cornea position in image `img`, which should be full color. The `thresh` parameter is the fraction of the max in the green channel intensity that a pixel should have to be considered part of the sclera. Returns a vector of zero to two tuples of (center_i,center_j,radius)
""" 
function initvals(X::AbstractMatrix{T} where T <: AbstractRGB,thresh=0.5)
	m,n = size(X)
	purk = findpurkinje(X,0.4)
	
	# vertically, use the purkinje
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
	guess = r_c > m/4 ? [(i_c,j_c,r_c)] : []

	# horizontal guess 2: use the purkinje heuristically
	if length(purk) > 50
		# purkinje is roughly 1/3 of the way across the cornea
		r_c = 0.75(jright - median(i[2] for i in purk))
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
		#dx = minimum([i-1,m-i])
		#dy = minimum([j-1,n-j])
		#return f(i/m,0.05,0.95) + f(j/n,0.05,0.95) + f(r/m,0.25,min(dx,dy)/m)
		return f(i/m,0.05,0.95) + f(j/n,0.05,0.95) + f(r/m,0.25,0.6)
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
	uinit = [i0,j0,r0]
	lower = [0.25m,0.25n,0.25m]
	upper = [0.75m,0.75n,0.6m]
	#res = optimize(objective,lower,upper,uinit,Fminbox(BFGS()))
	res = optimize(objective,uinit;options...)
	#println(res.minimum)
	if !Optim.converged(res)
		#@warn "Did not converge"
		#show(res)
	end
	return res.minimizer,res.minimum
end

# apply detection with some default choices
function fitcircle(img)
	m,n = size(img)
	X = imfilter(green.(img),KernelFactors.gaussian((m/80,m/80)))
	Z = interpimage(X)
	u_init = initvals(X)
	fitcircle(Z,m,n,u_init[1]...,π*[-2/3,2/3])
end
