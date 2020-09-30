function findpurkinje(img,thresh=0.25)
	B = blue.(img)
	idx = findall(B .> thresh*maximum(B))
	m,n = size(img)
	keep = i -> (5 < i[1] < 0.66m) && (5 < i[2] < 0.5n)
	return filter(keep,idx)
end

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
	initvals(img[,thresh])
Find the midpoint of the image portion containing all sufficiently large pixel values, as measured by the average values along horizontal and vertical lines. Default threshhold is 0.1.
""" 
function initvals(X::AbstractMatrix{T} where T <: AbstractRGB,thresh=0.5)
	m,n = size(X)
	
	# horizontally, use the sclera's brightness
	G = green.(X)

	## try to screen out the purkinje 
	purk = findpurkinje(X)
	G[purk] .= 0
	idx = findall(G.>thresh*maximum(G))

	jcen = jr = jleft = jright = 0
	try
		jleft,jright = findpeaks(idx,(m,n),2)
		@debug "jleft,jright = $jleft,$jright"
		jcen = (jleft+jright)/2
		jr = (jright-jleft)/2
	catch 
	end

	if jr < m/4
		jcen = n÷2
		jr = m÷3 
	end

	# vertically, use the purkinje
	icen = length(purk) > 50 ? median(i[1] for i in purk) : m/2

	# modification to use the purkinje horizontally as well
	if length(purk) > 50
		# purkinje is roughly 1/3 of the way across the cornea
		jr = 0.75(jright - median(i[2] for i in purk))
		jcen = jright - jr
	end

	return icen,jcen,jr
end

function initvals(X::AbstractMatrix{T} where T <: Gray{S} where S,thresh=0.75)
	m,n = size(X)
	
	# use the sclera's brightness if possible
	idx = findall(X.>thresh*maximum(X))
	jleft,jright = findpeaks(idx,n÷2,2)

	jcen = (jleft+jright)/2
	jr = (jright-jleft)/2

	try
		keep = idx -> jleft < idx[2] < jright
		itop,ibot = findpeaks(filter(keep,idx),m÷2,1)
		icen = (itop+ibot)/2
	catch
		icen = m/2
	end

	return icen,jcen,jr
end

"""
	totalgrad(icen,jcen,r,Z[,trange])
Computes the sum of the radial derivative in an image along the specified arc of a circle. `Z` must be a callable function that returns numerical values for real arguments (i.e., allow interpolation). If given, `trange` is a 2-vector defining the range of angles to be used, measured ccw from "straight down" in the usual image visualization (i.e., vertically flipped).
"""
function totalgrad(ic,jc,r,Z,θ)
    g = zeros(size(θ))
	p2c = (r,t) -> (ic+r*cos(t),jc+r*sin(t))
    for (k,t) in enumerate(θ)
		Zplus = Z(p2c(r+3,t)...) 
		Zminus = Z(p2c(r-3,t)...) 
		g[k] = Zplus - Zminus
    end
    return sum(g)
end

"""
	fitcircle(imgfun,i0,j0,r0[,t_range][;options=missing])
Find the circle that optimizes the total r-difference criterion, given initial values for the center and radius. Optionally specify the range of angles to be used (measured from straight down in the image) and options for the `optimize` solver.
"""
function fitcircle(imgfun,m,n,i0,j0,r0,θ=[-π,π];options=missing)
	# don't let the center get too close to the edges, or the radius get too small or too close to the nearest edge
	function penalty(u)
		f = (x,lo,hi) -> -log( max(1e-4,(x-lo)*(hi-x)) )
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
	fitcircle(Z,m,n,u_init...,π*[-2/3,2/3])
end
