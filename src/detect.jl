
ishidden = s -> startswith(basename(s),'.')
isimg = s -> !ishidden(s) && any(endswith.(lowercase(s),[".tif",".tiff",".png",".jpg",".jpeg"]))

"""
	interpimage (img)
Returns a callable function that performs interpolation/extrapolation of a filtered image. 
"""
function interpimage(img)
	X = interpolate(img, BSpline(Linear()))
	Z = extrapolate(X,0)  # extend by zero
	return (i,j) -> float(Z(i,j))
end

"Return data relevant for cornea detection to be used on the given image."
function detectiondata(img,m,n)
	# select the angles for optimization of the circle residual
	θ = 2π*(-180:180)/360
	select = falses(size(θ))
	for rng in DETECTION_ANGLES
		select[ rng[1] .<= θ .<= rng[2] ] .= true
	end
	θ = θ[select]
	
	# optimization methods
	options = ( 
		(method = NewtonTrustRegion(),x_tol = 5e-2,f_tol = 1e-6),
		(method = NelderMead(initial_simplex=Optim.AffineSimplexer(m÷5,0)),x_tol = 5e-2,f_tol = 1e-6) 
	)
	
	# optimization initializations
	u_init = [ (m/size(img,1)).*i for i in initvals(img) ]  # "smart"
	append!(u_init,INITIALIZATION(size(img)...))            # dumb

	X = imresize(img,m,n)
	G,B = green.(X),blue.(X)

	# try to screen out the purkinje and eyelid lines, using bright blue channel values as the indicator
	Bmax = mapwindow(maximum,B,(15,15))  # windowed max
	maxB = maximum(B)
	for i in 1:m, j in 1:n 
		if (Bmax[i,j] > 0.6*maxB) && G[i,j] > 0.25
			G[i,j] = 0.33
		end
	end
	
	# smooth and interpolate green channel
	ker = KernelFactors.gaussian((m/BLUR_WIDTH,m/BLUR_WIDTH))
	Z = interpimage(imfilter(G,ker))

	return Z,θ,u_init,options
end

function detect(sz,Z,θ,u_init,options)
	# Given data and multiple intitializations and options, try them all and find the best.
	u,fmin,best = [],Inf,[]
	for ui in u_init, opt in options
		unew,fnew = fitcircle(Z,sz...,ui...,θ,options=opt)
		@debug "latest: $unew,$fnew,$ui,$opt"
		if fnew < fmin 
			u,fmin = unew,fnew
			best = (ui,opt)
		end
	end
	return u,fmin,best
end

function detect(img::AbstractMatrix{T} where T <: Colorant,sz=size(img))
	detect(sz,detectiondata(img,sz...)...)
end

"""
	detect(files[,size])
	detect(img[,size])
Perform cornea detection on all the image files named in the vector `files`, or on the image `img`. If `size` is given, all images to be resized to match it. 

Returns a named tuple with the following fields:
* `cenrow` row location of the cornea center, normalized by the image height
* `cencol` column location of the cornea center, normalized by the image height
* `radius` radius of the cornea center, normalized by the image height
* `fmin` optimized value of the optimization objective function
* `init` initial condition that led to the optimal result
* `method` info about the optimizer that got the optimal result
"""
function detect(folder::AbstractVector{String},sz=[])
	result = (cenrow = Float32[], cencol = Float32[], radius = Float32[], fmin = Float32[], init = [], method = [] )
	if isempty(sz)
		sz = size(load(folder[1]))
	end
	@showprogress for fname in folder
		img = load(fname)
		Z,θ,u_init,options = detectiondata(img,sz...)
		if !isempty(result.radius)
			new_ui = sz[1].*(result.cenrow[end],result.cencol[end],result.radius[end])
			push!(u_init,new_ui)
		end
		u,fmin,best = detect(sz,Z,θ,u_init,options)

		push!(result.cenrow,u[1]/sz[1])
		push!(result.cencol,u[2]/sz[1])
		push!(result.radius,u[3]/sz[1])
		push!(result.fmin,fmin)
		push!(result.init,best[1])
		push!(result.method,best[2])		
	end
	return result
end

