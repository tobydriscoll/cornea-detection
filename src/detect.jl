
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
function detectiondata(img,m,n;options=get_defaults())
	# select the angles for optimization of the circle residual
	θ = 2π*(-180:180)/360
	select = falses(size(θ))
	for rng in options.angle_ranges
		select[ rng[1] .<= θ .<= rng[2] ] .= true
	end
	θ = θ[select]

	X = imresize(img,m,n)
	purk = findpurkinje(X)
	
	# optimization initializations
	u_init = [ (m/size(X,1)).*i for i in initvals(X,purk) ]  # "smart"
	append!(u_init,options.initializer(m,n))            # dumb

	G = green.(X)

	# try to screen out the eyelid lines, using bright blue values as the indicator.
	C = blue.(X)
	Cmax = mapwindow(maximum,C,(15,15))  # windowed max
	thresh = 0.75*maximum(C)
	G[@. (Cmax > thresh) & (G > 0.25) ] .= 0.33

	# Screen out purkinje, using a generous overestimate
	iran,jran = [ extrema(i[k] for i in purk) for k in 1:2 ]
	h,w = iran[end]-iran[1],jran[end]-jran[1]
	cen = round.(Int,(median(iran),median(jran)))
	rows = clamp.(cen[1]-2h:cen[1]+2h,1,m)
	cols = clamp.(cen[2]-2w:cen[2]+2w,1,n)
	G[rows,cols] .= median(G[rows,cols])

	# smooth and interpolate green channel
	k = m/options.blur_width
	ker = KernelFactors.gaussian((k,k))
	Z = interpimage(imfilter(G,ker))

	return Z,θ,u_init
end

function detect(sz,Z,θ,u_init;options=get_defaults())
	# Given data and multiple intitializations and options, try them all and find the best.
	u,fmin,best = [],Inf,[]
	for ui in u_init, method in options.optimizers
		unew,fnew = fitcircle(Z,sz...,ui...,method,θ;options)
		@debug "latest: $unew,$fnew,$ui,$opt"
		if fnew < fmin 
			u,fmin = unew,fnew
			best = (ui,method)
		end
	end
	return u,fmin,best
end

function detect(img::AbstractMatrix{T} where T <: Colorant,sz=size(img);options=get_defaults())
	detect(sz,detectiondata(img,sz...)...;options)
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
function detect(folder::AbstractVector{String},sz=[];options=get_defaults())
	result = (cenrow = Float32[], cencol = Float32[], radius = Float32[], fmin = Float32[], init = [], method = [] )
	if isempty(sz)
		sz = size(load(folder[1]))
	end
	@showprogress for fname in folder
		img = load(fname)
		Z,θ,u_init = detectiondata(img,sz...;options)
		if !isempty(result.radius)
			new_ui = sz[1].*(result.cenrow[end],result.cencol[end],result.radius[end])
			push!(u_init,new_ui)
		end
		u,fmin,best = detect(sz,Z,θ,u_init;options)

		push!(result.cenrow,u[1]/sz[1])
		push!(result.cencol,u[2]/sz[1])
		push!(result.radius,u[3]/sz[1])
		push!(result.fmin,fmin)
		push!(result.init,best[1])
		push!(result.method,best[2])		
	end
	return result
end

