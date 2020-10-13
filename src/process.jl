
ishidden = s -> startswith(basename(s),'.')
isimg = s -> !ishidden(s) && any(endswith.(lowercase(s),[".tif",".tiff",".png",".jpg",".jpeg"]))

"""
	interpimage (img;post=float)
Returns a callable function that performs interpolation/extrapolation of a filtered image. If `post` keyword is given, the provided value is a function applied to the result of the interpolation; default is `float`.
"""
function interpimage(img;post=float)
	X = interpolate(img, BSpline(Linear()))
	Z = extrapolate(X,0)
	return (i,j) -> post(Z(i,j))
end

"Return data relevant for cornea detection to be used on the given image."
function detectiondata(img,m,n)
	# select the angles for optimization of the circle residual
	θ = 2π*(-180:180)/360
	select = falses(size(θ))
	#for range in [ pi*[-2/3,-1/6], pi*[1/6,2/3] ]
	for range in [ pi*[-3/4,-1/4], pi*[1/4,3/4] ]
		select[ range[1] .<= θ .<= range[2] ] .= true
	end
	θ = θ[select]
	
	# optimization methods
	options = ( 
		(method = NewtonTrustRegion(),x_tol = 5e-2,f_tol = 1e-6),
		(method = NelderMead(initial_simplex=Optim.AffineSimplexer(m÷5,0)),x_tol = 5e-2,f_tol = 1e-6) 
	)
	
	# optimization initialization
	u_init = [ (m/size(img,1)).*i for i in initvals(img) ]
	push!(u_init, (m/2.,n/2.,m/2.7))

	# try to screen out the purkinje and eyelid lines
	X = imresize(img,m,n)
	G,B = green.(X),blue.(X)
	Bmax = mapwindow(maximum,B,(15,15))
	maxB = maximum(B)
	for i in 1:m, j in 1:n 
		if (Bmax[i,j] > 0.6*maxB) && G[i,j] > 0.25
			G[i,j] = 0.33
		end
	end
	
	# smooth and interpolate green channel
	ker = KernelFactors.gaussian((m/80,m/80))
	Z = interpimage(imfilter(G,ker))

	return Z,θ,u_init,options
end

"Detect the cornea in an image using default choices."
function detect(img::AbstractMatrix{T} where T <: AbstractRGB,sz=size(img))
	detect(sz,detectiondata(img,sz...)...)
end

"Detect the cornea in an image using provided choices."
function detect(sz,Z,θ,u_init,options)
	u,fmin,best = [],Inf,[]
	for ui in u_init, opt in options
		unew,fnew = fitcircle(Z,sz...,ui...,θ,options=opt)
		@debug "latest: $fnew,$ui,$opt"
		if fnew < fmin 
			u,fmin = unew,fnew
			best = (ui,opt)
		end
	end
	return u,fmin,best
end

function detectfolder(root,subj,vis,tri,sz=[];dosave=true)
	indir = joinpath(root,makedirname(subj,vis,tri))
	result = (fname = String[], cenrow = Float32[], cencol = Float32[], radius = Float32[], fmin = Float32[], init = [], method = [] )
	try
		append!(result.fname,filter(isimg,readdir(indir,join=false)))
	catch
		@warn "Unable to read anything for $subj/$vis/$tri."
	end
	if isempty(sz)
		sz = size(load(joinpath(indir,result.fname[1])))
	end
	@showprogress for fname in copy(result.fname)
		img = load(joinpath(indir,fname))

		# too dark to bother?
		if count(green.(img) .> 0.25 ) < 0.1*prod(size(img))
			@info "Skipping $fname."
			deleteat!(result.fname,length(result.radius)+1)
			continue
		end

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

	if dosave
		outfile = makefilename(subj,vis,tri)
		save("$(outfile).jld2","result",result)
		CSV.write("$(outfile).csv",result)
	end

	return result
end
