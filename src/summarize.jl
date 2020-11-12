"""
    resize(files,dest,size)
Resize all the images in the image files named in the string vector `files`, using 2-vector `size` for [height,width]. Rewrite them as PNG files in directory `dest`. 

Also computes the same data as `summary` and returns a named tuple of the result.
"""
function resize(folder::AbstractVector{String},dest,sz)
	N = length(folder)
	summary = ( 
		size = fill([0,0],N), 
		avgG = zeros(N), 
		purkrow = fill(NaN,N), 
		purkcol = fill(NaN,N), 
		freqR = fill([NaN],N), 
		freqG = fill([NaN],N),
		freqB = fill([NaN],N)
	)
	@showprogress for (idx,fname) in enumerate(folder)
		X = load(fname)
		X = imresize(X,sz)
		pngname = splitext(basename(fname))[1]*".png"
		save(joinpath(dest,pngname),X)
		update_summary!(summary,idx,X)
	end
	return summary
end

"""
    summarize(files)
Get summary data for all the images in the image files named in the string vector `files`. Returns named tuple with fields:

* `size` [height,width] in pixels
* `purkrow` row index of the purkinje (float), `NaN` if none found
* `purkcol` column index of the purkinje (float), `NaN` if none found
* `avgG` average intensity in the green channel
* `freqR` histogram bin values for intensity in the red channel
* `freqG` histogram bin values for intensity in the green channel
* `freqB` histogram bin values for intensity in the blue channel
"""
function summarize(folder::AbstractVector{String})
	N = length(folder)
	summary = ( 
		size = fill([0,0],N), 
		purkrow = fill(NaN,N), 
		purkcol = fill(NaN,N), 
		avgG = zeros(N), 
		freqR = fill([NaN],N), 
		freqG = fill([NaN],N),
		freqB = fill([NaN],N)
	)
	@showprogress for (idx,fname) in enumerate(folder)
		update_summary!(summary,idx,load(fname))
	end
	return summary
end

# utility for summary
function update_summary!(summary,idx,X)
	sz = size(X)
	summary.size[idx] = [sz...]
	
	purk = findpurkinje(X,0.65)
	if length(purk) > 40
		summary.purkrow[idx] = median(i[1] for i in purk)
		summary.purkcol[idx] = median(i[2] for i in purk)
	end
	edges = collect(LinRange(0,1,33))
	edges[1] = -eps()
	
	summary.avgG[idx] = mean(green.(X))
	
	w = fit(Histogram,vec(red.(X)),edges,closed=:right).weights
	summary.freqR[idx] = w/prod(sz)
	w = fit(Histogram,vec(green.(X)),edges,closed=:right).weights
	summary.freqG[idx] = w/prod(sz)
	w = fit(Histogram,vec(blue.(X)),edges,closed=:right).weights
	summary.freqB[idx] = w/prod(sz)

end

"""
	intensity(data)
Return a vector of an intensity measure in the green channel of each image of a movie. The input should be a DataFrame of the data returned by `summary`.
"""
function intensity(data::DataFrame)
	# Construct an average from the intensity histogram, but leaving out the last bin. 
	# The reason is that most of the pixels in that bin are in the sclera, and the 
	# average value becomes sensitive to small motions of the eyelids. 
	# Might be better to use a median...
	bins = (0.5:31.5)/33
	intensity = Float64[]
	for (i,distr) in enumerate(data.freqG)
		# data might be in string form, if it was loaded from CSV–so parse it
		if !(eltype(distr) <: Number)
			distr = parse.(Float64,split(replace(distr[2:end-1],","=>" ")))
		end
		push!(intensity,sum(bins[j]*distr[j] for j in 1:31))
	end
	return intensity
end

"""
	darkframes(data)
Return a boolean vector indicating which images in a folder are likely too dark to be used. The input should be a DataFrame of the data returned by `summary`.
"""
function darkframes(data::DataFrame)
	q = intensity(data)
	ip,jp = data.purkrow,data.purkcol
	nopurk = isnan.(ip)
	isdark = @. ( q < 0.15 ) | ( nopurk & (q < 0.25) )

	# Some "bright" frames might still need to get marked dark, if they are a good 
	# bit below the "typical bright" frame. Iterate until this stabilizes.
	ct = 3
	while ct < count(isdark)
		ct = count(isdark)
		μ = median(q[.!isdark])
		gap = 0.3*(μ - median(q[isdark]))
		isdark = @. isdark | (q < μ - gap)
	end
	return isdark
end

"""
	goodframes(data)
Return a boolean vector indicating which images in a folder appear to be in the first illuminated interblink period. The input should be a DataFrame of the data returned by `summary`.
"""
function goodframes(data::DataFrame)
	isdark = darkframes(data)
	N = length(isdark)
	# discard where at least 2 frames within each window of three are marked as dark
	w = [ sum(isdark[i:i+2]) for i in 1:N-2 ]
	fd = findall(w.>1)

	# Look for a jump in the indices of dark frames. The period within that jump is our 
	# target.
	idx = findfirst(diff(fd).>2)
	!isnothing(idx) && (isdark[1:fd[idx]] .= true)
	!isnothing(idx) && (isdark[fd[idx+1]:end] .= true)

	return .!isdark
end
