"""
	drawcircle!(img,i,j,r[,trange])
Modifies the image to superimpose a circle of specified center and radius. If `trange` is given, it's a 2-vector defining the range of angles to be used, measured ccw from "straight down" in the usual image visualization (i.e., vertically flipped).
"""
function drawcircle!(img,i,j,r,trange::AbstractVector=[-π,π];usecolor=RGB(1,0,0.8))
	t = π/2 .+ LinRange(trange[1],trange[2],5size(img,2))
	valid = i -> checkbounds(Bool,img,i)
	for r in 0.99r:1.01r
		idx = CartesianIndex.(zip(round.(Int,i.+r*sin.(t)),round.(Int,j.+r*cos.(t))))
		idx = filter(valid,idx)
		img[idx] .= usecolor 
	end
	return img
end

"""
	makemovie(files,result,size;numframes=Inf)
Make a movie to show the detected cornea on each image named in the string vector `files`, using the detection results in `result` and the given `size` as [height,width]. If a value is given for `numframes`, frames will be skipped to bring the total frame count to be no greater than that value. Returns a vector of images.
"""
function makemovie(folder::AbstractVector{String},result,sz;numframes=Inf)
	imgstack = []
	circ = sz[1].*[ result.cenrow result.cencol result.radius ]
	N = length(folder)
	select = 1:max(1,ceil(Int,(N-1)/numframes)):N
	@showprogress for k in select
		img = imresize(RGB.(load(folder[k])),sz...)
		drawcircle!(img,circ[k,:]...)
		push!(imgstack,img)
	end
	return imgstack
end
