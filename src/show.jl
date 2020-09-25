"""
	drawcircle!(img,i,j,r[,trange])
Modifies the image to superimpose a circle of specified center and radius. If `trange` is given, it's a 2-vector defining the range of angles to be used, measured ccw from "straight down" in the usual image visualization (i.e., vertically flipped).
"""
function drawcircle!(img,i,j,r,trange::AbstractVector=[-π,π])
    colr = ndims(channelview(img))>2 ? RGB(1,0,0.8) : Gray(1)
	t = π/2 .+ LinRange(trange[1],trange[2],size(img,2))
	for r in 0.99r:1.01r
		idx = CartesianIndex.(zip(round.(Int,i.+r*sin.(t)),round.(Int,j.+r*cos.(t))))
		img[idx] .= colr 
	end
	return img
end
