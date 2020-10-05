function makedirname(subject,visit,trial)
	if isa(subject,Number)
		subject = subject < 10 ? "0$(subject)_" : "$(subject)_"
	elseif !endswith(subject,'_')
		subject *= "_"
	end

	return "$subject/visit$visit/t$trial"
end

function makefilename(subject,visit,trial)
	if isa(subject,Number)
		subject = subject < 10 ? "0$(subject)" : "$(subject)"
	elseif endswith(subject,'_')
		subject = subject[1:end-1]
	end

	return "S$(subject)_V$(visit)_T$trial"
end

function getsummary(root,subject,visit,trial)
	fname = makefilename(subject,visit,trial)
	return load(joinpath(root,makedirname(subject,visit,trial),"summary.jld2"))[fname] |> DataFrame 
end

function getresults(root,subject,visit,trial)
	fname = makefilename(subject,visit,trial)
	return load(joinpath(root,"$(fname).jld2"))["result"] |> DataFrame 
end

"""
	drawcircle!(img,i,j,r[,trange])
Modifies the image to superimpose a circle of specified center and radius. If `trange` is given, it's a 2-vector defining the range of angles to be used, measured ccw from "straight down" in the usual image visualization (i.e., vertically flipped).
"""
function drawcircle!(img,i,j,r,trange::AbstractVector=[-π,π];usecolor=RGB(1,0,0.8))
    #colr = ndims(channelview(img))>2 ? RGB(1,0,0.8) : Gray(1)
	t = π/2 .+ LinRange(trange[1],trange[2],size(img,2))
	valid = i -> checkbounds(Bool,img,i)
	for r in 0.99r:1.01r
		idx = CartesianIndex.(zip(round.(Int,i.+r*sin.(t)),round.(Int,j.+r*cos.(t))))
		idx = filter(valid,idx)
		img[idx] .= usecolor 
	end
	return img
end

function makemovie(datadir::String,outname::String,result,sz=(470,706);numframes=Inf)
	imgstack = []
	circ = sz[1].*[ result.cenrow result.cencol result.radius ]
	N = length(result.fname)
	select = 1:max(1,ceil(Int,(N-1)/numframes)):N
	@showprogress for k in select
		fn = joinpath(datadir,result.fname[k])
		img = imresize(RGB.(load(fn)),sz...)
		drawcircle!(img,circ[k,:]...)
		push!(imgstack,img)
	end
	encodevideo(outname*".mp4",imgstack,framerate=10)
	return outname*".mp4"
end

function makemovie(dataroot::String,subject::Integer,visit::Integer,trial::Integer)
	fn = makefilename(subject,visit,trial)
	dn = joinpath(dataroot,makedirname(subject,visit,trial))
	makemovie(dn,fn,getresults(".",subject,visit,trial))
end

# function makemovie(dataroot,detectroot,subject,visit,trial,sz=(470,706);numframes=Inf)
# 	imgstack = []
# 	result = getresults(detectroot,subject,visit,trial)
# 	circ = sz[1].*[ result.cenrow result.cencol result.radius ]
# 	summary = getsummary(dataroot,subject,visit,trial)
# 	purk = (sz[1]/summary.size[1])*[summary.purkrow summary.purkcol]
# 	N = length(result.fname)
# 	select = 1:max(1,ceil(Int,(N-1)/numframes)):N
# 	@showprogress for k in select
# 		img = imresize(RGB.(load(result.fname[k])),sz...)
# 		drawcircle!(img,circ[k,:]...)

# 		push!(imgstack,img)
# 	end
# 	outname = makefilename(subject,visit,trial)
# 	encodevideo(outname*".mp4",imgstack,framerate=8)
# 	return outname*".mp4"
# end