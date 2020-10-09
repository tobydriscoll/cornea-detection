resize(root::String,subj::Number,vis::Number,tri::Number,args...) = resize(root,".",subj,vis,tri,args...)
function resize(rootin::String,rootout::String,subj,vis,tri,sz=(2824÷2,4240÷2))
	indir = joinpath(rootin,makedirname(subj,vis,tri))
	outdir = joinpath(rootout,makedirname(subj,vis,tri))
	ishidden = s -> startswith(basename(s),'.')
	isimg = s -> !ishidden(s) && any(endswith.(s,[".tif",".tiff",".png"]))
	summary = (fname = String[], size = [], purkrow = Float64[], purkcol = Float64[], freqR = Vector{Float64}[], freqG = Vector{Float64}[], freqB = Vector{Float64}[]  )

	try
		append!(summary.fname,filter(isimg,readdir(indir,join=false)))
	catch
		@warn "Unable to read anything for $subj/$vis/$tri."
	end
	@showprogress for fname in summary.fname
		X = load(joinpath(indir,fname))
		X = imresize(X,sz)
		push!(summary.size,sz)
		pngname = splitext(fname)[1]*".png"
		save(joinpath(outdir,pngname),X)
	end
end

function summarize(root,subj,vis,tri;resize=false)
	indir = joinpath(root,makedirname(subj,vis,tri))
	ishidden = s -> startswith(basename(s),'.')
	isimg = s -> !ishidden(s) && any(endswith.(s,[".tif",".tiff",".png"]))
	summary = (fname = String[], size = [], avgG = Float32[], purkrow = Float32[], purkcol = Float32[], freqR = [], freqG = [], freqB = []  )

	try
		append!(summary.fname,filter(isimg,readdir(indir,join=false)))
	catch
		@warn "Unable to read anything for $subj/$vis/$tri."
	end
	@showprogress for fname in summary.fname
		X = load(joinpath(indir,fname))
		sz = size(X)

		if resize 
			sz = div.(sz,2)
			X = imresize(X,sz)
			pngname = splitext(fname)[1]*".png"
			save(joinpath("/home/driscoll/tmp",makedirname(subj,vis,tri),pngname),X)
		end	

		push!(summary.size,[sz...])

		purk = findpurkinje(X,0.5)
		if length(purk) > 40
			push!(summary.purkrow, median(i[1] for i in purk))
			push!(summary.purkcol, median(i[2] for i in purk))
		else
			push!(summary.purkrow, NaN)
			push!(summary.purkcol, NaN)
		end
		edges = collect(LinRange(0,1,33))
		edges[1] = -eps()

		push!(summary.avgG,mean(green.(X)))

		w = fit(Histogram,vec(red.(X)),edges,closed=:right).weights
		push!(summary.freqR,w/prod(sz))
		w = fit(Histogram,vec(green.(X)),edges,closed=:right).weights
		push!(summary.freqG,w/prod(sz))
		w = fit(Histogram,vec(blue.(X)),edges,closed=:right).weights
		push!(summary.freqB,w/prod(sz))
	end
	csvname = "summary.csv"
	CSV.write(joinpath(indir,csvname),summary)
	tmp = joinpath(tempname(),"summary.jld2")
	save(tmp,makefilename(subj,vis,tri),summary)
	cp(tmp,joinpath(indir,"summary.jld2"),force=false)
	return summary
end

