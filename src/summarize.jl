using StatsBase,ImageTransformations,Images,ProgressMeter,CSV,JLD2

function resize(root,subj,vis,tri,sz=(2824÷2,4240÷2))
	indir = joinpath(root,makedirname(subj,vis,tri))
	outdir = joinpath("/home/driscoll/tmp",makedirname(subj,vis,tri))
	ishidden = s -> startswith(basename(s),'.')
	isimg = s -> !ishidden(s) && any(endswith.(s,[".tif",".tiff",".png"]))
	summary = (fname = [], size = [], purkrow = [], purkcol = [], freqR = [], freqG = [], freqB = []  )

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

function summarize(root,subj,vis,tri)

	indir = joinpath(root,makedirname(subj,vis,tri))
	ishidden = s -> startswith(basename(s),'.')
	isimg = s -> !ishidden(s) && any(endswith.(s,[".tif",".tiff",".png"]))
	summary = (fname = [], size = [], purkrow = [], purkcol = [], freqR = [], freqG = [], freqB = []  )

	try
		append!(summary.fname,filter(isimg,readdir(indir,join=false)))
	catch
		@warn "Unable to read anything for $subj/$vis/$tri."
	end
	@showprogress for fname in summary.fname
		X = load(joinpath(indir,fname))
		sz = size(X)
		push!(summary.size,sz)

		purk = findpurkinje(X,0.33)
		push!(summary.purkrow, median(i[1] for i in purk))
		push!(summary.purkcol, median(i[2] for i in purk))

		edges = collect(LinRange(0,1,33))
		edges[1] = -eps()
		w = fit(Histogram,vec(red.(X)),edges,closed=:right).weights
		push!(summary.freqR,w/prod(sz))
		w = fit(Histogram,vec(green.(X)),edges,closed=:right).weights
		push!(summary.freqG,w/prod(sz))
		w = fit(Histogram,vec(blue.(X)),edges,closed=:right).weights
		push!(summary.freqB,w/prod(sz))
	end
	csvname = "summary.csv"
	CSV.write(joinpath(indir,csvname),summary)
	save(joinpath(indir,"summary.jld2"),makefilename(subj,vis,tri),summary)
	return summary
end
