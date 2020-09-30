using StatsBase,ImageTransformations,Images,ProgressMeter,CSV,JLD2

function summarize(root,subj,vis,tri,sz=(2824÷2,4240÷2))
	if isa(subj,Number)
		subj = subj < 10 ? "0$(subj)_" : "$(subj)_"
	elseif !endswith(subj,'_')
		subj *= "_"
	end

	indir = joinpath(root,"$subj/visit$vis/t$tri")
	outdir = joinpath("/home/driscoll/tmp","$subj/visit$vis/t$tri")
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
	CSV.write(joinpath(outdir,csvname),summary)
	save(joinpath(outdir,"summary.jld2"),"s$(subj[1:2])-v$vis-t$tri",summary)
	return summary
end
