resize(root::String,subj::Number,vis::Number,tri::Number,args...) = resize(root,".",subj,vis,tri,args...)
function resize(rootin::String,rootout::String,subj,vis,tri,sz=(2824÷2,4240÷2))
	T = Trial(rootin,subj,vis,tri)
	outdir = joinpath(rootout,dirname(T))
	summary = (fname = filenames(T), size = [], purkrow = Float64[], purkcol = Float64[], freqR = Vector{Float64}[], freqG = Vector{Float64}[], freqB = Vector{Float64}[]  )
	@showprogress for fname in get(T)
		X = load(fname)
		X = imresize(X,sz)
		push!(summary.size,sz)
		pngname = splitext(basename(fname))[1]*".png"
		save(joinpath(outdir,pngname),X)
	end
end

function summarize(root,subj,vis,tri)
	trial = Trial(root,subj,vis,tri)
	summary = (fname = filenames(trial), size = [], avgG = Float32[], purkrow = Float32[], purkcol = Float32[], freqR = [], freqG = [], freqB = []  )
	@showprogress for fname in get(trial)
		X = load(fname)
		sz = size(X)
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
	CSV.write(joinpath(fullname(trial),"summary.csv"),summary)
	tmp = joinpath(tempname(),"summary.jld2")
	save(tmp,shortname(trial),summary)
	cp(tmp,joinpath(fullname(trial),"summary.jld2"),force=false)
	return summary
end

