resize(root::String,subj::Number,vis::Number,tri::Number,args...) = resize(root,".",subj,vis,tri,args...)
function resize(rootin::String,rootout::String,subj,vis,tri,sz=(2824÷2,4240÷2))
	T = Trial(rootin,subj,vis,tri)
	outdir = joinpath(rootout,dirname(T))
	summary = (fname = filenames(T), 
	size = [], 
	avgG = Float32[], 
	purkrow = Float32[], 
	purkcol = Float32[], 
	freqR = [], 
	freqG = [], 
	freqB = []  )
	@showprogress for (idx,fname) in enumerate(filenames(T,join=true))
		X = load(fname)
		X = imresize(X,sz)
		push!(summary.size,sz)
		pngname = splitext(basename(fname))[1]*".png"
		save(joinpath(outdir,pngname),X)
		summary.fname[idx] = pngname
		
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
	CSV.write(joinpath(outdir,"summary.csv"),summary)
	tmp = joinpath(tempname(),"summary.jld2")
	save(tmp,shortname(T),summary)
	cp(tmp,joinpath(outdir,"summary.jld2"),force=false)
	return summary
end

summarize(root,subj,vis,tri) = summarize(Trial(root,subj,vis,tri))
function summarize(trial::Trial)
	summary = (fname = filenames(trial), 
	size = [], 
	avgG = Float32[], 
	purkrow = Float32[], 
	purkcol = Float32[], 
	freqR = [], freqG = [], freqB = []  )
	@showprogress for fname in filenames(trial,join=true)
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

function cleanup(T::Trial)

end

backfill(V::Visit) = foreach(backfill,trials(V))
function backfill(T::Trial)
	s = summary(T)
	resdir = "/Users/driscoll/Dropbox/research/tearfilm/cornea/version4/data"
	r = rightjoin(results(T,resdir),s,on=:fname)
	folder = filenames(T,join=true)
	sz = size(load(folder[1]))
	@showprogress for (idx,fname) in enumerate(folder)
		if !ismissing(r.cenrow[idx])
			continue
		end
		img = load(fname)
		Z,θ,u_init,options = detectiondata(img,sz...)
		if idx > 1
			new_ui = sz[1].*(r.cenrow[idx-1],r.cencol[idx-1],r.radius[idx-1])
			push!(u_init,new_ui)
		end
		u,fmin,best = detect(sz,Z,θ,u_init,options)
		r.cenrow[idx] = u[1]/sz[1]
		r.cencol[idx] = u[2]/sz[1]
		r.radius[idx] = u[3]/sz[1]
		r.fmin[idx] = fmin
		r.init[idx] = best[1]
		r.method[idx] = best[2]
	end
	outfile = shortname(T)
	save("$(outfile).jld2","result",r)
	CSV.write("$(outfile).csv",r)
	return r
end