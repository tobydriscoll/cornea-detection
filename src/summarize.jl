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

function summarize(folder::AbstractVector{String})
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
		update_summary!(summary,idx,load(fname))
	end
	return summary
end

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


function intensity(T::Trial)
	s = summary(T)
	bins = (0.5:31.5)/33
	intensity = Float64[]
	for (i,distr) in enumerate(s.freqG)
		# data might be in string form if it was loaded from CSV
		if !(eltype(distr) <: Number)
			distr = parse.(Float64,split(replace(distr[2:end-1],","=>" ")))
		end
		push!(intensity,sum(bins[j]*distr[j] for j in 1:31))
	end
	intensity
end

function finddark(T::Trial)
	res = outerjoin(summary(T),results(T),on=:fname,makeunique=true)
	ip,jp = res.purkrow,res.purkcol
	nopurk = isnan.(ip)
	q = intensity(T)
	dark = @. ( q < 0.15 ) | ( nopurk & (q < 0.25) )
	ct = 3
	while ct < count(dark)
		ct = count(dark)
		bright = q[.!dark]
		μ = median(bright)
#		if count(dark) > 3
		σ = 0.15*(μ - median(q[dark]))
#		end
		dark = @. dark | (q < μ - 2σ)
	end

	# N = length(dark)
	# Δ = diff(intensity)
	# # for n in 1:N
	# # 	if n < N
	# # 		dark[n] |= (dark[n+1] & (intensity[n] < μ-σ))			
	# # 	end
	# # 	if n > 1
	# # 		dark[n] |= (dark[n-1] & (intensity[n] < μ-σ))
	# # 	end
	# # end 
	# dark[1] |= dark[2]
	# dark[N] |= dark[N-1]
	# # for n in 2:N-1
	# # 	s = intensity[n+1] - intensity[n-1]
	# # 	dark[n] |= ( dark[n+1] & (s < -2σ) ) || ( dark[n-1] & (s > 2σ) )
	# end
	return dark

end

function goodframes(T::Trial)
	isdark = CorneaDetection.finddark(T);

	N = length(isdark)
	w = [ sum(isdark[i:i+2]) for i in 1:N-2 ]
	#@. isdark[2:N-1] |= (isdark[1:N-2] & isdark[3:N])  # peer pressure
	#fd = findall(isdark)
	fd = findall(w.>1)
	idx = findfirst(diff(fd).>2)
	!isnothing(idx) && (isdark[1:fd[idx]] .= true)
	!isnothing(idx) && (isdark[fd[idx+1]:end] .= true)

	return .!isdark
end


backfill(V::Visit) = foreach(backfill,trials(V))
function backfill(T::Trial)
	s = summary(T)
	resdir = "/Users/driscoll/Dropbox/research/tearfilm/cornea/version4/data"
	r = rightjoin(results(T,resdir),s,on=:fname)
	sort!(r,:fname)
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