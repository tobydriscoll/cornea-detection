# All the types and methods specific to our data storage hierarchy are here. The root of the image 
# data is assumed to have subdirectories in the form "SS_/visitV/tT", where SS is a zero-padded 
# subject number, V is a visit number, and T is a trial number.

struct Subject
	dataroot::String
	number::Integer
end

struct Visit
	subject::Subject
	number::Integer
end
Visit(dataroot,s::Integer,number::Integer) = Visit(Subject(dataroot,s),number)

struct Trial 
	visit::Visit
	number::Integer 
	resultdir::String
end
Trial(v::Visit,number::Integer) = Trial(v,number,".")
Trial(dataroot,s::Integer,v::Integer,number::Integer,resultdir::String=".") = Trial(Visit(Subject(dataroot,s),v),number,resultdir)

#
# Subject methods
#

dataroot(s::Subject) = s.dataroot
number(S::Subject) = S.number
String(s::Subject) = s.number < 10 ? "0$(s.number)" : "$(s.number)"
dirname(s::Subject,full=false) = full ? joinpath(s.dataroot,String(s)*"_") : String(s)*"_"
fullname(s::Subject) = dirname(s,true)
function visits(s::Subject) 
	vis = Visit[]
	screen = t->startswith(t,"visit") 
	contents = filter(screen,readdir(fullname(s)))
	for v in 1:2
		newvis = Visit(s,v)
		if String(newvis) in contents
			push!(vis,newvis)
		end
	end
	return vis
end
numvisits(s::Subject) = length(visits(s))
trials(s::Subject) = reduce(vcat,trials(v) for v in visits(s))
numtrials(s::Subject) = sum( numtrials(v) for v in visits(s) )
length(s::Subject) = sum( numframes(t) for t in trials(s) )

#
# Visit methods
#

subject(v::Visit) = v.subject
dataroot(v::Visit) = dataroot(v.subject)
number(V::Visit) = [number(V.subject),V.number]
String(v::Visit) = "visit$(v.number)"
dirname(v::Visit,full=false) = joinpath(dirname(v.subject,full),String(v))
fullname(v::Visit) = dirname(v,true)
function trials(v::Visit) 
	tri = Trial[]
	screen = t->startswith(t,"t") 
	contents = filter(screen,readdir(fullname(v)))
	for t in 1:10
		newtri = Trial(v,t)
		if String(newtri) in contents 
			push!(tri,newtri)
		end
	end
	return tri
end
numtrials(v::Visit) = length(trials(v))
length(v::Visit) = sum( numframes(t) for t in trials(v) )

#
# Trial methods
#

visit(t::Trial) = t.visit
subject(t::Trial) = subject(visit(t))
dataroot(t::Trial) = dataroot(t.visit)
number(T::Trial) = append!(number(T.visit),T.number)
String(t::Trial) = "t$(t.number)"
dirname(t::Trial,full=false) = joinpath(dirname(t.visit,full),String(t))
fullname(t::Trial) = dirname(t,true)
function shortname(t::Trial) 
	v = visit(t)
	s = subject(v)
	"S$(String(s))_V$(v.number)_T$(t.number)"
end
numframes(t::Trial) = count(isimg,readdir(fullname(t)))
filenames(t::Trial;join=false) = filter(isimg,readdir(fullname(t),join=join))
length(t::Trial) = numframes(t)

#
# Versions of key methods that work with trials
#

function resize(T::Trial,destroot::String,sz=(2824÷2,4240÷2))
	summary = resize(filenames(T,join=true),joinpath(destroot,dirname(T)),sz)
	T2 = Trial(destroot,number(T)...)
	summary = (fname=filenames(T2),summary...)
	CSV.write(joinpath(fullname(T2),"summary.csv"),summary)
	tmp = joinpath(tempname(),"summary.jld2")
	save(tmp,shortname(T2),summary)
	cp(tmp,joinpath(fullname(T2),"summary.jld2"),force=false)
	return T
end

summarize(root,subj,vis,tri) = summarize(Trial(root,subj,vis,tri))
function summarize(T::Trial;dosave=true)
	summary = summarize(filenames(T,join=true))
	summary = (fname=filenames(T),summary...)
	if dosave
		outfile = shortname(T)
		CSV.write(joinpath(dataroot(T),"summaries","$(outfile).csv"),summary)
		tmp = joinpath(tempname(),"$(outfile).jld2")
		save(tmp,shortname(T),summary)
		cp(tmp,joinpath(dataroot(T),"summaries","$(outfile).jld2"),force=false)
	end
	return summary
end

function detect(T::Trial;dosave=true,size=[])
	result = detect(filenames(T,join=true),size)
	result = (fname=filenames(T),result...)
	if dosave
		outfile = shortname(T)
		save("$(outfile).jld2","result",result)
		CSV.write("$(outfile).csv",result)
	end
	return result
end

#
# Accessing results/summary data for a Trial
#

results(root,subject,visit,trial) = results(Trial(root,subject,visit,trial))
results(subject,visit,trial) = results(Trial(".",subject,visit,trial))
results(t::Trial,dir=t.resultdir) = load(joinpath(dir,shortname(t)*".jld2"))["result"] |> DataFrame 

summary(root,subject,visit,trial) = summary(Trial(root,subject,visit,trial))
function summary(T::Trial) 
	sn = shortname(T)
	load(joinpath(dataroot(T),"summaries","$(sn).jld2"))[sn] |> DataFrame
end

function makemovie(T::Trial,sz=(470,706);numframes=Inf)
	folder = filenames(T,join=true)
	imgstack = makemovie(folder,results(T),sz,numframes=numframes)
	encodevideo(shortname(T)*".mp4",imgstack,framerate=10)
	return shortname(T)*".mp4"
end

#
# Utilities
#

intensity(T::Trial) = intensity(summary(T))
darkframes(T::Trial) = darkframes(summary(T))
goodframes(T::Trial) = goodframes(summary(T))

#
# ImageFolder: iterator of images for a given trial
#
struct ImageFolder
	t::Trial
	file::AbstractVector
end

ImageFolder(t::Trial) = ImageFolder(t,filenames(t,join=true))
images(t::Trial) = ImageFolder(t)
ImageFolder(dataroot,s::Integer,v::Integer,t::Integer) = ImageFolder(Trial(dataroot,s,v,t))

iterate(f::ImageFolder) = isempty(f.file) ? nothing : load(f.file[1]),1
iterate(f::ImageFolder,state) = state==length(f.file) ? nothing : (load(f.file[state+1]),state+1)
length(f::ImageFolder) = length(f.file)
isempty(f::ImageFolder) = isempty(f.file)
IteratorEltype(ImageFolder) = EltypeUnknown()
getindex(f::ImageFolder,ind::Integer) = load(f.file[ind])
getindex(f::ImageFolder,inds) = [getindex(f,i) for i in inds]
