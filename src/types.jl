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
Trial(dataroot,s::Integer,v::Integer,number::Integer,resultdir::String=".") = Trial(Visit(Subject(dataroot,s),v),number,resultdir)

# Subject
dataroot(s::Subject) = s.dataroot
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

# Visit
subject(v::Visit) = v.subject
dataroot(v::Visit) = dataroot(v.subject)
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

# Trial
visit(t::Trial) = t.visit
subject(t::Trial) = subject(visit(t))
dataroot(t::Trial) = dataroot(t.visit)
String(t::Trial) = "t$(t.number)"
dirname(t::Trial,full=false) = joinpath(dirname(t.visit,full),String(t))
fullname(t::Trial) = dirname(t,true)
function shortname(t::Trial) 
	v = visit(t)
	s = subject(v)
	"S$(String(s))_V$(v.number)_T$(t.number)"
end
numframes(t::Trial) = count(isimg,readdir(fullname(t)))
filenames(t::Trial,join=false) = filter(isimg,readdir(fullname(t),join=join))
length(t::Trial) = numframes(t)

results(t::Trial,dir=t.resultdir) = load(joinpath(dir,shortname(t)*".jld2"))["result"] |> DataFrame 
summary(t::Trial) = CSV.File(joinpath(dirname(t,true),"summary.csv")) |> DataFrame

# ImageFolder: produce iterator of file names for a given trial
struct ImageFolder
	t::Trial
	file::AbstractVector
end

ImageFolder(t::Trial) = ImageFolder(t,filenames(t,true))
get(t::Trial) = ImageFolder(t)
ImageFolder(dataroot,s::Integer,v::Integer,t::Integer) = ImageFolder(Trial(dataroot,s,v,t))
iterate(f::ImageFolder) = isempty(f.file) ? nothing : f.file[1],1
iterate(f::ImageFolder,state) = state==length(f.file) ? nothing : (f.file[state+1],state+1)
length(f::ImageFolder) = length(f.file)
isempty(f::ImageFolder) = isempty(f.file)
IteratorEltype(ImageFolder) = EltypeUnknown()
getindex(f::ImageFolder,ind::Integer) = f.file[ind]
getindex(f::ImageFolder,inds) = [getindex(f,i) for i in inds]
