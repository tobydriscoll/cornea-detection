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
end
Trial(dataroot,s::Integer,v::Integer,number::Integer) = Trial(Visit(Subject(dataroot,s),v),number)

# Subject
dataroot(s::Subject) = s.dataroot
String(s::Subject) = s.number < 10 ? "0$(s.number)" : "$(s.number)"
dirname(s::Subject,full=false) = full ? joinpath(s.dataroot,String(s)*"_") : String(s)*"_"
fullname(s::Subject) = dirname(s,true)
function numvisits(s::Subject) 
	screen = t->startswith(t,"visit") 
	count(screen,readdir(fullname(s)))
end
numtrials(s::Subject) = sum( numtrials(Visit(s,k)) for k in 1:numvisits(s) )

# Visit
subject(v::Visit) = v.subject
dataroot(v::Visit) = dataroot(v.subject)
String(v::Visit) = "visit$(v.number)"
dirname(v::Visit,full=false) = joinpath(dirname(v.subject,full),String(v))
fullname(v::Visit) = dirname(v,true)
function numtrials(v::Visit) 
	screen = t->startswith(t,"t") 
	count(screen,readdir(fullname(v)))
end

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

results(t::Trial) = load(shortname(t)*".jld2")["result"] |> DataFrame 
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
