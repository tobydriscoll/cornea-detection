module CorneaDetection

using ImageIO,VideoIO
using Images,ImageFiltering,ImageTransformations,Interpolations
using Optim,Statistics,ProgressMeter,StatsBase
using JLD2,CSV,DataFrames

export findpurkinje,detect,detectfolder
export summarize,getsummary,getresults
export makedirname,makefilename,makemovie
import Base: dirname, length, isempty, getindex, iterate, get, fullname
export Subject,Visit,Trial
export dirname,fullname,filenames,shortname,numvisits,numtrials,numframes,summary,results,ImageFolder

include("types.jl")
include("summarize.jl")
include("optimize.jl")
include("process.jl")
include("show.jl")

end # module
