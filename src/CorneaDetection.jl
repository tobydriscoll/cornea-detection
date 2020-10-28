module CorneaDetection

using ImageIO,VideoIO
using Images,ImageFiltering,ImageTransformations,Interpolations
using Optim,Statistics,ProgressMeter,StatsBase
using JLD2,CSV,DataFrames

import Base: dirname, length, isempty, getindex, iterate, get, fullname, summary
export findpurkinje,detect,detectfolder
export summarize,summary,results,makemovie
export Subject,Visit,Trial,ImageFolder,images
export dataroot,dirname,fullname,filenames,shortname
export numvisits,numtrials,numframes,length,visits,trials,summary,results

include("types.jl")
include("summarize.jl")
include("optimize.jl")
include("process.jl")
include("show.jl")

end # module
