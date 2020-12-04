module CorneaDetection

using ImageIO,VideoIO
using Images,ImageFiltering,ImageTransformations,Interpolations
using Optim,Statistics,ProgressMeter,StatsBase
using JLD2,CSV,DataFrames

import Base: dirname, length, isempty, getindex, iterate, get, fullname, summary
export findpurkinje,detect,detectfolder
export resize,summarize,intensity,darkframes,goodframes,makemovie
export Subject,Visit,Trial,subject,visit,trial,ImageFolder,images
export dataroot,dirname,fullname,filenames,shortname
export numvisits,numtrials,numframes,length,visits,trials,summary,results

include("utils.jl")
include("summarize.jl")
include("optimize.jl")
include("detect.jl")
include("show.jl")
include("interface.jl")

end # module
