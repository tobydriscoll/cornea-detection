function largest_rect(X::AbstractArray{Bool})
	m,n = size(X)
	H = count.(X[1,:])
	col,height,areamax = largest_rect(H)
	rowrange = 1:1
	colrange = col .+ (0:height-1)
	for i in 2:m
		# in each column, heights of consecutive ones from row 1 to i
		H = [ X[i,j]*(1+H[j]) for j in 1:n ]
		newcol,newheight,newarea = largest_rect(H)
		if newarea > areamax
			areamax = newarea 
			rowrange = i-newheight+1:i
			colrange = newcol .+ (0:newarea÷newheight-1)
		end
	end
	return rowrange,colrange,areamax

end

function largest_rect(x::AbstractVector{T} where T <: Number)
    maxarea = 0
    x = [-1;x;-1]
    stack = [1]  # index into x

	idx = 0
	height = 0
    for i in 1:length(x)
        while x[i] < x[stack[end]]
			h = x[pop!(stack)]
			area = h*(i-stack[end]-1)
			if area > maxarea
				idx,height = stack[end],h 
				maxarea = area 
			end
		end
		push!(stack,i)
	end
	return idx,height,maxarea
end

"""
	findpurkinje(img,thresh=0.5)
Return indices of pixels believed to be in the purkinje of the image `img`. Decrease `thresh` toward zero to include more pixels, or increase to one to include fewer. Returns a vector of `CartesianIndex`.
"""
function findpurkinje(img::AbstractMatrix{T} where T<:AbsrtactRGB,thresh=0.5)
	B = blue.(img)
	idx = findall(B .> thresh*maximum(B))
	m,n = size(img)
	keep = i -> (5 < i[1] < 0.66m) && (5 < i[2] < 0.66n)
	return filter(keep,idx)
end

function findpurkinje(X::AbstractMatrix{T} where T<:Number,thresh=0.5)
	idx = findall(X .> thresh*maximum(X))
	m,n = size(img)
	keep = i -> (5 < i[1] < 0.66m) && (5 < i[2] < 0.66n)
	return filter(keep,idx)
end

"""
	findpeaks(index,sz,dim)
Look for peaks in the distribution of a vector `index` of indexes. Tuple `sz` is the size of the image that the indexes are taken from, and `dim` is the dimension (1=rows,2=columns) to operate in. Returns a pair of real values at the peaks of the distributions in the two halves of the range of the index.
"""
function findpeaks(idx,sz,dim)
	# look for the peaks in a distribution of indices 
	vals = getindex.(idx,dim)

	edges = LinRange(1,sz[dim]+1,101)
	hist = fit(Histogram,vals,edges)
	w = hist.weights

	## search left/right sides separately 
	k = argmax(w[1:50])
	jleft = mean(edges[k:k+1])
	k = argmax(w[51:end])
	jright = mean(edges[50+k:50+k+1])

	return jleft,jright
end

