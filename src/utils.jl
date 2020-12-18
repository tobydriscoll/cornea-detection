# Find the largest possible axes-aligned rectangle of true values within an array.
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

# Find the largest rectangle underneath a histogram-count vector.
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

# Grow a rectangular region to include all the adjacent pixels with similar value. Smaller threshold means stricter check.
function grow_rectangle(X::AbstractMatrix,irange,jrange,thresh=0.1)
	m,n = size(X)
	rect = Array(CartesianIndices((irange,jrange)))
	isempty(rect) && (return rect)
	region = vec(rect)
	ifirst,ilast = irange[[1,end]]
	jfirst,jlast = jrange[[1,end]]
	avg = mean(X[region])  
	for kk = 1:100
		ifirst,ilast = max(1,ifirst-1),min(m,ilast+1)
		jfirst,jlast = max(1,jfirst-1),min(n,jlast+1)
		newpts = setdiff(CartesianIndices((ifirst:ilast,jfirst:jlast)),region)
		add = @. abs(X[newpts]-avg) < thresh
		#@debug "number added = $(count(add))"
		if count(add)==0
			break
		end
		union!(region,newpts[add])
		avg = mean(X[region])
	end
	return region
end

"""
	findpurkinje(img,thresh=0.5)
Return indices of pixels believed to be in the purkinje of the image `img`. Decrease `thresh` toward zero to include more pixels, or increase to one to include fewer. Returns a vector of `CartesianIndex`.
"""
function findpurkinje(img::AbstractMatrix{T} where T<:Colorant,thresh=0.5;options=get_defaults(),rectangle=false)
	X = options.purkinje_channel.(img)
	idx = findpurkinje(X,thresh;options,rectangle=rectangle)
	return idx
end

function findpurkinje(X::AbstractMatrix{T} where T<:Number,thresh=0.5;options=get_defaults(),rectangle=false)
	m,n = size(X)
	θ,area = 0.95,0
	iran = jran = cen = NaN
	B = X.>θ
	while θ >= thresh
		iran,jran,area = largest_rect(B)
		@debug "iran = $iran, jran = $jran"
		# If area is too small, there aren't enough pixels at this brightness.
		if area ≤ m*n*options.purkinje_minarea
			θ -= 0.05
			B = X.>θ
			@debug "θ = $θ"
		else
			h,w = (iran[end]-iran[1]+1,jran[end]-jran[1]+1)
			# Must have appropriate aspect ratio. Bright lines at eyelids can give short, skinny rectangles.
			if h > 0.8w 
				# Must have contrast with neighboring pixels.
				cen = round.(Int,(median(iran),median(jran)))
				# Expand to larger box
				rows = clamp.(cen[1]-h:cen[1]+h,1,m)
				cols = clamp.(cen[2]-w:cen[2]+w,1,n)
				inner = sum(X[iran,jran])
				outer = sum(X[rows,cols]) - inner 
				@debug "inner = $inner, outer = $outer"
				if inner/area - outer/((2h+1)*(2w+1)-area) > 0.25
					# we've got it
					break
				end
			end
			# Remove this region from consideration.
			B[iran,jran] .= false
		end
	end
	return rectangle ? CartesianIndices((iran,jran)) : grow_rectangle(X,iran,jran,.2)
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

