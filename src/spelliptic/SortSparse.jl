include("./MutHeap.jl")
using StaticArrays
using LinearAlgebra
#############################################################################
#Introducing the "daycare", which keeps track of the descendants of every node
#The struct is essentially a buffered lower triangular CSC sparse matrix,
#together with an of the degrees of freedom
#############################################################################
mutable struct Daycare{Tv,Ti<:Integer}
  NParents::Ti
  NChildren::Ti
  NBuffer::Ti

  #This array gives for contains the ordering. The i-th parent in the 
  #daycare has id P[i]
  P::Array{Ti,1}
  #This array contains as the i-th element the number that the ith parent 
  #has with respect to the multiresolution ordering.
  revP::Array{Ti,1}
  
  #The array that contains the first "child" for every parent
  colptr::Array{Ti,1}

  #The array that contains the global id-s of the children 
  rowval::Array{Node{Tv,Ti},1}
end

#Function that begins a new parent aka column in daycare
function newParent( dc::Daycare{Tv,Ti}, IdParent::Ti ) where 
          {Ti <: Integer,Tv <: Real}
  dc.NParents += 1 
  dc.P[ dc.NParents ] = IdParent
  dc.colptr[ dc.NParents ] = dc.NChildren + 1
  dc.colptr[ dc.NParents + 1 ] = dc.NChildren + 1
  dc.revP[ IdParent ] = dc.NParents
end

function newChild( dc::Daycare{Tv,Ti}, newChild::Node{Tv,Ti}) where 
          {Ti <: Integer, Tv  <: Real}
  if dc.NChildren >= dc.NBuffer
    dc.NBuffer = 2 * dc.NBuffer 
    resize!( dc.rowval, dc.NBuffer )
  end
  dc.NChildren += 1
  dc.colptr[ dc.NParents + 1 ] += 1
  dc.rowval[ dc.NChildren ] = newChild 
end

function newChildren( dc::Daycare{Tv,Ti}, 
                      newChildren::SubArray{Node{Tv,Ti}}) where 
                        {Ti <: Integer, Tv <: Real}
  while dc.NChildren + size(newChildren,1) >= dc.NBuffer - 1
    dc.NBuffer = 2 * dc.NBuffer 
    resize!( dc.rowval, dc.NBuffer )
  end
  dc.NChildren += size(newChildren,1)
  dc.colptr[dc.NParents + 1] += size(newChildren,1)
  dc.rowval[dc.NChildren - size(newChildren,1) + 1 : dc.NChildren] .= newChildren
end

function _determineChildren!(h::MutHeap{Tv,Ti},
                             dc::Daycare{Tv,Ti},
                             parents::Array{Node{Tv,Ti},1},
                             pivot::Node{Tv,Ti},
                             buffer::Array{Node{Tv,Ti},1},
                             rho::Tv,
                             dist2Func) where {Tv<:Real,Ti<:Integer}
#Function to determine the children of a node in the ordering and sparsity
#pattern.
#TODO Update description
#Inputs:
# h:
#   Heap to keep track of the elements and their distance to the points added
#   already.
# dc: 
#   The "daycare" keeping track of the nodes and their children.
# parents:
#   Array that in its i-th position contains the a node with the id of the 
#   preferred parent of the i-th node and the distance of the i-th node to 
#   its preferred parent.
# Id: 
#   The Id of the point, the children of which are being determined.
# rho:
#   The oversampling parameter determining the size of the sparsity pattern.
# dist:
#   The function dist(i,j) gives the distance between the points with ids 
#   i and j.

  #adding the new parent
  distToParent::Tv = parents[pivot.id].val
  lengthscale = pivot.val
  iterBuffer::Ti = zero(Ti)
  for index = dc.colptr[dc.revP[parents[pivot.id].id]] : (dc.colptr[
                    dc.revP[parents[pivot.id].id] + one(Ti)] - one(Ti))
    #The candidate point for the pivots children
    candidate::Node{Tv,Ti} = dc.rowval[ index ]
    #Distance of the candidate to the pivot:
    dist2::Tv = dist2Func( candidate.id, pivot.id )
    #Check whether the candidate is close enough to be added as a child
    if (dc.revP[candidate.id] == zero(Ti)) && (dist2 <= (lengthscale * rho)^2)
      dist = sqrt(max(dist2,zero(Tv)))
      #Step 1: We add a new child to the pivot:
      #Increase the count of added children by one
      iterBuffer += 1
      #Add node representing the new child to buffer
      buffer[iterBuffer] = Node{Tv,Ti}(dist, candidate.id )
      #Step 2: We possibly update the parent update the distance of the point 
      newDist = update!( h, candidate.id, dist )
      #Step 3: If possible and useful, update the preferred parent:
      if (dist + rho * newDist <= rho * lengthscale) && 
         (dist < parents[candidate.id].val)  
         parents[candidate.id] = Node{Tv,Ti}(dist,pivot.id)
      end
    #If the candidate is far enough away from the pivots parent, such that it can
    #not possibly be within reach, break:
    elseif candidate.val > distToParent + lengthscale * rho
      break
    end
  end
  viewBuffer = view(buffer, 1 : iterBuffer)
  sort!(viewBuffer)
  newParent( dc, pivot.id )
  newChildren(dc, viewBuffer)
end

#Function that constructs the ordering and sparsity pattern from the
#a function evaluating the squared distance
function sortSparse( N::Ti, rho::Tv, dist2Func, initInd = one(Ti) ) where 
                    {Ti<:Integer, Tv<:Real}
  #Constructing the heap and initialising all variables have maximal distance
  h::MutHeap{Tv,Ti} = MutHeap{Tv,Ti}( Array{Node{Tv,Ti},1}(undef, N), 
                                      one(Ti) : N )
  for i = one(Ti) : N
    h.nodes[i] = Node(typemax(Tv), i)
  end
  #Constructing the Daycare The permutation matrices are initialized to the 
  #maximal element to force errors in case an index is being used before it
  #has been initialized.
  dc::Daycare{Tv,Ti} = Daycare{Tv,Ti}(zero(Ti), 
                        zero(Ti),
                        N,
                        zeros(Ti,N),
                        zeros(Ti,N),
                        zeros(Ti, N + one(Ti)),
                        fill(Node{Tv,Ti}(zero(Tv),zero(Ti)), N))

  #Initializing the Buffer used by _determineChildren!
  nodeBuffer::Array{Node{Tv,Ti},1} = Array{Node{Tv,Ti},1}(undef, N)

  #Initializing the array that will hold the distances measured in the end
  distances::Array{Tv,1} = - ones( Tv, N )

  #Performing the first step of the algorithm:
  #Adding the note initInd as the first parent and making all other nodes its
  #children, as well as updating their distance:
  newParent(dc, initInd)
  distances[1] = typemax(Tv)
  for i = one(Ti) : N
    #adds a new Child and updates the corresponding distance of the child
    nodeBuffer[i] = Node{Tv,Ti}(update!(h,i,sqrt(max(dist2Func(i,initInd),zero(Tv)))),i ) 
  end
  viewBuffer = view(nodeBuffer, 1 : N)
  sort!(viewBuffer)
  newChildren(dc, viewBuffer)

  parents::Array{Node{Tv,Ti},1} = Array{Node{Tv,Ti},1}(undef, N)
  for i = one(Ti) : N
    parents[i] = Node{Tv,Ti}(sqrt(max(dist2Func(initInd,i),zero(Tv))),initInd)
  end

  for i = (one(Ti) + one(Ti) ) : N 
    distances[i] = topNode(h).val
    _determineChildren!(h,dc,parents,topNode(h),nodeBuffer,rho,dist2Func)
  end

  dc.rowval = dc.rowval[1:dc.colptr[end - one(Ti)]]
  
  for k = one(Ti) : size( dc.rowval,1 ) 
    dc.rowval[k] = Node{Tv,Ti}( dc.rowval[k].val, dc.revP[dc.rowval[k].id] )
  end

  return dc.colptr, dc.rowval, dc.P, dc.revP, distances
end


#Function constructing the ordering and sparsity pattern for the 
#euclidian distance from a function and $d \times N$ matrix containing 
#N points in d-dimensional space.
function sortSparse(x::Array{Tv,2},
                    rho::Tv,
                    initInd::Ti) where {Tv<:Real,Ti<:Integer}
  N = size(x,2)
  d = size(x,1)
  #Recast as static arrays to use fast methods provided by StaticArrays.jl
  #Possibly remove, doesn't seem to yield a lot.


  function dist2Func( i::Ti, j::Ti )
    out = zero(Tv)
    @fastmath @inbounds @simd for k = 1 : d
      out += ( x[k,i] - x[k,j] )^2
    end
    return out
  end


  colptr, rowval, P, revP, distances = 
    sortSparse( N, rho, dist2Func, initInd )
  return colptr, rowval, P, revP, distances
end

