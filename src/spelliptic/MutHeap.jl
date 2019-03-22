#############################################################################
#Implementation of a mutable maximal Heap
#############################################################################
#Struct representing a node of the heap
struct Node{Tv, Ti} 
  val::Tv
  id::Ti
end

#Mutable Heap (maximal value first)
struct MutHeap{Tv, Ti}
  nodes::Array{Node{Tv,Ti},1}
  lookup::Array{Ti,1}
end

#A function to swap two heapNodes in 
function _swap!(h::MutHeap{Tv,Ti}, 
                a::Ti, 
                b::Ti) where {Tv,Ti}
  #Assining the new values to the lookup table 
  h.lookup[ h.nodes[a].id ] = b
  h.lookup[ h.nodes[b].id ] = a
  tempNode::Node{Tv,Ti} = h.nodes[a]
  h.nodes[a] = h.nodes[b]
  h.nodes[b] = tempNode
end

#Node comparisons
import Base.isless
function isless( a::Node{Tv,Ti}, b::Node{Tv,Ti} ) where {Tv,Ti}
  isless(a.val, b.val)
end

import Base.>=
function >=( a::Node{Tv,Ti}, b::Node{Tv,Ti} ) where {Tv,Ti}
  a.val >= b.val
end

function >=( a::Node{Tv,Ti}, b::Tv ) where {Tv,Ti}
  a.val >= b
end


import Base.>
function >( a::Node{Tv,Ti}, b::Node{Tv,Ti} ) where {Tv,Ti}
  a.val > b.val
end

function >( a::Node{Tv,Ti}, b::Tv ) where {Tv,Ti}
  a.val > b
end


#Function that looks at element h.nodes[hInd] and moves it down the tree 
#if it is sufficiently small. Returns the new index if a move took place, 
#and lastindex(h.nodes), else
function _moveDown!( h::MutHeap{Tv,Ti}, hInd::Ti ) where {Tv,Ti}
  val::Tv = h.nodes[hInd].val
  #If both children exist:
  if 2 * hInd + 1 <= lastindex( h.nodes )
    #If the left child is larger:
    if h.nodes[2 * hInd] >= h.nodes[ 2 * hInd + 1]
      #Check if the child is larger than the parent:
      if h.nodes[2 * hInd] >= val
        _swap!( h, hInd, 2 * hInd )
        return 2 * hInd
      else
        #No swap occuring:
        return lastindex( h.nodes )
      end
    #If the left child is larger:
    else
      #Check if the Child is larger than the paren:
      if h.nodes[2 * hInd + 1] >= val
        _swap!( h, hInd, 2 * hInd + 1 )
        return  2 * hInd + 1
      else
        #No swap occuring:
        return lastindex( h.nodes )
      end
    end
    #If only one child exists:
  elseif 2 * hInd <= lastindex( h.nodes )
    if h.nodes[2 * hInd] > val
      _swap!( h, hInd, 2 * hInd )
      return 2 * hInd 
    end
  end
  #No swap occuring:
  return lastindex( h.nodes )
end

#Get the leading node
function topNode( h::MutHeap )
  return first( h.nodes )
end

#Updates (decreases) an element of the heap and restores the heap property
function update!( h::MutHeap{Tv,Ti}, id::Ti, val::Tv ) where {Tv,Ti}
  tempInd::Ti = h.lookup[ id ]
  if h.nodes[tempInd].val > val
    h.nodes[tempInd] = Node{Tv,Ti}(val,id)
    while ( tempInd < lastindex( h.nodes ) )
      tempInd = _moveDown!( h, tempInd )
    end
    return val
  else
    return h.nodes[id].val
  end
end




