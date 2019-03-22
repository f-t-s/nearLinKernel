#This file contains the incomplete factorization routine

#This function computes the inner procuct of two sparse vectors that are 
#stores using the same nzind and nzval, but different offsets iter1:u1 resp.
#iter2:u2
function _innerprod( iter1::Int64,
                     u1::Int64,
                     iter2::Int64,
                     u2::Int64,
                     nzind::Array{Int64,1},
                     nzval::Array{Float64,1} )
  @inbounds begin
    res = zero(Float64)
    while ( iter1 <= u1 ) && ( iter2 <= u2 )
      while (nzind[iter1] == nzind[iter2])  &&( iter1 <= u1 ) && ( iter2 <= u2 )
        res += nzval[iter1] * nzval[iter2]
        iter1 += 1
        iter2 += 1
      end
      if nzind[iter1] < nzind[iter2]
        iter1 += 1
      else
        iter2 += 1
      end
    end
    return res

  end
end


##This function computes the in place incomplete Cholesky factorisation using 
##an inner product wrapped in a function.
#function _el_icholRWHL!( s::Array{Float64},
#                          ind::Array{Int64, 1}, 
#                          jnd::Array{Int64, 1}, 
#                          )
#  @inbounds begin
#    #piv_i goes from 1 to N and designates the number, which is being treated 
#    for piv_i = 1 : ( size( ind, 1 ) - 1 ) 
#      #piv_j ranges over pointers and designates where to find the present column
#      for piv_j = ind[ piv_i ] : ( ind[ piv_i + 1 ] - 1 )
#        #iter designates the pointer, at which the i-iterator starts
#        iter::Int64 = ind[ piv_i ]
#        #jter designates the pointer, at which the j-iterator starts
#        jter::Int64 = ind[ jnd[ piv_j ] ]
#        #The condition makes sure, that only the columns up to th pivot 
#        #are being summed.
#  
#        s[ piv_j ] -= _innerprod( iter, ind[ piv_i + 1 ] - 2,
#                                  jter, ind[ jnd[ piv_j ] + 1 ] - 2,
#                                  jnd,
#                                  s)
#   
#        if jnd[ piv_j ] < piv_i
#          s[ piv_j ] /= ( s[ ind[ jnd[ piv_j ] + 1 ] - 1 ] )
#        elseif jnd[ piv_j ] == piv_i
#          if ( s[ ind[ piv_i + 1 ] - 1 ] ) <= 0.
#            s[ ind[ piv_i + 1 ] - 1 ] = 0
#            #Treating the case where the matrix is near-low rank and our 
#            #factorization returns a low-rank matrix.
#            maxj = piv_i
#            for newpiv_i = ( piv_i + 1 ) : ( size( ind, 1 ) - 1 )               
#              newpiv_j::Int64 = ind[ newpiv_i ] 
#              while jnd[ newpiv_j ] < maxj
#                #iter designates the pointer, at which the i-iterator starts
#                iter = ind[ newpiv_i ]
#                #jter designates the pointer, at which the j-iterator starts
#                jter = ind[ jnd[ newpiv_j ] ]
#                #The condition makes sure, that only the columns up to th pivot 
#                #are being summed.
#          
#                s[ newpiv_j ] -= _innerprod( iter, ind[ newpiv_i + 1 ] - 2,
#                                          jter, ind[ jnd[ newpiv_j ] + 1 ] - 2,
#                                          jnd,
#                                          s)
#
#                s[ newpiv_j ] /= ( s[ ind[ jnd[ newpiv_j ] + 1 ] - 1 ] )
#                newpiv_j += 1
#              end
#              while newpiv_j <= ( ind[ newpiv_i + 1 ] - 1 ) 
#                s[ newpiv_j ] = 0.
#                newpiv_j += 1
#              end
#            end
#            println( "Pivot in i = ", piv_i, " is nonpositive", 
#                     " returning low rank approximation" )
#            return 0
#          else
#            s[ piv_j ] = sqrt( s[ piv_j ] ) 
#          end
#        else
#          println( "Warning: jnd[ piv_j ] > piv_i" )
#        end
#      end
#    end
#  end
#end

#This function computes the in place incomplete Cholesky factorisation using 
#an inner product wrapped in a function.
function _el_icholRWHL!( s::Array{Float64},
                          ind::Array{Int64, 1}, 
                          jnd::Array{Int64, 1}, 
                          )
  @inbounds begin
    #piv_i goes from 1 to N and designates the number, which is being treated 
    for piv_i = 1 : ( size( ind, 1 ) - 1 ) 
      #piv_j ranges over pointers and designates where to find the present column
      for piv_j = ind[ piv_i ] : ( ind[ piv_i + 1 ] - 1 )
        #iter designates the pointer, at which the i-iterator starts
        iter::Int64 = ind[ piv_i ]
        #jter designates the pointer, at which the j-iterator starts
        jter::Int64 = ind[ jnd[ piv_j ] ]
        #The condition makes sure, that only the columns up to th pivot 
        #are being summed.
  
        s[ piv_j ] -= _innerprod( iter, ind[ piv_i + 1 ] - 2,
                                  jter, ind[ jnd[ piv_j ] + 1 ] - 2,
                                  jnd,
                                  s)
   
        if ( s[ ind[ jnd[ piv_j ] + 1 ] - 1 ] ) > 0.0
          if jnd[ piv_j ] < piv_i
            s[ piv_j ] /= ( s[ ind[ jnd[ piv_j ] + 1 ] - 1 ] )
            #debugging:
            if isnan( s[ piv_j ] ) || isinf( s[ piv_j ] )
              println( "ALERT! Not a number" )
            end
            #end debugging

          else 
            #debugging:
            if ind[ piv_i + 1 ] - 1 != piv_j
              println( "ALERT! Not selecting diagonal element" )
            end
            #end debugging
            s[ piv_j ] = sqrt( s[ piv_j ] ) 
            #debugging:
            if isnan( s[ piv_j ] ) || isinf( s[ piv_j ] )
              println( "ALERT! Not a number" )
            end
            #end debugging

          end
        else
          s[ piv_j ] = 0.0
        end
      end
    end
  end
end


#In-place Cholesky factorization of an uppper triangular matrix
function icholU_high_level!( U::SparseMatrixCSC{Float64, Int64} )
  @fastmath _el_icholRWHL!( U.nzval, U.colptr, U.rowval )
end

