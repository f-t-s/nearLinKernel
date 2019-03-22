#Helper function used to create some of the numerical experiments
using LinearAlgebra

function generateLaplaceBounds2D( q::Ti, dcond::Ti ) where Ti<:Integer 
  n = 2^q + 1
  N = n^2 

  function n_k( k::Ti ) 
    return ( 2 ^ k + 1 )
  end

  function N_k( k::Ti )
    if k == 0
      return 0
    else
      return (2 ^ k + 1) ^ 2 
    end
  end
  
  sizeP = zero(Ti)
  for k = 1 : q
    sizeP += N_k(k)
  end
  sizeP += 4 * n * dcond
  #Adding 4 times the n boundary points to either side
  P = zeros(Ti, sizeP ) 

  function linInd( i::Ti, j::Ti )
    return i + ( j - 1 ) * n
  end

  function ind( linInd::Ti )
    return mod( linInd - 1, n ) + 1 
  end

  function jnd( linInd::Ti )
    return div( linInd - 1, n ) + 1
  end

  iterP::Ti = 1 
  for i = 1 : n
    for j = 1 : dcond 
      P[iterP] = linInd( i, j) 
      iterP += 1
    end
  end

  for i = 1 : n
    for j = n : -1 : ( n -  dcond + 1 )
      P[iterP] = linInd( i, j) 
      iterP += 1
    end
  end

  for j = 1 : n
    for i = 1 : dcond 
      P[iterP] = linInd( i, j) 
      iterP += 1
    end
  end

  for j = 1 : n
    for i = n : -1 : ( n -  dcond + 1 )
      P[iterP] = linInd( i, j) 
      iterP += 1
    end
  end


  NBound = dcond * 4 * n - dcond^2 * 4

  for klev = 1 : q
    for i = 1 : 2^(q-klev) : n
      for j = 1 : 2^(q-klev) : n
        P[iterP] = linInd( i, j )
        iterP += 1
      end
    end
  end

  P = unique( P ) 
  revPRed = sortperm( P[(NBound + 1):end ])
  revP = Array{Ti,1}(undef, N)
  PRed = Array{Ti,1}(undef, N-NBound)
  revP[P] = 1:N
  PRed[revPRed] = 1:(N - NBound)

  h = Float64( 1/(n-1) )
  x::Array{Float64, 2} = zeros( 2, N ) 
  for i = 1 : n
    for j = 1 : n
      x[ 1, linInd(i,j)] = i * h
      x[ 2, linInd(i,j)] = j * h
    end
  end


  Lap = zeros( N, N )
  for piv = 1 : N
    x_i = ind( piv )
    x_j = jnd( piv )
    for off_i = -1 : 1
      for off_j = -1 : 1 
        y_i = min( max( 1, x_i + off_i ), n )
        y_j = min( max( 1, x_j + off_j ), n )

        if x_i == y_i
          if x_j == y_j 
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = 4.0
          else
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = -1.0
          end
        else
          if x_j == y_j 
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = -1.0
          else
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = 0.0
          end
        end
      end
    end
  end
  return x[:,sort(P[(1+NBound ) : end])],
         Lap[sort(P[(1+NBound ) : end ]),sort(P[(1+NBound) : end ]) ],
         PRed,
         revPRed,
         x,
         P,
         revP,
         NBound
end


function generateLaplace2DNoBounds( q::Ti ) where Ti<:Integer 
  n = 2^q + 1
  N = n^2 

  function n_k( k::Ti ) 
    return ( 2 ^ k + 1 )
  end

  function N_k( k::Ti )
    if k == 0
      return 0
    else
      return (2 ^ k + 1) ^ 2 
    end
  end
  
  sizeP = zero(Ti)
  for k = 1 : q
    sizeP += N_k(k)
  end

  P = zeros(Ti, sizeP )

  function linInd( i::Ti, j::Ti )
    return i + ( j - 1 ) * n
  end

  function ind( linInd::Ti )
    return mod( linInd - 1, n ) + 1 
  end

  function jnd( linInd::Ti )
    return div( linInd - 1, n ) + 1
  end

  iterP::Ti = 1 
  for klev = 1 : q
    for i = 1 : 2^(q-klev) : n
      for j = 1 : 2^(q-klev) : n
        P[iterP] = linInd( i, j )
        iterP += 1
      end
    end
  end

  P = unique( P ) 
  revP = Array{Ti,1}(undef, N)
  revP[P] = 1:N

  h = Float64( 1/(n+1) )
  x::Array{Float64, 2} = zeros( 2, N ) 
  for i = 1 : n
    for j = 1 : n
      x[ 1, linInd(i,j)] = i * h
      x[ 2, linInd(i,j)] = j * h
    end
  end

  Lap = zeros( N, N )
  for piv = 1 : N
    x_i = ind( piv )
    x_j = jnd( piv )
    for off_i = -1 : 1
      for off_j = -1 : 1 
        y_i = min( max( 1, x_i + off_i ), n )
        y_j = min( max( 1, x_j + off_j ), n )

        if x_i == y_i
          if x_j == y_j 
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = 4.0
          else 
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = -1.0
          end
        else 
          if x_j == y_j 
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = -1.0
          else 
            Lap[ linInd( x_i, x_j ), linInd( y_i, y_j ) ] = 0.0
          end
        end
      end
    end
  end
  return x , Lap, P
end
