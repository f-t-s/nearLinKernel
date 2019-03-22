#These functions are used to generate the datasets displayed in Figure 2.7
#of the preprint, used to create Table 2.8
using Rotations
using Distributions
function generateEllipse( N::Int64, 
                          jitter::Float64 = 0.1, 
                          ambientDim::Int64 = 20,
                          l1 = 0.5,
                          l2 = 0.5,
                          l3 = 1.0)

  #Uniform sampling from Sphere:
  x = zeros( 3, N )
  x[1, :] .= 1.
  for k = 1 : N 
    x[:,k] = rand( RotMatrix{3} ) * x[:,k]
  end
  
  #Dilation by factor two to create ellipse:
  x[1,:] .= x[1,:] .* l1
  x[2,:] .= x[2,:] .* l2
  x[3,:] .= x[3,:] .* l3

  #Adding random jitter:
  x = x + rand( Uniform(-jitter,jitter), 3, N )

  #Random translation of the Ellipse in threedimensional space
  x = x + repeat( rand( Uniform( -0.50, 0.50), 3 ), 1, N )

  #Random rotation in three-dimensional space:
  x = rand(RotMatrix{3}) * x

  #Random rotation in higher ambient dimension
  #Generate random orthogonal matrix:
  Q = randn(ambientDim,ambientDim)
  Q,~ = qr(Q)
  xHighD = Q * vcat(x, zeros(ambientDim - 3, N) )

  return x, xHighD
end

function generateSpiral( N::Int64,
                         jitter::Float64 = 0.05,
                         ambientDim::Int64 = 20,
                         radius::Float64 = 1.0)
  r = rand( Uniform(-1,1), N )
  rotfact::Int64 = 3
  x = hcat( sin.( rotfact * 2 * pi * r ) * radius,
            cos.( rotfact * 2 * pi * r ) * radius , r )'

  #Adding random jitter:
  x = x + rand( Uniform(-jitter,jitter), 3, N )


  #Random translation of the Spiral in threedimensional space
  x = x + repeat( rand( Uniform( -0.25, 0.25), 3 ), 1, N )


  #Random rotation in three-dimensional space:
  x = rand(RotMatrix{3}) * x

  #Random rotation in higher ambient dimension
  #Generate random orthogonal matrix:
  Q = randn(ambientDim,ambientDim)
  Q,~ = qr(Q)
  xHighD = Q * vcat(x, zeros(ambientDim - 3, N) )

  return x, xHighD
end

function generateSet( N::Int64,
                      jitter::Float64 = 0.05,
                      ambientDim = 20 )
  N1 = round( Int64, N/4 )
  N2 = round( Int64, N/4 )
  N3 = round( Int64, N/4 )
  N4 = N - N1 - N2 - N3 

  x1, x1HighD = generateSpiral( N1, jitter, ambientDim, 1.0)
  x2, x2HighD = generateSpiral( N2, jitter, ambientDim, 0.5)
  x3, x3HighD = generateEllipse( N3, jitter, ambientDim, 1.0, 1.0, 0.5 )
  x4, x4HighD = generateEllipse( N4, jitter, ambientDim, 0.3, 0.4, 1.0 )

  return  hcat( x1HighD, x2HighD, x3HighD, x4HighD ), x1, x2, x3, x4
end

