#This file implements a number of popular covariance functions
using SpecialFunctions
function matern( r::Float64, l::Float64, nu::Float64 )
  if r == 0
    return 1.0
  else
    return 2.0^(1.0 - nu) / gamma(nu) * (sqrt( 2.0 * nu ) * r / l )^nu * besselk( nu, sqrt( 2.0 * nu ) * r / l ) 
  end
end

function gaussian( r::Float64, l::Float64 )
  return exp( - ( r/ l )^2 )
end

function exponential( r::Float64, l::Float64 )
  return exp( -r/l )
end

function invMultiquadratic( r::Float64, l::Float64, c::Float64 = 1. )
  return 1/sqrt( c + ( r / l )^2 )
end

function ratQuadratic( r::Float64, l::Float64, c::Float64 = 1. )
  return 1 - (r/l)^2/( (r/l)^2 + c )
end

function cauchy( r::Float64, l::Float64, alpha::Float64, beta::Float64 )
  return ( 1 + (r/l)^alpha )^(-beta/alpha)
end
