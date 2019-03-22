include("./Spelliptic.jl")
using .Spelliptic
using Distributions
using Random
using LinearAlgebra
using Base.Threads

function testAccuracyStat( K::AbstractSpelMat{Tv,Ti}, 
                           x::Array{Tv,2}, 
                           covFunc, 
                           eta::Float64 = 1.,
                           sigma::Float64 = 0.,
                           NTest::Ti = 1000 ) where{Tv<:AbstractFloat,Ti<:Integer}
  testI = rand( 1 : K.N, NTest )
  testJ = rand( 1 : K.N, NTest )
  errorVec = zeros( NTest )
  distVec = zeros( NTest )
  trueVec = zeros( NTest )
  appVec = zeros( NTest )

  for i = 1 : NTest
    trueVec[ i ] = eta * covFunc( norm( x[ :, testI[ i ] ] -
           x[ :, testJ[ i ] ] ) )
    if ( testI[i] == testJ[i] )
      trueVec[ i ] += sigma
    end
    appVec[ i ] = K[ testI[ i ], testJ[ i ]]
    distVec[ i ] = norm( x[ :, testI[ i ]] - x[ :, testJ[ i ]] )
  end
  errorVec = abs.( trueVec - appVec )
  return trueVec, appVec, errorVec, distVec
end

function testAccuracyStatInt( K::AbstractSpelMat{Tv,Ti}, 
                              x::Array{Tv,2}, 
                              covFunc, 
                              eta::Float64 = 1.,
                              sigma::Float64 = 0.,
                              NTest::Ti = 1000 ) where{Tv<:AbstractFloat,Ti<:Integer}

  distBound = 0.05
  testI = rand( 1 : K.N, NTest )
  testJ = rand( 1 : K.N, NTest )
  errorVec = zeros( NTest )
  distVec = zeros( NTest )
  trueVec = zeros( NTest )
  appVec = zeros( NTest )

  for i = 1 : NTest
    trueVec[ i ] = eta * covFunc( norm( x[ :, testI[ i ] ] -
           x[ :, testJ[ i ] ] ) )
    if ( testI[i] == testJ[i] )
      trueVec[ i ] += sigma
    end
    appVec[ i ] = K[ testI[ i ], testJ[ i ]]
    distVec[ i ] = norm( x[ :, testI[ i ]] - x[ :, testJ[ i ]] )
  end
  errorVec = abs.( trueVec - appVec )

  indsI = maximum( abs.( x[:,testI] .- 0.5 ), dims = 1) .< ( 0.5 - distBound )

  indsJ = maximum( abs.( x[:,testJ] .- 0.5 ), dims = 1) .< ( 0.5 - distBound )

  inds = findall( vec(indsI .& indsJ) )

  return trueVec[inds], appVec[inds], errorVec[inds], distVec[inds]
end


function createScatterVals( K::AbstractSpelMat{Tv,Ti}, 
                            x::Array{Tv,2}, 
                            NTest::Ti = 1000 ) where{Tv<:AbstractFloat,Ti<:Integer}
  NTotalSamples = 10000000
  testI = rand( 1 : K.N, NTotalSamples )
  testJ = rand( 1 : K.N, NTotalSamples )
  distVec = zeros( NTotalSamples )
  indices = zeros( Ti, NTest)
  for i = 1 : NTotalSamples
    distVec[ i ] = norm( x[ :, testI[ i ]] - x[ :, testJ[ i ]] )
  end

  #Sorting by distance
  perm = sortperm( distVec ) 
  distVec = distVec[ perm ]
  testI = testI[ perm ]
  testJ = testJ[ perm ]

  distSamp = rand( NTest )
  for i = 1 : NTest
    indices[i] = searchsortedfirst( distVec, distSamp[i])
  end

  indices = unique( indices )

  distVec = distVec[ indices ]
  testI = testI[ indices ]
  testJ = testJ[ indices ]


  appVec = zeros( size(indices,1) )

  for i = 1 : size(indices,1)
    appVec[ i ] = K[ testI[ i ], testJ[ i ]]
  end

  return distVec, appVec
end


function estErrorFrobeniusStat( K::AbstractSpelMat{Tv,Ti}, 
                                x::Array{Tv,2}, 
                                covFunc, 
                                eta::Tv,
                                sigma::Tv,
                                NSamples = 10,
                                NTest=10000 ) where{Tv<:AbstractFloat,Ti<:Integer}
  estFrob = zeros(Tv, NSamples)
  for k = 1 : NSamples
    trueMat, appMat, errorMat, distMat =  testAccuracyStat( K, x, covFunc, eta, sigma, NTest )
    estFrob[k] = sqrt( K.N * K.N * mean( errorMat.^2 ) )
  end

  return mean( estFrob ), std( estFrob )
end

function estRelErrorFrobeniusStat( K::AbstractSpelMat{Tv,Ti}, 
                                   x::Array{Tv,2}, 
                                   covFunc, 
                                   eta::Tv,
                                   sigma::Tv,
                                   NSamples = 10,
                                   NTest=10000 ) where{Tv<:AbstractFloat,
                                                       Ti<:Integer}
  estFrob = zeros(Tv, NSamples)
  for k = 1 : NSamples
    trueMat, appMat, errorMat, distMat =  testAccuracyStat( K, x, covFunc, eta, sigma, NTest )
    estFrob[k] = sqrt( K.N * K.N * mean( errorMat.^2 ) ) / 
                  sqrt( K.N * K.N * mean( trueMat.^2 ) )
  end

  return mean( estFrob ), std( estFrob )
end


import LinearAlgebra.rank
function rank( K::SpelMat{Tv,Ti} ) where{Tv<:AbstractFloat,
                                         Ti<:Integer}
  return count( diag(K.U) .> zero(Tv) )
end

function estErrorFrobeniusStatInt( K::AbstractSpelMat{Tv,Ti}, 
                                   x::Array{Tv,2}, 
                                   covFunc, 
                                   eta::Tv,
                                   sigma::Tv,
                                   NSamples = 10,
                                   NTest=10000 ) where{Tv<:AbstractFloat,Ti<:Integer}
  estFrob = zeros(Tv, NSamples)
  for k = 1 : NSamples
    trueMat, appMat, errorMat, distMat =  testAccuracyStatInt( K, x, covFunc, eta, sigma, NTest )
    estFrob[k] = sqrt( K.N * K.N * mean( errorMat.^2 ) )
  end

  return mean( estFrob ), std( estFrob )
end

function estRelErrorFrobeniusStatInt( K::AbstractSpelMat{Tv,Ti}, 
                                      x::Array{Tv,2}, 
                                      covFunc, 
                                      eta::Tv,
                                      sigma::Tv,
                                      NSamples = 10,
                                      NTest=10000 ) where{Tv<:AbstractFloat,
                                                       Ti<:Integer}
  estFrob = zeros(Tv, NSamples)
  for k = 1 : NSamples
    trueMat, appMat, errorMat, distMat =  testAccuracyStatInt( K, x, covFunc, eta, sigma, NTest )
    estFrob[k] = sqrt( K.N * K.N * mean( errorMat.^2 ) ) / 
                  sqrt( K.N * K.N * mean( trueMat.^2 ) )
  end

  return mean( estFrob ), std( estFrob )
end


import LinearAlgebra.rank
function rank( K::SpelMat{Tv,Ti} ) where{Tv<:AbstractFloat,
                                         Ti<:Integer}
  return count( diag(K.U) .> zero(Tv) )
end


