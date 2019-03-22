#This module defines the Spelliptic module and the corresponding operator 
#overloading. 
module Spelliptic
using IterativeSolvers
using SparseArrays
using Distributions
using Random
using LinearAlgebra
using Base.Threads

include("SortSparse.jl")
include("Factorization.jl")

export SpelMat
export AbstractSpelMat
export setCovStat!
export SpelMatIter
export computeConditionNumber
export computeMaxEigenvalue
#export testAccuracy

#An abstract operator represented by the sparse Cholesky factorization.
abstract type AbstractSpelMat{Tv<:AbstractFloat,Ti<:Integer} end

struct SpelMat{Tv<:AbstractFloat, Ti<:Integer} <: AbstractSpelMat{Tv,Ti}
  #Number of degrees of freedom
  N::Ti 
  #Number of nonzeros
  nnz::Ti
  #An integer array containing the IDs of the dofs in the order of elimination
  P::Array{Ti,1}
  #An integer array that contains as ith entry the position that the dof with 
  #Id i has in the elimination ordering
  revP::Array{Ti,1}
  #The upper-triangular matrix containing the sparse factor 
  U::SparseMatrixCSC{Tv, Ti}
end


#Constructors
#Constructor to generate a new SpelMat from a d \times N array of N points in 
#d-dimensional space and an accuracy factor \rho.
function SpelMat{Tv,Ti}( x::Array{Tv,2}, rho::Tv ) where{Tv<:AbstractFloat, 
                                                         Ti<:Integer}
  N::Ti = size(x,2)
  ind, jnd, P, revP = sortSparse(x, rho, one(Ti))
  jnd = getfield.(jnd,:id)
  #sorting the jnd vectors to comply with expectation of Julia's sparse matrix
  #format
  for k = 1 : N
    sort!(view( jnd, ind[k] : ( ind[k+1] - 1 ) ) )
  end
  nnz = size( jnd, 1 )
  U = SparseMatrixCSC{Tv,Ti}( N,
                               N,
                               ind,
                               jnd,
                               ones(Tv,nnz) )
  U = SparseMatrixCSC( U' )
  GC.gc()
  return SpelMat( N,
                  nnz,
                  P,
                  revP,
                  U )
end

import LinearAlgebra.Matrix
function Matrix( K::SpelMat ) 
  return ( Matrix(K.U)' * Matrix(K.U) )[K.revP, K.revP]
end


#Operator and function overloading:
#compute the log determinant
import LinearAlgebra.logdet
function logdet( K::SpelMat{Tv,Ti} ) where{Tv<:AbstractFloat,
                                           Ti<:Integer}
  return 2 * sum( log.( diag( K.U )))
end

#Operator overloading:
import Base.*
function *(K::SpelMat{Tv,Ti}, v::AbstractVector ) where{Tv<:AbstractFloat,
                                                            Ti<:Integer}
  return ( K.U' * ( K.U * v[K.P] ) )[K.revP]
end

function *(K::SpelMat{Tv,Ti}, v::AbstractMatrix ) where{Tv<:AbstractFloat,
                                                            Ti<:Integer}
  return ( K.U' * ( K.U * v[K.P,:] ) )[K.revP,:]
end

function *( v::AbstractMatrix, K::SpelMat{Tv,Ti} ) where{Tv<:AbstractFloat,
                                                                Ti<:Integer}
  return ( ( v[:,K.P] * K.U' ) * K.U )[:,K.revP]
end

import Base.\

function \( K::SpelMat{Tv,Ti}, v::AbstractVector{Tv} ) where{Tv<:AbstractFloat,
                                                             Ti<:Integer}
  return ( K.U \ ( K.U' \ v[K.P] ) )[K.revP]
end

function \( K::SpelMat{Tv,Ti}, v::AbstractMatrix{Tv} ) where{Tv<:AbstractFloat,
                                                             Ti<:Integer}
  return ( K.U \ ( K.U' \ v[K.P,:] ) )[K.revP,:]
end



import Base./
function /( v::AbstractMatrix{Tv}, K::SpelMat{Tv,Ti} ) where{Tv<:AbstractFloat,
                                                             Ti<:Integer}
  return ( ( v[:,K.P] / K.U ) / K.U' )[:,K.revP]
end

#function /( v::AbstractVector{Tv}, K::SpelMat{Tv,Ti} ) where{Tv<:AbstractFloat,
#                                                             Ti<:Integer}
#  return (( K.U \ ( v[K.P] / K.U )' )[K.revP ] )'
#end


import Base.getindex
function getindex( K::SpelMat, i::Int64, j::Int64 )
  return dot( K.U[ :, K.revP[i]], K.U[ :, K.revP[j]])
end

function getindex( K::SpelMat, ivec, jvec )
  return K.U[ :, K.revP[ivec]]' * K.U[ :, K.revP[jvec]]
end


#Function to set the content of a SpelMat struct to the those given by a 
#stationary covariance function. 
#The variable sigma allows to add an optional nugget term to the covariance.
#Note that for larger nugget terms the approximation property will detoriate,
#Use SpelMatIter when working with large nugget terms.
function setCovStat!( K::SpelMat{Tv,Ti},
                      x::Array{Tv,2}, 
                      covFunc, 
                      eta::Tv = 1.,
                      sigma::Tv = 0. ) where{Tv<:AbstractFloat,Ti<:Integer}
  for i = 1 : K.N
    for j = K.U.colptr[ i ] : ( K.U.colptr[ i + 1 ] - 1 )
      K.U.nzval[ j ] = norm( view( x, :, K.P[ i ] ) - 
                             view( x, :, K.P[ K.U.rowval[ j ] ] ) )
    end
  end
  K.U.nzval .= map( covFunc, K.U.nzval )
  K.U.nzval .*= eta
  if sigma != 0.
    for i = 1 : K.N
      for j = K.U.colptr[ i ] : ( K.U.colptr[ i + 1 ] - 1 )
        if K.U.rowval[ j ] == i
          K.U.nzval[ j ] += sigma 
        end
      end
    end
  end
  #returns the time taken by the Cholesky factorization
  (@timed icholU_high_level!( K.U ))[2]
end



###############################################################################
#An extension of SpelMat that treats the "nugget" by using iterative methods.
#The compression is only applied to the part that comes from a Green's function,
#then, conjugate gradient is used to invert the sum of this compressed matrix
#and its nugget.
mutable struct SpelMatIter{Tv<:AbstractFloat, Ti<:Integer} <: AbstractSpelMat{Tv,Ti}
  #Number of degrees of freedom
  N::Ti 
  #Number of nonzeros
  nnz::Ti
  #An integer array containing the IDs of the dofs in the order of elimination
  P::Array{Ti,1}
  #An integer array that contains as ith entry the position that the dof with 
  #Id i has in the elimination ordering
  revP::Array{Ti,1}
  #The upper-triangular matrix containing the sparse factor 
  U::SparseMatrixCSC{Tv,Ti}
  #The upper-triangular matrix containing the sparse approximatte factor that is used 
  #as a preconditioner for the inversion
  UPrec::SparseMatrixCSC{Tv,Ti}
  sigma::Tv
  tol::Tv
end

#Constructors
#Constructor to generate a new SpelMatIter from a d \times N array of N points in 
#d-dimensional space and an accuracy factor \rho.
function SpelMatIter{Tv,Ti}( x::Array{Tv,2}, rho::Tv ) where{Tv<:AbstractFloat, 
                                                         Ti<:Integer}
  N::Ti = size(x,2)
  ind, jnd, P, revP = sortSparse(x, rho, one(Ti))
  jnd = getfield.(jnd,:id)
  #sorting the jnd vectors to comply with expectation of Julia's sparse matrix
  #format
  for k = 1 : N
    sort!(view( jnd, ind[k] : ( ind[k+1] - 1 ) ) )
  end
  nnz = size( jnd, 1 )
  U = SparseMatrixCSC{Tv, Ti}( N,
                               N,
                               ind,
                               jnd,
                               ones(Tv,nnz) )
  U = SparseMatrixCSC( U' )
  GC.gc()
  UPrec = SparseMatrixCSC{Tv,Ti}( N,
                                  N,
                                  U.colptr,
                                  U.rowval,
                                  ones(Tv,nnz) )
                                  
  return SpelMatIter{Tv,Ti}( N,
                             nnz,
                             P,
                             revP,
                             U,
                             UPrec,
                             zero(Tv),
                             1e-7 )
end

import LinearAlgebra.Matrix
function Matrix( K::SpelMatIter )
  return ( Matrix(K.U)' * Matrix(K.U) + K.sigma * I )[K.revP, K.revP]
end




function setCovStat!( K::SpelMatIter{Tv,Ti},
                      x::Array{Tv,2}, 
                      covFunc, 
                      eta::Tv, 
                      sigma::Tv,
                      sigmaimp::Tv = zero(Tv) ) where{Tv<:AbstractFloat,Ti<:Integer}
  for i = 1 : K.N
    for j = K.U.colptr[ i ] : ( K.U.colptr[ i + 1 ] - 1 )
      K.U.nzval[ j ] = norm( view( x, :, K.P[ i ] ) - 
                             view( x, :, K.P[ K.U.rowval[ j ] ] ) )
    end
  end
  K.U.nzval .= map( covFunc, K.U.nzval )
  K.U.nzval .*= eta
  if sigmaimp != zero(Tv)
    for i = 1 : K.N
      for j = K.U.colptr[ i ] : ( K.U.colptr[ i + 1 ] - 1 )
        if K.U.rowval[ j ] == i
          K.U.nzval[ j ] += sigmaimp
        end
      end
    end
  end

  K.sigma = sigma
  K.UPrec.nzval .= K.U.nzval
  K.UPrec .+= sigma * sparse( I, K.N, K.N ) 
  @time icholU_high_level!( K.U )
  @time icholU_high_level!( K.UPrec )
end


#Operator overloading:
import Base.*
function *(K::SpelMatIter{Tv,Ti}, v::AbstractVector ) where{Tv<:AbstractFloat,
                                                                Ti<:Integer}
  return ( K.U' * ( K.U * v[K.P] ) )[K.revP] + K.sigma * v
end

function *(K::SpelMatIter{Tv,Ti}, v::AbstractMatrix ) where{Tv<:AbstractFloat,
                                                                Ti<:Integer}
  return ( K.U' * ( K.U * v[K.P,:] ) )[K.revP,:] + K.sigma * v
end


import Base.*
function *( v::AbstractMatrix, K::SpelMatIter{Tv,Ti} ) where{Tv<:AbstractFloat,
                                                                 Ti<:Integer}
  return ( ( v[:,K.P] * K.U' ) * K.U )[:,K.revP]  + v * K.sigma
end


#Wrapper for the preconditioned K::SpelMatIter, which is obtained as
# UPrec^{-\top} U^{\top} U UPrec^{-1}
struct _SpelMatIterPrec{Tv<:AbstractFloat,Ti<:Integer}
  K::SpelMatIter{Tv,Ti}
end

#Implementing left matrix multiply for the preconditioned part
import Base.*
function *( P::_SpelMatIterPrec{Tv,Ti}, v::AbstractVector ) where{Tv<:AbstractFloat,
                                                                      Ti<:Integer}
  v = P.K.UPrec \ v
  v = P.K.UPrec' \ ( P.K.U' * ( P.K.U * v ) + P.K.sigma * v )
  return v
end

#Implementing right matrix multiply for the preconditioned part
import Base.*
function *( v, P::_SpelMatIterPrec{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  v = v / P.K.UPrec'
  v = ( v * P.K.sigma + ( v * P.K.U' ) * P.K.U ) / P.K.UPrec
  return v
end

#Implement left inverse of the preconditioned matrix, used only for diagnostic purposes
#(computing condition numbers)
import Base.\
function \( P::_SpelMatIterPrec{Tv,Ti}, v ) where{Tv<:AbstractFloat,Ti<:Integer}
    cg( P,
        v; 
        verbose = false,
        tol = 1e-5 )
end


import Base.\
function \( K::SpelMatIter{Tv,Ti}, v::AbstractVector{Tv} ) where{Tv<:AbstractFloat,
                                                                 Ti<:Integer}
  return ( K.UPrec \ cg( _SpelMatIterPrec{Tv,Ti}(K) ,
                         K.UPrec' \ v[K.P]  ; 
                         verbose = true,
                         tol = K.tol  ) )[K.revP]
end

function \( K::SpelMatIter{Tv,Ti}, v::AbstractMatrix{Tv} ) where{Tv<:AbstractFloat,
                                                                 Ti<:Integer}
  return ( K.UPrec \ cg( _SpelMatIterPrec{Tv,Ti}(K) ,
                        K.UPrec' \ v[K.P,:]  ; 
                         verbose = true,
                         tol = K.tol ) )[K.revP,:]
end


import Base./
function /( v, K::SpelMatIter{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  return ( cg( _SpelMatIterPrec{Tv,Ti}(K),
               v[K.P] / K.UPrec;
               verbose = true,
               tol = K.tol ) / K.UPrec' )[K.revP] 
end

import LinearAlgebra.mul!
function mul!( y , P::_SpelMatIterPrec, v::AbstractVector )
y .= P * v
end


import Base.getindex
function getindex( K::SpelMatIter{Tv,Ti}, i, j ) where{Tv<:AbstractFloat,Ti<:Integer}
  return dot( K.U[ :, K.K.revP[i]], K.U[ :, K.revP[j]]) + 
  K.sigma * convert( Tv, i==j )
end

function getindex( K::SpelMatIter{Tv,Ti}, ivec, jvec ) where{Tv<:AbstractFloat,
                                                             Ti<:Integer}
  return K.U[ :, K.revP[ivec]]' * K.U[ :, K.revP[jvec]] +
  K.sigma * speye( K.N )[ivec, jvec]
end

import Base.eltype
function eltype( K::SpelMatIter{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  return eltype( K.U )
end

function eltype( P::_SpelMatIterPrec{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  return eltype( P.K.U )
end

import Base.size
function size( K::SpelMatIter{Tv,Ti}, d ) where{Tv<:AbstractFloat,Ti<:Integer}
  return size( K.U, d )
end

import Base.size
function size( P::_SpelMatIterPrec{Tv,Ti}, d ) where{Tv<:AbstractFloat,Ti<:Integer}
  return size( P.K.U, d )
end

function size( P::_SpelMatIterPrec{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  return size( P.K.U )
end

import LinearAlgebra.issymmetric
function issymmetric( P::_SpelMatIterPrec{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  return true
end

struct _SpelMatIterPrecInv{Tv<:AbstractFloat,Ti<:Integer}
  P::_SpelMatIterPrec{Tv,Ti}
end

function *( iP::_SpelMatIterPrecInv{Tv,Ti}, v ) where{Tv<:AbstractFloat,Ti<:Integer}
  return iP.P \ v
end

function eltype( iP::_SpelMatIterPrecInv{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  return eltype( iP.P )
end

function size( iP::_SpelMatIterPrecInv{Tv,Ti} ) where{Tv<:AbstractFloat,Ti<:Integer}
  return size( iP.P )
end

function size( iP::_SpelMatIterPrecInv{Tv,Ti}, d ) where{Tv<:AbstractFloat,Ti<:Integer}
  return size( iP.P, d )
end

import LinearAlgebra.mul!
function mul!( y , P::_SpelMatIterPrecInv, v::AbstractVector )
y .= P * v
end

function mul!( y , K::SpelMatIter, v::AbstractVector )
y .= K * v
end



function computeConditionNumber( P::_SpelMatIterPrec{Tv,Ti} ) where{Tv<:AbstractFloat, Ti<:Integer}
  @show λmax = powm( P; tol=10e-3, verbose=true )[1] 
  @show λmin = one(Tv)/powm(_SpelMatIterPrecInv{Tv,Ti}(P); tol = 10e-3, verbose=true)[1]
  return λmax/λmin
end

function computeMaxEigenvalue( K::SpelMatIter{Tv,Ti}) where{Tv<:AbstractFloat,
                                                           Ti<:Integer}
  return powm( K; tol = 10e-3 )[1]
end

function computeConditionNumber( K::SpelMatIter{Tv,Ti} ) where{Tv<:AbstractFloat,
                                                               Ti<:Integer}
  return computeConditionNumber( _SpelMatIterPrec{Tv,Ti}(K) ) 
end


#the end of Module Spelliptic
end
