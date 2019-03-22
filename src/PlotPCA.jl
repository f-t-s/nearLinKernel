#This script creates Figure 1.9 of the article
include("./spelliptic/Spelliptic.jl")
include("./spelliptic/CovFuncs.jl")
include("./spelliptic/Diagnostics.jl")

using .Spelliptic
using SparseArrays
using Plots 
using LinearAlgebra
using Distances
#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.5)
using LaTeXStrings

NVec = [ 10, 100, 1000, 2000, 4000, 6000, 8000, 10000]
errChol = zeros( size(NVec) ) 
errEig = zeros( size(NVec) ) 

N = 2 * NVec[end]
# create points and order them in maximin ordering, by constructing an 
# auxillary SpelMat
x = rand(2,N)
K = SpelMat{Float64,Int}(x,5.0)
x = x[:,K.P]

global l = 0.5

cov1 = r -> matern( r, l, 1.0)

KTrue = cov1.( pairwise( Euclidean(), x, x ) )
eigFact = eigen( KTrue ) 
cholFact = cholesky( KTrue )

for k = 1 : length( NVec ) 
  n = NVec[k]
  errChol[k] = opnorm( KTrue - cholFact.L[:,1:n] * cholFact.U[1:n,:])
  errEig[k] = opnorm( KTrue - eigFact.vectors[:, (N-n+1) : end] * 
                              diagm( 0 => eigFact.values[(N-n+1):end] ) *
                              eigFact.vectors[:, (N-n+1) : end]' )
end

errChol = log10.( errChol ) 
errEig = log10.( errEig ) 

plot1 = plot( xlabel = "Rank",
             ylabel = L"\log_{10}(\| \Theta - \tilde{\Theta}\|)")
plot!( plot1, NVec, errEig, linestyle = :dash, markershape = :circle,
       label = "PCA"  )
plot!( plot1, NVec, errChol, linestyle = :dashdot, markershape = :star4,
       label = "Cholesky" )

savefig( plot1, "./out/Plot_PCA_nu_1.pdf")

cov2 = r -> matern( r, l, 2.0)

KTrue = cov2.( pairwise( Euclidean(), x, x ) )
GC.gc()
eigFact = eigen( KTrue ) 
cholFact = cholesky( KTrue )

for k = 1 : length( NVec ) 
  n = NVec[k]
  errChol[k] = opnorm( KTrue - cholFact.L[:,1:n] * cholFact.U[1:n,:])
  errEig[k] = opnorm( KTrue - eigFact.vectors[:, (N-n+1) : end] * 
                              diagm( 0 => eigFact.values[(N-n+1):end] ) *
                              eigFact.vectors[:, (N-n+1) : end]' )
end

errChol = log10.( errChol ) 
errEig = log10.( errEig ) 

plot2 = plot( xlabel = "Rank",
             ylabel = L"\log_{10}(\| \Theta - \tilde{\Theta}\|)")
plot!( plot2, NVec, errEig, linestyle = :dash, markershape = :circle,
       label = "PCA"  )
plot!( plot2, NVec, errChol, linestyle = :dashdot, markershape = :star4,
       label = "Cholesky" )

savefig( plot2, "./out/Plot_PCA_nu_2.pdf")







