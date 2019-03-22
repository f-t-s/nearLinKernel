#This script creates Figure 2.5 of the article
include("./spelliptic/CovFuncs.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using Distances
Plots.scalefontsizes()
Plots.scalefontsizes(2.5)
xlimit = [-1.5, 1.5]

nu1 = 0.5
nu2 = 1.0
nu3 = 1.5
nu4 = 2.0

cov1 = r -> matern( abs(r), 1., nu1)
cov2 = r -> matern( abs(r), 1., nu2)
cov3 = r -> matern( abs(r), 1., nu3)
cov4 = r -> matern( abs(r), 1., nu4)

plot( cov1, xlim = xlimit, label = "\$\\nu = $nu1\$", linestyle  = :solid )
plot!( cov2, xlim = xlimit, label = "\$\\nu = $nu2\$", linestyle = :dash )
plot!( cov3, xlim = xlimit, label = "\$\\nu = $nu3\$", linestyle = :dot )
plot!( cov4, xlim = xlimit, label = "\$\\nu = $nu4\$", linestyle = :dashdot )

savefig( "../figures/maternKernels.pdf")

#Number of points used to compute spectrum
N = 8000
x = rand( 2, N )
maxNumEig = N

eigvals1 = sort( eigvals( cov1.( pairwise( Euclidean(), x, x) ./ 0.3 ) ),
                 order=Base.Reverse )
eigvals1 = log10.(eigvals1)


Plots.scalefontsizes()
Plots.scalefontsizes(2.5)

plot( eigvals1[1:maxNumEig], label = "\$\\nu = $nu1\$", linestyle  = :solid, 
     xlabel = "\$k\$",
     ylabel = "\$\\log_{10}(\\lambda_{k})\$",
     )

eigvals2 = sort( eigvals( cov2.( pairwise( Euclidean(), x, x) ./ 0.3 ) ),
                 order=Base.Reverse )
eigvals2 = log10.(eigvals2)
plot!( eigvals2[1:maxNumEig], label = "\$\\nu = $nu2\$", linestyle  = :dash )

eigvals3 = sort( eigvals( cov3.( pairwise( Euclidean(), x, x) ./ 0.3 ) ),
                 order=Base.Reverse )
eigvals3 = log10.(eigvals3)
plot!( eigvals3[1:maxNumEig], label = "\$\\nu = $nu3\$", linestyle  = :dot )

eigvals4 = sort( eigvals( cov4.( pairwise( Euclidean(), x, x) ./ 0.3 ) ),
                 order=Base.Reverse )
eigvals4 = log10.(eigvals4)
plot!( eigvals4[1:maxNumEig], label = "\$\\nu = $nu4\$", linestyle  = :dashdot )


savefig( "../figures/maternKernelsSpec.pdf")
