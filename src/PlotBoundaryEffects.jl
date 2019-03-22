#This script computes Figure 2.2 of the article
using LinearAlgebra
import PyPlot 
using LaTeXStrings
using Plots
#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.3)
using Distances
include("./Utils.jl")
include("./spelliptic/Spelliptic.jl")
using .Spelliptic 

#Setting the random seed, for reproducibility
import Random.seed!
seed!(1234)

#We only use this function to create a boundary layer, on which we can 
#condition on
x, Lap, P, revP, xFull, PFull, revPFull, NBound = generateLaplaceBounds2D(7, 10) 

#Rescaling the points so that they form an external boundary layer around
#[0,1]^2
xFull = ( xFull .- minimum( x ) ) ./ ( maximum( x ) - minimum( x ) )
xBound = xFull[ :, P[ 1 : NBound ] ]

#Constructing inner points
N = 30000
x = rand( 2, N )


l = 0.5
cov = r -> exp(-r/l)
kmin = floor(Int, N/5)

KNoBounds = cov.( pairwise( Euclidean(), x, x) )
KTotal = cov.( pairwise( Euclidean(), hcat(xBound, x), hcat(xBound, x) ) )
KBounds = Hermitian( KTotal[ (NBound + 1) : end, (NBound + 1) : end] - 
                     KTotal[ (NBound + 1) : end, 1 : NBound ] *
                     inv( KTotal[ 1 : NBound, 1 : NBound ] ) * 
                     KTotal[ 1 : NBound, (NBound + 1) : end ] )

KSP = SpelMat{Float64,Int}(x, 3.0 )
setCovStat!( KSP, x, cov, 1.0, 0.0)
KApprox = Matrix( KSP.U )' * Matrix( KSP.U )

P = KSP.P 
revp = KSP.revP

#Reordering everything to avoid confusion regarding the different orderings
x = x[ :, P ]
KNoBounds = KNoBounds[P,P]
KBounds = KBounds[P,P]

maxEl = findmax( x[ 1, kmin : end ] )[2] + ( kmin - 1 )

centEl = findmin( sum( (x[:, maxEl : (maxEl + 10) ] .- 0.5).^2, dims=1) )[2][2] + ( maxEl - 1 )

LNoBounds = cholesky( KNoBounds ).L
LBounds = cholesky( KBounds ).L
LApprox = Matrix( KSP.U )'

SNoBounds = ( LNoBounds[:, maxEl : end] * LNoBounds[:, maxEl : end]' )

SNoBoundsMaxEl = SNoBounds[maxEl, maxEl : end] ./ SNoBounds[maxEl,maxEl]
SNoBoundsCentEl = SNoBounds[centEl, centEl : end] ./ SNoBounds[centEl,centEl]

SBounds = ( LBounds[:, maxEl : end] * LBounds[:, maxEl : end]' )

SBoundsMaxEl = SBounds[maxEl, maxEl : end] ./ SBounds[maxEl,maxEl]
SBoundsCentEl = SBounds[centEl, centEl : end] ./ SBounds[centEl,centEl]

SApprox = (LApprox[:, maxEl : end] * LApprox[:, maxEl : end]')
SApproxMaxEl = SApprox[maxEl, maxEl : end] ./ SApprox[maxEl,maxEl]
SApproxCentEl = SApprox[centEl, centEl : end] ./ SApprox[centEl,centEl]

climits = ( -16.00, 0.0 )
scatter( x[1, maxEl:end], x[2, maxEl:end], legend = false, colorbar = true,
        marker_z = log10.(abs.(SNoBoundsMaxEl)), clims = climits, color = :viridis)

savefig("../figures/scatterSchurMaxNoBoundVal")

scatter( x[1, maxEl:end], x[2, maxEl:end], legend = false, colorbar = true,
        marker_z = log10.(abs.(SBoundsMaxEl)), clims = climits, color = :viridis)

savefig("../figures/scatterSchurMaxBoundVal")

scatter( x[1, centEl:end], x[2, centEl:end], legend = false, colorbar = true,
        marker_z = log10.(abs.(SNoBoundsCentEl)), clims = climits, color = :viridis)

savefig("../figures/scatterSchurCentNoBoundVal")

scatter( x[1, centEl:end], x[2, centEl:end], legend = false, colorbar = true,
        marker_z = log10.(abs.(SBoundsCentEl)), clims = climits, color = :viridis)

savefig("../figures/scatterSchurCentBoundVal")


LNoBounds = 0
LBounds = 0
SNoBounds = 0
SBounds = 0 
GC.gc()

scatterByDist = plot( r -> cov(abs(r)),
                      xlabel = L"r = |x_k - x_j|",
                      label = L"\exp(-2r)",
                      xlim = [0.0,sqrt(2.0)],
                      ylim = [0.0,1.0] )

rvals = pairwise( Euclidean(), x, x)[:, maxEl]

scatter!( scatterByDist, 
          rvals, 
          KApprox[:,maxEl],
          xlim = [0.0, sqrt(2)],
          label = L"\left( L^{\rho} L^{\rho, \top} \right)_{kj}",
          markershape = :circle,
          markersize = 0.05 )

savefig(scatterByDist, "../figures/scatterByDistNoBoundVal")

