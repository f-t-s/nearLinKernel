#This script creates Figure 2.3 of the article
using LinearAlgebra
using PyPlot
using Distances
include("./Utils.jl")

x, Lap, P, revP, xFull, PFull, revPFull, NBound = generateLaplaceBounds2D(7, 5) 

#Need to first include boundary points:
l = 0.5
ExpK = exp.( - pairwise( Euclidean(), xFull, xFull)/l )
#Conditioning the boundary to be zero:
ExpKOrd = ExpK[PFull, PFull]
ExpKOrd = ExpKOrd[(NBound+1) : end, (NBound+1) : end] -
  ExpKOrd[(NBound+1) : end, 1 : NBound] * inv( ExpKOrd[ 1 : NBound, 1 : NBound] ) *
    ExpKOrd[1 : NBound, (NBound+1) : end]

ExpK = ExpKOrd[ revP, revP ]
ExpKI = Hermitian(inv(ExpK))
ExpKIOrd = ExpKI[ P[end:-1:1], P[end:-1:1]] 

PyPlot.set_cmap(PyPlot.ColorMap("viridis") ) 
maxVal= 1.0
minVal= -16.0
for (sigma,k) in [ (0.0, 1), (0.1,2), (1.0,3), (10.0,4) ]
  PyPlot.matshow(log10.(abs.(cholesky(Hermitian(ExpKOrd + sigma * I)).L )),
          vmax = maxVal, vmin = minVal)

  PyPlot.axis("off")

  cb = colorbar()
  cb[:ax][:tick_params](labelsize = 20 )
  PyPlot.savefig("../figures/CholKSigma$k.jpg")
end

for (sigma,k) in [ (0.0, 1), (0.1,2), (1.0,3), (10.0,4) ]
  PyPlot.matshow(log10.(abs.(cholesky(Hermitian(ExpKIOrd + sigma * I)).L )),
          vmax = maxVal, vmin = minVal)

  PyPlot.axis("off")
  cb = PyPlot.colorbar()
  cb[:ax][:tick_params](labelsize = 20 )
  PyPlot.savefig("../figures/CholKISigma$k.jpg")
end
