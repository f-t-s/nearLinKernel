using LinearAlgebra
import PyPlot
using Plots
#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.5)

using Distances
include("./Utils.jl")

x, Lap, P, revP, xFull, PFull, revPFull, NBound = generateLaplaceBounds2D(7, 5) 

N = size(x,2)  
log4 = r -> log(4.0,r)



BiLap = Lap * Lap'

K = Hermitian( inv( Lap ) )
KOrd = K[ P, P ]

BiK = Hermitian( inv(BiLap) )
BiKOrd = BiK[P,P] 

LapOrd = Lap[ P[end:-1:1], P[end:-1:1] ]
BiLapOrd = BiLap[ P[end:-1:1], P[end:-1:1] ]

colorchoice = "blue"
scatter( x[1,:], x[2,:], 
         markersize = 15 * ( 1.0 ./ (1.0 : 1.0 : Float64(N)) ).^(1/3), 
         malpha = 0.8 * ( 1.0 ./ (1.0 : 1.0 : Float64(N)) ).^(1/6),
         markerstrokealpha = 0.1,
         markercolor = colorchoice,
         legend = false,
         size = ( 525, 500 ) )

savefig("../figures/OrdLex")

scatter( x[1,P], x[2,P], 
         markersize = 30 * ( 1.0 ./ (1.0 : 1.0 : Float64(N)) ).^(1/2), 
         malpha = ( 1.0 ./ (1.0 : 1.0 : Float64(N)) ).^(1/4),
         markerstrokealpha = 0.2,
         markercolor = colorchoice,
         legend = false,
         size = ( 525, 500 ) )

savefig("../figures/OrdMaxMin")


PyPlot.set_cmap(PyPlot.ColorMap("viridis") ) 
maxVal= 1.0
minVal= -16.0
PyPlot.matshow(log10.(abs.(cholesky(BiK).L )),
vmax = maxVal, vmin = minVal)
PyPlot.axis("off")

cb = PyPlot.colorbar()
cb[:ax][:tick_params](labelsize = 20 )
PyPlot.savefig("../figures/IBiLap.jpg")

PyPlot.matshow(log10.(abs.(cholesky(BiKOrd).L )),
        vmax = maxVal, vmin = minVal)

PyPlot.axis("off")
cb = PyPlot.colorbar()
cb[:ax][:tick_params](labelsize = 20 )
PyPlot.savefig("../figures/IBiLapOrd.jpg")
