#This script computes Table 2.8 of the article
include("./spelliptic/Spelliptic.jl")
include("./spelliptic/CovFuncs.jl")
include("./spelliptic/Diagnostics.jl")
include("./GenerateHighDimPoints.jl")

using .Spelliptic 
using SparseArrays

using TexTables
using LaTeXStrings

using Plots
Plots.backend(:pyplot)
#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.5)

N = 1000000
NPlot = 10000
rhoList = [ 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 ]
NDataPoints = length( rhoList )
m_reps = 50
m_samples = 500000

x, x1, x2, x3, x4 = generateSet( N )

plotx1 = scatter( x1[1,1:NPlot], x1[2,1:NPlot], x1[3,1:NPlot],
                  legend = false,
                  grid = false,
                  markerstrokewidth = 0.1,
                  markersize = 1.1,
                  zlim = [-1.0,1.0],
                  xlim = [-1.0,1.0],
                  ylim = [-1.0,1.0]
                )

plotx2 = scatter( x2[1,1:NPlot], x2[2,1:NPlot], x2[3,1:NPlot],
                  legend = false,
                  grid = false,
                  markerstrokewidth = 0.1,
                  markersize = 1.1,
                  zlim = [-1.0,1.0],
                  xlim = [-1.0,1.0],
                  ylim = [-1.0,1.0]
                )

plotx3 = scatter( x3[1,1:NPlot], x3[2,1:NPlot], x3[3,1:NPlot],
                  legend = false,
                  grid = false,
                  markerstrokewidth = 0.1,
                  markersize = 1.1,
                  zlim = [-1.0,1.0],
                  xlim = [-1.0,1.0],
                  ylim = [-1.0,1.0]
                )

plotx4 = scatter( x4[1,1:NPlot], x4[2,1:NPlot], x4[3,1:NPlot],
                  legend = false,
                  grid = false,
                  markerstrokewidth = 0.1,
                  markersize = 1.1,
                  zlim = [-1.0,1.0],
                  xlim = [-1.0,1.0],
                  ylim = [-1.0,1.0]
                )

savefig( plotx1, "./out/points_x1_N_$N.pdf")
savefig( plotx2, "./out/points_x2_N_$N.pdf")
savefig( plotx3, "./out/points_x3_N_$N.pdf")
savefig( plotx4, "./out/points_x4_N_$N.pdf")



global l = 0.5
global nu = 0.5
global cov = r -> matern( r, l, nu )
η = 1.0
sigma = 0.0

timingCholesky = zeros( NDataPoints )
timingSortSparse = zeros( NDataPoints )
timingSetCov = zeros( NDataPoints )
timingComputeEntries = zeros( NDataPoints )
timingTotal = zeros( NDataPoints )
nnzCholesky = zeros( NDataPoints )
rankCholesky = zeros( Integer, NDataPoints )
meanFrobeniusError = zeros( NDataPoints )
stdFrobeniusError = zeros( NDataPoints )
meanFrobeniusRelError = zeros( NDataPoints )
stdFrobeniusRelError = zeros( NDataPoints )


for k = 1 : NDataPoints
  rho = rhoList[k]
  @show k
  K, timingSortSparse[k] = @timed SpelMat{Float64,Int}(x, rho)
  #Computing the entries on the sparsity pattern and computing the Cholesky
  #factorization
  timingCholesky[k], timingSetCov[k] = 
    @timed setCovStat!( K, x, cov, η, sigma)
  #adjusting for the double counting of the Cholesky time 
  timingComputeEntries[k] = timingSetCov[k] - timingCholesky[k]
  timingTotal[k] = timingSetCov[k] + timingSortSparse[k]
  meanFrobeniusError[k], stdFrobeniusError[k] = 
    estErrorFrobeniusStat( K, x, cov, η, sigma, m_reps, m_samples) 
  meanFrobeniusRelError[k], stdFrobeniusRelError[k] = 
    estRelErrorFrobeniusStat( K, x, cov, η, sigma, m_reps, m_samples) 

  nnzCholesky[k] = nnz( K.U )/(N^2)
  rankCholesky[k] = rank( K )

  K = 0 
  GC.gc()
end

keyVals = [ Symbol( :"\$\\rho= ", k, "\$" ) for k in rhoList ] 

t_timingCholesky = TableCol( "\$\\tich\$", keyVals, timingCholesky)
setfield!.(t_timingCholesky.data.vals, :format, "{:,2f}")
t_timingSortSparse = TableCol( "\$\\tsos\$", keyVals, timingSortSparse)
setfield!.(t_timingSortSparse.data.vals, :format, "{:,2f}")
t_timingComputeEntries = TableCol( "\$\\tent\$", keyVals, timingComputeEntries)
setfield!.(t_timingComputeEntries.data.vals, :format, "{:,2f}")

t_timingTotal = TableCol( "\$\\ttot\$", keyVals, timingTotal)
setfield!.(t_timingTotal.data.vals, :format, "{:,2f}")
t_nnzCholesky = TableCol( "\$\\nnz(L)/N^2\$", keyVals, nnzCholesky)
setfield!.(t_nnzCholesky.data.vals, :format, "{:,2e}")
t_rankCholesky = TableCol( "\$\\rank(L)\$", keyVals, rankCholesky)
t_FrobeniusRelError = TableCol( "\$ E \$", 
                                keyVals, 
                                meanFrobeniusRelError,
                                stdFrobeniusRelError )




setfield!.(t_FrobeniusRelError.data.vals, :format, "{:,2e}")
setfield!.(t_FrobeniusRelError.data.vals, :format_se, "{:,2e}")

tab = hcat( t_nnzCholesky,
            t_rankCholesky,
            t_timingSortSparse, 
            t_timingComputeEntries,
            t_timingCholesky, 
            t_FrobeniusRelError )

f = open("./out/vary_rho_highDim_l_$(l)_nu_$(nu)_N$(N).tex", "w")
write(f, to_tex( tab, se_pos = :inline ) )
close(f)

