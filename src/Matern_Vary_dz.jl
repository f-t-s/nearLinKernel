#This script computes Table 2.7 in the article
include("./spelliptic/Spelliptic.jl")
include("./spelliptic/CovFuncs.jl")
include("./spelliptic/Diagnostics.jl")

using .Spelliptic 
using SparseArrays


using Plots
Plots.backend(:pyplot)
#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.5)

using TexTables
using LaTeXStrings

plotRhoErr = plot( xlabel = L"$\rho$",
                  ylabel = L"$\log_{10}(E,\bar{E})$" )


plotCompErr = plot( xlabel = L"$ t_{\texttt{ichol(0)}}$",
                    ylabel = L"$\log_{10}(E)$" )



#This experiment runs the factorization for different values of $N$ and
#a fixed value of $\rho$, for an exponential kernel. 
dzList = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 ]
NDataPoints = length( dzList )

N = 1000000
m_reps = 50
m_samples = 500000
m_scatterSamples = 10000

NPlot = min( 20000, N )

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
meanFrobeniusRelErrorInt = zeros( NDataPoints )
stdFrobeniusRelErrorInt = zeros( NDataPoints )

l = 0.20
η = 1.0
sigma = 0.0
rho = 3.0
nu = 0.5
cov = r -> matern( r, l, nu )

x = rand( 3, N ) 


outputTables = Array{IndexedTable{1,1},1}(undef,NDataPoints)
keyVals = 
[Symbol(:"\$\\rank(L)\$"), 
 Symbol(:"\$ E \\defeq \\frac{\\left\\|LL^T - \\Theta \\right\\|_{\\FRO}}{\\left\\| \\Theta \\right\\|_{\\FRO}}\$" ) ]


for k = 1 : NDataPoints
  dz = dzList[k]
  x[3,:] = -dz * sin.( 6 * x[1,:] ) .* cos.( 2 * ( 1 .- x[2,:] ) ) .+ 
    dz * randn(N) * 0.001
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
  meanFrobeniusRelErrorInt[k], stdFrobeniusRelErrorInt[k] = 
    estRelErrorFrobeniusStatInt( K, x, cov, η, sigma, m_reps, m_samples) 
    nnzCholesky[k] = nnz( K.U )/(N^2)
    rankCholesky[k] = rank( K )

    rankTab = TableCol( "\$\\delta_{z} = $dz\$", [Symbol(:"\$\\rank(L)\$")], 
                       rankCholesky[k:k] )
    errTab = TableCol( "\$\\delta_{z} = $dz\$",[Symbol(:"\$ E \$")],
                      meanFrobeniusRelError[k:k],stdFrobeniusRelError[k:k] ) 

    errTabInt = TableCol( "\$\\delta_{z} = $dz\$",[Symbol(:"\$ \\bar{E} \$")],
                      meanFrobeniusRelErrorInt[k:k],stdFrobeniusRelErrorInt[k:k] ) 

    
    spaceTab = TableCol( "\$\\delta_{z} = $dz\$",[Symbol(:"\$\\nnz(L)/N^2\$")],
                      nnzCholesky[k:k] ) 

    timeTab = TableCol( "\$\\delta_{z} = $dz\$",[Symbol(:"\$ \\tich \$")],
                      timingCholesky[k:k] ) 



    setfield!.(spaceTab.data.vals, :format, "{:,2e}")
    setfield!.(timeTab.data.vals, :format, "{:,2f}")

    setfield!.(errTab.data.vals, :format, "{:,2e}")
    setfield!.(errTab.data.vals, :format_se, "{:,2e}")

    setfield!.(errTabInt.data.vals, :format, "{:,2e}")
    setfield!.(errTabInt.data.vals, :format_se, "{:,2e}")
    
    outputTables[k] = vcat( spaceTab, timeTab,  rankTab, errTab)

    if dz <= 0.5
      fig = scatter( x[1, 1:NPlot], x[2,1:NPlot], x[3,1:NPlot], 
                        legend = false,
                        grid = false,
                        markerstrokewidth = 0.1,
                        markersize = 1.1,
                        zlim = [-0.5, 0.5],
                        xlim = [0.0,1.0],
                        ylim = [0.0,1.0],
                        camera = (30, 45) ) 
      savefig( fig, "./out/points_dz_$(dz)_N_$N.pdf" )
    end
  K = 0
  GC.gc()
end


tab = hcat( outputTables[:]... )

f = open("./out/vary_dz_matern_l$(l)_rho$(rho)_N$(N).tex", "w")
write(f, to_tex( tab, se_pos = :below ) )
close(f)
