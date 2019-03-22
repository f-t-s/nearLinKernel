#This script computes Table 2.5 of the article
include("./spelliptic/Spelliptic.jl")
include("./spelliptic/CovFuncs.jl")
include("./spelliptic/Diagnostics.jl")

using .Spelliptic 
using SparseArrays
using Plots
#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.5)

using TexTables
using LaTeXStrings


nuList = [ 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7 ]
d = 3
rho = 5.0

NDataPoints = length( nuList )

N = 1000000
m_reps = 50
m_samples = 500000
m_scatterSamples = 10000

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

outputTables = Array{IndexedTable{1,1},1}(undef,NDataPoints)

l = 0.20
η = 1.0
sigma = 0.0

keyVals = 
[Symbol(:"\$\\rank(L)\$"), 
 Symbol(:"\$ E \\defeq \\frac{\\left\\|LL^T - \\Theta \\right\\|_{\\FRO}}{\\left\\| \\Theta \\right\\|_{\\FRO}}\$" ) ]

for k = 1 : NDataPoints
  nu = nuList[k]
  cov = r -> matern( r, l, nu )
  @show k
  x = rand( d, N ) 
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

    rankTab = TableCol( "\$\\nu = $nu\$", [Symbol(:"\$\\rank(L)\$")], 
                       rankCholesky[k:k] )
    errTab = TableCol( "\$\\nu = $nu\$",[Symbol(:"\$ E \$")],
                      meanFrobeniusRelError[k:k],stdFrobeniusRelError[k:k] ) 

    errTabInt = TableCol( "\$\\nu = $nu\$",[Symbol(:"\$ \\bar{E} \$")],
                      meanFrobeniusRelErrorInt[k:k],stdFrobeniusRelErrorInt[k:k] ) 

    setfield!.(errTab.data.vals, :format, "{:,2e}")
    setfield!.(errTab.data.vals, :format_se, "{:,2e}")

    setfield!.(errTabInt.data.vals, :format, "{:,2e}")
    setfield!.(errTabInt.data.vals, :format_se, "{:,2e}")
    
    outputTables[k] = vcat( rankTab, errTab, errTabInt ) 
  K = 0
  GC.gc()
end


tab = hcat( outputTables[:]... )

f = open("./out/vary_nu_matern_l$(l)_rho$(rho)_N$(N)_d$(d).tex", "w")
write(f, to_tex( tab, se_pos = :below ) )
close(f)

