#This script produces Table 2.6 of the article
include("./spelliptic/Spelliptic.jl")
include("./spelliptic/CovFuncs.jl")
include("./spelliptic/Diagnostics.jl")

using .Spelliptic 
using SparseArrays
using Plots
Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

using TexTables
using LaTeXStrings

rhoList = [ 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 ]
d = 2

NDataPoints = length( rhoList )

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

l = 0.40
α = 0.5
β = 0.025
η = 1.0
sigma = 0.0
cov = r -> cauchy( r, l, α, β )

keyVals = 
[Symbol(:"\$\\rank(L)\$"), 
 Symbol(:"\$ E \\defeq \\frac{\\left\\|LL^T - \\Theta \\right\\|_{\\FRO}}{\\left\\| \\Theta \\right\\|_{\\FRO}}\$" ) ]

for k = 1 : NDataPoints
  rho = rhoList[k]
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

    rankTab = TableCol( "\$\\rho = $rho\$", [Symbol(:"\$\\rank(L)\$")], 
                       rankCholesky[k:k] )
    errTab = TableCol( "\$\\rho = $rho\$",[Symbol(:"\$ E \$")],
                      meanFrobeniusRelError[k:k],stdFrobeniusRelError[k:k] ) 

    errIntTab = TableCol( "\$\\rho = $rho\$",[Symbol(:"\$ \\bar{E} \$")],
                      meanFrobeniusRelErrorInt[k:k],stdFrobeniusRelErrorInt[k:k] ) 


    setfield!.(errTab.data.vals, :format, "{:,2e}")
    setfield!.(errTab.data.vals, :format_se, "{:,2e}")

    setfield!.(errIntTab.data.vals, :format, "{:,2e}")
    setfield!.(errIntTab.data.vals, :format_se, "{:,2e}")
    
    outputTables[k] = vcat( rankTab, errTab, errIntTab ) 
  K = 0
  GC.gc()
end
tab = hcat( outputTables[:]... )
f = open("./out/vary_rho_cauchy_l$(l)_alpha$(α)_beta$(β)_N$(N)_d$(d).tex", "w")
write(f, to_tex( tab, se_pos = :below ) )
close(f)


outputTables = Array{IndexedTable{1,1},1}(undef,NDataPoints)

l = 0.20
α = 1.0
β = 0.2
η = 1.0
sigma = 0.0
cov = r -> cauchy( r, l, α, β )

keyVals = 
[Symbol(:"\$\\rank(L)\$"), 
 Symbol(:"\$ E \\defeq \\frac{\\left\\|LL^T - \\Theta \\right\\|_{\\FRO}}{\\left\\| \\Theta \\right\\|_{\\FRO}}\$" ) ]

for k = 1 : NDataPoints
  rho = rhoList[k]
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

    rankTab = TableCol( "\$\\rho = $rho\$", [Symbol(:"\$\\rank(L)\$")], 
                       rankCholesky[k:k] )

    errTab = TableCol( "\$\\rho = $rho\$",[Symbol(:"\$ E \$")],
                      meanFrobeniusRelError[k:k],stdFrobeniusRelError[k:k] ) 

    errIntTab = TableCol( "\$\\rho = $rho\$",[Symbol(:"\$ \\bar{E} \$")],
                      meanFrobeniusRelErrorInt[k:k],stdFrobeniusRelErrorInt[k:k] ) 


    setfield!.(errTab.data.vals, :format, "{:,2e}")
    setfield!.(errTab.data.vals, :format_se, "{:,2e}")

    setfield!.(errIntTab.data.vals, :format, "{:,2e}")
    setfield!.(errIntTab.data.vals, :format_se, "{:,2e}")
    
    outputTables[k] = vcat( rankTab, errTab, errIntTab ) 
  K = 0
  GC.gc()
end
tab = hcat( outputTables[:]... )
f = open("./out/vary_rho_cauchy_l$(l)_alpha$(α)_beta$(β)_N$(N)_d$(d).tex", "w")
write(f, to_tex( tab, se_pos = :below ) )
close(f)

