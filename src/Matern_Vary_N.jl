#This script computes Tables 2.1 and 2.2 and the upper two panels of Figure 2.4
#in the article
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

plotNT = plot( xlabel = L"$N$",
               ylabel = L"$t$" )





#This experiment runs the factorization for different values of $N$ and
#a fixed value of $\rho$, for an exponential kernel. 
NList = [ 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000  ]
N = 1000000
m_reps = 50
m_samples = 500000
m_scatterSamples = 10000

NDataPoints = length( NList )

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
d = 2
nu = 0.5
rho = 3.0
cov = r -> matern( r, l, nu )

for k = 1 : NDataPoints
  N = NList[k]
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

  K = 0 
  GC.gc()
end

keyVals = [ Symbol( :"\$N = ", k, "\$" ) for k in NList ] 

t_timingCholesky = TableCol( "\$\\tich\$", keyVals, timingCholesky)
setfield!.(t_timingCholesky.data.vals, :format, "{:,2f}")
t_timingSortSparse = TableCol( "\$\\tsos\$", keyVals, timingSortSparse)
setfield!.(t_timingSortSparse.data.vals, :format, "{:,2f}")
#timingSetCov = zeros( NDataPoints )
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

t_FrobeniusRelErrorInt = TableCol( "\$ \\bar{E} \$", 
                                   keyVals, 
                                   meanFrobeniusRelErrorInt,
                                   stdFrobeniusRelErrorInt )



setfield!.(t_FrobeniusRelError.data.vals, :format, "{:,2e}")
setfield!.(t_FrobeniusRelError.data.vals, :format_se, "{:,2e}")
setfield!.(t_FrobeniusRelErrorInt.data.vals, :format, "{:,2e}")
setfield!.(t_FrobeniusRelErrorInt.data.vals, :format_se, "{:,2e}")

tab = hcat( t_nnzCholesky,
            t_rankCholesky,
            t_timingSortSparse, 
            t_timingComputeEntries,
            t_timingCholesky, 
            t_FrobeniusRelError,
            t_FrobeniusRelErrorInt )

f = open("./out/vary_N_matern_l$(l)_nu$(nu)_rho$(rho)_d$(d).tex", "w")
write(f, to_tex( tab, se_pos = :inline ) )
close(f)

plot!( plotNT, 
       NList, 
       timingCholesky, 
       label = "\$d = 2\$", 
       markershape=:circle,
       linestyle = :dash)


###############################################################################
#This experiment runs the factorization for different values of $N$ and
#a fixed value of $\rho$, for a matern kernel 
NDataPoints = length( NList )

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

d = 3
nu = 0.5
cov = r -> matern( r, l, nu )

for k = 1 : NDataPoints
  N = NList[k]
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

  K = 0 
  GC.gc()
end

keyVals = [ Symbol( :"\$N = ", k, "\$" ) for k in NList ] 

t_timingCholesky = TableCol( "\$\\tich\$", keyVals, timingCholesky)
setfield!.(t_timingCholesky.data.vals, :format, "{:,2f}")
t_timingSortSparse = TableCol( "\$\\tsos\$", keyVals, timingSortSparse)
setfield!.(t_timingSortSparse.data.vals, :format, "{:,2f}")
#timingSetCov = zeros( NDataPoints )
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

t_FrobeniusRelErrorInt = TableCol( "\$ \\bar{E} \$", 
                                   keyVals, 
                                   meanFrobeniusRelErrorInt,
                                   stdFrobeniusRelErrorInt )



setfield!.(t_FrobeniusRelError.data.vals, :format, "{:,2e}")
setfield!.(t_FrobeniusRelError.data.vals, :format_se, "{:,2e}")
setfield!.(t_FrobeniusRelErrorInt.data.vals, :format, "{:,2e}")
setfield!.(t_FrobeniusRelErrorInt.data.vals, :format_se, "{:,2e}")

tab = hcat( t_nnzCholesky,
            t_rankCholesky,
            t_timingSortSparse, 
            t_timingComputeEntries,
            t_timingCholesky, 
            t_FrobeniusRelError,
            t_FrobeniusRelErrorInt )

f = open("./out/vary_N_matern_l$(l)_nu$(nu)_rho$(rho)_d$(d).tex", "w")
write(f, to_tex( tab, se_pos = :inline ) )
close(f)

plot!( plotNT, 
       NList, 
       timingCholesky, 
       label = "\$ d = 3 \$", 
       markershape=:star4, 
       linestyle=:dashdotdot   )




savefig( plotNT, "./out/Matern_Vary_N_plotNT.pdf")
