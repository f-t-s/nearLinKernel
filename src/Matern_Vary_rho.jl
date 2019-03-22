#This script computes Tables 2.3 and 2.4, as well as the lower two panels of 
#Figure 2.4
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

plotRhoErr = plot( xlabel = L"$\rho$",
                  ylabel = L"$\log_{10}(E,\bar{E})$" )


plotCompErr = plot( xlabel = L"$ t_{\texttt{ichol(0)}}$",
                    ylabel = L"$\log_{10}(E)$" )



rhoList = [ 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 ]
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

l = 0.20
η = 1.0
sigma = 0.0
d = 2
nu = 1.0
cov = r -> matern( r, l, nu )
plotScatter2d = plot( cov, 
                      xlabel = L"|x_i - x_j|",
                      label = L"G(x_i,x_j)",
                      xlim = [0.0,1.0],
                      ylim = [0.0,1.0])

plot!( plotScatter2d, 
       inset = bbox( 0.20, 0.00, 0.50, 0.35, :center),
       subplot = 2,
       xlim = [0.35, 0.4],
       ylim = [0.12, 0.20],
       legend = false )

plot!( plotScatter2d,
       r -> cov(abs(r)),
       subplot = 2 )


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

  if rho in [2.0, 3.0]
    scatterVals = createScatterVals( K, x, m_scatterSamples)
    if rho == 2.0
      malphaL = 1.0
      malpha = 0.6
      mshape = :diamond
    else
      mshape = :circle
      malphaL = 0.15
      malpha = 0.30
    end
    scatter!( plotScatter2d, 
              scatterVals...,
              label = latexstring("\\rho = $rho"),
              markersize = 0.50,
              markershape = mshape,
              markeralpha = malphaL, 
              markerstrokealpha = 0.0 )

    scatter!( plotScatter2d, 
              scatterVals...,
              subplot = 2,
              markersize = 2.00,
              markershape = mshape,
              markeralpha = malpha,
              markerstrokealpha = 0.0 )
  end
  K = 0 
  GC.gc()
end

keyVals = [ Symbol( :"\$\\rho= ", k, "\$" ) for k in rhoList ] 

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
#meanFrobeniusError = zeros( NDataPoints )
#stdFrobeniusError = zeros( NDataPoints )
#t_FrobeniusRelError = TableCol( "\$ E \\defeq \\frac{\\left\\|LL^T - \\Theta \\right\\|_{\\FRO}}{\\left\\| \\Theta \\right\\|_{\\FRO}}\$", 
#                                keyVals, 
#                                meanFrobeniusRelError,
#                                stdFrobeniusRelError )
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
            #t_timingTotal,
            t_FrobeniusRelError,
            t_FrobeniusRelErrorInt )

f = open("./out/vary_rho_matern_l$(l)_nu$(nu)_N$(N)_d$(d).tex", "w")
write(f, to_tex( tab, se_pos = :inline ) )
close(f)

plot!( plotRhoErr, 
       rhoList, 
       log10.(meanFrobeniusRelError), 
       label = "\$E: d = 2, \\nu = 1.0\$", 
       markershape=:circle,
       linestyle = :dash)

plot!( plotRhoErr, 
       rhoList, 
       log10.(meanFrobeniusRelErrorInt), 
       label = "\$\\bar{E}: d = 2, \\nu = 1.0\$", 
       markershape=:pentagon,
       linestyle = :dashdot)


plot!( plotCompErr, 
       timingCholesky, 
       log10.(meanFrobeniusRelError), 
       label = "\$d = 2, \\nu = 1.0\$", 
       markershape=:circle,
       linestyle=:dot)

###############################################################################
#This experiment runs the factorization for different values of $N$ and
#a fixed value of $\rho$, for a matern kernel 
rhoList = [ 2.0, 3.0, 4.0, 5.0, 6.0 ] 
NDataPoints = length( rhoList )

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
plotScatter3d = plot( cov, 
                      xlabel = L"|x_i - x_j|",
                      label = L"G(x_i,x_j)",
                      xlim = [0.0, 1.0],
                      ylim = [0.0, 1.0])
          
plot!( plotScatter3d, 
       inset = bbox( 0.20, 0.00, 0.50, 0.35, :center),
       subplot = 2,
       xlim = [0.35, 0.4],
       ylim = [0.12, 0.20],
       legend = false )

plot!( plotScatter3d,
       r -> cov(abs(r)),
       subplot = 2 )


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
  if rho in [2.0, 3.0]
    scatterVals = createScatterVals( K, x, m_scatterSamples)
    if rho == 2.0
      malphaL = 1.0
      malpha = 0.6
      mshape = :diamond
    else
      mshape = :circle
      malphaL = 0.15
      malpha = 0.3
    end
    scatter!( plotScatter3d, 
              scatterVals...,
              label = latexstring("\\rho = $rho"),
              markersize = 0.50,
              markershape = mshape,
              markeralpha = malphaL, 
              markerstrokealpha = 0.0 )

    scatter!( plotScatter3d, 
              scatterVals...,
              subplot = 2,
              markersize = 2.00,
              markershape = mshape,
              markeralpha = malpha,
              markerstrokealpha = 0.0 )
  end
  K = 0 
  GC.gc()
end

keyVals = [ Symbol( :"\$\\rho= ", k, "\$" ) for k in rhoList ] 

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
            #t_timingTotal,
            t_FrobeniusRelError,
            t_FrobeniusRelErrorInt )

f = open("./out/vary_rho_matern_l$(l)_nu$(nu)_N$(N)_d$(d).tex", "w")
write(f, to_tex( tab, se_pos = :inline ) )
close(f)

plot!( plotRhoErr, 
       rhoList, 
       log10.(meanFrobeniusRelError), 
       label = "\$ E: d = 3, \\nu = 0.5\$", 
       markershape=:star4, 
       linestyle=:dashdotdot   )

plot!( plotRhoErr, 
       rhoList, 
       log10.(meanFrobeniusRelErrorInt), 
       label = "\$ \\bar{E}: d = 3, \\nu = 0.5\$", 
       markershape=:cross, 
       linestyle=:dot   )





plot!( plotCompErr, 
       timingCholesky, 
       log10.(meanFrobeniusRelError), 
       label = "\$d = 3, \\nu = 0.5\$", 
       markershape=:circle,
       linestyle=:dashdot)


savefig( plotRhoErr, "./out/Matern_Vary_rho_plotRhoErr.pdf")
savefig( plotCompErr, "./out/Matern_Vary_rho_plotCompErr.pdf")
#output as png, since pdf output seems to produce bug with subplot
savefig( plotScatter2d, "./out/Matern_Vary_rho_plotScatter2d.png")
savefig( plotScatter3d, "./out/Matern_Vary_rho_plotScatter3d.png")
