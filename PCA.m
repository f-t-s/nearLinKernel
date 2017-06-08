%This Matlab file runs the numerical experiments corresponding to the subsection "Sparse approximate PCA"
clear all
close all
q  = 7 ;
%experiment with $q = 7 , l = 0.2 , nu = 1,  rho = 6 , d_x = 0.4 , dz = 0
errorPCA1 = zeros ( q - 1 , 1 ) ;
errorChol1 = zeros ( q - 1 , 1 ) ;
dims = 1 : q - 1 ;
dims = (2.^dims + 1 ).^2 ;

data = matChol ( q , 0.2 , 1 , 6 , 0.2 , 0 ) ;
[ data.eigV data.eigD ] = eig ( data.K ) ;
[ ~ , P ] = sort ( diag ( data.eigD ) , 'descend' ) ;
data.eigV = data.eigV( : , P ) ;
data.eigD = data.eigD( P , P ) ;

for k = 1 : q - 1 
  appK = data.Lichol ( : , 1 : dims ( k )  ) * data.Lichol ( : , 1 : dims ( k ) )' ;
  appKPCA = data.eigV( : , 1 : dims ( k ) ) * data.eigD ( 1 : dims ( k ) , 1 : dims ( k ) ) * data.eigV( : , 1 : dims ( k ) )' ;
  errorPCA1 ( k )  = norm ( appKPCA - data.K ) ;
  errorChol1 ( k )  = norm ( appK - data.K ) ;
end



%experiment with $q = 7 , l = 0.2 , nu = 2 ,  rho = 6 , d_x = 0.4 , dz = 0
errorPCA2 = zeros ( q - 1 , 1 ) ;
errorChol2 = zeros ( q - 1 , 1 ) ;

data = matChol ( q , 0.2 , 2 , 8 , 0.2 , 0 ) ;
[ data.eigV data.eigD ] = eig ( data.K ) ;
[ ~ , P ] = sort ( diag ( data.eigD ) , 'descend' ) ;
data.eigV = data.eigV( : , P ) ;
data.eigD = data.eigD( P , P ) ;


for k = 1 : q - 1
  appK = data.Lichol ( : , 1 : dims ( k ) ) * data.Lichol ( : , 1 : dims ( k ) )' ;
  appKPCA = data.eigV( : , 1 : dims ( k ) ) * data.eigD ( 1 : dims ( k ) , 1 : dims ( k ) ) * data.eigV( : , 1 : dims ( k ) )' ;
  errorPCA2 ( k ) = norm ( appKPCA - data.K ) ;
  errorChol2 ( k ) = norm ( appK - data.K ) ;

end




figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
plot ( dims  , log10 ( errorChol1 ) , '-r'  , dims , log10 ( errorPCA1 ) , '--b' ) ;
legend (  'Cholesky'  , 'PCA' ) ;
xlabel ( 'rank' ) ;
ylabel ( 'log_{10} ( Error )' ) ;
saveas ( figPlot , './figures/PCA1' , 'jpg' ) ;
close all


figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
plot ( dims  , log10 ( errorChol2 ) , '-r'  , dims , log10 ( errorPCA2 ) , '--b' ) ;
legend (  'Cholesky'  , 'PCA' ) ;
xlabel ( 'rank' ) ;
ylabel ( 'log_{10} ( Error )' ) ;
saveas ( figPlot , './figures/PCA2' , 'jpg' ) ;
close all



clear all
close all




   
  

