%This Matlab-file contains the code producing the numerical examples in subsection "Detoriation of exponential decay at the boundary"
clear all
close all
q = 7 ;
data = matChol ( q , 0.4 , 1 , 2 , 0.2 , 0 ) ;
data.appK = data.Lichol * data.Lichol' ;
data.L = chol ( data.K )' ;
distx1 = sqrt ( sum ( ( data.x - repmat ( data.x ( : , data.n ) , 1 , data.n ) ).^2 ) ) ; 

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( distx1 , data.K ( : , data.n ) , 1 , 'filled') ; 
hold on 
scatter ( distx1( 1 : 5 : end ) , data.appK ( 1 : 5 : end , data.n ) , 10 , '+' ) ; 
legend ( 'K', 'K_{approx}' ) ;
xlabel ( '||x_{N} - x_{j}||' ) ;
ylabel ( 'K_{N,j}' ) ;
saveas ( figPlot , './figures/boundaryScatter' , 'jpg' ) ;
close all

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( data.x(1 , : ) , data.x ( 2 , : ) , 100 ,  data.appK ( data.n , : )  , 'filled') ; 
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
colorbar;
saveas ( figPlot , './figures/boundaryScatterSpatial' , 'jpg' ) ;
close all

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( data.x(1 , : ) , data.x ( 2 , : ) , 100 ,  log10 ( abs ( data.L ( : , (2^(q-1)+1)^2 + round(( (2^(q-0)+1)^2 - (2^(q-1)+1)^2)/2 ) ) ) )  , 'filled') ; 
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
legend ( 'log_{10} ( L_{:,i} ) ' )
colorbar;
saveas ( figPlot , './figures/interiorScatterL' , 'jpg' ) ;
close all


figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( data.x(1 , : ) , data.x ( 2 , : ) , 100 ,  log10 ( abs ( data.L ( : , (2^(q-1)+1)^2 + 1 ) ) )  , 'filled') ; 
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
colorbar;
legend ( 'log_{10} ( L_{:,j} ) ' )
saveas ( figPlot , './figures/boundaryScatterL' , 'jpg' ) ;
close all


clear all
close all
