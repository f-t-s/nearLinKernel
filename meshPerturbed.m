%This matlab file runs the experiments of the subsection "Irregularity of the grid and points on a submanifold"
close all ;
clear all ;
q  = 7 ;
%experiment with $q = 7 , l = 0.2 , nu = 0.1 , rho = 4 , d_x = {0.2 , 0.4 , 2 , 4 } , dz = 0
errorAbsOp1 = zeros ( 4 , 1 ) ;
errorRelOp1 = zeros ( 4 , 1 ) ;
errorAbsFro1 = zeros ( 4 , 1 ) ;
errorRelFro1 = zeros ( 4 , 1 ) ;
absnnz1 = zeros ( 4 , 1 ) ;
relnnz1 = zeros ( 4 , 1 ) ;
deltaxVec1 = [0.2 , 0.4 , 2 , 4]' ;
x1 = cell ( 4 , 1 ) ;

input.dataFormat = {'%1.1f' , 1 , '%1.3e' , 6 } ;

for k = 1 : 4
  deltax = deltaxVec1 ( k ) ;
  data = matChol ( q , 0.2 , 1 , 5 , deltax , 0 ) ;
  data.appK = data.Lichol * data.Lichol';
  normOpK = norm ( data.K ) ;
  normFroK = norm ( data.K , 'fro' ) ;
  errorAbsOp1 ( k )  = norm ( data.appK - data.K ) ;
  errorAbsFro1 ( k ) = norm ( data.appK - data.K , 'fro' ) ;
  errorRelOp1 ( k )  = errorAbsOp1 ( k ) / normOpK ; 
  errorRelFro1 ( k )  = errorAbsOp1 ( k ) / normFroK ; 
  absnnz1 ( k ) = nnz ( data.mask ) ;
  relnnz1 ( k ) = nnz ( data.mask )/ (data.n)^2 ;
  x1{k} = data.x ;
end
input.tableCaption = ( 'Compression and Accuracy for $q = 7, \ l = 0.2 , \ \rho = 5, \ \nu  = 1$  and different values of $\delta_x$' )
input.tableLabel = ( 'perturbedX' ) ;
input.data =  horzcat ( deltaxVec1 , errorAbsOp1 , errorRelOp1 , errorAbsFro1 , errorRelFro1 , absnnz1 , relnnz1 ) ; 
input.tableColLabels = {'$\delta_x$' , '$\| \KMC^{\rho} - \KMC\|$' ,  '$\| \KMC^{\rho} - \KMC\| / \| \KMC\|$' , '$\| \KMC^{\rho} - \KMC\|_{\FRO}$', '$\| \KMC^{\rho} - \KMC\|_{\FRO} / \| \KMC\|_{\FRO}$' , '$\# S$' , '$\#S / N^2$' } ;
textCellArray = latexTable ( input ) ;

delete ./figures/table_perturbedX.tex ;
fid = fopen ( './figures/table_perturbedX.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( x1{1}( 1 , : ) , x1{1}( 2 , : ) ,  3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
saveas ( figPlot , './figures/scatterPlotdx02' , 'jpg' ) ;

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( x1{2}( 1 , : ) , x1{2}( 2 , : ) ,  3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
saveas ( figPlot , './figures/scatterPlotdx04' , 'jpg' ) ;

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( x1{3}( 1 , : ) , x1{3}( 2 , : ) ,  3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
saveas ( figPlot , './figures/scatterPlotdx20' , 'jpg' ) ;

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter ( x1{4}( 1 , : ) , x1{4}( 2 , : ) ,  3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
saveas ( figPlot , './figures/scatterPlotdx40' , 'jpg' ) ;


%experiment with $q = 7 , l = 0.2 , nu = 1 ,  rho = 4 , d_x = 0.2 , dz = { 0 , 0.1 , 0.2 , 0.4 }
errorAbsOp2 = zeros ( 4 , 1 ) ;
errorRelOp2 = zeros ( 4 , 1 ) ;
errorAbsFro2 = zeros ( 4 , 1 ) ;
errorRelFro2 = zeros ( 4 , 1 ) ;
absnnz2 = zeros ( 4 , 1 ) ;
relnnz2 = zeros ( 4 , 1 ) ;
deltazVec2 = [0.0 , 0.1 , 0.2 , 0.4]' ;
x2 = cell ( 4 , 1 ) ;

input.dataFormat = {'%1.1f' , 1 , '%1.3e' , 6 } ;
inputTrunc.dataFormat = {'%1.1f' , 1 , '%1.3e' , 3 } ;
inputTrunc.transposeTable = 1 ;

for k = 1 : 4
  deltaz = deltazVec2 ( k ) ;
  data = matChol ( q , 0.2 , 1 , 5 , 2.0 , deltaz ) ;
  data.appK = data.Lichol * data.Lichol';
  normOpK = norm ( data.K ) ;
  normFroK = norm ( data.K , 'fro' ) ;
  errorAbsOp2 ( k )  = norm ( data.appK - data.K ) ;
  errorAbsFro2 ( k ) = norm ( data.appK - data.K , 'fro' ) ;
  errorRelOp2 ( k )  = errorAbsOp1 ( k ) / normOpK ; 
  errorRelFro2 ( k )  = errorAbsOp1 ( k ) / normFroK ; 
  absnnz2 ( k ) = nnz ( data.mask ) ;
  relnnz2 ( k ) = nnz ( data.mask )/ (data.n)^2 ;
  x2{k} = data.x ;
end
input.tableCaption = ( 'Compression and Accuracy for $q = 7, \ l = 0.2 , \ \rho = 5, \ \nu  = 1 , \ \delta_x = 2$  and different values of $\delta_z$' )
input.tableLabel = ( 'perturbedZ' ) ;
input.data =  horzcat ( deltazVec2 , errorAbsOp2 , errorRelOp2 , errorAbsFro2 , errorRelFro2 , absnnz2 , relnnz2 ) ; 
input.tableColLabels = {'$\delta_z$' , '$\| \KMC^{\rho} - \KMC\|$' ,  '$\| \KMC^{\rho} - \KMC\| / \| \KMC\|$' , '$\| \KMC^{\rho} - \KMC\|_{\FRO}$', '$\| \KMC^{\rho} - \KMC\|_{\FRO} / \| \KMC\|_{\FRO}$' , '$\# S$' , '$\#S / N^2$' } ;
textCellArray = latexTable ( input ) ;

delete ./figures/table_perturbedZ.tex ;
fid = fopen ( './figures/table_perturbedZ.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter3 ( x2{1}( 1 , : ) , x2{1}( 2 , : ) , x2{1}( 3 , : ),   3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
zlim ( [ -0.5 , 0.5 ] )
saveas ( figPlot , './figures/scatterPlotdz00' , 'jpg' ) ;

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter3 ( x2{2}( 1 , : ) , x2{2}( 2 , : ) , x2{2}( 3 , : ) ,  3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
zlim ( [ -0.5 , 0.5 ] )
saveas ( figPlot , './figures/scatterPlotdz01' , 'jpg' ) ;

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter3 ( x2{3}( 1 , : ) , x2{3}( 2 , : ) , x2{3}( 3 , : ) ,  3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
zlim ( [ -0.5 , 0.5 ] )
saveas ( figPlot , './figures/scatterPlotdz02' , 'jpg' ) ;
close all

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
scatter3 ( x2{4}( 1 , : ) , x2{4}( 2 , : ) , x2{4}( 3 , : ) ,  3 , 'filled'  ) ;
xlim ( [ -0.05 , 1.05 ] ) ;
ylim ( [ -0.05 , 1.05 ] ) ;
zlim ( [ -0.5 , 0.5 ] )
saveas ( figPlot , './figures/scatterPlotdz04' , 'jpg' ) ;
close all









clear all ;
close all ;
