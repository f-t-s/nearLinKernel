%This file contains the code creating the numerical experiments in subsection "Fractional operators".
clear all
close all
q  = 7 ;
%experiment with $q = 7 , l = 0.2 , nu = 0.4 ,  rho = 2 : 2 : 10 , d_x = 0.2 , dz = 0
errorAbsOp1 = zeros ( 7 , 1 ) ;
errorRelOp1 = zeros ( 7 , 1 ) ;
errorAbsFro1 = zeros ( 7 , 1 ) ;
errorRelFro1 = zeros ( 7 , 1 ) ;
absnnz1 = zeros ( 7 , 1 ) ;
relnnz1 = zeros ( 5 , 1 ) ;
nuVec1 = [ 1 , 1.15 , 1.33 , 1.5 , 1.62, 1.83, 2 ]' ;


input.dataFormat = {'%1.1f' , 1 , '%1.3e' , 6 } ;

for k = 1 : 7
  nu = nuVec1 ( k ) ;
  data = matChol ( q , 0.2 , nu , 6 , 0.2 , 0 ) ;
  data.appK = data.Lichol * data.Lichol';
  normOpK = norm ( data.K ) ;
  normFroK = norm ( data.K , 'fro' ) ;
  errorAbsOp1 ( k )  = norm ( data.appK - data.K ) ;
  errorAbsFro1 ( k ) = norm ( data.appK - data.K , 'fro' ) ;
  errorRelOp1 ( k )  = errorAbsOp1 ( k ) / normOpK ; 
  errorRelFro1 ( k )  = errorAbsOp1 ( k ) / normFroK ; 
  absnnz1 ( k ) = nnz ( data.mask ) ;
  relnnz1 ( k ) = nnz ( data.mask )/ (data.n)^2 ;
end
input.tableCaption = ( 'Compression and Accuracy for $q = 7, \ l = 0.2 , \ \rho  = 6 , \ \delta_x = 0.2$ and different values of $\nu$' )
input.tableLabel = ( 'fractional' ) ;
input.data =  horzcat ( nuVec1 , errorAbsOp1 , errorRelOp1 , errorAbsFro1 , errorRelFro1 , absnnz1 , relnnz1 ) ; 
input.tableColLabels = {'$\nu$' , '$\| \KMC^{\rho} - \KMC\|$' ,  '$\| \KMC^{\rho} - \KMC\| / \| \KMC\|$' , '$\| \KMC^{\rho} - \KMC\|_{\FRO}$', '$\| \KMC^{\rho} - \KMC\|_{\FRO} / \| \KMC\|_{\FRO}$' , '$\# S$' , '$\#S / N^2$' } ;
textCellArray = latexTable ( input ) ;

delete ./figures/table_fractional.tex ;
fid = fopen ( './figures/table_fractional.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;

clear all
close all
