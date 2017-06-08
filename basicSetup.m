%This Matlab-file creates the result of subsection "The basic setup"
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
rhoVec1 = [2:2:14]' ;
errorRelOpRest1 = [7:1]' ;

errorFroTruncCommutator1 = zeros ( 7 ,1 ) ;
errorFroTrunc1 = zeros ( 7 ,1 ) ;
errorFroTruncRel1 = zeros ( 7 ,1 ) ;

input.dataFormat = {'%1.1f' , 1 , '%1.3e' , 6 } ;
inputTrunc.dataFormat = {'%1.1f' , 1 , '%1.3e' , 3 } ;
inputTrunc.transposeTable = 1 ;

for k = 1 : 7
  rho = rhoVec1 ( k ) ;
  data = matChol ( q , 0.2 , 1 , rho , 0.2 , 0 ) ;
  data.appK = data.Lichol * data.Lichol';
  data.L = chol ( data.K )' ;
  data.LTrunc = data.mask .* data.L ;
  data.appKTrunc = data.LTrunc * data.LTrunc' ;
  errorFroTruncCommutator1 ( k )  = norm ( data.appK - data.appKTrunc , 'fro' ) ;
  errorFroTrunc1 ( k )  = norm ( data.K - data.appKTrunc , 'fro' ) ;
  errorFroTruncRel1 ( k ) = errorFroTruncCommutator1 ( k ) / errorFroTrunc1 ( k ) ;
  normOpK = norm ( data.K ) ;
  normFroK = norm ( data.K , 'fro' ) ;
  errorAbsOp1 ( k )  = norm ( data.appK - data.K ) ;
  errorAbsFro1 ( k ) = norm ( data.appK - data.K , 'fro' ) ;
  errorRelOp1 ( k )  = errorAbsOp1 ( k ) / normOpK ; 
  errorRelFro1 ( k )  = errorAbsOp1 ( k ) / normFroK ; 
  absnnz1 ( k ) = nnz ( data.mask ) ;
  relnnz1 ( k ) = nnz ( data.mask )/ (data.n)^2 ;

  data.restInd = find ( ( abs (data.x(1 , : ) - 0.5 ) < 0.4) .* ( abs (data.x(2 , : ) - 0.5 ) < 0.4 ) ) ;
  errorRelOpRest1 ( k ) = norm ( data.appK ( data.restInd , data.restInd ) - data.K ( data.restInd , data.restInd ) ) / norm ( data.K ( data.restInd , data.restInd ) ) ;

end
input.tableCaption = ( 'Compression and Accuracy for $q = 7, \ l = 0.2 , \ \nu  = 1 , \ \delta_x = 0.2$ and different values of $\rho$' )
input.tableLabel = ( 'basicSetup1' ) ;
input.data =  horzcat ( rhoVec1 , errorAbsOp1 , errorRelOp1 , errorAbsFro1 , errorRelFro1 , absnnz1 , relnnz1 ) ; 
input.tableColLabels = {'$\rho$' , '$\| \KMC^{\rho} - \KMC\|$' ,  '$\| \KMC^{\rho} - \KMC\| / \| \KMC\|$' , '$\| \KMC^{\rho} - \KMC\|_{\FRO}$', '$\| \KMC^{\rho} - \KMC\|_{\FRO} / \| \KMC\|_{\FRO}$' , '$\# S$' , '$\#S / N^2$' } ;
textCellArray = latexTable ( input ) ;

delete ./figures/table_basicSetup1.tex ;
fid = fopen ( './figures/table_basicSetup1.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;

inputTrunc.tableCaption = ( 'Error induced by the incomplete factorisation for $q = 7, \ l = 0.2 , \ \nu  = 1 , \ \delta_x = 0.2$ and different values of $\rho$' )
inputTrunc.tableLabel = ( 'trunc1' ) ;
inputTrunc.data =  horzcat ( rhoVec1 , errorFroTruncCommutator1 , errorFroTrunc1 , errorFroTruncRel1) ; 
inputTrunc.tableColLabels = {'$\rho$' , '$\| L^{\rho} L^{\rho,T} - (L|_S) (L|_S)^T\|$' ,  '$\| (L|_S)(L|_S)^T - \KMC\|$' , '$\frac{\| L^{\rho} L^{\rho,T} - (L|_S) (L|_S)^T\|}{\| (L|_S)(L|_S)^T - \KMC\|}$' } ;
textCellArray = latexTable ( inputTrunc ) ;

delete ./figures/table_trunc1.tex ;
fid = fopen ( './figures/table_trunc1.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;












%experiment with $q = 7 , l = 0.5 , nu = 0.4 ,  rho = 2 : 2 : 10 , d_x = 0.2 , dz = 0
errorAbsOp2 = zeros ( 7 , 1 ) ;
errorRelOp2 = zeros ( 7 , 1 ) ;
errorAbsFro2 = zeros ( 7 , 1 ) ;
errorRelFro2 = zeros ( 7 , 1 ) ;
absnnz2 = zeros ( 7 , 1 ) ;
relnnz2 = zeros ( 7 , 1 ) ;
rhoVec2 = [2:2:14]' ;
errorRelOpRest2 = [7:1]' ;

errorFroTruncCommutator2 = zeros ( 7 ,1 ) ;
errorFroTrunc2 = zeros ( 7 ,1 ) ;
errorFroTruncRel2 = zeros ( 7 ,1 ) ;



for k = 1 : 7
  rho = rhoVec2 ( k ) ;
  data = matChol ( q , 0.4 , 1 , rho , 0.2 , 0 ) ;
  data.appK = data.Lichol * data.Lichol';
  data.L = chol ( data.K )' ;
  data.LTrunc = data.mask .* data.L ;
  data.appKTrunc = data.LTrunc * data.LTrunc' ;
  errorFroTruncCommutator2 ( k )  = norm ( data.appK - data.appKTrunc , 'fro' ) ;
  errorFroTrunc2 ( k )  = norm ( data.K - data.appKTrunc , 'fro' ) ;
  errorFroTruncRel2 ( k ) = errorFroTruncCommutator2 ( k ) / errorFroTrunc2 ( k ) ;
  normOpK = norm ( data.K ) ;
  normFroK = norm ( data.K , 'fro' ) ;
  errorAbsOp2 ( k )  = norm ( data.appK - data.K ) ;
  errorAbsFro2 ( k ) = norm ( data.appK - data.K , 'fro' ) ;
  errorRelOp2 ( k )  = errorAbsOp2 ( k ) / normOpK ; 
  errorRelFro2 ( k )  = errorAbsOp2 ( k ) / normFroK ; 
  absnnz2 ( k ) = nnz ( data.mask ) ;
  relnnz2 ( k ) = nnz ( data.mask )/ (data.n)^2 ;

  data.restInd = find ( ( abs (data.x(1 , : ) - 0.5 ) < 0.4) .* ( abs (data.x(2 , : ) - 0.5 ) < 0.4 ) ) ;
  errorRelOpRest2 ( k ) = norm ( data.appK ( data.restInd , data.restInd ) - data.K ( data.restInd , data.restInd ) ) / norm ( data.K ( data.restInd , data.restInd ) ) ;
end
input.tableCaption = ( 'Compression and Accuracy for $q = 7, \ l = 0.4 , \ \nu  = 1 , \ \delta_x = 0.2$ and different values of $\rho$' )
input.tableLabel = ( 'basicSetup2' ) ;
input.data =  horzcat ( rhoVec2 , errorAbsOp2 , errorRelOp2 , errorAbsFro2 , errorRelFro2 , absnnz2 , relnnz2 ) ; 
input.tableColLabels = {'$\rho$' , '$\| \KMC^{\rho} - \KMC\|$' ,  '$\| \KMC^{\rho} - \KMC\| / \| \KMC\|$' , '$\| \KMC^{\rho} - \KMC\|_{\FRO}$', '$\| \KMC^{\rho} - \KMC\|_{\FRO} / \| \KMC\|_{\FRO}$' , '$\# S$' , '$\#S / N^2$' } ;
textCellArray = latexTable ( input ) ;

delete ./figures/table_basicSetup2.tex ;
fid = fopen ( './figures/table_basicSetup2.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;

inputTrunc.tableCaption = ( 'Error induced by the incomplete factorisation for $q = 7, \ l = 0.4 , \ \nu  = 1 , \ \delta_x = 0.2$ and different values of $\rho$' )
inputTrunc.tableLabel = ( 'trunc2' ) ;
inputTrunc.data =  horzcat ( rhoVec2 , errorFroTruncCommutator2 , errorFroTrunc2 , errorFroTruncRel2) ; 
inputTrunc.tableColLabels = {'$\rho$' , '$\| L^{\rho} L^{\rho,T} - (L|_S) (L|_S)^T\|$' ,  '$\| (L|_S)(L|_S)^T - \KMC\|$' , '$\frac{\| L^{\rho} L^{\rho,T} - (L|_S) (L|_S)^T\|}{\| (L|_S)(L|_S)^T - \KMC\|}$' } ;
textCellArray = latexTable ( inputTrunc ) ;

delete ./figures/table_trunc2.tex ;
fid = fopen ( './figures/table_trunc2.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;





%experiment with $q = 7 , l = 0.5 , nu = 0.4 ,  rho = 2 : 2 : 10 , d_x = 0.2 , dz = 0
errorAbsOp3 = zeros ( 5 , 1 ) ;
errorRelOp3 = zeros ( 5 , 1 ) ;
errorAbsFro3 = zeros ( 5 , 1 ) ;
errorRelFro3 = zeros ( 5 , 1 ) ;
absnnz3 = zeros ( 5 , 1 ) ;
relnnz3 = zeros ( 5 , 1 ) ;
rhoVec3 = [6:2:14]' ;
errorRelOpRest3 = [5:1]' ;


errorFroTruncCommutator3 = zeros ( 5 ,1 ) ;
errorFroTrunc3 = zeros ( 5 ,1 ) ;
errorFroTruncRel3 = zeros ( 5 ,1 ) ;

for k = 1 : 5
  rho = rhoVec3 ( k ) ;
  data = matChol ( q , 0.2 , 2 , rho , 0.2 , 0 ) ;
  data.appK = data.Lichol * data.Lichol';
  data.L = chol ( data.K )' ;
  data.LTrunc = data.mask .* data.L ;
  data.appKTrunc = data.LTrunc * data.LTrunc' ;
  errorFroTruncCommutator3 ( k )  = norm ( data.appK - data.appKTrunc , 'fro' ) ;
  errorFroTrunc3 ( k )  = norm ( data.K - data.appKTrunc , 'fro' ) ;
  errorFroTruncRel3 ( k ) = errorFroTruncCommutator3 ( k ) / errorFroTrunc3 ( k ) ;
  normOpK = norm ( data.K ) ;
  normFroK = norm ( data.K , 'fro' ) ;
  errorAbsOp3 ( k )  = norm ( data.appK - data.K ) ;
  errorAbsFro3 ( k ) = norm ( data.appK - data.K , 'fro' ) ;
  errorRelOp3 ( k )  = errorAbsOp3 ( k ) / normOpK ; 
  errorRelFro3 ( k )  = errorAbsOp3 ( k ) / normFroK ; 
  absnnz3 ( k ) = nnz ( data.mask ) ;
  relnnz3 ( k ) = nnz ( data.mask )/ (data.n)^2 ;

  data.restInd = find ( ( abs (data.x(1 , : ) - 0.5 ) < 0.4) .* ( abs (data.x(2 , : ) - 0.5 ) < 0.4 ) ) ;
  errorRelOpRest3 ( k ) = norm ( data.appK ( data.restInd , data.restInd ) - data.K ( data.restInd , data.restInd ) ) / norm ( data.K ( data.restInd , data.restInd ) ) ;
end
input.tableCaption = ( 'Compression and Accuracy for $q = 7, \ l = 0.2 , \ \nu  = 2 , \ \delta_x = 0.2$ and different values of $\rho$' )
input.tableLabel = ( 'basicSetup3' ) ;
input.data =  horzcat ( rhoVec3 , errorAbsOp3 , errorRelOp3 , errorAbsFro3 , errorRelFro3 , absnnz3 , relnnz3 ) ; 
input.tableColLabels = {'$\rho$' , '$\| \KMC^{\rho} - \KMC\|$' ,  '$\| \KMC^{\rho} - \KMC\| / \| \KMC\|$' , '$\| \KMC^{\rho} - \KMC\|_{\FRO}$', '$\| \KMC^{\rho} - \KMC\|_{\FRO} / \| \KMC\|_{\FRO}$' , '$\# S$' , '$\#S / N^2$' } ;
textCellArray = latexTable ( input ) ;

delete ./figures/table_basicSetup3.tex ;
fid = fopen ( './figures/table_basicSetup3.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;

inputTrunc.tableCaption = ( 'Error induced by the incomplete factorisation for $q = 7, \ l = 0.2 , \ \nu  = 2 , \ \delta_x = 0.2$ and different values of $\rho$' )
inputTrunc.tableLabel = ( 'trunc3' ) ;
inputTrunc.data =  horzcat ( rhoVec3 , errorFroTruncCommutator3 , errorFroTrunc3 , errorFroTruncRel3) ; 
inputTrunc.tableColLabels = {'$\rho$' , '$\| L^{\rho} L^{\rho,T} - (L|_S) (L|_S)^T\|$' ,  '$\| (L|_S)(L|_S)^T - \KMC\|$' , '$\frac{\| L^{\rho} L^{\rho,T} - (L|_S) (L|_S)^T\|}{\| (L|_S)(L|_S)^T - \KMC\|}$' } ;
textCellArray = latexTable ( inputTrunc ) ;

delete ./figures/table_trunc3.tex ;
fid = fopen ( './figures/table_trunc3.tex' , 'wt' ) ;
for i = 1 : size ( textCellArray , 1 )
  fprintf(fid, '%s\t\n', textCellArray{i,:}); 
end
fclose ( fid ) ;





figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
plot ( rhoVec1 , log10 ( errorRelOp1 ) , '-r'  , rhoVec2 , log10 ( errorRelOp2 ) , '--b' ,  rhoVec3 , log10 ( errorRelOp3 ) , ':k' ) ;
legend (  '\nu = 1, l = 0.2'  , '\nu = 1, l = 0.4' , '\nu = 2, l = 0.2' ) ;
xlabel ( '\rho' ) ;
ylabel ( 'log_{10} ( E_{rel} )' ) ;
saveas ( figPlot , './figures/convergencePlotBasic' , 'jpg' ) ;
close all


figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
plot ( absnnz1 , log10 ( errorRelOp1 ) , '-r'  , absnnz2 , log10 ( errorRelOp2 ) , '--b' ,  absnnz3 , log10 ( errorRelOp3 ) , ':k' ) ;
legend (  '\nu = 1, l = 0.2'  , '\nu = 1, l = 0.4' , '\nu = 2, l = 0.2' ) ;
xlabel ( '# S' ) ;
ylabel ( 'log_{10} ( E_{rel} )' ) ;
saveas ( figPlot , './figures/convergencePlotVsSBasic' , 'jpg' ) ;
close all

figPlot = figure ( 'DefaultAxesFontsize' , 18 ) ;
plot ( rhoVec1 , log10 ( errorRelOpRest1 ) , '-r'  , rhoVec2 , log10 ( errorRelOpRest2 ) , '--b' ,  rhoVec3 , log10 ( errorRelOpRest3 ) , ':k' ) ;
legend (  '\nu = 1, l = 0.2'  , '\nu = 1, l = 0.4' , '\nu = 2, l = 0.2' ) ;
xlabel ( '\rho' ) ;
ylabel ( 'log_{10} ( E_{rel} )' ) ;
saveas ( figPlot , './figures/convergencePlotRest' , 'jpg' ) ;
close all








   
  

clear all
close all
