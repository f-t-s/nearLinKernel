%This function creates the numerical experiment of the sparifying Cholesky factorisation on an artificially created grid
function data = matChol ( q , l , nu , rho , dx , dz ) 

%Create mesh of evaluation points
h = 1/2 ;
xAll = cell ( q , 1 ) ;
nAll = zeros ( q ,1 ) ;
for k = 1 : q
  xTemp1d = 0 : h^k : 1 ;
  nTemp = length ( xTemp1d ) ;
  nAll ( k ) = nTemp^2 ;
  xTemp2d = repmat ( xTemp1d , nTemp  , 1 ) ;
  xAll { k } = vertcat ( reshape ( xTemp2d , 1 , nTemp^2 ) , reshape ( xTemp2d' , 1 , nTemp^2 ) ) ; 
end
%lTracker tracks the level of each point
lTracker = zeros ( sum ( nAll , 1 ) , 1 ) ;
for k = 1 : q 
  lTracker ( sum ( nAll ( 1 : k - 1 ) , 1 ) + 1 : sum ( nAll ( 1 : k ) , 1 ) ) =  k * ones ( nAll ( k ) , 1 ) ;
end
x = cat ( 2 , xAll { : } ) ;
[ x , ind , ~ ] = unique ( x' , 'rows' , 'stable' ) ; 
x = x' ;
lTracker = lTracker ( ind , 1 ) ;
lTracker = ( 1 / h ).^lTracker * h ;
x = x + 2 * h^q * dx * ( rand ( 2 , length ( x ) ) - 0.5 * ones ( 2 , length ( x ) ) ) ;
x = vertcat ( x , - dz.* sin( 6 * x(1,:) ).*cos(2 * x(2 , : ) ) + 2 * 2 * dz * h^q * dx * ( rand ( 1 , length ( x ) ) - 0.5 * ones ( 1 , length ( x ) ) )) ;
data.x = x ;
data.xAll = xAll ;
data.n = length ( x ) ;

%Create Kernel Matrix
maternFunction = @(r) 2.^(1-nu)./gamma(nu)  .* ( sqrt(2 * nu ) .* ( r + eps )  ).^nu .* besselk(nu , sqrt( 2 * nu ) .* ( r + eps  ) ) .*  ( r ~= 0 ) + 1 * (r ==0 ) ;
data.kernelFunction = maternFunction;
data.distK =  squareform ( pdist ( x' )  ) ; 
K = maternFunction ( squareform ( pdist ( x' )  ) / l ) ;
K = K + ( K == 0 ).* eps ; 
data.K = K ;

%Create Mask for sparsification
minTracker = sparse( ( repmat (  lTracker , 1 , data.n )  + ( repmat (  lTracker , 1 , data.n ) )' - abs ( repmat(lTracker , 1 , data.n )  - ( repmat (  lTracker , 1 , data.n ) )')  )/2 );
data.mask = ( ( squareform ( pdist ( x' ) ) .* minTracker ) < rho ) ;
data.sparseK = sparse ( data.mask .* K ) ;

%computing sparse Cholesky decomposition
opts.type = 'nofill' ;
data.Lichol = ichol ( data.sparseK , opts ) ;
end



