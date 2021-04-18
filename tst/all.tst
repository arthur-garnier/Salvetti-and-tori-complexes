gap> XX := "A";; n := 3;; k := 3;;
gap> D := SalvettiComplex( XX, n, k );; 
gap> W := D.coxeter_group;; QW := D.group_algebra;; D := D.differentials;;
gap> B := Basis( QW, List( Elements( W ), w -> w * One( QW ) ) );;
gap> D0 := List( D, d -> TO_INTEGER_MATRIX( W, ListWithIdenticalEntries( DimensionsMat( d )[1], TrivialSubgroup( W ) ), ListWithIdenticalEntries( DimensionsMat( d )[2], TrivialSubgroup( W ) ), d ) );;
gap> F0 := List( D0, d -> LeftModuleHomomorphismByMatrix( Basis( Rationals ^ ( DimensionsMat( d )[1] ) ), d, Basis( Rationals ^ ( DimensionsMat( d )[2] ) ) ) );;
gap> List( [ 1 .. Size( F0 ) - 1], i -> Dimension( Kernel( F0[i] ) / Image( F0[i+1] ) ) );
[ 0, 0 ]
gap> H1 := HOMOLOGY_REP( W, D0[1], [ TrivialSubgroup( W ) ], NullMat( Size( W ), 1, Rationals ) );;
gap> List( ConjugacyClasses( W ), c -> TraceMat( Image( H1[2], Representative( c ) ) ) );
[ 1, 1, 1, 1, 1 ]

gap> XX := "A";; n := 2;;
gap> W := WeylGroup( RootSystem( SimpleLieAlgebra( XX, n, Rationals ) ) );;
gap> Rad := CellularComplexT( XX, n, [ 1 .. n ] ).homology_representations;; 
gap> Rsc := ComplexForFiniteCoxeterGroup( XX, n ).homology_representations;;
gap> CHad := List( Rad, r -> List( ConjugacyClasses( W ), c -> TraceMat( Image( r[2], Representative( c ) ) ) ) );;
gap> CHsc := List( Rsc, r -> List( ConjugacyClasses( W ), c -> TraceMat( Image( r[2], Representative( c ) ) ) ) );;
gap> chi := First( Irr( W ), c -> ValuesOfClassFunction( c ) = CHad[2] );;
gap> ForAll( [ 2 .. Size( CHad ) ], i -> ValuesOfClassFunction( AntiSymmetricParts( CharacterTable( W ), [ chi ], i - 1 )[1] ) = CHad[i]); 
true
gap> CHad=CHsc;
true

gap> R := ComplexForFiniteCoxeterGroup( "I2", 5 );;
gap> R.homology_characters[ 2 ];
[ 4, 0, -1, -1 ]
gap> R.coeffs_in_tbl[ 2 ];
[ 0, 0, 1, 1 ]

gap> R := ComplexForFiniteCoxeterGroup( "H", 3 );;
gap> R.homology_characters;
[ [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ], [ 11, 3, -1, -1, 1, -1, 1, -1, -1, -1 ], 
  [ 11, -3, -1, -1, 1, 1, 1, 1, 1, 1 ], [ 1, -1, 1, 1, 1, -1, 1, -1, -1, -1 ] ]
gap> R.coeffs_in_tbl;
[ [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 1, 0, 0, 1, 0 ], [ 0, 0, 1, 1, 0, 0, 0, 0, 0, 1 ], 
  [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ] ]
gap> STOP_TEST("freeintegralmodules02.tst", 1 );