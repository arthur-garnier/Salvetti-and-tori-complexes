#! @Chapter Examples and Tests


LoadPackage( "Salvetti-and-tori-complexes" );

#! @Example
LoadPackage( "Salvetti-and-tori-complexes" );
#! true
#! First, we test the function SalvettiComplex on a simple example in type $A_3$, 
#! with the four first differentials. 
#! We check that there is no homology in positive degrees and that the zero
#! homology space is the trivial $\mathbb{Q}[W]$-module.
XX := "A";; n := 3;; k := 4;;
D := SalvettiComplex( XX, n, k );; 
W := D.coxeter_group;; QW := D.group_algebra;; D := D.differentials;;
B := Basis( QW, List( Elements( W ), w -> w * One( QW ) ) );;
D0 := List( D, d -> TO_INTEGER_MATRIX( W, 
ListWithIdenticalEntries( DimensionsMat( d )[1], TrivialSubgroup( W ) ), 
ListWithIdenticalEntries( DimensionsMat( d )[2], TrivialSubgroup( W ) ), d ) );;
F0 := List( D0, d -> LeftModuleHomomorphismByMatrix( 
Basis( Rationals ^ ( DimensionsMat( d )[1] ) ), d, 
Basis( Rationals ^ ( DimensionsMat( d )[2] ) ) ) );;
List( [ 1 .. Size( F0 ) - 1], i -> 
Dimension( Kernel( F0[i] ) / Image( F0[i+1] ) ) );
#! [ 0, 0, 0 ]
H1 := HOMOLOGY_REP( W, D0[1], [ TrivialSubgroup( W ) ], 
NullMat( Size( W ), 1, Rationals ) );;
List( ConjugacyClasses( W ), c -> TraceMat( Image( H1[2], 
Representative( c ) ) ) );
#! [ 1, 1, 1, 1, 1 ]
# @EndExample

#! @Example
#! Next, we compute the representations associated to the $W$-equivariant cellular 
#! homology chain complexes of the adjoint and simply-connected tori in type $A_2$ 
#! using the functions "CellularComplexT" and "CellularComplexTsc".
#! We check that the obtained representations are indeed the exterior powers of 
#! the geometric representation of $W$.
XX := "A";; n := 2;;
W := WeylGroup( RootSystem( SimpleLieAlgebra( XX, n, Rationals ) ) );;
Rad := CellularComplexT( XX, n, [ 1 .. n ] ).homology_representations;; 
Rsc := CellularComplexTsc( XX, n ).homology_representations;;
CHad := List( Rad, r -> List( ConjugacyClasses( W ), c -> 
TraceMat( Image( r[2], Representative( c ) ) ) ) );;
CHsc := List( Rsc, r -> List( ConjugacyClasses( W ), c -> 
TraceMat( Image( r[2], Representative( c ) ) ) ) );;
chi := First( Irr( W ), c -> ValuesOfClassFunction( c ) = CHad[2] );;
ForAll( [ 2 .. Size( CHad ) ], i -> ValuesOfClassFunction( 
AntiSymmetricParts( CharacterTable( W ), [ chi ], i - 1 )[1] ) = CHad[i]); 
#! true
CHad=CHsc;
#! true
# @EndExample

#! @Example
#! Now, we compute the representations associated to the complex of $I_2(5)$,
#! which is the dihedral group of order 10, using the same combinatorics as 
#! in the Weyl case, using "ComplexForFiniteCoxeterGroup".
#! We notice that the homology representation $H_2$ is no longer irreducible.
R := ComplexForFiniteCoxeterGroup( "I2", 5 );;
R.homology_characters[ 2 ];
#! [ 4, 0, -1, -1 ]
R.coeffs_in_tbl[ 2 ];
#! [ 0, 0, 1, 1 ]
#! @EndExample

#! @Example
#! Finally, we compute the characters of the homology representations of
#! the complex associated to the group $H_3$. Here again, the
#! representations are no longer irreducible.
R := ComplexForFiniteCoxeterGroup( "H", 3 );;
R.homology_characters;
#! [ [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ], [ 11, 3, -1, -1, 1, -1, 1, -1, -1, -1 ], 
#!   [ 11, -3, -1, -1, 1, 1, 1, 1, 1, 1 ], [ 1, -1, 1, 1, 1, -1, 1, -1, -1, -1 ] ]
R.coeffs_in_tbl;
#! [ [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 1, 0, 0, 1, 0 ], 
#!   [ 0, 0, 1, 1, 0, 0, 0, 0, 0, 1 ], [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ] ]
#! @EndExample

