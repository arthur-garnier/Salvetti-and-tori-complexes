#
# Salvetti-and-tori-complexes: Routine to build the Salvetti complex of a finite Coxeter group and the Weyl-equivariant cellular homology chain complex of the torus of an irreducible root datum.
#
# Implementations
#


InstallGlobalFunction( LEFT_COSETS,
		[ IsGroup, IsGroup ],
	function ( G, H )
		return Set( List( RightCosets( G, H ), r -> Set( List( r, x -> x ^ -1 ) ) ) );
end );



InstallGlobalFunction( LEFT_COSET,
		[ IsMultiplicativeElementWithInverse, IsGroup ],
	function ( g, H )
    	return Set( List( H, h -> g * h ) );
end );



InstallGlobalFunction( ON_LEFT_COSET,
		[ IsGroup, IsMultiplicativeElementWithInverse, IsListOrCollection ],
	function ( G, x, gH )
    	return LEFT_COSET( x * gH[1], Subgroup( G, gH[1] ^ -1 * gH ) );
end );



InstallGlobalFunction( MATRIX_TRANSITIVE_ACTION,
		[ IsGroup, IsGroup, IsMultiplicativeElementWithInverse ],
	function ( G, H, g )
    	return 
     	TransposedMat( PermutationMat( PermListList( LEFT_COSETS( G, H ), List( LEFT_COSETS( G, H ), function ( r )
        	          return ON_LEFT_COSET( G, g, r );
            	  end ) ), Index( G, H ) ) );
end );



InstallGlobalFunction( MATRIX_RIGHT_TRANSLATION,
		[ IsGroup, IsGroup, IsGroup, IsMultiplicativeElementWithInverse ],
	function ( G, K, H, g )
	    local M, LCK, LCH, R, p, L, gH;
	    M := [  ];
	    LCK := LEFT_COSETS( G, K );
 	   LCH := LEFT_COSETS( G, H );
	    gH := LEFT_COSET( g, H );
    	if IsSubgroup( H, K ^ g ) then
        	for R in LCK do
            	p := Position( LCH, LEFT_COSET( R[1] * g, H ) );
            	L := ListWithIdenticalEntries( Index( G, H ), 0 );
            	L[p] := 1;
            	Add( M, L );
        	od;
    	else
        	M := fail;
    	fi;
    	return M;
end );



InstallGlobalFunction( MATRIX_RIGHT_TRANSLATION_ZG, 
		[ IsGroup, IsGroup, IsGroup, IsMultiplicativeElement ],
	function ( G, K, H, x )
	    local X0, i, C;
		C := CoefficientsAndMagmaElements( x );
    	X0 := NullMat( Index( G, K ), Index( G, H ), Rationals );
		for i in [ 1 .. Size( C ) / 2 ] do
			X0 := X0 + C[ 2 * i ] * MATRIX_RIGHT_TRANSLATION( G, K, H, C[ 2 * i - 1 ] );
		od;
    	return X0;
end );



InstallGlobalFunction( TO_INTEGER_MATRIX,
		[ IsGroup, IsList, IsList, IsMatrix ],
	function ( G, LH, LK, d )
    	local D, SH, SK, h, k, H, K;
    	D := NullMat( Sum( List( LH, function ( H )
        	        return Index( G, H );
            	end ) ), Sum( List( LK, function ( K )
                	return Index( G, K );
            	end ) ), Rationals );
    	SH := 0;
    	for h in [ 1 .. Size( LH ) ] do
        	H := LH[h];
        	SK := 0;
        	for k in [ 1 .. Size( LK ) ] do
            	K := LK[k];
            	D{[ SH + 1 .. SH + Index( G, H ) ]}{[ SK + 1 .. SK + Index( G, K ) ]} := MATRIX_RIGHT_TRANSLATION_ZG( G, H, K, d[h, k] );
            	SK := SK + Index( G, K );
      	od;
      	SH := SH + Index( G, H );
  	od;
  	return D;
end );



InstallGlobalFunction( HOMOLOGY_REP,
		[ IsGroup, IsMatrix, IsList, IsMatrix ],
	function ( G, D, LH, Dp )
  	local QQ, Dm, Dpm, V, pi, B0, MM, M, g, b, b0, H, h, bH;
  	QQ := Rationals;
  	Dm := LeftModuleHomomorphismByMatrix( Basis( QQ ^ DimensionsMat( D )[1] ), D, Basis( QQ ^ DimensionsMat( D )[2] ) );
  	Dpm := LeftModuleHomomorphismByMatrix( Basis( QQ ^ DimensionsMat( Dp )[1] ), Dp, Basis( QQ ^ DimensionsMat( Dp )[2] ) );
  	V := Kernel( Dpm ) / Image( Dm );
  	pi := NaturalHomomorphismBySubspace( Kernel( Dpm ), Image( Dm ) );
  	if Dimension( V ) = 0 then
    	  MM := ListWithIdenticalEntries( Size( GeneratorsOfGroup( G ) ), [ [ 1 ] ] );
  	else
    	  B0 := Basis( Image( pi ) );
      	MM := [  ];
      	for g in GeneratorsOfGroup( G ) do
        	  M := [  ];
          	for b in B0 do
            	  b := ShallowCopy( PreImagesRepresentative( pi, b ) );
            	  b0 := [  ];
            	  for H in LH do
                	  bH := [  ];
                	  for h in [ 1 .. Index( G, H ) ] do
                    	  Add( bH, b[1] );
                    	  Remove( b, 1 );
                  	od;
                  	b0 := Concatenation( b0, TransposedMat( MATRIX_TRANSITIVE_ACTION( G, H, g ) * TransposedMat( [ bH ] ) )[1] );
              	od;
              	Add( M, Coefficients( B0, Image( pi, b0 ) ) );
          	od;
          	Add( MM, TransposedMat( M ) );
      	od;
  	fi;
  	return [ Dimension( V ), GroupHomomorphismByImages( G, Group( MM ), GeneratorsOfGroup( G ), MM ) ];
end );



InstallMethod( REFLECTIONS,
		[ IsGroup ],
	function ( W )
  	local R, s;
  	R := [  ];
  	for s in GeneratorsOfGroup( W ) do
    	  R := Concatenation( R, List( ConjugateSubgroups( W, Subgroup( W, [ s ] ) ), function ( H )
        	        return First( H, function ( h )
            	            return h <> One( H );
                	    end );
            	end ) );
  	od;
  	return Set( R );
end );



InstallMethod( LONGEST_ELEMENT,
		[ IsGroup ],
	function ( W )
  	return First( W, function ( w )
    	      return Length( Factorization( W, w ) ) = Size( PositiveRoots( RootSystem( W ) ) );
      	end );
end );



InstallMethod( PROJ_OMEGA,
		[ IsGroup, IsList ],
	function ( W, Jinp )
  	local J, Omp, Om0, str, Ws, n, w0, j, Wj, wj, Om, a0, S, R;
 	 R := RootSystem( W );
  	S := SimpleSystem( R );
  	a0 := 0;
  	str := SemiSimpleType( UnderlyingLieAlgebra( R ) );
  	n := IntChar( str[2] ) - 48;
  	Ws := GeneratorsOfGroup( W );
  	Om := [ One( W ) ];
  	if str[1] = 'A' then
    	  w0 := LONGEST_ELEMENT( W );
      	Om0 := [ One( W ) ];
      	for j in [ 1 .. n ] do
        	  Wj := ShallowCopy( Ws );
        	  Remove( Wj, j );
        	  wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = j * (j - 1) / 2 + (n - j + 1) * (n - j) / 2;
            	  end );
          	Add( Om0, wj * w0 );
      	od;
      	a0 := Sum( S );
      	Omp := Elements( Subgroup( Group( Om0 ), Filtered( Om0, function ( o )
        	          return Position( Om0, o ) - 1 in Jinp;
            	  end ) ) );
      	J := List( Filtered( [ 2 .. n + 1 ], function ( i )
        	        return Om0[i] in Omp;
            	end ), function ( j )
             	 return j - 1;
          	end );
      	Om := Concatenation( [ One( W ) ], List( J, function ( j )
       	         return Om0[j + 1];
        	    end ) );
  	elif str[1] = 'B' then
    	  if Jinp = [  ] then
        	  J := [  ];
      	else
        	  J := [ 1 ];
         	 w0 := LONGEST_ELEMENT( W );
         	 Wj := ShallowCopy( Ws );
          	Remove( Wj, 1 );
         	 wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = (n - 1) ^ 2;
            	  end );
          	Add( Om, wj * w0 );
      	fi;
      	a0 := 2 * Sum( S ) - S[1];
  	elif str[1] = 'C' then
    	  if Jinp = [  ] then
        	  J := [  ];
      	else
        	  J := [ n ];
        	  w0 := LONGEST_ELEMENT( W );
        	  Wj := ShallowCopy( Ws );
        	  Remove( Wj, n );
        	  wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = n * (n - 1) / 2;
              	end );
         	 Add( Om, wj * w0 );
      	fi;
      	a0 := 2 * Sum( S ) - S[n];
  	elif str[1] = 'D' then
    	  if Size( Jinp ) > 1 then
        	  J := [ 1, n - 1, n ];
      	elif Size( Jinp ) = 0 then
        	  J := [  ];
      	else
        	  if ( IsOddInt( n ) and (Jinp = [ n - 1 ] or Jinp = [ n ]) ) then
            	  J := [ 1, n - 1, n ];
          	else
            	  J := Jinp;
          	fi;
      	fi;
      	if J <> [  ] then
        	  w0 := LONGEST_ELEMENT( W );
      	fi;
      	if 1 in J then
        	  Wj := ShallowCopy( Ws );
         	 Remove( Wj, 1 );
         	 wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = (n - 1) ^ 2 - n + 1;
            	  end );
         	 Add( Om, wj * w0 );
      	fi;
      	if n - 1 in J then
        	  Wj := ShallowCopy( Ws );
        	  Remove( Wj, n - 1 );
        	  wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = n * (n - 1) / 2;
            	  end );
          	Add( Om, wj * w0 );
      	fi;
      	if n in J then
        	  Wj := ShallowCopy( Ws );
          	Remove( Wj, n );
          	wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = n * (n - 1) / 2;
              	end );
          	Add( Om, wj * w0 );
      	fi;
      	a0 := 2 * Sum( S ) - S[1] - S[(n - 1)] - S[n];
  	elif str[1] = 'E' and n = 6 then
    	  if Jinp = [  ] then
        	  J := [  ];
      	else
        	  J := [ 1, 6 ];
        	  w0 := LONGEST_ELEMENT( W );
      	fi;
      	if 1 in J then
        	  Wj := ShallowCopy( Ws );
          	Remove( Wj, 1 );
          	wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = 20;
              	end );
          	Add( Om, wj * w0 );
      	fi;
      	if 6 in J then
        	  Wj := ShallowCopy( Ws );
         	 Remove( Wj, 6 );
         	 wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = 20;
             	 end );
          	Add( Om, wj * w0 );
      	fi;
      	a0 := S[1] + 2 * S[2] + 2 * S[3] + 3 * S[4] + 2 * S[5] + S[6];
  	elif str[1] = 'E' and n = 7 then
    	  if Jinp = [  ] then
        	  J := [  ];
      	else
        	  J := [ 7 ];
          	w0 := LONGEST_ELEMENT( W );
      	fi;
      	if 7 in J then
        	  Wj := ShallowCopy( Ws );
          	Remove( Wj, 7 );
          	wj := First( Subgroup( W, Wj ), function ( w )
            	      return Length( Factorization( W, w ) ) = 36;
              	end );
          	Add( Om, wj * w0 );
      	fi;
      	a0 := 2 * S[1] + 2 * S[2] + 3 * S[3] + 4 * S[4] + 3 * S[5] + 2 * S[6] + S[7];
  	elif str[1] = 'E' and n = 8 then
    	  J := [  ];
      	a0 := 2 * S[1] + 3 * S[2] + 4 * S[3] + 6 * S[4] + 5 * S[5] + 4 * S[6] + 3 * S[7] + 2 * S[8];
  	elif str[1] = 'F' then
    	  J := [  ];
      	a0 := 2 * S[2] + 3 * S[4] + 4 * S[3] + 2 * S[1];
  	elif str[1] = 'G' then
    	  J := [  ];
      	a0 := 3 * S[1] + 2 * S[2];
  	fi;
  	return [ Om, Concatenation( [ First( REFLECTIONS( W ), function ( w )
    	                return a0 * w = - a0;
        	        end ) ], Ws ), J ];
end );



InstallMethod( LATTICE_FROM_ROOT_SYSTEM,
		[ IsString, IsPosInt ],
	function ( XX, n )
  	local I, S, x;
  	I := [ 0 .. n ];
  	S := [  ];
  	for x in IteratorOfCombinations( [ 1 .. Size( I ) ] ) do
    	  Add( S, I{x} );
  	od;
  	Remove( S, 1 );
  	return S;
end );



InstallMethod( PERMS_OF_OMEGA,
		[ IsString, IsPosInt, IsList],
	function ( XX, n, Jp )
  	local J, Om, per, perm1, perm, i, k;
	Om := [ () ]; J := [ ];
  	if XX = "A" then
    	  Om := Elements( Group( List( Union( [ 0 ], Intersection( Jp, [ 1 .. n ] ) ), function ( k )
        	      return CycleFromList( [ 1 .. (n + 1) ] ) ^ k;
          	end ) ) );
		  J := List( Om, function ( i ) return Position( List( [ 0 .. n ], function ( k ) return CycleFromList( [ 1 .. (n + 1) ] ) ^ k; end ), i ) - 1; end);
		  Remove( J, 1 );
  	elif XX = "B" then
    	if Jp <> [  ] then
       	   Om := [ (), (1,2) ];
		   J := [ 1 ];
      	fi;
  	elif XX = "C" then
    	if Jp <> [  ] then
       	   per := ();
        	  for i in [ 1 .. Int( (n + 2) / 2 ) ] do
            	  if i <> n + 2 - i then
                	  per := per * (i,(n + 2 - i));
              	fi;
          	od;
          	Om := [ (), per ];
			J := [ n ];
      	fi;
  	elif XX = "D" then
    	  perm1 := ();
      	perm := ();
      	if IsEvenInt( n ) then
        	  for i in [ 2 .. n / 2 - 1 ] do
            	  if i <> n - i then
                	  perm1 := perm1 * ((i + 1),(n - i + 1));
                  	perm := perm * ((i + 1),(n - i + 1));
              	fi;
          	od;
          	perm1 := perm1 * (1,n) * (2,(n + 1));
          	perm := perm * (1,(n + 1)) * (n,2);
      	else
        	  for i in [ 2 .. (n - 1) / 2 ] do
            	  if i <> n - i then
                	  perm1 := perm1 * ((i + 1),(n - i + 1));
                 	 perm := perm * ((i + 1),(n - i + 1));
              	fi;
          	od;
          	perm1 := perm1 * (n,2,(n + 1),1);
          	perm := perm * ((n + 1),2,n,1);
      	fi;
      	if 1 in Jp then
        	  Add( Om, (1,2) * (n,(n + 1)) );
		fi;
      	if n - 1 in Jp then
        	  Add( Om, perm1 );
		fi;
      	if n in Jp then
        	  Add( Om, perm );
		fi;
      	Om := Elements( Group( Om ) );
		if (1,2) * (n,(n + 1)) in Om then
			Add( J, 1 );
		fi;
		if perm1 in Om then
			Add( J, n - 1 );
		fi;
		if perm in Om then
			Add(J, n );
		fi;
  	elif XX = "E" and n = 6 then
    	  if Jp <> [  ] then
        	  Om := [ (), (1,2,7) * (3,4,6), (1,7,2) * (4,3,6) ];
			  J := [ 1, 6 ];
      	fi;
  	elif XX = "E" and n = 7 then
    	  if Jp <> [  ] then
        	  Om := [ (), (1,8) * (2,7) * (4,6) ];
			  J := [ 7 ];
      	fi;
  	elif XX = "E" and n = 8 or XX = "F" and n = 4 or XX = "G" and n = 2 then
    	  Om := [ () ];
  	else
    	  Om := fail;
  	fi;
  	return [ Om, J ];
end );



InstallGlobalFunction( PERM_ON_LATTICE,
		[ IsPerm, IsListOrCollection ],
	function ( s, M )
  	local s1;
  	s1 := List( OnTuples( List( M, function ( m )
    	          return m + 1;
       	   end ), s ), function ( t )
       	   return t - 1;
     	 end );
  	Sort( s1 );
  	return s1;
end );



InstallGlobalFunction( PREORBS_OF_LATTICE,
		[ IsString, IsPosInt, IsList ],
	function ( XX, n, Jp )
  	local S, Om, Orbs, s, t, Oo, o;
  	S := ShallowCopy( LATTICE_FROM_ROOT_SYSTEM( XX, n ) );
  	Om := PERMS_OF_OMEGA( XX, n, Jp )[1];
  	Orbs := [  ];
  	while Size( S ) > 0 do
    	  t := S[1];
      	Remove( S, 1 );
      	Oo := [ [ 1, t ] ];
      	for o in Om{[ 2 .. Size( Om ) ]} do
        	s := PERM_ON_LATTICE( o, t );
          	Add( Oo, [ Position( Om, o ), s ] );
            if s in S then
                Remove( S, Position( S, s ) );
            fi;
      	od;
      	Add( Orbs, Oo );
  	od;
  	return Orbs;
end );



InstallGlobalFunction( STABS_IN_LATTICE,
		[ IsString, IsPosInt, IsList ],
	function ( XX, n, Jp )
  	local Om, Orbs, Stabs, ore, o1, sta, o;
  	Om := PERMS_OF_OMEGA( XX, n, Jp )[1];
  	Orbs := PREORBS_OF_LATTICE( XX, n, Jp );
  	Stabs := [  ];
  	for ore in Orbs do
  	    o1 := ore[1][2];
  	    sta := [  ];
  	    for o in Om do
  	        if PERM_ON_LATTICE( o, o1 ) = o1 then
  	            Add( sta, Position( Om, o ) );
  	        fi;
  	    od;
  	    Add( Stabs, [ ore, sta ] );
  	od;
  	return Stabs;
end );



InstallGlobalFunction( CHAINS_OF_POSET,
		[ IsList, IsFunction ],
	function ( L, ord )
  	local Ch, x, L0, s, l;
  	if Size( L ) = 1 then
  	    Ch := [ [ L[1] ] ];
 	 else
 	     Ch := [  ];
 	     for x in L do
 	         if ForAll( L, function ( s )
 	                   return ord( s, x ) = false or x = s;
 	               end ) then
 	             L0 := CHAINS_OF_POSET( Filtered( L, function ( s )
 	                       return ord( x, s ) and x <> s;
 	                   end ), ord );
 	             for l in L0 do
 	                 if ( ord( x, l[1] ) and x <> l[1] ) then
 	                     Add( L0, Concatenation( [ x ], l ) );
 	                 fi;
 	             od;
 	             Add( L0, [ x ] );
 	         fi;
 	         Ch := Concatenation( Ch, L0 );
 	     od;
  	fi;
  	return DuplicateFreeList( Ch );
end );



InstallGlobalFunction( CHAINS_BY_LENGTH,
		[ IsList, IsFunction ],
	function ( L, ord )
  	local CCk, Ck, Ch, k, c;
  	CCk := [  ];
  	Ch := CHAINS_OF_POSET( L, ord );
  	for k in [ 1 .. Maximum( List( Ch, Size ) ) ] do
  	    Ck := [  ];
  	    for c in Ch do
  	        if Size( c ) = k then
  	            Add( Ck, c );
  	        fi;
  	    od;
  	    Add( CCk, [ Ck, k ] );
  	od;
  	return CCk;
end );



InstallGlobalFunction( PERM_ACTS_ON_CHAIN,
		[ IsPerm, IsList ],
	function ( s, C )
  	local U;
  	U := List( C, function ( c )
  	        return PERM_ON_LATTICE( s, c );
  	    end );
  	return U;
end );



InstallGlobalFunction( PREORBS_OF_CHAINS,
		[ IsString, IsPosInt, IsList ],
	function ( XX, n, Jp )
  	local C, Om, Orbs, k, c, Ok, T, t, s, Oo, o;
  	C := CHAINS_BY_LENGTH( LATTICE_FROM_ROOT_SYSTEM( XX, n ), function( a, b ) return IsSubset( b, a ); end );
  	Om := PERMS_OF_OMEGA( XX, n, Jp )[1];
  	Orbs := [ [ PREORBS_OF_LATTICE( XX, n, Jp ), 1 ] ];
  	for k in [ 2 .. Maximum( List( C, function ( c )
  	                return c[2];
  	            end ) ) ] do
  	    Ok := [  ];
  	    T := C[k][1];
  	    while Size( T ) > 0 do
  	        t := T[1];
  	        Remove( T, 1 );
  	        Oo := [ [ 1, t ] ];
  	        for o in Om{[ 2 .. Size( Om ) ]} do
  	            s := PERM_ACTS_ON_CHAIN( o, t );
  	            if s in C[k][1] then
  	                Add( Oo, [ Position( Om, o ), s ] );
  	                if s in T then
  	                    Remove( T, Position( T, s ) );
  	                fi;
  	            fi;
  	        od;
  	        Add( Ok, Oo );
  	    od;
  	    Add( Orbs, [ Ok, k ] );
  	od;
  	return Orbs;
end );



InstallGlobalFunction( STABS_OF_CHAINS,
		[ IsString, IsPosInt, IsList ],
	function ( XX, n, Jp )
  	local Om, Stabs, P, k, Stk, ore, sta, o;
  	Om := PERMS_OF_OMEGA( XX, n, Jp )[1];
  	Stabs := List( STABS_IN_LATTICE( XX, n, Jp ), function ( k )
  	        return [ List( k[1], function ( o )
  	                    return [ o[1], [ o[2] ] ];
  	                end ), k[2] ];
  	    end );
  	Stabs := [ [ Stabs, 1 ] ];
  	P := PREORBS_OF_CHAINS( XX, n, Jp );
  	for k in [ 2 .. Size( P ) ] do
  	    Stk := [  ];
  	    for ore in P[k][1] do
  	        sta := [  ];
  	        for o in Om do
  	            if PERM_ACTS_ON_CHAIN( o, ore[1][2] ) = ore[1][2] then
  	                Add( sta, Position( Om, o ) );
  	            fi;
  	        od;
  	        Add( Stk, [ ore, sta ] );
  	    od;
  	    Add( Stabs, [ Stk, k ] );
  	od;
  	return Stabs;
end );



InstallGlobalFunction( LIST_OF_STABS,
		[ IsString, IsPosInt, IsList ],
	function ( XX, n, Jp )
  	local St, SS, k, Sk, q, x;
  	St := STABS_OF_CHAINS( XX, n, Jp );
  	SS := [  ];
  	for k in [ 1 .. Size( St ) ] do
  	    Sk := [  ];
  	    for q in St[k][1] do
  	        Add( Sk, [ Minimum( List( q[1], function ( x )
  	                      return x[2];
  	                  end ) ), q[2] ] );
  	    od;
  	    q := ShallowCopy( Sk );
  	    Sort( q );
  	    Add( SS, [ q, k ] );
  	od;
  	return SS;
end );



InstallGlobalFunction( BOUNDARY_SIMPLEX,
		[ IsList ],
	function ( L )
  	local L0, i, Bou;
  	Bou := [  ];
  	for i in [ 0 .. Size( L ) - 1 ] do
  	    L0 := ShallowCopy( [ 1 .. Size( L ) ] );
  	    Remove( L0, i + 1 );
  	    Add( Bou, [ (-1) ^ i, L{L0} ] );
  	od;
  	return Bou;
end );



InstallGlobalFunction( SEARCH_FOR_ORBIT,
		[ IsList, IsList, IsList ],
	function ( C, Ls, Om )
  	local L, l0, o0, o, l;
  	L := Ls[Size( C )][1];
  	l0 := First( L, function ( l )
  	        return ForAny( Om, function ( o )
  	                return PERM_ACTS_ON_CHAIN( o, C ) = l[1];
  	            end );
  	    end );
  	o0 := First( Om, function ( o )
  	        return PERM_ACTS_ON_CHAIN( o, C ) = l0[1];
  	    end );
  	return [ Position( L, l0 ), o0 ];
end );



InstallGlobalFunction( ComplexForOmega,
#		[ IsString, IsPosInt, IsList, IsList ],
	function ( arg )
  	local XX, n, J, degs, Ls, Om, QO, Hs, L, l, x, i, d, dl, b, D, O, H0;
	XX := arg[1];; n := arg[2];; J := arg[3];; 
	if Length( arg ) < 4 then
		degs := [ 0 .. n ];;
		else
		degs := arg[4];;
	fi;
  	Om := PERMS_OF_OMEGA( XX, n, J );
	J := Om[2]; Om := Om[1];
  	O := Group( Om );
  	QO := GroupRing( Rationals, O );
  	Ls := LIST_OF_STABS( XX, n, J );
  	D := [ NullMat( Size( Ls[1][1] ), 1, QO ) ];
  	Hs := List( Ls, function ( L )
  	        return List( L[1], function ( l )
  	                return List( l[2], function ( x )
  	                        return Om[x];
  	                    end );
  	            end );
  	    end );
  	for i in [ 2 .. Size( Ls ) ] do
  	    d := [  ];
  	    for l in Ls[i][1] do
  	        dl := ListWithIdenticalEntries( Size( Ls[i - 1][1] ), Zero( QO ) );
  	        for b in BOUNDARY_SIMPLEX( l[1] ) do
  	            x := SEARCH_FOR_ORBIT( b[2], Ls, Om );
  	            dl[x[1]] := dl[x[1]] + b[1] * One( QO ) * ( x[2] ^ -1 );
  	        od;
  	        Add( d, dl );
  	    od;
  	    Add( D, d );
  	od;
  	Add( D, NullMat( 1, Size( Ls[Size( Ls )][1] ), QO ) );
  	H0 := [  ];
  	for l in Hs do
  	    Add( H0, List( l, function ( x )
  	            return Subgroup( O, x );
  	        end ) );
  	od;
  	return rec( fundamental_group_permutations := Om, group_algebra := QO, stabilizers := H0, differentials := D{[ 2 .. n + 1 ]} );
end );



InstallGlobalFunction( CellularComplexT,
#		[ IsString, IsPosInt, IsList, IsList ],
	function ( arg )
  	local cm, XX, n, J, degs, Jp, W, Wf, Om, QW, QWf, O, Ws, Ls, D, Hs, L, l, x, H0, i, j, k, d, dl, b, D0, Reps, q, CH;
	XX := arg[1];; n := arg[2];; J := arg[3];; 
	if Length( arg ) < 4 then
		degs := [ 0 .. n ];;
		else
		degs := arg[4];;
	fi;
  	Wf := WeylGroup( RootSystem( SimpleLieAlgebra( XX, n, Rationals ) ) );
  	O := PROJ_OMEGA( Wf, J );
  	Ws := O[2];
  	Jp := O[3];
  	O := O[1];
  	Om := PERMS_OF_OMEGA( XX, n, Jp )[1];
	q := IsomorphismPermGroup( Wf ); Ws := List( Ws, function( i ) return Image( q, i ); end ); O := List( O, function( i ) return Image( q, i ); end ); W := Image( q );
  	QW := GroupRing( Rationals, W );
  	Ls := LIST_OF_STABS( XX, n, Jp );
  	D := [ NullMat( Size( Ls[1][1] ), 1, QW ) ];
  	Hs := List( Ls, function ( L )
  	        return List( L[1], function ( l )
  	                return Concatenation( Ws{List( Difference( [ 0 .. n ], l[1][Size( l[1] )] ), function ( i )
  	                           return i + 1;
  	                       end )}, List( l[2], function ( x )
  	                          return O[x];
  	                      end ) );
  	            end );
  	    end );
  	H0 := [  ];
  	for l in Hs do
  	    Add( H0, List( l, function ( x )
  	            return Subgroup( W, x );
  	        end ) );
  	od;
  	for i in [ 2 .. Size( Ls ) ] do
  	    d := [  ];
  	    for l in Ls[i][1] do
  	        dl := ListWithIdenticalEntries( Size( Ls[i - 1][1] ), Zero( QW ) );
  	        for b in BOUNDARY_SIMPLEX( l[1] ) do
  	            x := SEARCH_FOR_ORBIT( b[2], Ls, Om );
  	            dl[x[1]] := dl[x[1]] + b[1] * One( QW ) * O[Position( Om, x[2] ^ -1 )];
  	        od;
  	        Add( d, dl );
  	    od;
  	    Add( D, d );
  	od;
  	Add( D, NullMat( 1, Size( Ls[Size( Ls )][1] ), QW ) );
  	D0 := [ TO_INTEGER_MATRIX( W, H0[1], [ W ], D[1] ) ];
  	D0 := Concatenation( D0, List( [ 2 .. n + 1 ], function ( i )
            if ( i - 1 in degs ) or ( i - 2 in degs ) then
  	            return TO_INTEGER_MATRIX( W, H0[i], H0[i - 1], D[i] );
            else
                return TO_INTEGER_MATRIX( W, H0[i], H0[i - 1], Zero( QW ) * D[i] );
            fi;
  	    end ) );
  	Add( D0, TO_INTEGER_MATRIX( W, [ W ], H0[Size( Hs )], D[Size( D )] ) );
  	Reps := List( degs, function ( i )
  	        return HOMOLOGY_REP( W, D0[i + 2], H0[i + 1], D0[i + 1] );
  	    end );
	QWf := GroupRing( Rationals, Wf );
	D := List( D, function( d )
		k := ShallowCopy( NullMat( DimensionsMat( d )[1], DimensionsMat( d )[2], QWf ) );
		for i in [ 1 .. DimensionsMat( d )[1] ] do
			for j in [ 1.. DimensionsMat( d )[2] ] do
				if d[ i, j ] = Zero( QW ) then
					k[ i, j ] := Zero( QWf );
					else
					cm := CoefficientsAndMagmaElements( d[ i, j ] );
					k[ i, j ] := One( QWf ) * cm[ 2 ] * PreImagesRepresentative(q, cm[ 1 ]);
				fi;
			od; 
		od;
		return k;
		end );
	H0 := List( H0, function( d ) 
		return List( d, function( i ) 
				return Subgroup( Wf, List( GeneratorsOfGroup( i ), function( j ) 
						return PreImagesRepresentative( q, j ); 
						end ) );
				end ); 
		end );
	Reps := List( Reps, function( i )
		return [ i[1], GroupHomomorphismByImages( Wf, Group( List( GeneratorsOfGroup( W ), function( k ) 
				return Image( i[2], k ); end ) ), GeneratorsOfGroup( Wf ), List( GeneratorsOfGroup( Wf ), function( d ) 
						return Image( i[2], Image( q, d ) ); end ) ) ];
			end );
    CH := List( Reps, function ( r )
            return List( ConjugacyClasses( Wf ), function ( c )
                    return TraceMat( Image( r[2], Representative( c ) ) );
                end );
        end );
  	return rec( homology_representations := Reps, homology_characters := CH, weyl_group := Wf, group_algebra := QWf, stabilizers := H0, differentials := D{[ 2 .. n + 1 ]} );
end );




InstallGlobalFunction( ComplexForFiniteCoxeterGroup,
function ( arg )
    local XX, n, degs, Wf, W, s0, s1, s2, s3, s4, Ws, R, r, Subs, QW, QWf, x, Hs, k, S, D, H0, i, d, L, dl, I0L, D0, Reps, CH, Bch, q, j, a0;
    XX := arg[1];
    n := arg[2];
    if XX = "I2" then
        W := FreeGroup( [ "s1", "s2" ] );
        s1 := W.1; s2 := W.2;
        Wf := W / [ s1 ^ 2, s2 ^ 2, (s1 * s2) ^ n ];
        s1 := Wf.1; s2 := Wf.2;
        R := REFLECTIONS( Wf );
		s0 := ( (s1 * s2) ^ Int( (n - 1) / 2 ) ) * s1;
        Ws := ShallowCopy( [ s0, s1, s2 ] );
        r := 2;
    elif XX = "H" and n = 3 then
        W := FreeGroup( [ "s1", "s2", "s3" ] );
        s1 := W.1; s2 := W.2; s3 := W.3;
        Wf := W / [ s1 ^ 2, s2 ^ 2, s3 ^ 2, (s1 * s2) ^ 5, (s2 * s3) ^ 3, (s1 * s3) ^ 2 ];
        s1 := Wf.1; s2 := Wf.2; s3 := Wf.3;
        #s0 := (s2 * s1) ^ 2 * s3 * s2 * s1 * s2 * s3 * (s1 * s2) ^ 2;
		s0 := s3 ^ ( (s2 * s1) ^ 2 );
        Ws := ShallowCopy( [ s0, s1, s2, s3 ] );
        r := 3;
    elif XX = "H" and n = 4 then
        W := FreeGroup( [ "s1", "s2", "s3", "s4" ] );
        s1 := W.1; s2 := W.2; s3 := W.3; s4 := W.4;
        Wf := W / [ s1 ^ 2, s2 ^ 2, s3 ^ 2, s4 ^ 2, (s1 * s2) ^ 5, (s2 * s3) ^ 3, (s3 * s4) ^ 3, (s1 * s3) ^ 2, (s1 * s4) ^ 2, (s2 * s4) ^ 2 ];
        s1 := Wf.1; s2 := Wf.2; s3 := Wf.3; s4 := Wf.4;
        #s0 := (s4 * s3 * (s2 * s1) ^ 2 * s3 * s2 * s1 * s2 * s3) ^ 2 * s4 * (s3 * s2 * s1 * s2 * s3 * (s1 * s2) ^ 2 * s3 * s4) ^ 2;
		s0 := s4 ^ ( (s3 * s2 * s1 * s2 * s3 * (s1 * s2) ^ 2 * s3 * s4) ^ 2 );
        Ws := ShallowCopy( [ s0, s1, s2, s3, s4 ] );
        r := 4;
    else
        Wf := WeylGroup( RootSystem( SimpleLieAlgebra( XX, n, Rationals ) ) );
		S := SimpleSystem( RootSystem( Wf ) );
		if XX = "A" then
			a0 := Sum( S );
		elif XX = "B" then
			a0 := 2 * Sum( S ) - S[1];
		elif XX = "C" then
			a0 := 2 * Sum( S ) - S[n];
		elif XX = "D" then
			a0 := 2 * Sum( S ) - S[1] - S[(n - 1)] - S[n];
		elif XX = "E" and n = 6 then
			a0 := S[1] + 2 * S[2] + 2 * S[3] + 3 * S[4] + 2 * S[5] + S[6];
		elif XX = "E" and n = 7 then
			a0 := 2 * S[1] + 2 * S[2] + 3 * S[3] + 4 * S[4] + 3 * S[5] + 2 * S[6] + S[7];
		elif XX = "E" and n = 8 then
			a0 := 2 * S[1] + 3 * S[2] + 4 * S[3] + 6 * S[4] + 5 * S[5] + 4 * S[6] + 3 * S[7] + 2 * S[8];
		elif XX = "F" then
			a0 := 2 * S[2] + 3 * S[4] + 4 * S[3] + 2 * S[1];
		elif XX = "G" then
			a0 := 3 * S[1] + 2 * S[2];
		fi;
        Ws := ShallowCopy( Concatenation( [ First( REFLECTIONS( Wf ), function( k ) return a0 * k = - a0; end ) ], GeneratorsOfGroup( Wf ) ) );
        r := n;
    fi;
	q := IsomorphismPermGroup( Wf ); Ws := List( Ws, function( i ) return Image( q, i ); end ); W := Image( q );
    if Length( arg ) < 3 then
        degs := [ 0 .. r ];
    else
        degs := arg[3];
    fi;
    Subs := [  ];
    QW := GroupRing( Rationals, W );
    for x in IteratorOfCombinations( [ 0 .. r ] ) do
        Add( Subs, AsSet( x ) );
    od;
    Hs := List( [ 0 .. r ], function ( k )
            return Filtered( Subs, function ( S )
                    return Size( S ) = r - k;
                end );
        end );
    D := [ NullMat( r + 1, 1, QW ) ];
    H0 := List( Hs, function ( k )
            return List( k, function ( x )
                    return Subgroup( W, Ws{List( x, function ( i )
                               return i + 1;
                           end )} );
                end );
        end );
    for i in [ 2 .. r + 1 ] do
        d := [  ];
        for L in Hs[i] do
            dl := ListWithIdenticalEntries( Size( Hs[i - 1] ), Zero( QW ) );
            I0L := SortedList( Difference( [ 0 .. r ], L ) );
            for x in I0L do
                dl[Position( Hs[i - 1], AsSet( Union( L, [ x ] ) ) )] := (- One( QW )) ^ (Position( I0L, x ) + 1);
            od;
            Add( d, dl );
        od;
        Add( D, d );
    od;
    Add( D, NullMat( 1, 1, QW ) );
    D0 := [ TO_INTEGER_MATRIX( W, H0[1], [ W ], D[1] ) ];
    D0 := Concatenation( D0, List( [ 2 .. r + 1 ], function ( i )
              if i - 1 in degs or i - 2 in degs then
                  return TO_INTEGER_MATRIX( W, H0[i], H0[i - 1], D[i] );
              else
                  return TO_INTEGER_MATRIX( W, H0[i], H0[i - 1], Zero( QW ) * D[i] );
              fi;
              return;
          end ) );
    Add( D0, TO_INTEGER_MATRIX( W, [ W ], H0[r + 1], D[r + 2] ) );
    Reps := List( degs, function ( i )
            return HOMOLOGY_REP( W, D0[i + 2], H0[i + 1], D0[i + 1] );
        end );
	QWf := GroupRing( Rationals, Wf );
	D := List( D, function( d )
		k := ShallowCopy( NullMat( DimensionsMat( d )[1], DimensionsMat( d )[2], QWf ) );
		for i in [ 1 .. DimensionsMat( d )[1] ] do
			for j in [ 1.. DimensionsMat( d )[2] ] do
				if d[ i, j ] = One( QW ) then
					k[ i, j ] := One( QWf );
					elif d[ i, j ] = - One( QW ) then
					k[ i, j ] := - One( QWf );
					else
					k[ i, j ] := Zero( QWf );
				fi;
			od; 
		od;
		return k;
		end );
	H0 := List( H0, function( d ) 
		return List( d, function( i ) 
				return Subgroup( Wf, List( GeneratorsOfGroup( i ), function( j ) 
						return PreImagesRepresentative( q, j ); 
						end ) );
				end ); 
		end );
	Reps := List( Reps, function( i )
		return [ i[1], GroupHomomorphismByImages( Wf, Group( List( GeneratorsOfGroup( W ), function( k ) 
				return Image( i[2], k ); end ) ), GeneratorsOfGroup( Wf ), List( GeneratorsOfGroup( Wf ), function( d ) 
						return Image( i[2], Image( q, d ) ); end ) ) ];
			end );
    CH := List( Reps, function ( r )
            return List( ConjugacyClasses( Wf ), function ( c )
                    return TraceMat( Image( r[2], Representative( c ) ) );
                end );
        end );
	Bch := List( CH, function( i ) return List( Irr( Wf ), function( k ) return ScalarProduct( CharacterTable( Wf ), Character( Wf, i ), k ); end ); end );
    return rec(
        coxeter_group := Wf,
        group_algebra := QWf,
        stabilizers := H0,
        differentials := D{[ 2 .. r + 1 ]},
        homology_representations := Reps,
        homology_characters := CH,
		coeffs_in_tbl := Bch );
end );









InstallGlobalFunction( GAMMA,
		[ IsList, IsInt ],
	function ( S, k )
  	  local G, i, Gm, C, g, c;
  	  G := [  ];
  	  if k = 0 then
  	      G := [ [  ] ];
  	  else
  	      for i in [ 1 .. k - 1 ] do
  	          Gm := GAMMA( S, k - i );
  	          C := Combinations( S, i );
  	          for c in C do
  	              for g in Gm do
  	                  if IsSubsetSet( g[Size( g )], c ) then
  	                      Add( G, Concatenation( g, [ c ] ) );
  	                  fi;
  	              od;
  	          od;
  	      od;
  	      for c in Combinations( S, k ) do
  	          Add( G, [ c ] );
  	      od;
  	  fi;
  	  return SortedList( Filtered( G, function ( g )
  	            return Sum( List( g, Size ) ) = k;
  	        end ) );
end );



InstallGlobalFunction( MINIMAL_LEFT_COSETS,
		[ IsList, IsList ],
	function ( S, T )
  	  local W, R, r, u;
  	  W := Group( S );
  	  R := List( RightTransversal( W, Subgroup( W, T ) ), function ( r )
  	          return List( Union( T, [ One( W ) ] ), function ( u )
  	                  return u * r;
  	              end );
  	      end );
  	  return List( R, function ( r )
  	          return r[PositionMinimum( List( r, function ( u )
  	                       return Length( Factorization( W, u ) );
  	                   end ) )] ^ -1;
  	      end );
end );



InstallGlobalFunction( PARITY_INVERSIONS,
		[ IsList, IsList ],
	function ( S, T )
  	  local s, t;
  	  if Size( S ) > 1 then
  	      s := PARITY_INVERSIONS( S{[ 2 .. Size( S ) ]}, T{[ 2 .. Size( S ) ]} ) * (-1) ^ Number( T{[ 2 .. Size( S ) ]}, function ( x )
  	                  return (x < T[1]);
  	              end );
  	  else
  	      s := 1;
  	  fi;
  	  return s;
end );


InstallGlobalFunction( SalvettiComplex,
		[ IsString, IsPosInt, IsInt ],
	function ( XX, n, k )
  	  local s1, s2, s3, s4, F, S, W, QW, o, D, Gkm, Gk, kk, d, p, g, i, tau, gti, R, beta, gp, q, mualpha;
	  	if XX = "I2" then
	  		F := FreeGroup( [ "s1", "s2" ] );; 
			s1 := F.1;; s2 := F.2;;
			S := ShallowCopy( GeneratorsOfGroup( F / [ s1^2, s2^2, (s1*s2)^n ] ) );;
	  	elif XX = "H" and n = 3 then
	  		F := FreeGroup( [ "s1", "s2", "s3"] );;
			s1 := F.1;; s2 := F.2;; s3 := F.3;;
			S := ShallowCopy( GeneratorsOfGroup( F / [s1^2, s2^2, s3^2, (s1*s2)^5, (s2*s3)^3, (s1*s3)^2 ] ) );;
	  	elif XX = "H" and n = 4 then
	  		F := FreeGroup( [ "s1", "s2", "s3", "s4" ] );;
			s1 := F.1;; s2 := F.2;; s3 := F.3;; s4 := F.4;;
			S := ShallowCopy( GeneratorsOfGroup( F / [ s1^2, s2^2, s3^2, s4^2, (s1*s2)^5, (s2*s3)^3, (s3*s4)^3, (s1*s3)^2, (s1*s4)^2, (s2*s4)^2 ] ) );;
	  	else
  	  		S := ShallowCopy( GeneratorsOfGroup( WeylGroup( RootSystem( SimpleLieAlgebra( XX, n, Rationals ) ) ) ) );
		fi;
  	  W := Group( S );
  	  QW := GroupRing( Rationals, W );
  	  o := One( QW );
  	  D := [  ];
  	  Gkm := ShallowCopy( GAMMA( S, 0 ) );
  	  for kk in [ 1 .. k ] do
  	      Gk := ShallowCopy( GAMMA( S, kk ) );
  	      d := NullMat( Size( Gk ), Size( Gkm ), QW );
  	      for p in [ 1 .. Size( Gk ) ] do
  	          g := ShallowCopy( Gk[p] );
  	          Add( g, [  ] );
  	          for i in [ 1 .. Size( g ) - 1 ] do
  	              if Size( g[i] ) > Size( g[i + 1] ) then
  	                  for tau in g[i] do
  	                      gti := Filtered( g[i], function ( x )
  	                              return x <> tau;
  	                          end );
  	                      R := MINIMAL_LEFT_COSETS( g[i], gti );
  	                      for beta in R do
  	                          if IsSubsetSet( gti, List( g[i + 1], function ( x )
  	                                      return x ^ beta;
  	                                  end ) ) then
  	                              gp := Iterated( [ g{[ 1 .. i - 1 ]}, [ gti ], List( g{[ i + 1 .. Size( g ) ]}, function ( x )
  	                                            return List( x, function ( y )
  	                                                    return y ^ beta;
  	                                                end );
  	                                        end ) ], Concatenation );
  	                              while [  ] in gp do
  	                                  Remove( gp, Position( gp, [  ] ) );
  	                              od;
  	                              q := Position( Gkm, List( gp, AsSet ) );
  	                              mualpha := (-1) ^ (Sum( List( g{[ 1 .. i - 1 ]}, Size ) ) + Size( Filtered( g[i], function ( x )
  	                                                return (x <= tau);
  	                                            end ) ) + i * Length( Factorization( W, beta ) )) * Product( List( g{[ (i + 1) .. Size( g ) ]}, function ( x )
  	                                          return PARITY_INVERSIONS( x, List( x, function ( y )
  	                                                    return y ^ beta;
  	                                                end ) );
  	                                      end ) );
  	                              d[p, q] := d[p, q] + mualpha * o * beta;
  	                        fi;
  	                    od;
  	                od;
  	            fi;
  	        od;
  	    od;
  	    Add( D, ShallowCopy( d ) );
  	    Gkm := Gk;
  	od;
  	return rec( coxeter_group := W, group_algebra := QW, differentials := D );
end );