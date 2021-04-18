#
# Salvetti-and-tori-complexes: Routine to build the Salvetti complex of a finite Coxeter group and the Weyl-equivariant cellular homology chain complex of the torus of an irreducible root datum.
#
#! @Chapter Introduction
#!
#! This package provides several functions to build the Salvetti complex of a finite 
#! irreducible Coxeter group and the W-equivariant cellular homology chain complex 
#! of any maximal torus of a simple compact Lie group, described by its root datum.
#! It also allows to recover the homology of the torus as a rational representation
#! of the Weyl group W, given as a map into a matrix group.
#! Moreover, this combinatorics is extended to finite irreducible Coxeter groups,
#! using the longest reflection as the element playing the role of the additional reflection appearing in the Weyl case.
#! As we also provide functions to translate the differentials of the complex
#! as matrices with integer coefficients, one can also study these complexes
#! through an adapted package, such as "Homalg" or "CAP" for instance.
#! The Salvetti complex is constructed as in 
#! <Cite Key="deconcini-salvetti"/> and the one for
#! the tori is as described in 
#! <Cite Key="torus21"/>.
#!
#! @Chapter Functionality
#!
#!
#! This chapter describes the
#! methods of Salvetti-and-tori-complexes



#! @Section Main functions

#! @Description 
#! The input is a string $X$, the type of the Coxeter group, an integer $n$, its rank
#! and a positive integer $k$, which is the maximal degree of the differential to be computed.
#! In the crystallographic case, the input of the rank is as for the command $\texttt{SimpleLieAlgebra}$.
#! The non-crystallographic cases are given by ["H",3], ["H",4] and ["I2",p].
#! It returns a record with the following entries:
#!
#! $\texttt{coxeter}$_$\texttt{group}$: the Coxeter group $W$ of type $X_n$ (or $I_2(n)$),
#!
#! $\texttt{group}$_$\texttt{algebra}$: the rational group algebra of the Coxeter group of type $X_n$,
#!
#! $\texttt{differentials}$: the $k$ first differentials (as right multiplications by matrices with entries in $\mathbb{Z}[W]$) of the Salvetti complex.
#!
#! @Arguments X,n,k
#! @Returns a record
DeclareGlobalFunction( "SalvettiComplex", [ IsString, IsPosInt, IsInt ] );



#! @Description

#! The input is a string $X$, the type of the root datum, a positive integer $n$, its rank (except in the case of $I_2$),
#! and a subset $L$ of $\{0,\dotsc,n\}$ indexing the minuscule coweights generating the lattice of the given root datum. 
#! The set $\{0,\dotsc,n\}$ is indexing the simple roots, together with minus the highest root (which corresponds to 0 in the set, the other roots being indexed as usual, see
#! <Cite Key="humlie"/> for instance).
#! The optional argument $\texttt{degs}$ is a list of degrees at which the homology representations are to be computed.
#! The default value of $\texttt{degs}$ is $[0..n]$ that is, all the homology representations are computed.
#! It returns a record with the following entries:
#!
#! $\texttt{homology}$_$\texttt{representations}$: the list of homology representations $H_i(T,\mathbb{Q})$ (in homogeneous degrees $\texttt{degs}$ if this argument is specified), each one of which being given as a matrix representation,
#!
#! $\texttt{weyl}$_$\texttt{group}$: the Weyl group of the root system of type $X_n$,
#!
#! $\texttt{group}$_$\texttt{algebra}$: the rational group algebra of the Weyl group of type $X_n$,
#!
#! $\texttt{stabilizers}$: the (ordered) list $L$ of stabilizers of simplicies in the barycentric subdivision of the fundamental alcove $\overline{\mathcal{A}_0}$, i.e. the $k^\text{th}$ homogeneous component of the complex is $\bigoplus_{H\in L[k]}\mathbb{Z}[W/H]$
#!
#! $\texttt{differentials}$: the differentials (as right multiplications by matrices with entries in $\mathbb{Z}[W]$) of the complex.
#!
#! @Arguments X,n,L[,degs]
#! @Returns a record
DeclareGlobalFunction( "CellularComplexT", [ IsString, IsPosInt, IsList, IsList ] );

#! @Description
#! The input is a string $X$, the type of the Coxeter group and a positive integer $n$, its rank (except in the case of $I_2$).
#! In the crystallographic case, the input of the rank is as for the command $\texttt{SimpleLieAlgebra}$.
#! The non-crystallographic cases are given by ["H",3], ["H",4] and ["I2",p]. In these cases the \emph{distinguished reflection} $s_0$, which is associated to the highest root in the Weyl case.
#! The computed complex is denoted by $C$.
#! The optional argument $\texttt{degs}$ is a list of degrees at which the homology representation is to be computed.
#! The default value is $[0..n]$ (or $[0..2]$ in the case of $I_2$).
#! The output is a record with the following entries:
#!
#! $\texttt{homology}$_$\texttt{representations}$: the list of homology representations $H_i(C)$ (in the homogeneous degrees $\texttt{degs}$ if this argument is specified), each one of which being given as a matrix representation,
#!
#! $\texttt{homology}$_$\texttt{characters}$: the list of characters of $\texttt{homology}$_$\texttt{representations}$ (note that in the crystallographic case, the representations are irreducible),
#!
#! $\texttt{coeffs}$_$\texttt{in}$_$\texttt{tbl}$: the coefficients of the characters $\texttt{homology}$_$\texttt{characters}$, as linear combinations of irreducible characters of $W$. The indexing of the irreducibles is the same as for $\texttt{Irr}(W)$,
#!
#! $\texttt{coxeter}$_$\texttt{group}$: the Coxeter group of type $X_n$ (or $I_2(n)$),
#!
#! $\texttt{group}$_$\texttt{algebra}$: the rational group algebra of the Coxeter group of type $X_n$,
#!
#! $\texttt{stabilizers}$: the (ordered) list $L$ of stabilizers of ''simplices" (represented by ordered sublists of $[0..n]$ or $[0..2]$), i.e. the $k^\text{th}$ homogeneous component of the complex is $\bigoplus_{H\in L[k]}\mathbb{Z}[W/H]$,
#!
#! $\texttt{differentials}$: the differentials (as right multiplication by matrices with entries in $\mathbb{Z}[W]$) of the complex $C$.
#!
#! @Arguments X,n,[,degs]
#! @Returns a record
DeclareGlobalFunction( "ComplexForFiniteCoxeterGroup", [ IsString, IsPosInt, IsList ] );

#! @Description
#! The input is a string $X$, the type of the root datum, a positive integer $n$, its rank,
#! and a subset $L$ of $\{0,\dotsc,n\}$ indexing the minuscule weights belonging to the lattice of the given root datum.
#! The output is a record with the following entries:
#!
#! $\texttt{fundamental}$_$\texttt{group}$_$\texttt{permutations}$: the elements of the fundamental group $\Omega_L$ of the root datum $(L^\vee,\Phi,L,\Phi^\vee)$, with $\Phi$ the root system of type $X_n$; given as a group of permutations on the set $\{0,\dotsc,n\}$,
#!
#! $\texttt{group}$_$\texttt{algebra}$: the rational group algebra of $\Omega_L$ (in the permutation form),
#!
#! $\texttt{stabilizers}$: the (ordered) list $L$ of stabilizers under $\Omega_L$ (in the permutation form again) of simplicies in the barycentric subdivision of the fundamental alcove $\overline{\mathcal{A}_0}$, i.e. the $k^\text{th}$ homogeneous component of the complex is $\bigoplus_{H\in L[k]}\mathbb{Z}[\Omega_L/H]$
#!
#! $\texttt{differentials}$: the differentials (as right multiplications by matrices with entries in $\mathbb{Z}[\Omega_L]$) of the complex.
#! 
#! @Arguments X,n,L
#! @Returns a record
DeclareGlobalFunction( "ComplexForOmega", [ IsString, IsPosInt, IsList, IsList ] );




#! @Section Side functions


#! @Description
#! The input is a group $G$ ans a subgroup $H$ of $G$.
#! The output is the set of left cosets of $G$ modulo $H$, given as sets partitioning $G$.
#! @Arguments G,H
#! @Returns a set
DeclareGlobalFunction( "LEFT_COSETS", [ IsGroup, IsGroup ] );

#! @Description
#! The input is a element $g$ of an ambient group $G$ and $H$ is a subgroup of the ambient group $G$.
#! The output is the left coset $gH$, given a set of elements of $G$.
#! @Arguments g,H
#! @Returns a set
DeclareGlobalFunction( "LEFT_COSET", [ IsMultiplicativeElementWithInverse, IsGroup ] );

#! @Description
#! The input is a group $G$, an element $x$ of $G$ and a list (or set) $gH$, representing a left coset in $G$.
#! The output is the left coset $xgH$, given by the natural left action of $G$.
#! @Arguments G,x,gH
#! @Returns a set
DeclareGlobalFunction( "ON_LEFT_COSET", [ IsGroup, IsMultiplicativeElementWithInverse, IsListOrCollection ] );

#! @Description
#! The input is a group $G$, a subgroup $H$ of $G$ and an element $g$ of $G$. 
#! The output is a permutation matrix corresponding to the $xH \mapsto gxH$ which permutes the left cosets $G/H$.
#! (the labelling of the permutation is given by a suitable ordering that GAP puts on elements of a given group)
#! @Arguments G,H,g
#! @Returns a permutation matrix
DeclareGlobalFunction( "MATRIX_TRANSITIVE_ACTION", [ IsGroup, IsGroup, IsMultiplicativeElementWithInverse ] );

#! @Description
#! The input is a group $G$, two subgroups $H$ and $K$ of $G$ and an element $g$ of $G$.
#! The output is the matrix of the map $xH \mapsto xgK$ (here again, with labelling given by a suitable ordering on 
#! the elements of $G$). $\textbf{Caution:}$ the element $g$ must conjugate $H$ into $K$, i.e. $H^g \le K$ in order
#! for this map to be well-defined. If this is not the case, the function returns 'fail'. 
#! @Arguments G,H,K,g
#! @Returns an integer matrix
DeclareGlobalFunction( "MATRIX_RIGHT_TRANSLATION", [ IsGroup, IsGroup, IsGroup, IsMultiplicativeElementWithInverse ] );

#! @Description
#! The input is a group $G$, two subgroups $H$ and $K$ of $G$ and an element $x$ of the group algebra $\mathbb{Q}[G]$. 
#! The output is a matrix corresponding to the map $\mathbb{Z}[G/H] \to \mathbb{Z}[G/K]$, linearly extending the map $xH \mapsto xgK$.
#! Here again, for this to be well-defined, one needs that all the non-zero basis vectors appearing in $x$ to conjugate $H$ into $K$.
#! @Arguments G,H,K,x
#! @Returns an integer matrix
DeclareGlobalFunction( "MATRIX_RIGHT_TRANSLATION_ZG", [ IsGroup, IsGroup, IsGroup, IsMultiplicativeElement ] );

#! @Description
#! The input is a group $G$, two lists $L_H$ and $L_K$ of subgroups of $G$ and a matrix $d$ with coefficients in $\mathbb{Z}[G]$. 
#! The output is an integer matrix. The matrix $d$ represents a $G$-homomorphism given by right multiplication by the matrix 
#! $d : \bigoplus_{H\in L_H}\mathbb{Z}[G/H] \to \bigoplus_{K\in L_K}\mathbb{Z}[G/K]$ and the function returns the corresponding
#! homomorphism of abelian groups $\bigoplus_{H\in L_H}\mathbb{Z}^{[G:H]}\to\bigoplus_{K\in L_K}\mathbb{Z}^{[G:K]}$, given by
#! right multiplication by an integer matrix. 
#! @Arguments G,LH,LK,d
#! @Returns an integer matrix
DeclareGlobalFunction( "TO_INTEGER_MATRIX", [ IsGroup, IsList, IsList, IsMatrix ] );

#! @Description
#! The input is a group $G$, an integer matrix $D$, a list $L_H$ of subgroups of $G$ and another matrix $Dp$, where $D$ and $Dp$
#! are the integer versions of $G$-homomorphisms $\stackrel{\tiny{d}}\longrightarrow \bigoplus_{H\in L_H}\mathbb{Z}[G/H] \stackrel{\tiny{dp}}\longrightarrow$.
#! The matrices $D$ and $Dp$ are supposed to verify $Dp*D=0$.
#! The output is the matrix representation of $G$ on the vector space $\ker(dp)/\mathrm{im}(d)$, given as pairs $[dim,rep]$, where $dim$ is the dimension
#! of the representation, and $rep$ is the homomorphism $G\to GL_{dim}(\mathbb{Q})$.
#! @Arguments G,D,LH,Dp
#! @Returns a liste
DeclareGlobalFunction( "HOMOLOGY_REP", [ IsGroup, IsMatrix, IsList, IsMatrix ] );

#! @Description
#! The input is a Weyl group $W$, with GeneratorsOfGroup(W) known to be the simple reflections of $W$.
#! The output is the list of the reflections of $W$.
#! @Arguments W
#! @Returns a list
DeclareAttribute( "REFLECTIONS", IsGroup );

#! @Description
#! The input is a Weyl group $W$.
#! The output is its longest element.
#! @Arguments W
#! @Returns an element of W
DeclareAttribute( "LONGEST_ELEMENT", IsGroup );

#! @Description
#! The input is a Weyl group $W$ and a subset $J$ of $\{0,\dotsc,n\}$, with $n=\mathrm{rk}(W)$.
#! The output is a list $[O,Ws,Jp]$ where $O$ is the image in $W$ of the subgroup of the group $\Omega$
#! of the symmetries of the fundamental alcove generated by those $\omega_i\in \Omega$ such that $i$ is in $J$,
#! $Ws$ is the set of simple reflections of $W$, together with the reflection associated to the highest root and
#! $Jp$ is the subset of $\{0,\dotsc,n\}$ indexing the elements of $O$.
#! @Arguments W,J
#! @Returns a list
DeclareOperation( "PROJ_OMEGA", [ IsGroup, IsList ] );

#! @Description
#! The input is a string $X$, and a positive integer $n$.
#! Thet output is the list $S$ of subsets of $\{0,\dotsc,n\}$ indexing the faces of the polytope (fundamental alcove) $\overline{\mathcal{A}_0}$
#! associated to the root system of type $X_n$.
#! @Arguments X,n
#! @Returns a list
DeclareOperation( "LATTICE_FROM_ROOT_SYSTEM", [ IsString, IsPosInt ] );

#! @Description
#! The input is a string $X$, a positive integer $n$ and a subset $J$ of $\{0,\dotsc,n\}$.
#! The output is a list $[Om,J']$ with $Om$ the elements of the subgroup of $\Omega$ generated by $J$ (with same ordering as in 
#! <Ref Oper="PROJ_OMEGA" Label="for IsGroup, IsList"/>), given as a group of permutations of $\{0,\dotsc,n\}$ and $J'$ is the list of fundamental coweights in the lattice spanned by $J$.
#! This description is convenient in order to make $\Omega$ act on our face lattice.
#! @Arguments X,n,J
#! @Returns a list
DeclareOperation( "PERMS_OF_OMEGA", [ IsString, IsPosInt, IsList] );

#! @Description
#! The input is a permutation $s$ and a subset $A$ of $\{0,\dotsc,n\}$ corresponding to a face of $\overline{\mathcal{A}_0}$.
#! The output is the subset $s\cdot A=s(A)$.
#! @Arguments s,A
#! @Returns a pair
DeclareGlobalFunction( "PERM_ON_LATTICE", [ IsPerm, IsListOrCollection ] );

#! @Description
#! The input is a string $X$, a positive integer $n$ and a list $J$.
#! The output is the list of the $\Omega_J$-orbits in $\mathcal{P}(\{0,\dotsc,n\})$
#! @Arguments X,n,J
#! @Returns a list
DeclareGlobalFunction( "PREORBS_OF_LATTICE", [ IsString, IsPosInt, IsList ] );

#! @Description
#! The input is a string $X$, a positive integer $n$ and a list $J$.
#! The output is a list indexed by 
#! <Ref Func="PREORBS_OF_LATTICE"/>, with the stabilizer of a representative
#! of each orbit.
#! @Arguments X,n,J
#! @Returns a list
DeclareGlobalFunction( "STABS_IN_LATTICE", [ IsString, IsPosInt, IsList ] );

#! @Description
#! The input is a set $Z$ with a partial order given as a boolean function $ord$.
#! The output is the set of all the increasing chains of $(Z,ord)$.
#! @Arguments Z,ord
#! @Returns a list
DeclareGlobalFunction( "CHAINS_OF_POSET", [ IsList, IsFunction ] );

#! @Description
#! The input is a set $Z$ with a partial order given as a boolean function $ord$.
#! The output is a list having as $k^\text{th}$ entry the set of increasing chains of $(Z,ord)$ of length $k$.
#! @Arguments Z,ord
#! @Returns a list 
DeclareGlobalFunction( "CHAINS_BY_LENGTH", [ IsList, IsFunction ] );

#! @Description
#! The input is a permutation $s$ and a chain $C=\left(Z_1\lneq Z_2\lneq \cdots\lneq Z_k\right)$ of some poset $Z$
#! on which $s$ acts.
#! The output is the chain $s(C):=(s(Z_1)\lneq s(Z_2)\lneq \cdots\lneq s(Z_k))$.
#! @Arguments s,C
#! @Returns a list
DeclareGlobalFunction( "PERM_ACTS_ON_CHAIN", [ IsPerm, IsList ] );

#! @Description
#! The input is a string $X$, a positive integer $n$ and a list $J$.
#! The output is the list of orbits of chains of the lattice $\mathcal{P}(\{0,\dotsc,n\})$ under the
#! action of the subgroup $\Omega_J$ of $\Omega$ corresponding to $J$.
#! @Arguments X,n,J
#! @Returns a list
DeclareGlobalFunction( "PREORBS_OF_CHAINS", [ IsString, IsPosInt, IsList ] );

#! @Description
#! The input is a string $X$, a positive integer $n$ and a list $J$.
#! The output is a list indexed by 
#! <Ref Func="PREORBS_OF_CHAINS"/>, with the stabilizer of a representative
#! of each orbit.
#! @Arguments X,n,J
#! @Returns a list
DeclareGlobalFunction( "STABS_OF_CHAINS", [ IsString, IsPosInt, IsList ] );

#! @Description
#! The input is a string $X$, a positive integer $n$ and a list $J$.
#! The output is a list indexed by lengths of chains of $\mathcal{P}(\{0,\dotsc,n\})$ and each entry
#! is itself a list containing the stabilizer of every orbit of chains.
#! @Arguments X,n,J
#! @Returns a list
DeclareGlobalFunction( "LIST_OF_STABS", [ IsString, IsPosInt, IsList ] );

#! @Description
#! The input is a sorted list $A$.
#! The output is a list containing $A\setminus\{x\}$ as entries (for $x\in A$), together with an alternate sign.
#! The "boundary" would then by the sum of this list, if this would make any sense here.
#! @Arguments A
#! @Returns a list
DeclareGlobalFunction( "BOUNDARY_SIMPLEX", [ IsList ] );

#! @Description
#! The input is a chain represented by a list $C$, a list of stabilizers $Ls$ (as for 
#! <Ref Func="LIST_OF_STABS"/> above)
#! and a list $Om$ (representing elements of the group $\Omega_J$ here).
#! The output is a pair consisting of the position of the orbit $\Omega_J\cdot C$ in the list of stabilizers
#! and the elements of $Om$ sending $C$ on the representative of that orbit.
#! @Arguments C,Ls,Om
#! @Returns a pair
DeclareGlobalFunction( "SEARCH_FOR_ORBIT", [ IsList, IsList, IsList ] );

#! @Description
#! The input is a set $S$ and a non-negative integer $n$.
#! The output is the list of flags of $S$ of cardinality $k$.
#! @Arguments S,k
#! @Returns a list
DeclareGlobalFunction( "GAMMA", [ IsList, IsInt ] );

#! @Description
#! The input is a list $S$ (indexing the simple reflections of an irreducible Weyl group) and a sublist $T$ of $S$.
#! The output is the list of minimal length coset representatives modulo the parabolic subgroup generated by $T$.
#! @Arguments S,T
#! @Returns a list
DeclareGlobalFunction( "MINIMAL_LEFT_COSETS", [ IsList, IsList ] );

#! @Description
#! The input is two lists $A$ and $B$ of same cardinality and a common (implicite) ordering.
#! The output is the parity of the number of inversions between $A$ and $B$.
#! @Arguments A,B
#! @Returns +1 or -1
DeclareGlobalFunction( "PARITY_INVERSIONS", [ IsList, IsList ] );
