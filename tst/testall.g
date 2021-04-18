#
# Salvetti-and-tori-complexes: Routine to build the Salvetti complex of a finite Coxeter group and the Weyl-equivariant cellular homology chain complex of the torus of an irreducible root datum.
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#

LoadPackage( "Salvetti-and-tori-complexes" );

options := rec(
    exitGAP := true,
    testOptions := rec(
        compareFunction := "uptowhitespace"
    ),
);

TestDirectory( DirectoriesPackageLibrary( "Salvetti-and-tori-complexes", "tst" ), options );

FORCE_QUIT_GAP( 1 ); # if we ever get here, there was an error
