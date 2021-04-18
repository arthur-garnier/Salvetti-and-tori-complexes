#
# Salvetti-and-tori-complexes: Routine to build the Salvetti complex of a finite Coxeter group and the Weyl-equivariant cellular homology chain complex of the torus of an irreducible root datum.
#
# This file is a script which compiles the package manual.
#
if fail = LoadPackage("AutoDoc", "2018.02.14") then
    Error("AutoDoc version 2018.02.14 or newer is required.");
fi;
LoadPackage("AutoDoc");
AutoDoc( "Salvetti-and-tori-complexes" : 
        scaffold :=
          rec(
            gapdoc_latex_options := rec(
            LateExtraPreamble := Concatenation( "\\usepackage{tikz}\n\\usetikzlibrary{arrows}\n\\usepackage{amsmath}",
                                  "\\pgfarrowsdeclarecombine{twohead}{twohead}{latex}{latex}{latex}{latex}" ) ),
            bib := ("bibliography.bib")
          ),
        autodoc :=
         rec( files := [ "doc/Intros.autodoc" ],
         scan_dirs := [ "gap", "examples", "doc" ] ),
         maketest := rec( folder := "tst",
                          commands :=
                            [ "LoadPackage( \"Salvetti-and-tori-complexes\" );",
                             ]
                           )
);

QUIT;


## bib := ("Salvetti-and-tori-complexes.bib")
