#
# Salvetti-and-tori-complexes: Routine to build the Salvetti complex of a finite Coxeter group and the Weyl-equivariant cellular homology chain complex of the torus of an irreducible root datum.
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "Salvetti-and-tori-complexes",
Subtitle := "Routine to build the Salvetti complex of a finite Coxeter group and the Weyl-equivariant cellular homology chain complex of the torus of an irreducible root datum.",
Version := "1.0",
Date := "28/12/2020", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Arthur",
    LastName := "Garnier",
    #WWWHome := TODO,
    Email := "arthur.garnier@u-picardie.fr",
    IsAuthor := true,
    IsMaintainer := true,
    PostalAddress := Concatenation(
               "LAMFA\n",
               "Universite de Picardie Jules Verne\n",
               "33 Rue Saint-Leu\n",
               "80000 Amiens\n",
               "France" ),
    Place := "Amiens",
    Institution := "UPJV",
  ),
],

#SourceRepository := rec( Type := "TODO", URL := "URL" ),
#IssueTrackerURL := "TODO",
PackageWWWHome := "https://github.com/arthur-garnier/",
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL     := Concatenation( ~.PackageWWWHome,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "Salvetti-and-tori-complexes",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Routine to build the Salvetti complex of a finite Coxeter group and Weyl-equivariant cellular homology chain complex of a torus of an irreducible root datum.",
),

Dependencies := rec(
  GAP := ">= 4.9",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


