#
# nofoma: Normal forms of matrices
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "nofoma",
Subtitle := "Normal forms of matrices",
Version := "1.0",
Date := "01/06/2022", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Meinolf",
    LastName := "Geck",
    WWWHome := "https://pnp.mathematik.uni-stuttgart.de/idsr/idsr1/geckmf/",
    Email := "meinolf.geck@mathematik.uni-stuttgart.de",
    IsAuthor := true,
    IsMaintainer := true,
    PostalAddress := "Fachbereich Mathematik, Pfaffenwaldring 57, 70569 Stuttgart, Germany",
    Place := "Stuttgart",
    Institution := "University of Stuttgart",
  ),
  rec(
    FirstNames := "Alia",
    LastName := "Bonnet",
    Email := "alia.bonnet@rwth.aachen.de",
    IsAuthor := false;
    IsMaintainer := true;
    Place := "Aachen",
    Institution := "RWTH Aachen"
  )
],

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/gap-packages/", ~.PackageName ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", ~.PackageName ),
README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),
ArchiveFormats := ".tar.gz .tar.bz2",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "This package computes maximal vectors, the Frobenius normal form and 
                    the Jordan-Chevalley decomposition of a (square) matrix over any field 
                    that is available in GAP. In particular, it extends the already 
                    existing GAP function 'MininmalPolynomial' by computing also a vector 
                    whose local minimal polynomial is equal to the minimal polynomial of a 
                    given matrix.",

BannerString := Concatenation(
"──────────────────────────────────────────────────────────────────────────\n",
"Loading  nofoma ", ~.Version, " (Normal forms of matrices), \n",
"by Meinolf Geck (https://pnp.mathematik.uni-stuttgart.de/idsr/idsr1/geckmf/)\n",
"Help about the main functions in this package is obtained by typing:\n",
"    ?MaximalVectorMat     ?FrobeniusNormalForm     ?JordanChevalleyDecMat\n",
"(or type:    ?nofoma: The nofoma package)\n",
"──────────────────────────────────────────────────────────────────────────\n"),

PackageDoc := rec(
  BookName  := "nofoma",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Normal forms of matrices",
),

Dependencies := rec(
  GAP := ">= 4.11",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));

