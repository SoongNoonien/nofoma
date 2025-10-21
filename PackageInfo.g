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
],

#SourceRepository := rec( Type := "TODO", URL := "URL" ),
#IssueTrackerURL := "TODO",
PackageWWWHome := "https://pnp.mathematik.uni-stuttgart.de/idsr/idsr1/geckmf/frobenius.g/",
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

BannerString := Concatenation(
"──────────────────────────────────────────────────────────────────────────\n",
"Loading  nofoma 1.0 (Normal forms of matrices), \n",
"by Meinolf Geck (https://pnp.mathematik.uni-stuttgart.de/idsr/idsr1/geckmf/)\n",
"Help about the main functions in this package is obtained by typing:\n",
"    ?MaximalVectorMat     ?FrobeniusNormalForm     ?JordanChevalleyDecMat\n",
"(or type:    ?nofoma: The nofoma package)\n",
"──────────────────────────────────────────────────────────────────────────\n"),

PackageDoc := rec(
  BookName  := "nofoma",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Normal forms of matrices",
  Autoload := true
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

