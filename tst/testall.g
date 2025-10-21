#
# nofoma: Normal forms of matrices
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "nofoma" );

TestDirectory(DirectoriesPackageLibrary( "nofoma", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
