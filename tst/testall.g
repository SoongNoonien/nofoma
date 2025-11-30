LoadPackage("nofoma");

TestDirectory(DirectoriesPackageLibrary("nofoma","tst"),
              rec(exitGAP     := true,
                  testOptions := rec(compareFunction := "uptowhitespace",
                                     transformFunction := "removenl") ) );