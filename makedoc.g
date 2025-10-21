#
# nofoma: Normal forms of matrices
#
# This file is a script which compiles the package manual.
#
if fail = LoadPackage("AutoDoc", "2018.02.14") then
    Error("AutoDoc version 2018.02.14 or newer is required.");
fi;

#AutoDoc( rec( scaffold := true, autodoc := true ) );
AutoDoc( rec( scaffold := rec( TitlePage := false,
                includes := [ "_Chapter_The_nofoma_package.xml" ], 
                index := false, bib := false ) #, autodoc := true 
));

# Thomas: Setze TitlePage=false und in 'doc/title.xml' wird in der 
# '<Title>'-Zeile 'nofoma' durch nofoma<Index>nofoma</Index>' ersetzt.
