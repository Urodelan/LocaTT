.onAttach<-function(libname,pkgname) {
  packageStartupMessage(
    paste0("Welcome to the Local Taxa Tool (LocaTT)! This is version ",
           utils::packageVersion(pkgname),".\n",
           "   ______  ______  ______   Geographically-conscious\n",
           "    |  | \\/ |  | \\/ |  |    taxonomic assignment\n",
           "   _|__|_/\\_|__|_/\\_|__|_   for DNA metabarcoding.\n",
           "For a tutorial, see: https://urodelan.github.io/LocaTT_Tutorial/\n",
           "Command line BLAST is required, see tutorial for details.\n",
           "See citation('LocaTT') if used in publications."))
}