// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifdef BOINC
#ifdef WIN32
#ifndef __CYGWIN__
#undef read
// this redifintion is done by boinc headers.
// if you compile your .cc file without this change,
// but later use the .hh file after the boinc header is included somewhere up in the include tree,
// you'll get linking problems like
// some_class::_read unresolved symbol
#endif
#endif
#endif
