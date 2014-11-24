// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file StrandAssemblyCommon.hh
/// @brief Has most basic inclusions and namespace usings
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_features_strand_assembly_StrandAssemblyCommon_hh
#define INCLUDED_protocols_features_strand_assembly_StrandAssemblyCommon_hh

//Basic rosetta
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/strand_assembly.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

//C library
#include <math.h> // for round, floor, ceil, trunc, sqrt

// c++
#include <algorithm>	// for avg,min,max
#include <cmath>	// for std::abs				// reference:	http://www.cplusplus.com/reference/cmath/abs/
#include <fstream>
#include <iostream>
#include <numeric>
//#include <stdio.h>     //for remove( ) and rename( )
#include <stdlib.h> // for std::abs()
#include <vector>

//Core
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh> // for dssp application
#include <core/pose/PDBInfo.hh> // maybe for PDBInfoCOP
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh> // ScoreFunction.hh seems required for compilation of InterfaceAnalyzerMover.hh
#include <core/scoring/ScoreFunctionFactory.hh> // maybe needed for "get_score_function" ?
#include <core/types.hh>

//DSSP
#include <core/scoring/dssp/Dssp.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <cppdb/frontend.h>

// for get_sw_can_by_sh_id, get_central_residues_in_each_of_two_edge_strands
#include <vector>

// for parse_my_tag
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMap.hh>

// for string return
#include <string>

//for vector
#include <numeric/xyzVector.hh>
#include <core/id/NamedAtomID.hh>

//Others
#include <numeric/xyz.functions.hh> // for torsion calculations
#include <protocols/analysis/InterfaceAnalyzerMover.hh> // for SASA

//Protocols
#include <protocols/features/FeaturesReporter.hh>

// Utility
#include <utility/numbers.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh> // for utility::vector1<Column> primary_key_columns;

// Utility: exception handling
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>


#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


template <typename T, size_t N> const T* mybegin(const T (&a)[N]) { return a; }
template <typename T, size_t N> const T* myend  (const T (&a)[N]) { return a+N; }
// reference:	http://stackoverflow.com/questions/9874802/how-can-i-get-the-max-or-min-value-in-a-vector-c


namespace protocols {
namespace features {
namespace strand_assembly {

using core::id::NamedAtomID;
using numeric::xyzVector;


} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
