// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers

// Package Headers
#include <protocols/noesy_assign/CovalentCompliance.hh>


// Project Headers
#include <core/id/NamedAtomID.hh>

// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>

// Utility headers
//#include <ObjexxFCL/format.hh>
//#include <ObjexxFCL/string.functions.hh>

//#include <utility/string_util.hh>
// #include <utility/excn/Exceptions.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <core/util/prof.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

#include <utility/vector1.hh>

// #include <core/options/option.hh>
// #include <core/options/keys/abinitio.OptionKeys.gen.hh>
// #include <core/options/keys/run.OptionKeys.gen.hh>
//#include <core/options/keys/templates.OptionKeys.gen.hh>

// Utility headers
//#include <core/options/option_macros.hh>


//// C++ headers
//#include <cstdlib>
//#include <string>


namespace protocols {
namespace noesy_assign {

using namespace core;


/// WARNING WARNING WARNING THE COVALENTCOMPLIANCE CLASS IS THREAD UNSAFE
bool covalent_compliance( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2 ) {
	return CovalentCompliance::get_instance()->is_compliant( atom1, atom2 );
}
// core::Real compute_RPF_score( CrossPeakList const& cpl, pose::Pose const& pose, core::Real  ) {

//   Real model_sum_dist( 0 );
//   for ( ResonanceList::const_iterator it1 = resonances.begin(); it1 != resonances.end(); ++it1 ) {
// 		if ( it1->second.atom().atom()[0]=='C' ) continue;
// 		ResonanceList::const_iterator it2 = it1;
// 		++it2;
//     for ( ; it2 != resonances.end(); ++it2 ) {
//       //cross peaks have also symmetry.. so double count all things later...

//       //abuse the constraint code to parse atom names...
// 			if ( it2->second.atom().atom()[0]=='C' ) continue;
// 			core::scoring::constraints::AmbiguousNMRDistanceConstraint a_cst( it1->second.atom(), it2->second.atom(), pose, NULL );
//       bool proton_at_it1( pose.residue_type( a_cst.atom( 1 ).rsd() ).atom_is_hydrogen( a_cst.atom( 1 ).atomno() ) );
//       Size second_atom_index( a_cst.natoms( 1 ) + 1 );
//       bool proton_at_it2( pose.residue_type( a_cst.atom( second_atom_index ).rsd() ).atom_is_hydrogen( a_cst.atom( second_atom_index ).atomno() ) );
//       if ( proton_at_it1 && proton_at_it2 ) {
// 				Real dist( a_cst.dist( pose ) );
// 				Real invd = 1.0/dist;
// 				Real invd3 = invd*invd*invd;
// 				Real invd6 = invd3*invd3;
// 				model_sum_dist+=2*invd6; //times 2 due to symmetry
//       }
// 		}
// 	}

// 	Real peak_sum_dist( 0 );
// 	Size nr_true_peaks( 0 );
// 	for ( CrossPeakList::const_iterator it = cpl.begin(); it != cpl.end(); ++it ) {
// 		//		core::scoring::constraints::AmbiguousNMRConstraintOP peak_cst( (*it)->create_constraint( cpl.resonances(), pose ) );
// 		//Real dist( peak_cst->dist( pose ) );
// 		//		peak_sum_dist += pow( dist, -6 );
// 		//		nr_true_peaks += dist <= dcut ;
//   }

// 	Real RPF_precision(  model_sum_dist/peak_sum_dist );
// 	Real RPF_recall( nr_true_peaks/cpl.size() );
// 	return 2*RPF_recall*RPF_precision/( RPF_recall + RPF_precision );

// }

}
}
