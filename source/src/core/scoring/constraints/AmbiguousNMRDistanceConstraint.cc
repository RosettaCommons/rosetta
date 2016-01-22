// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/AmbiguousNMRDistanceConstraint.cc
///
/// @brief
/// @author Oliver Lange

// Unit Headers
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>

// Package Headers
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AA.hh>
// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/prof.hh>
// Numeric Headers
#include <numeric/deriv/distance_deriv.hh>

//Auto Headers
#include <core/id/SequenceMapping.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "core.io.constraints" );


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

inline bool is_aromatic( pose::Pose const& pose, core::Size res ) {
	using namespace core::chemical;
	return pose.residue_type( res ).aa() == aa_phe
		|| pose.residue_type( res ).aa() == aa_tyr
		|| pose.residue_type( res ).aa() == aa_trp ;
}

bool is_aromatic( core::chemical::AA aa ) {
	using core::chemical::AA;
	using namespace core::chemical;
	return ( aa == aa_trp || aa == aa_phe || aa == aa_tyr || aa == aa_his );
}

void parse_NMR_name( std::string name, core::Size res, core::chemical::AA aa, NamedAtoms& atoms ) {
	using core::chemical::AA;
	using namespace core::chemical;
	using core::id::NamedAtomID;
	///now fix all the problems with atomnames:
	if ( (name.substr(0,2) == "QA" || name == "HA") && ( aa == aa_gly ) ) {
		atoms.push_back( NamedAtomID( "1HA", res ) );
		atoms.push_back( NamedAtomID( "2HA", res ) );
	} else if ( ( name == "QB" || name == "HB" ) && (
			aa == aa_ala
			|| aa == aa_leu
			|| aa == aa_ser
			|| aa == aa_asn
			|| aa == aa_gln
			|| aa == aa_pro
			|| aa == aa_lys
			|| aa == aa_cys
			|| aa == aa_asp
			|| aa == aa_glu
			|| aa == aa_arg
			|| aa == aa_tyr
			|| aa == aa_phe
			|| aa == aa_trp
			|| aa == aa_his
			|| aa == aa_met ) ) {
		atoms.push_back( NamedAtomID( "1HB", res ) );
		atoms.push_back( NamedAtomID( "2HB", res ) );
		if ( aa == aa_ala ) atoms.push_back( NamedAtomID( "3HB", res ) );
	} else if ( ( name == "QD1" ||  name == "HD1" ) && (
			aa == aa_ile
			|| aa == aa_leu ) ) {
		atoms.push_back( NamedAtomID( "1HD1", res ) );
		atoms.push_back( NamedAtomID( "2HD1", res ) );
		atoms.push_back( NamedAtomID( "3HD1", res ) );
	} else if ( ( name == "QD2" || name == "HD2" ) && (
			aa == aa_leu
			|| aa == aa_asn ) ) {
		atoms.push_back( NamedAtomID( "1HD2", res ) );
		atoms.push_back( NamedAtomID( "2HD2", res ) );
		if ( aa == aa_leu ) atoms.push_back( NamedAtomID( "3HD2", res ) );
	} else if ( ( name == "QQD" || name == "QD" ) && aa == aa_leu ) {
		atoms.push_back( NamedAtomID( "1HD2", res ) );
		atoms.push_back( NamedAtomID( "2HD2", res ) );
		atoms.push_back( NamedAtomID( "3HD2", res ) );
		atoms.push_back( NamedAtomID( "1HD1", res ) );
		atoms.push_back( NamedAtomID( "2HD1", res ) );
		atoms.push_back( NamedAtomID( "3HD1", res ) );
	} else if ( ( name == "QD" || name == "HD" ) && (
			aa == aa_pro
			|| aa == aa_lys
			|| aa == aa_arg
			|| aa == aa_tyr
			|| aa == aa_phe
			|| aa == aa_his ) ) {
		if ( aa == aa_arg || aa == aa_lys || aa == aa_pro ) {
			atoms.push_back( NamedAtomID( "1HD", res ) );
			atoms.push_back( NamedAtomID( "2HD", res ) );
		} else { // tyr, phe, his (his doesn't get down here... )
			atoms.push_back( NamedAtomID( "HD1", res ) );
			atoms.push_back( NamedAtomID( "HD2", res ) );
		}
	} else if ( ( name == "QE" || name == "HE" ) && (
			aa == aa_tyr
			|| aa == aa_phe
			|| aa == aa_trp
			|| aa == aa_met
			|| aa == aa_lys ) ) {
		if ( aa == aa_phe || aa == aa_tyr ) {
			atoms.push_back( NamedAtomID( "HE1", res ) );
			atoms.push_back( NamedAtomID( "HE2", res ) );
		} else if ( aa == aa_trp ) {
			atoms.push_back( NamedAtomID( "HE1", res ) );
			atoms.push_back( NamedAtomID( "HE3", res ) );
		} else if ( aa == aa_met || aa == aa_lys ) { //MET LYS
			atoms.push_back( NamedAtomID( "1HE", res ) );
			atoms.push_back( NamedAtomID( "2HE", res ) );
			if ( aa == aa_met ) {
				atoms.push_back( NamedAtomID( "3HE", res ) );
			}
		} //QE MET LYS
	}  else if ( ( name == "QG" || name == "HG" ) && (
			aa == aa_met
			|| aa == aa_gln
			|| aa == aa_pro
			|| aa == aa_lys
			|| aa == aa_glu
			|| aa == aa_arg ) ) {
		atoms.push_back( NamedAtomID( "1HG", res ) );
		atoms.push_back( NamedAtomID( "2HG", res ) );
		//  atoms.push_back( NamedAtomID( "3HE", res ) );
	} else if ( ( name == "QG2" || name == "HG2" ) && (
			aa == aa_ile
			|| aa == aa_thr
			|| aa == aa_val ) ) {
		atoms.push_back( NamedAtomID( "1HG2", res ) );
		atoms.push_back( NamedAtomID( "2HG2", res ) );
		atoms.push_back( NamedAtomID( "3HG2", res ) );
	} else if ( ( name == "QQG" || name == "HG" ) && (
			aa == aa_ile
			|| aa == aa_val ) ) {
		atoms.push_back( NamedAtomID( "1HG2", res ) );
		atoms.push_back( NamedAtomID( "2HG2", res ) );
		atoms.push_back( NamedAtomID( "3HG2", res ) );
		atoms.push_back( NamedAtomID( "1HG1", res ) );
		atoms.push_back( NamedAtomID( "2HG1", res ) );
		if ( aa != aa_ile ) atoms.push_back( NamedAtomID( "3HG1", res ) );
	} else if ( ( name == "QG1" || name == "HG1" ) && (
			aa == aa_ile
			|| aa == aa_val ) ) {
		atoms.push_back( NamedAtomID( "1HG1", res ) );
		atoms.push_back( NamedAtomID( "2HG1", res ) );
		if ( aa != aa_ile ) atoms.push_back( NamedAtomID( "3HG1", res ) );
	} else if ( ( name == "QE2" || name == "HE2" || name =="HE" ) && aa == aa_gln ) {
		atoms.push_back( NamedAtomID( "1HE2", res ) );
		atoms.push_back( NamedAtomID( "2HE2", res ) );
	} else if ( ( name == "HZ" || name == "QZ" ) && ( aa == aa_trp || aa == aa_lys ) ) {
		if ( aa == aa_lys ) {
			atoms.push_back( NamedAtomID( "1HZ", res ) );
			atoms.push_back( NamedAtomID( "2HZ", res ) );
			atoms.push_back( NamedAtomID( "3HZ", res ) );
		} else { //aa==trp
			atoms.push_back( NamedAtomID( "HZ2", res ) );
			atoms.push_back( NamedAtomID( "HZ3", res ) );
		}
	} else  if ( name == "HB1" ) {
		atoms.push_back( NamedAtomID( "1HB", res ) );
	} else  if ( name == "HB2" ) {
		atoms.push_back( NamedAtomID( "2HB", res ) );
	} else  if ( name == "HB3" ) {
		if (  aa != aa_ala ) {
			atoms.push_back( NamedAtomID( "1HB", res ) ); //yeah they call it 2HB and 3HB...
		} else {
			atoms.push_back( NamedAtomID( "3HB", res ) );
		}
	} else  if ( name == "HD1" && !is_aromatic( aa ) ) {
		atoms.push_back( NamedAtomID( "1HD", res ) );
	} else  if ( name == "HD2" && !is_aromatic( aa ) ) {
		atoms.push_back( NamedAtomID( "2HD", res ) );
	} else  if ( name == "HD3" ) { //LYS, PRO, ARG  no other has HD3
		atoms.push_back( NamedAtomID( "1HD", res ) );

	} else  if ( name == "HG1" && aa != aa_thr ) {
		atoms.push_back( NamedAtomID( "1HG", res ) );
	} else  if ( name == "HG2" ) {
		atoms.push_back( NamedAtomID( "2HG", res ) );
	} else  if ( name == "HG3" ) { //GLU, ARG, GLN, MET
		atoms.push_back( NamedAtomID( "1HG", res ) );

	} else  if ( name == "HA1" ) {
		atoms.push_back( NamedAtomID( "1HA", res ) );
	} else  if ( name == "HA2" ) {
		atoms.push_back( NamedAtomID( "2HA", res ) );
	} else  if ( name == "HA3" ) { //GLY
		atoms.push_back( NamedAtomID( "1HA", res ) );

	} else  if ( name == "HZ1" && aa == aa_lys ) {
		atoms.push_back( NamedAtomID( "1HZ", res ) );
	} else  if ( name == "HZ2" && aa == aa_lys ) {
		atoms.push_back( NamedAtomID( "2HZ", res ) );
	} else  if ( name == "HZ3" && aa == aa_lys ) {
		atoms.push_back( NamedAtomID( "3HZ", res ) );
	} else if ( ( name == "HE1" || name == "HE2" || name=="HE3" ) && !is_aromatic( aa ) ) { //trp and similar is already done
		if ( ( name == "HE3" && aa != aa_met ) || name == "HE1" ) {
			atoms.push_back( NamedAtomID( "1HE", res ) ); //e.g. LYS
		} else if ( name == "HE2" ) {
			atoms.push_back( NamedAtomID( "2HE", res ) );
			// atoms.push_back( id::AtomID( pose.residue_type(res).atom_index(name.substr(2,1)+name.substr(0,2)), res ) );
		} else if ( name == "HE3" ) {
			atoms.push_back( NamedAtomID( "3HE", res ) );
		}
	} else if ( name == "HD11" ) {
		atoms.push_back( NamedAtomID( "1HD1", res ) );
	} else if ( name == "HD12" ) {
		atoms.push_back( NamedAtomID( "2HD1", res ) );
	} else  if ( name == "HD13" ) {
		atoms.push_back( NamedAtomID( "3HD1", res ) );
	} else  if ( name == "HD21" ) {
		atoms.push_back( NamedAtomID( "1HD2", res ) );
	} else if ( name == "HD22" ) {
		atoms.push_back( NamedAtomID( "2HD2", res ) );
	} else  if ( name == "HD23" ) {
		atoms.push_back( NamedAtomID( "3HD2", res ) );
	} else  if ( name == "HG11" ) {
		atoms.push_back( NamedAtomID( "1HG1", res ) );
	} else  if ( name == "HG12" ) {
		atoms.push_back( NamedAtomID( "2HG1", res ) );
	} else if ( name == "HG13" ) {
		if ( aa == aa_ile ) {
			atoms.push_back( NamedAtomID( "1HG1", res ) );
		} else {
			atoms.push_back( NamedAtomID( "3HG1", res ) );
		}
	} else  if ( name == "HG21" ) {
		atoms.push_back( NamedAtomID( "1HG2", res ) );
	} else  if ( name == "HG22" ) {
		atoms.push_back( NamedAtomID( "2HG2", res ) );
	} else if ( name == "HG23" ) {
		atoms.push_back( NamedAtomID( "3HG2", res ) );
	} else if ( name == "HE11" ) {
		atoms.push_back( NamedAtomID( "1HE1", res ) );
	} else if ( name == "HE12" ) {
		atoms.push_back( NamedAtomID( "2HE1", res ) );
	} else if ( name == "HE13" ) {
		atoms.push_back( NamedAtomID( "3HE1", res ) );
	} else if ( name == "HE21" ) {
		atoms.push_back( NamedAtomID( "1HE2", res ) );
	} else if ( name == "HE22" ) {
		atoms.push_back( NamedAtomID( "2HE2", res ) );
	} else if ( name == "HH11" ) {
		atoms.push_back( NamedAtomID( "1HH1", res ) );
	} else if ( name == "HH12" ) {
		atoms.push_back( NamedAtomID( "2HH1", res ) );
	} else if ( name == "HH21" ) {
		atoms.push_back( NamedAtomID( "1HH2", res ) );
	} else if ( name == "HH22" ) {
		atoms.push_back( NamedAtomID( "2HH2", res ) );
	} else if ( (name == "HH1" || name == "QH1" ) && aa == aa_arg ) {
		atoms.push_back( NamedAtomID( "1HH1", res ) );
		atoms.push_back( NamedAtomID( "2HH1", res ) );
	} else if ( (name == "HH2" || name == "QH2" ) && aa == aa_arg ) {
		atoms.push_back( NamedAtomID( "1HH2", res ) );
		atoms.push_back( NamedAtomID( "2HH2", res ) );
	} else if ( (name == "HH" || name == "QQH" ) && aa == aa_arg ) {
		atoms.push_back( NamedAtomID( "1HH1", res ) );
		atoms.push_back( NamedAtomID( "2HH1", res ) );
		atoms.push_back( NamedAtomID( "1HH2", res ) );
		atoms.push_back( NamedAtomID( "2HH2", res ) );
	} else {
		atoms.push_back( NamedAtomID( name, res ) );
	}
}


void parse_NMR_name( std::string name, core::Size res, AmbiguousNMRDistanceConstraint::Atoms& atoms, core::pose::Pose const& pose ) {
	using core::chemical::AA;
	using namespace core::chemical;
	using core::id::NamedAtomID;
	using core::pose::named_atom_id_to_atom_id;
	runtime_assert( res >= 1 || res <= pose.total_residue() );
	AA const aa( pose.residue_type( res ).aa() );

	if ( ( name.substr(0,2) == "HD" ) && aa == aa_his ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "CD2", res ), pose ) );
		return;
	}
	if ( ( name.substr(0,2) == "HE" ) && aa == aa_his ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "CE1", res ), pose ) );
		return;
	}

	///now fix terminus problems
	if ( name == "H" && res == 1 && pose.is_fullatom() ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1H", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2H", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3H", res ), pose ) );
		return;
	}

	NamedAtoms named_atoms;
	parse_NMR_name( name, res, aa, named_atoms );

	for ( NamedAtoms::const_iterator it= named_atoms.begin(); it!=named_atoms.end(); ++it ) {
		atoms.push_back( named_atom_id_to_atom_id( *it, pose ) );
	}
	while ( atoms.size() && !atoms.back().valid() ) atoms.pop_back();
}


void parse_NMR_name_old( std::string name, core::Size res, AmbiguousNMRDistanceConstraint::Atoms& atoms, core::pose::Pose const& pose ) {
	using core::chemical::AA;
	using namespace core::chemical;
	using core::id::NamedAtomID;
	using core::pose::named_atom_id_to_atom_id;
	runtime_assert( res >= 1 || res <= pose.total_residue() );
	AA const aa( pose.residue_type( res ).aa() );
	//tr.Debug << "[ERROR]: name is " << name << " res " << res << " and name " << pose.residue( res ).name3() << std::endl;
	//use named_atom_id_to_atom_id because it can throw an exception if atoms are missing instead of hard-exit...

	//put histidine protons onto the the respective Carbon atom... it is just unpredictable if the respective H is present...
	if ( ( name.substr(0,2) == "HD" ) && aa == aa_his ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "CD2", res ), pose ) );
		return;
	}
	if ( ( name.substr(0,2) == "HE" ) && aa == aa_his ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "CE1", res ), pose ) );
		return;
	}

	///now fix all the problems with atomnames:
	if ( name == "H" && res == 1 && pose.is_fullatom() ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1H", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2H", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3H", res ), pose ) );
	} else if ( (name.substr(0,2) == "QA" || name == "HA") && ( aa == aa_gly ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HA", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HA", res ), pose ) );
	} else if ( ( name == "QB" || name == "HB" ) && (
			aa == aa_ala
			|| aa == aa_leu
			|| aa == aa_ser
			|| aa == aa_asn
			|| aa == aa_gln
			|| aa == aa_pro
			|| aa == aa_lys
			|| aa == aa_cys
			|| aa == aa_asp
			|| aa == aa_glu
			|| aa == aa_arg
			|| aa == aa_tyr
			|| aa == aa_phe
			|| aa == aa_trp
			|| aa == aa_his
			|| aa == aa_met ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HB", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HB", res ), pose ) );
		if ( aa == aa_ala ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HB", res ), pose ) );
	} else if ( ( name == "QD1" ||  name == "HD1" ) && (
			aa == aa_ile
			|| aa == aa_leu ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD1", res ), pose ) );
	} else if ( ( name == "QD2" || name == "HD2" ) && (
			aa == aa_leu
			|| aa == aa_asn ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD2", res ), pose ) );
		if ( aa == aa_leu ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD2", res ), pose ) );
	} else if ( ( name == "QQD" || name == "QD" ) && aa == aa_leu ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD1", res ), pose ) );
	} else if ( ( name == "QD" || name == "HD" ) && (
			aa == aa_pro
			|| aa == aa_lys
			|| aa == aa_arg
			|| aa == aa_tyr
			|| aa == aa_phe
			|| aa == aa_his ) ) {
		if ( aa == aa_arg || aa == aa_lys || aa == aa_pro ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD", res ), pose ) );
		} else { // tyr, phe, his (his doesn't get down here... )
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HD1", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HD2", res ), pose ) );
		}
	} else if ( ( name == "QE" || name == "HE" ) && (
			aa == aa_tyr
			|| aa == aa_phe
			|| aa == aa_trp
			|| aa == aa_met
			|| aa == aa_lys ) ) {
		if ( aa == aa_phe || aa == aa_tyr ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE1", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE2", res ), pose ) );
		} else if ( aa == aa_trp ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE1", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE3", res ), pose ) );
		} else if ( aa == aa_met || aa == aa_lys ) { //MET LYS
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE", res ), pose ) );
			if ( aa == aa_met ) {
				atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HE", res ), pose ) );
			}
		} //QE MET LYS
	}  else if ( ( name == "QG" || name == "HG" ) && (
			aa == aa_met
			|| aa == aa_gln
			|| aa == aa_pro
			|| aa == aa_lys
			|| aa == aa_glu
			|| aa == aa_arg ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG", res ), pose ) );
		//  atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HE", res ), pose ) );
	} else if ( ( name == "QG2" || name == "HG2" ) && (
			aa == aa_ile
			|| aa == aa_thr
			|| aa == aa_val ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG2", res ), pose ) );
	} else if ( ( name == "QQG" || name == "HG" ) && (
			aa == aa_ile
			|| aa == aa_val ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG1", res ), pose ) );
		if ( aa != aa_ile ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG1", res ), pose ) );
	} else if ( ( name == "QG1" || name == "HG1" ) && (
			aa == aa_ile
			|| aa == aa_val ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG1", res ), pose ) );
		if ( aa != aa_ile ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG1", res ), pose ) );
	} else if ( ( name == "QE2" || name == "HE2" || name =="HE" ) && aa == aa_gln ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE2", res ), pose ) );
	} else if ( ( name == "HZ" || name == "QZ" ) && ( aa == aa_trp || aa == aa_lys ) ) {
		if ( aa == aa_lys ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HZ", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HZ", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HZ", res ), pose ) );
		} else { //aa==trp
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HZ2", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HZ3", res ), pose ) );
		}
	} else  if ( name == "HB1" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HB", res ), pose ) );
	} else  if ( name == "HB2" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HB", res ), pose ) );
	} else  if ( name == "HB3" ) {
		if (  aa != aa_ala ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HB", res ), pose ) ); //yeah they call it 2HB and 3HB...
		} else {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HB", res ), pose ) );
		}
	} else  if ( name == "HD1" && !is_aromatic( pose, res ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD", res ), pose ) );
	} else  if ( name == "HD2" && !is_aromatic( pose, res ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD", res ), pose ) );
	} else  if ( name == "HD3" ) { //LYS, PRO, ARG  no other has HD3
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD", res ), pose ) );

	} else  if ( name == "HG1" && aa != aa_thr ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG", res ), pose ) );
	} else  if ( name == "HG2" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG", res ), pose ) );
	} else  if ( name == "HG3" ) { //GLU, ARG, GLN, MET
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG", res ), pose ) );

	} else  if ( name == "HA1" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HA", res ), pose ) );
	} else  if ( name == "HA2" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HA", res ), pose ) );
	} else  if ( name == "HA3" ) { //GLY
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HA", res ), pose ) );
	} else if ( ( name == "HE1" || name == "HE2" || name=="HE3" ) && !is_aromatic( pose, res ) ) {
		if ( name == "HE3" && aa != aa_met ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE", res ), pose ) ); //e.g. LYS
		} else {
			atoms.push_back( id::AtomID( pose.residue_type(res).atom_index(name.substr(2,1)+name.substr(0,2)), res ) );
		}
	} else if ( name == "HD11" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD1", res ), pose ) );
	} else if ( name == "HD12" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD1", res ), pose ) );
	} else  if ( name == "HD13" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD1", res ), pose ) );
	} else  if ( name == "HD21" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD2", res ), pose ) );
	} else if ( name == "HD22" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD2", res ), pose ) );
	} else  if ( name == "HD23" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD2", res ), pose ) );
	} else  if ( name == "HG11" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
	} else  if ( name == "HG12" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG1", res ), pose ) );
	} else if ( name == "HG13" ) {
		if ( aa == aa_ile ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
		} else {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG1", res ), pose ) );
		}
	} else  if ( name == "HG21" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG2", res ), pose ) );
	} else  if ( name == "HG22" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG2", res ), pose ) );
	} else if ( name == "HG23" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG2", res ), pose ) );
	} else if ( name == "HE11" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE1", res ), pose ) );
	} else if ( name == "HE12" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE1", res ), pose ) );
	} else if ( name == "HE13" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HE1", res ), pose ) );
	} else if ( name == "HE21" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE2", res ), pose ) );
	} else if ( name == "HE22" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE2", res ), pose ) );
	} else if ( name == "HH11" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HH1", res ), pose ) );
	} else if ( name == "HH12" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HH1", res ), pose ) );
	} else if ( name == "HH21" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HH2", res ), pose ) );
	} else if ( name == "HH22" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HH2", res ), pose ) );
	} else if ( (name == "HH1" || name == "QH1" ) && aa == aa_arg ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HH1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HH1", res ), pose ) );
	} else if ( (name == "HH2" || name == "QH2" ) && aa == aa_arg ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HH2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HH2", res ), pose ) );
	} else if ( (name == "HH" || name == "QQH" ) && aa == aa_arg ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HH1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HH1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HH2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HH2", res ), pose ) );
	} else {
		tr.Trace << "adding " << id::NamedAtomID( name, res ) << " "
			<< pose.residue_type( res ).name3() << std::endl;
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( name, res ), pose ) );
		tr.Trace << "as atom: " << atoms.back();
	}

	while ( atoms.size() && !atoms.back().valid() ) atoms.pop_back();
}

bool requires_CB_mapping( AmbiguousNMRDistanceConstraint::Atoms atoms, pose::Pose const& pose ) {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_REQUIRES_CB_MAPPING );
	using core::chemical::AA;
	using namespace core::chemical;
	Size ct_backbone( 0 );
	for ( Size i=1; i<=atoms.size(); ++i ) {
		chemical::ResidueType const& res_type( pose.residue_type( atoms[ i ].rsd() ) );
		ct_backbone += res_type.atom_is_backbone( atoms[ i ].atomno() );
	}
	chemical::ResidueType const& res_type( pose.residue_type( atoms[ 1 ].rsd() ) );
	ct_backbone -= ( ct_backbone == 1 && res_type.aa() != aa_gly && res_type.atom_index( "CB" ) == atoms[ 1 ].atomno() );
	if ( ct_backbone == 0 ) return true;
	if ( ct_backbone == atoms.size() ) return false;
	runtime_assert( 0 ); //this should not happen, as there should only be methyl-Hs combined into one constraint
	return false;
}


void combine_NMR_atom_string( AmbiguousNMRDistanceConstraint::Atoms atoms, std::string &atom_str, pose::Pose const& pose) {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NMR_STRING );
	atom_str = " BOGUS ";
	if ( atoms.size() == 1 ) {
		atom_str = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() );
	} else if ( atoms.size() == 2 && atoms[ 1 ].rsd() == atoms[ 2 ].rsd()
			&& pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 2 ) == "HB"
			&& pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 2 ) == "HB" ) {
		atom_str = "QB";
	} else if ( atoms.size() == 3
			&& atoms[ 1 ].rsd() == atoms[ 2 ].rsd() && atoms[ 1 ].rsd() == atoms[ 3 ].rsd() ) {
		std::string methyl = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 3 );
		if ( pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 3 ) == methyl
				&& pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ).substr( 1, 3 ) == methyl ) {
			if ( atoms[ 1 ].rsd() == 1 ) { //1H, 2H, 3H Nterminus
				atom_str = "H";
			} else {
				atom_str = "Q"+methyl.substr(1);
			}
		} else {
			atom_str = " untranslatable-3-atom-combi";
			tr.Trace << "could not translate: "
				<< pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ) << " "
				<< pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ) << " "
				<< pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ) << std::endl;
		}
	} else if ( atoms.size() == 2 //QE or QE2 or QD/PHE
			&& atoms[ 1 ].rsd() == atoms[ 2 ].rsd() ) {
		if ( pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr(0,2)==" H" ) {
			atom_str = "Q"+pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr(2,1);
		} else {
			std::string methyl = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 3 );
			if ( pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 3 ) == methyl ) {
				atom_str = "Q"+methyl.substr(1);
			} else {
				atom_str = " untranslatable-2-methyl-atom-combi ";
				tr.Trace << "could not translate: " << pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ) << " " << pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ) << std::endl;
			}
		}
	} else if ( atoms.size() == 6
			&& atoms[ 1 ].rsd() == atoms[ 2 ].rsd() && atoms[ 1 ].rsd() == atoms[ 3 ].rsd()
			&& atoms[ 1 ].rsd() == atoms[ 4 ].rsd() && atoms[ 1 ].rsd() == atoms[ 5 ].rsd()
			&& atoms[ 1 ].rsd() == atoms[ 6 ].rsd() ) {
		std::string methyl = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 2 );
		if ( pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 2 ) == methyl
				&& pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ).substr( 1, 2 ) == methyl ) {
			atom_str = "QQ"+methyl.substr(1);
		} else {
			atom_str = " untranslatable-6-atom-combi ";
			tr.Trace << "could not translate: "
				<< pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ) << " "
				<< pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ) << " "
				<< pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ) << " "
				<< pose.residue(  atoms[ 4 ].rsd() ).atom_name(  atoms[ 4 ].atomno() ) << " "
				<< pose.residue(  atoms[ 5 ].rsd() ).atom_name(  atoms[ 5 ].atomno() ) << " "
				<< pose.residue(  atoms[ 6 ].rsd() ).atom_name(  atoms[ 6 ].atomno() ) << std::endl;
		}
	}
	std::istringstream eat_white_space( atom_str );
	eat_white_space >> atom_str;
}


///c-tor
AmbiguousNMRDistanceConstraint::AmbiguousNMRDistanceConstraint(
	Atoms const & a1,
	Atoms const & a2,
	func::FuncOP func,
	ScoreType scoretype /* = atom_pair_constraint */
):
	Constraint( scoretype ),
	atoms1_(a1),
	atoms2_(a2),
	func_( func )
{}

AmbiguousNMRDistanceConstraint::AmbiguousNMRDistanceConstraint(
	id::NamedAtomID const & a1, //digests names like "QG1"
	id::NamedAtomID const & a2,
	core::pose::Pose const& pose,
	func::FuncOP func,
	ScoreType scoretype
) : Constraint( scoretype ),
	func_( func )
{
	parse_NMR_name( a1.atom(), a1.rsd(), atoms1_, pose );
	parse_NMR_name( a2.atom(), a2.rsd(), atoms2_, pose );

	if ( atoms1_.size() == 0 || atoms2_.size() == 0 ) {
		tr.Warning << "Error constructing from atoms: read in atom names("
			<< a1.atom() << "," << a2.atom() << "), " << std::endl;
	}
}

AmbiguousNMRDistanceConstraint::AmbiguousNMRDistanceConstraint() :
	Constraint( atom_pair_constraint ),
	func_( /* NULL */ )
{}

ConstraintOP AmbiguousNMRDistanceConstraint::clone() const {
	return ConstraintOP( new AmbiguousNMRDistanceConstraint( *this ));
}

///
ConstraintOP AmbiguousNMRDistanceConstraint::clone( func::FuncOP func ) const {
	return ConstraintOP( new AmbiguousNMRDistanceConstraint( atoms1_, atoms2_, func, score_type() ) );
}

bool AmbiguousNMRDistanceConstraint::operator == ( Constraint const & other ) const {
	if ( !       same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	AmbiguousNMRDistanceConstraint const & other_amb( static_cast< AmbiguousNMRDistanceConstraint const & > (other));
	if ( atoms1_ != other_amb.atoms1_ ) return false;
	if ( atoms2_ != other_amb.atoms2_ ) return false;

	return func_ == other_amb.func_ || ( func_ && other_amb.func_ && *func_ == *other_amb.func_ );
}

bool AmbiguousNMRDistanceConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< AmbiguousNMRDistanceConstraint const * > (&other);
}


/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP AmbiguousNMRDistanceConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	Atoms ids1, ids2;
	for ( Atoms::const_iterator it = atoms1_.begin(); it != atoms1_.end(); ++it ) {
		id::NamedAtomID atom(pose::atom_id_to_named_atom_id( *it, src ) );
		if ( smap ) {
			atom.rsd() = (*smap)[ it->rsd() ];
		}
		id::AtomID id1( core::pose::named_atom_id_to_atom_id(atom, dest ));
		if ( !id1.valid() ) return NULL;
		ids1.push_back( id1 );
	}
	for ( Atoms::const_iterator it = atoms2_.begin(); it != atoms2_.end(); ++it ) {
		id::NamedAtomID atom(atom_id_to_named_atom_id( *it, src ) );
		if ( smap ) {
			atom.rsd() = (*smap)[ it->rsd() ];
		}
		id::AtomID id2( core::pose::named_atom_id_to_atom_id( atom, dest ));
		if ( !id2.valid() ) return NULL;
		ids2.push_back( id2 );
	}
	return ConstraintOP( new AmbiguousNMRDistanceConstraint( ids1, ids2, func_, score_type() ) );
}

/// @details one line definition "AmbiguousNMRDistance atom1 res1 atom2 res2 function_type function_definition"
void
AmbiguousNMRDistanceConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & func_factory
) {
	Size res1, res2;
	std::string tempres1, tempres2;
	std::string name1, name2;
	std::string func_type;
	std::string type;

	data
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );

	tr.Debug << "read: " << name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
	if ( res1 > pose.total_residue() || res2 > pose.total_residue() ) {
		tr.Warning  << "ignored constraint (residue number to high for pose: " << pose.total_residue() << " !)"
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	parse_NMR_name( name1, res1, atoms1_, pose );
	parse_NMR_name( name2, res2, atoms2_, pose );

	if ( atoms1_.size() == 0 || atoms2_.size() == 0 ) {
		tr.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "), " << std::endl;
		//   << "and found AtomIDs (" << atom1_ << "," << atom2_ << ")" << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( data );

	if ( data.good() ) {
		//chu skip the rest of line since this is a single line defintion.
		while ( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}

	if ( tr.Debug.visible() ) {
		func_->show_definition( tr.Debug );
		tr.Debug << std::endl;
	}
} // read_def

void AmbiguousNMRDistanceConstraint::show( std::ostream& out ) const {
	out << "AmbiguousNMRDistanceConstraint (";
	for ( Atoms::const_iterator it = atoms1_.begin(); it != atoms1_.end(); ++it ) {
		out << it->atomno() << "," << it->rsd() << "||";
	}
	out << " - ";
	for ( Atoms::const_iterator it = atoms2_.begin(); it != atoms2_.end(); ++it ) {
		out << it->atomno() << "," << it->rsd();
	}
	func_->show( out );
}



void AmbiguousNMRDistanceConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	std::string def_atoms1;
	std::string def_atoms2;
	combine_NMR_atom_string( atoms1_, def_atoms1, pose);
	combine_NMR_atom_string( atoms2_, def_atoms2, pose);
	out << type() << " "
		<< def_atoms1 << " " << atoms1_[ 1 ].rsd() << " "
		<< def_atoms2 << " " << atoms2_[ 1 ].rsd() << " ";
	if ( func_ ) func_->show_definition( out );
	else out << std::endl;
}

ConstraintOP AmbiguousNMRDistanceConstraint::map_to_CEN( pose::Pose const& fa_pose, pose::Pose const& centroid, core::Size &mapped, std::string const& map_atom ) const {
	bool still_ambiguous( false );
	std::string atom1,atom2;
	basic::ProfileThis doit( basic::NOESY_ASSIGN_MAP2CB );
	using core::id::NamedAtomID;
	using core::pose::named_atom_id_to_atom_id;
	mapped = 2;
	if ( requires_CB_mapping( atoms1_, fa_pose) ) atom1 = map_atom;
	else {
		--mapped;
		combine_NMR_atom_string( atoms1_, atom1, fa_pose );
		still_ambiguous |= ( atoms1_.size() > 1 );
	}
	if ( requires_CB_mapping( atoms2_, fa_pose) ) atom2 = map_atom;
	else {
		--mapped;
		combine_NMR_atom_string( atoms2_, atom2, fa_pose );
		still_ambiguous |= ( atoms2_.size() > 1 );
	}
	tr.Debug << "map_to_CEN:  " << atom1 << " " << resid( 1 ) << " --> " << atom2 << " " << resid( 2 ) << (still_ambiguous ? " ambiguous " : " straight ") << std::endl;
	{ //scope for profile
		basic::ProfileThis doit( basic::NOESY_ASSIGN_MAP2CB_NEW );
		if ( still_ambiguous ) {
			return ConstraintOP( new AmbiguousNMRDistanceConstraint( id::NamedAtomID( atom1, resid( 1 ) ), id::NamedAtomID( atom2, resid( 2 ) ), centroid, func_, score_type() ) );
		} else {
			return ConstraintOP( new AtomPairConstraint( named_atom_id_to_atom_id( NamedAtomID( atom1, resid( 1 ) ), centroid ),
				named_atom_id_to_atom_id( NamedAtomID( atom2, resid( 2 ) ), centroid ), func_, score_type() ) );
		}
	} //scope
	return NULL; // cannot be reached
}

Real
AmbiguousNMRDistanceConstraint::dist( pose::Pose const & pose ) const {
	return dist( func::ConformationXYZ( pose.conformation() ) );
}

Real
AmbiguousNMRDistanceConstraint::dist(
	func::XYZ_Func const & xyz
) const
{
	return pow( inv_dist6( xyz ), -1.0/6 );
}

Real
AmbiguousNMRDistanceConstraint::inv_dist6(
	func::XYZ_Func const & xyz
) const
{
	Real cum_dist( 0.0 );
	for ( Atoms::const_iterator it1 = atoms1_.begin(); it1 != atoms1_.end(); ++it1 ) {
		for ( Atoms::const_iterator it2 = atoms2_.begin(); it2 != atoms2_.end(); ++it2 ) {
			Vector const & xyz1( xyz( *it1 ) ), xyz2( xyz( *it2 ) );
			Vector const f2( xyz1 - xyz2 );
			Real const dist( f2.length() );
			Real const inv_dist( 1.0/dist );
			Real const inv_dist2( inv_dist*inv_dist );
			Real const inv_dist6( inv_dist2 * inv_dist2 * inv_dist2 );
			cum_dist += inv_dist6;
			//   tr.Trace << *it1 << " " << *it2 << " " << dist << " " << pow(cum_dist, -1.0/6) << std::endl;
		}
	}
	// tr.Trace << "finished distance" << std::endl;
	return cum_dist;
}

void AmbiguousNMRDistanceConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const {
	Real eff_dist = pow( inv_dist6( xyz ), -1.0/6 );
	emap[ this->score_type() ] += func( eff_dist );
}

Size AmbiguousNMRDistanceConstraint::show_violations(
	std::ostream& out,
	pose::Pose const& pose,
	Size verbose_level,
	Real threshold
) const {

	if ( verbose_level > 80 ) {
		out << "\nAmbiguousNMRDistanceConstraint ( "
			<< pose.residue_type(atoms1_.front().rsd() ).atom_name( atoms1_.front().atomno() ) << " : "
			<< atoms1_.front().rsd() << " - "
			<< pose.residue_type(atoms2_.front().rsd() ).atom_name( atoms2_.front().atomno() ) << " : "
			<< atoms2_.front().rsd() << " ) ";
	}

	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}

// atom deriv
void
AmbiguousNMRDistanceConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{

	bool not_methyl_1 = std::find( atoms1_.begin() , atoms1_.end() , atom ) == atoms1_.end();
	bool not_methyl_2 = std::find( atoms2_.begin() , atoms2_.end() , atom ) == atoms2_.end();
	if ( not_methyl_1 && not_methyl_2 ) {
		//  tr.Trace << "atom " << atom << " not found " << std::endl;
		return;
	}

	Real eff_dist = dist( xyz );
	Real out_wderiv( weights[ this->score_type() ] * dfunc( eff_dist ));
	Real in_deriv = -1.0/6.0 * pow( eff_dist, 7.0 );
	// tr.Trace << "deriv for atom " << atom << eff_dist << " " << out_wderiv << " " << in_deriv << std::endl;

	Atoms const& the_other_atoms( not_methyl_1 ? atoms1_ : atoms2_ );
	// tr.Trace << "the_other_atoms: " << the_other_atoms.size() << " " << the_other_atoms.front() << std::endl;
	for ( Atoms::const_iterator it = the_other_atoms.begin(); it != the_other_atoms.end(); ++it ) {
		AtomID other_atom = *it;
		//  tr.Trace << "contribution from " << other_atom << " to " << atom << std::endl;
		Real rdist(0.0);
		Vector f1(0.0), f2(0.0);
		numeric::deriv::distance_f1_f2_deriv( xyz( atom ), xyz( other_atom ), rdist, f1, f2 );
		Real wderiv = -6.0*pow(rdist,-7.0) * in_deriv * out_wderiv ;
		//  tr.Trace << "wderiv " << wderiv << std::endl;
		F1 += wderiv * f1;
		F2 += wderiv * f2;

	}

}

ConstraintOP
AmbiguousNMRDistanceConstraint::remap_resid( core::id::SequenceMapping const &smap ) const
{
	// runtime_assert( 0 );
	Atoms ids1, ids2;
	for ( Atoms::const_iterator it = atoms1_.begin(); it != atoms1_.end(); ++it ) {
		id::AtomID atom( it->atomno(), smap[ it->rsd() ] );
		if ( !atom.valid() ) return NULL;
		ids1.push_back( atom );
	}
	for ( Atoms::const_iterator it = atoms2_.begin(); it != atoms2_.end(); ++it ) {
		id::AtomID atom( it->atomno(), smap[ it->rsd() ] );
		if ( !atom.valid() ) return NULL;
		ids2.push_back( atom );
	}
	return ConstraintOP( new AmbiguousNMRDistanceConstraint( ids1, ids2, func_, score_type() ) );
}

core::Size
AmbiguousNMRDistanceConstraint::effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp ) const {
	return sp.dist( resid( 1 ), resid( 2 ) );
}

void AmbiguousNMRDistanceConstraint::atoms1( Atoms const & setting ) {
	atoms1_ = setting;
}

void AmbiguousNMRDistanceConstraint::atoms2( Atoms const & setting ) {
	atoms2_ = setting;
}

AmbiguousNMRDistanceConstraint::AmbiguousNMRDistanceConstraint( AmbiguousNMRDistanceConstraint const & src ) :
	Constraint( src ),
	atoms1_( src.atoms1_ ),
	atoms2_( src.atoms2_ ),
	func_( src.func_ ? src.func_->clone() : src.func_ )
{}


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::AmbiguousNMRDistanceConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( atoms1_ ) ); // Atoms
	arc( CEREAL_NVP( atoms2_ ) ); // Atoms
	arc( CEREAL_NVP( func_ ) ); // func::FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::AmbiguousNMRDistanceConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( atoms1_ ); // Atoms
	arc( atoms2_ ); // Atoms
	arc( func_ ); // func::FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::AmbiguousNMRDistanceConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::AmbiguousNMRDistanceConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_AmbiguousNMRDistanceConstraint )
#endif // SERIALIZATION
