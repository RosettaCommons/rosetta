// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/GenericBondedPotential.cc
/// @brief  A "generic" (atom-type-only-based) torsional potential
/// @author Hahnbeom Park and Frank DiMaio


// Unit headers
#include <core/scoring/GenericBondedPotential.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/ResidueConnection.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/Residue.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <math.h>
#include <iostream>
#include <core/kinematics/Jump.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// boost
#include <boost/unordered_set.hpp>

using namespace ObjexxFCL::format;

namespace core {
namespace scoring {

static basic::Tracer TR( "core.scoring.GenericBondedPotential" );

// enum to help in parsing
enum ReadMode {
	rmNONE,
	rmATOM,
	rmBOND,
	rmANGLE,
	rmTORSION,
	rmIMPROPER
};

// pre allocated "null" params
const SpringParams GenericBondedPotential::null_sp_param;
const GenTorsionParams GenericBondedPotential::null_tors_param;

// convert a "concrete" bond type to a bin index
// "concrete" means not unknown "name" or "ringness"
core::Size
bin_from_bond(core::chemical::BondName bn, core::chemical::BondRingness br) {
	int bn_index = (int)bn;
	runtime_assert (bn_index != 0); // "unknown"

	if ( br == core::chemical::UnknownRingness ) {
		br = core::chemical::BondNotInRing;
	}
	int br_index = (int)br;
	//runtime_assert (br_index != 0); // "unknown"

	int idx = ((bn_index-1)<<1) + (br_index);
	return ((core::Size)idx);
}

core::Size
bin_from_bond(core::chemical::Bond const &b) {
	return bin_from_bond(b.bond_name(), b.ringness());
}

// helper class to convert string specification of bondorders to indices
//   subset of SMARTS chemical description language
// everything is remapped via Rosetta's "bond" type
//
// string key:
//    − single bond
//    = double bond
//    # triple bond
//    ∶ aromatic bond
//    @ any ring bond
//    ∼ any bond: wild card
//
// bools
//    ! Not
//    & high-priority and
//    , or
//    ; low-priority and
//
// implement (and use) a simple recursive descent parser
// EBNF:
//   expr ::= oterm [';' oterm]*
//   oterm ::= aterm [',' aterm]*
//   aterm ::= nterm ['&' nterm]*
//   nterm ::= ['!'] factor
//   factor ::= '-' | '=' | '#' | ':' | '@' | '~'
//
// note: this is more complex than it needs to be
// however, that is because the way i have coded it makes it easy to
// allow much more complex specifications based on additional properties
class BondOrderParser {
private:
	std::string toparse_;
	int ptr_;
	char sym_[4];

	// "lexer"
	void nextsym(void) {
		ptr_++;
		sym_[0] = 0;
		if ( ptr_ < (int)toparse_.length() ) {
			const char* temp = toparse_.c_str();
			sym_[0] = (char) temp[ptr_];
		}
	}

	bool accept(char s) {
		if ( sym_[0] == s ) {
			nextsym();
			return true;
		}
		return false;
	}

	bool expect(char s) {
		if ( accept(s) ) return true;
		utility_exit_with_message("BondOrderParser: encountered unexpected symbol parsing "+toparse_);
		return false;
	}

	utility::vector1<core::Size>
	factor() {
		using namespace core::chemical;
		utility::vector1<core::Size> retval;
		if ( accept('-') ) {
			retval.push_back( bin_from_bond(SingleBond,BondInRing) );
			retval.push_back( bin_from_bond(SingleBond,BondNotInRing) );
		} else if ( accept('=') ) {
			retval.push_back( bin_from_bond(DoubleBond,BondInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondNotInRing) );
		} else if ( accept('#') ) {
			retval.push_back( bin_from_bond(TripleBond,BondInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondNotInRing) );
		} else if ( accept(':') ) {
			retval.push_back( bin_from_bond(AromaticBond,BondInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondNotInRing) );
		} else if ( accept('@') ) {
			retval.push_back( bin_from_bond(SingleBond,BondInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondInRing) );
		} else if ( accept('~') ) {
			retval.push_back( bin_from_bond(SingleBond,BondInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondInRing) );
			retval.push_back( bin_from_bond(SingleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondNotInRing) );
		} else {
			utility_exit_with_message("BondOrderParser: encountered unexpected symbol "+std::string(1,sym_[0])+"parsing "+toparse_);
		}
		return retval; // never reached but compiler complains
	}

	utility::vector1<core::Size>
	invfactor() {
		using namespace core::chemical;
		utility::vector1<core::Size> retval;
		if ( accept('-') ) {
			retval.push_back( bin_from_bond(DoubleBond,BondInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondNotInRing) );
		} else if ( accept('=') ) {
			retval.push_back( bin_from_bond(SingleBond,BondInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondInRing) );
			retval.push_back( bin_from_bond(SingleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondNotInRing) );
		} else if ( accept('#') ) {
			retval.push_back( bin_from_bond(SingleBond,BondInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondInRing) );
			retval.push_back( bin_from_bond(SingleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondNotInRing) );
		} else if ( accept(':') ) {
			retval.push_back( bin_from_bond(SingleBond,BondInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondInRing) );
			retval.push_back( bin_from_bond(SingleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondNotInRing) );
		} else if ( accept('@') ) {
			retval.push_back( bin_from_bond(SingleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(DoubleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(TripleBond,BondNotInRing) );
			retval.push_back( bin_from_bond(AromaticBond,BondNotInRing) );
		} else if ( accept('~') ) {
			return retval;
		} else {
			utility_exit_with_message("BondOrderParser: encountered unexpected symbol parsing "+toparse_);
		}
		return retval; // never reached but compiler complains
	}

	utility::vector1<core::Size>
	vector_intersect(
		utility::vector1<core::Size> const &i1,
		utility::vector1<core::Size> const &i2
	) {
		utility::vector1<core::Size> i1s=i1, i2s=i2, i1and2;
		std::sort(i1s.begin(), i1s.end());
		std::sort(i2s.begin(), i2s.end());
		std::set_intersection(
			i1s.begin(), i1s.end(), i2s.begin(), i2s.end(),
			std::back_inserter(i1and2));
		return i1and2;
	}

	utility::vector1<core::Size>
	vector_union(
		utility::vector1<core::Size> const &i1,
		utility::vector1<core::Size> const &i2
	) {
		utility::vector1<core::Size> i1s=i1, i2s=i2, i1and2;
		std::sort(i1s.begin(), i1s.end());
		std::sort(i2s.begin(), i2s.end());
		std::set_union(
			i1s.begin(), i1s.end(), i2s.begin(), i2s.end(),
			std::back_inserter(i1and2));
		return i1and2;
	}

	utility::vector1<core::Size> nterm() {
		if ( accept('!') ) {
			return invfactor();
		} else {
			return factor();
		}
	}

	utility::vector1<core::Size> aterm() {
		utility::vector1<core::Size> retval = nterm();
		while ( accept('&') ) {
			utility::vector1<core::Size> term2 = nterm();
			retval = vector_intersect(retval,term2);
		}
		return retval;
	}

	utility::vector1<core::Size>  oterm() {
		utility::vector1<core::Size> retval = aterm();
		while ( accept(',') ) {
			utility::vector1<core::Size> term2 = aterm();
			retval = vector_union(retval,term2);
		}
		return retval;
	}

	utility::vector1<core::Size> expr() {
		nextsym();
		utility::vector1<core::Size> retval = oterm();
		while ( accept(';') ) {
			utility::vector1<core::Size> term2 = oterm();
			retval = vector_intersect(retval,term2);
		}
		expect('\0'); // null term
		return retval;
	}

public:
	utility::vector1<core::Size> parse (std::string input) {
		toparse_ = input;
		ptr_ = -1;
		sym_[0] = sym_[1] = sym_[2] = sym_[3] = 0;
		return expr();
	}
};

utility::vector1< core::Size >
bondorders_map(std::string bt) {
	BondOrderParser bop;
	return (bop.parse(bt));
}

//
//

Real
SpringParams::energy ( core::Real value ) const {
	Real d = (value-delta_);
	return (k_*d*d);
}

Real
SpringParams::deriv ( core::Real value ) const {
	Real d = (value-delta_);
	return (2*k_*d);
}

Real
GenTorsionParams::energy ( core::Real value ) const {
	Real arg = k1_ * (1+cos( 1*value - f1_ ));
	arg +=  k2_ * (1+cos( 2*value - f2_ ));
	arg +=  k3_ * (1+cos( 3*value - f3_ ));
	arg +=  k4_ * (1+cos( 4*value - f4_ ));
	arg +=  k6_ * (1+cos( 6*value ));
	if ( k1_ < 0 ) arg += -2.0*k1_;
	if ( k2_ < 0 ) arg += -2.0*k2_;
	if ( k3_ < 0 ) arg += -2.0*k3_;
	if ( k4_ < 0 ) arg += -2.0*k4_;
	if ( k6_ < 0 ) arg += -2.0*k6_;

	//TR << "E("<<value<<") = " << arg << " ("<<k1_<<","<<k2_<<","<<k3_<<","<<f1_<<","<<f2_<<","<<f3_<<")"<<std::endl;

	return arg + offset_;
}

void
GenTorsionParams::calculate_offset() {
	core::Size nk_nonzero( 0 );
	if ( std::abs(k1_) > 1e-6 ) nk_nonzero++;
	if ( std::abs(k2_) > 1e-6 ) nk_nonzero++;
	if ( std::abs(k3_) > 1e-6 ) nk_nonzero++;

	core::Real const deg2rad( numeric::constants::d::pi/180.0 );

	if ( nk_nonzero > 1 ) {
		// get it in numeric way...
		std::vector< core::Real > values;
		core::Real angle( -180.0 );
		while ( angle < 180.0 ) {
			values.push_back( energy( angle*deg2rad ) );
			angle += 0.25;
		}
		std::sort( values.begin(), values.end() );
		offset_ = -values[0]; // take minimum as offset
	}

	return;
}

Real
GenTorsionParams::deriv ( core::Real value ) const {
	Real arg = -k1_ * (sin( 1*value - f1_ ));
	arg +=  -2*k2_ * (sin( 2*value - f2_ ));
	arg +=  -3*k3_ * (sin( 3*value - f3_ ));
	arg +=  -4*k4_ * (sin( 4*value - f4_ ));
	arg +=  -6*k6_ * (sin( 6*value ));
	return arg;
}

Real
GenTorsionParams::get_params(std::string keyword) const {
	if ( keyword == "k1" ) {
		return k1_;
	} else if ( keyword == "k2" ) {
		return k2_;
	} else if ( keyword == "k3" ) {
		return k3_;
	} else if ( keyword == "k4" ) {
		return k4_;
	} else if ( keyword == "k6" ) {
		return k6_;
	} else if ( keyword == "f1" ) {
		return f1_;
	} else if ( keyword == "f2" ) {
		return f2_;
	} else if ( keyword == "f3" ) {
		return f3_;
	} else if ( keyword == "f4" ) {
		return f4_;
	} else {
		return 0;//TR << keyword << " is invalid!" << std::endl;
	}
}


SpringParams const &
GenericBondedPotential::lookup_bond_params( Size type1, Size type2 ) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] ) return null_sp_param;

	auto it = bond_lookup_.find( get_parameter_hash(0, type1, type2) );
	if ( it == bond_lookup_.end() ) it = bond_lookup_.find( get_parameter_hash(0, type1, 0) );
	if ( it == bond_lookup_.end() ) it = bond_lookup_.find( get_parameter_hash(0, 0, type2) );
	if ( it == bond_lookup_.end() ) it = bond_lookup_.find( get_parameter_hash(0, 0, 0) );

	return (it == bond_lookup_.end()) ? null_sp_param : bond_pot_[it->second];
}

SpringParams const &
GenericBondedPotential::lookup_angle_params(
	Size type1, Size type2, Size type3
) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3] ) return null_sp_param;

	auto it = angle_lookup_.find( get_parameter_hash(0, type1, type2, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, 0, type2, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, type1, type2, 0) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, type1, 0, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, 0, type2, 0) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, type1, 0, 0) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, 0, 0, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, 0, 0, 0) );

	return (it == angle_lookup_.end()) ? null_sp_param : angle_pot_[it->second];
}


GenTorsionParams const &
GenericBondedPotential::lookup_tors_params(
	core::chemical::Bond const &bt, Size type1, Size type2, Size type3, Size type4
) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3]|| !defined_atom_types_[type4] ) return null_tors_param;

	core::Size btidx = bin_from_bond(bt);

	//fd look up tgt with best multiplicity
	auto it = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, type3, type4) );

	auto it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, type3, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, type3, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, type3, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_pot_[it2->second].multiplicity() < tors_pot_[it->second].multiplicity()) ) it = it2;

	// Final sanity check... (this should never get triggered)
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, 0, 0, 0, 0) );

	return tors_pot_[it->second];
}

// note!  improper torsions are order-specific, me might need to shuffle
SpringParams const &
GenericBondedPotential::lookup_improper_params(
	Size type1, Size type2, Size type3, Size type4, int &idx2, int &idx3, int &idx4
) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3]|| !defined_atom_types_[type4] ) return null_sp_param;

	Size i2in=idx2, i3in=idx3, i4in=idx4;

	auto it = improper_lookup_.find( get_parameter_hash(0, type1, type2, type3, type4) );
	if ( it != improper_lookup_.end() ) {
		return improper_pot_[it->second];
	}
	it = improper_lookup_.find( get_parameter_hash(0, type1, type2, type4, type3) );
	if ( it != improper_lookup_.end() ) {
		idx3 = i4in; idx4 = i3in;
		return improper_pot_[it->second];
	}


	it = improper_lookup_.find( get_parameter_hash(0, type1, type3, type2, type4) );
	if ( it != improper_lookup_.end() ) {
		idx2 = i3in; idx3 = i2in;
		return improper_pot_[it->second];
	}
	it = improper_lookup_.find( get_parameter_hash(0, type1, type3, type4, type2) );
	if ( it != improper_lookup_.end() ) {
		idx2 = i3in; idx3 = i4in; idx4 = i2in;
		return improper_pot_[it->second];
	}

	it = improper_lookup_.find( get_parameter_hash(0, type1, type4, type2, type3) );
	if ( it != improper_lookup_.end() ) {
		idx2 = i4in; idx3 = i2in; idx4 = i3in;
		return improper_pot_[it->second];
	}
	it = improper_lookup_.find( get_parameter_hash(0, type1, type4, type3, type2) );
	if ( it != improper_lookup_.end() ) {
		idx2 = i4in; idx4 = i2in;
		return improper_pot_[it->second];
	}

	return null_sp_param;
}

// load database on creation
GenericBondedPotential::GenericBondedPotential() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string db_file = option[ score::gen_bonded_params_file ]();
	read_database( db_file );
	if ( option[ corrections::genpotential::set_torsion_params ].user() ) {
		modify_torsion_params_from_cmd_line();
	}
}

void
GenericBondedPotential::read_database(
	std::string filename
) {
	utility::io::izstream instream;
	basic::database::open( instream, filename );

	chemical::AtomTypeSetCOP ats = chemical::ChemicalManager::get_instance()->atom_type_set( "fa_standard" );
	core::Size ntypes = ats->n_atomtypes();
	defined_atom_types_.resize( ntypes, false );
	//runtime_assert( ntypes < 65536 ); // we pack 4 atom indices into 64 bits
	runtime_assert( ntypes < 4096 ); // should be safe?

	// map atom names to ATS indices - handle wildcards
	ReadMode read_mode = rmNONE;
	Size linenum( 0 );
	natomtypes = 0;
	std::string fileline, tag;
	// special treatment of "fallbacks", e.g. 'X' in the params file
	utility::vector1< core::Size > wildcard_vec(1, 0);

	while ( instream ) {
		getline(instream, fileline);
		std::istringstream linestream(fileline);
		linenum ++;

		if ( fileline.length() < 2 ) continue;

		linestream >> tag;
		if ( tag.compare("ATOM") == 0  ) {
			read_mode = rmATOM;
			continue;
		} else if ( tag.compare("BOND") == 0 ) {
			read_mode = rmBOND;
			continue;
		} else if ( tag.compare("ANGLE") == 0 ) {
			read_mode = rmANGLE;
			continue;
		} else if ( tag.compare("TORSION") == 0 ) {
			read_mode = rmTORSION;
			continue;
		} else if ( tag.compare("IMPROPER") == 0 ) {
			read_mode = rmIMPROPER;
			continue;
		} else if ( tag.substr(0,1) == "#" ) {
			continue; // comment or empty line
		}

		if ( read_mode == rmATOM ) {       // in "ATOM" block
			core::Size atm_idx = ats->atom_type_index(tag);
			defined_atom_types_[ atm_idx ] = true;
			TR.Debug << "AtomicIndex: " << atm_idx << " " << tag << std::endl;
			name_index_map[ tag ].push_back( atm_idx );
			std::string alt_tag;
			while ( linestream >> alt_tag ) {
				if ( alt_tag != tag ) name_index_map[ alt_tag ].push_back( atm_idx );
			}
			natomtypes++;
		} else if ( read_mode == rmBOND ) {  // in "BOND" block
			std::string atm1 = tag,  atm2;
			core::Real k1, delta;
			linestream >> atm2 >> k1 >> delta;

			utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
			utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
			core::Size multiplicity = indices1.size() * indices2.size();

			utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
			utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;

			bond_pot_.push_back( SpringParams( k1, delta, multiplicity ) );
			core::Size tgt_pot_idx = bond_pot_.size();

			for ( auto i1 : indices1_loop ) {
				for ( auto i2 : indices2_loop ) {
					int64_t hashval = get_parameter_hash( 0, i1, i2 );
					auto it = bond_lookup_.find( hashval );
					// use the MOST SPECIFIC potential
					if ( it == bond_lookup_.end() ) {
						bond_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
						hashval = get_parameter_hash( 0, i2, i1 );
						bond_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
					} else if ( multiplicity < bond_pot_[ it->second ].multiplicity() ) {
						it->second = tgt_pot_idx;
						hashval = get_parameter_hash( 0, i2, i1 );
						it = bond_lookup_.find( hashval );
						it->second = tgt_pot_idx;
						//bond_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );???
					}
				}
			}

		} else if ( read_mode == rmANGLE ) {  // in "ANGLE" block
			std::string atm1 = tag, atm2, atm3;
			core::Real k1, delta;
			linestream >> atm2 >> atm3 >> k1 >> delta;

			utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
			utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
			utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
			core::Size multiplicity = indices1.size() * indices2.size() * indices3.size();

			utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
			utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
			utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;

			angle_pot_.push_back( SpringParams( k1, delta, multiplicity ) );
			core::Size tgt_pot_idx = angle_pot_.size();

			for ( auto i1 : indices1_loop ) {
				for ( auto i2 : indices2_loop ) {
					for ( auto i3 : indices3_loop ) {
						int64_t hashval = get_parameter_hash( 0, i1, i2, i3 );
						auto it = angle_lookup_.find( hashval );
						// use the MOST SPECIFIC potential
						if ( it == angle_lookup_.end() ) {
							angle_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
							hashval = get_parameter_hash( 0, i3, i2, i1 );
							angle_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
						} else if ( multiplicity < angle_pot_[ it->second ].multiplicity() ) {
							it->second = tgt_pot_idx;
							hashval = get_parameter_hash( 0, i3, i2, i1 );
							it = angle_lookup_.find( hashval );
							it->second = tgt_pot_idx;
						}
					}
				}
			}

		} else if ( read_mode == rmTORSION ) { // in "TORSION" block
			std::string atm1=tag, atm2, bt, atm3, atm4;
			core::Real k1, k2, k3, k4, k6=0.0, f1=0.0, f2=0.0, f3=0.0, f4=0.0;
			std::string dummy, extra;

			linestream >> atm2 >> bt >> atm3 >> atm4 >> dummy >> k1 >> k2 >> k3 >> k4 >> extra;
			if ( extra == "k6" ) {
				linestream >> k6 >> extra;
			}
			if ( extra == "phase" ) {
				linestream >> f1 >> f2 >> f3 >> f4 >> extra;
			}

			utility::vector1< core::Size > indicesBT = bondorders_map(bt);
			utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
			utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
			utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
			utility::vector1< core::Size > const & indices4 = name_index_map[atm4];

			// multiplicity is used for priority ...
			//   the applicable torsion with the lowest multiplicity is selcted
			// we want bondtype to dominate ... if there is a more-specific bt, always prefer that
			//   even if atomtyping is less specific
			// however we need to ensure "fallbacks" aren't favored too much
			utility::vector1< core::Size > const & indicesX = name_index_map["X"];
			core::Size const multi_max( indicesX.size()/2 );
			core::Size multBT = indicesBT.size()*multi_max*multi_max*multi_max*multi_max;

			core::Size multiplicity = multBT + indices1.size() * indices2.size() * indices3.size() * indices4.size();

			utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
			utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
			utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;
			utility::vector1< core::Size > const &indices4_loop = (indices4.size() == natomtypes) ? wildcard_vec : indices4;

			std::string torsion_type_in = atm1+","+atm2+","+bt+","+atm3+","+atm4;
			GenTorsionParams params( k1, k2, k3, k4, f1, f2, f3, f4, multiplicity, torsion_type_in );
			params.k6( k6 );

			// flipping atom order reverses phase -- generate a new constraint if necessary
			tors_pot_.push_back( params );

			core::Size tgt_pot_idx1 = tors_pot_.size();
			core::Size tgt_pot_idx2 = tgt_pot_idx1;

			if ( f1 != 0 || f2 != 0 || f3 != 0 || f4 != 0 ) {
				GenTorsionParams params2( k1, k2, k3, k4, -f1, -f2, -f3, -f4, multiplicity, torsion_type_in);
				params2.k6( k6 );
				tors_pot_.push_back( params2 );
				tgt_pot_idx2 = tors_pot_.size();
			}

			for ( auto i0 : indicesBT ) {
				for ( auto i1 : indices1_loop ) {
					for ( auto i2 : indices2_loop ) {
						for ( auto i3 : indices3_loop ) {
							for ( auto i4 : indices4_loop ) {

								int64_t hashval = get_parameter_hash( i0, i1, i2, i3, i4 );
								auto it = tors_lookup_.find( hashval );

								// use the MOST SPECIFIC potential
								if ( it == tors_lookup_.end()  ) {
									tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx1) );
									hashval = get_parameter_hash( i0, i4, i3, i2, i1 );
									tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx2) );
								} else if ( multiplicity < tors_pot_[ it->second ].multiplicity() ) {
									it->second = tgt_pot_idx1;
									hashval = get_parameter_hash( i0, i4, i3, i2, i1 );
									it = tors_lookup_.find( hashval );
									it->second = tgt_pot_idx2;
								}
							}
						}
					}
				}
			}
		} else if ( read_mode == rmIMPROPER ) { // in "IMPROPER" block
			std::string atm1=tag, atm2, atm3, atm4;
			core::Real k1, delta;
			linestream >> atm2 >> atm3 >> atm4 >> k1 >> delta;

			utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
			utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
			utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
			utility::vector1< core::Size > const & indices4 = name_index_map[atm4];
			core::Size multiplicity = indices1.size() * indices2.size() * indices3.size() * indices4.size();

			improper_pot_.push_back( SpringParams( k1, delta, multiplicity ) );
			core::Size tgt_pot_idx = improper_pot_.size();

			for ( auto i1 : indices1 ) {
				for ( auto i2 : indices2 ) {
					for ( auto i3 : indices3 ) {
						for ( auto i4 : indices4 ) {
							int64_t hashval = get_parameter_hash( 0, i1, i2, i3, i4 );

							auto it = improper_lookup_.find( hashval );

							// use the MOST SPECIFIC potential
							if ( it == improper_lookup_.end()  ) {
								improper_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
							} else if ( multiplicity < improper_pot_[ it->second ].multiplicity() ) {
								it->second = tgt_pot_idx;
							}
						}
					}
				}
			}

		}
	} //end while

	TR << "Added total " << bond_pot_.size() << " bond parameters corresponding to "
		<< bond_lookup_.size() << " unique bonds." << std::endl;
	TR << "Added total " << angle_pot_.size() << " angle parameters corresponding to "
		<< angle_lookup_.size() << " unique angles." << std::endl;
	TR << "Added total " << tors_pot_.size() << " torsion parameters corresponding to "
		<< tors_lookup_.size() << " unique torsions." << std::endl;
	TR << "Added total " << improper_pot_.size() << " improper parameters corresponding to "
		<< improper_lookup_.size() << " unique impropers." << std::endl;

}

void
GenericBondedPotential::modify_torsion_params_from_cmd_line() {
	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;

	utility::vector1< std::string > const & mods( option[ corrections::genpotential::set_torsion_params ] );

	string const errmsg( "-corrections:genpotential:set_torsion_params format should be::"
		"-corrections:genpotential:set_torsion_params <set1>:<atom1>:<param1>:"
		"<setting1> <set2>:<atom2>:<param2>:<setting2> ...; for example: "
		"'-corrections:genpotential:set_torsion_params fa_standard/C*/CS/~/CS/C*/k1/0.0/k2/0.0/k3/0.077"
		"fa_standard/CD/CS/@/CS/CD/k1/0.435/k2/0.039/k3/0.070/k4/0.000' ");

	for ( string const & mod : mods ) {
		utility::vector1< std::string > tags = utility::string_split( mod, '/' );

		// ensure at least one tag provided
		if ( tags.size() < 8 ) utility_exit_with_message( errmsg );
		if ( tags[1] != "fa_standard" ) continue; //? might want to exit with error here as well...

		// ensure all four atoms exist
		for ( core::Size i=2; i<=6; ++i ) {
			if ( i==4 ) continue;
			if ( name_index_map.find( tags[i] ) == name_index_map.end() ) {
				utility_exit_with_message( errmsg + ". Nonexistent atomname: " + tags[i] );
			}
		}

		core::Real k1=0,k2=0,k3=0,k4=0;
		core::Real f1=0,f2=0,f3=0,f4=0;
		for ( core::Size i=7; i<tags.size(); i+=2 ) {
			if ( tags[i] == "k1" ) {
				k1 = double_of( tags[i+1] );
			} else if ( tags[i] == "k2" ) {
				k2 = double_of( tags[i+1] );
			} else if ( tags[i] == "k3" ) {
				k3 = double_of( tags[i+1] );
			} else if ( tags[i] == "k4" ) {
				k4 = double_of( tags[i+1] );
			} else if ( tags[i] == "f1" ) {
				f1 = double_of( tags[i+1] );
			} else if ( tags[i] == "f2" ) {
				f2 = double_of( tags[i+1] );
			} else if ( tags[i] == "f3" ) {
				f3 = double_of( tags[i+1] );
			} else if ( tags[i] == "f4" ) {
				f4 = double_of( tags[i+1] );
			} else {
				utility_exit_with_message( "Error parsing tag "+mod );
			}
		}

		modify_tors_params(
			tags[2], tags[3], tags[4], tags[5], tags[6],
			k1,k2,k3,k4,f1,f2,f3,f4
		);

		TR << "modify_torsion_params_from_cmd_line: setting "
			<< "fa_standard" << ' ' << tags[2] << " " << tags[3] << " " << tags[4] << " " << tags[5] << " " << tags[6]
			<< " k=" << k1 << '/' << k2 << '/' << k3 << '/' << k4
			<< " f=" << f1 << '/' << f2 << '/' << f3 << '/' << f4 << endl;
	}// end for
}


void
GenericBondedPotential::setup_for_scoring(
	pose::Pose & /*pose*/,
	scoring::ScoreFunction const &/*sfxn*/
) const
{
}

//////////////////////////////////////////////////////////////////////

// energies (1b)
void
GenericBondedPotential::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & /*pose*/,
	EnergyMap & emap
) const {
	utility::vector1< DerivVectorPair > dummy; // I wanna get rid of this...

	TR.Debug << "Score residue " << rsd.seqpos() << "." << rsd.name() << std::endl;

	// TODO: also add logic for excl_improper...
	//auto const & excl_info
	// ( static_cast< GenBondedExclInfo const & >( pose.data().get( pose::datacache::CacheableDataType::GEN_BONDED_EXCL_INFO ) ) );
	//ResidueExclParamsCOP rsd_excl_info = excl_info.get_residue_data( rsd.seqpos() );

	Real E_bond  = eval_residue_energy_and_derivative_bond( rsd, dummy );
	Real E_angle   = eval_residue_energy_and_derivative_angle( rsd, dummy );
	Real E_torsion = eval_residue_energy_and_derivative_torsion( rsd, dummy); //,  rsd_excl_info );
	Real E_improper= eval_residue_energy_and_derivative_improper( rsd, dummy );

	TR.Debug << "SUMM: seqpos/bond/angle/tors/improp: " << I(3,int(rsd.seqpos()))
		<< " " << F(8,3,E_bond) << " " << F(8,3,E_angle)
		<< " " << F(8,3,E_torsion) << " " << F(8,3,E_improper) << std::endl;
	emap[ gen_bonded ]  += E_bond + E_angle + E_torsion + E_improper;
	emap[ gen_bonded_bond ]  += E_bond;
	emap[ gen_bonded_angle ]   += E_angle;
	emap[ gen_bonded_torsion ] += E_torsion;
	emap[ gen_bonded_improper ] += E_improper;
}

// energies (2b)
void
GenericBondedPotential::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	EnergyMap & emap
) const {
	// sanity checks
	//if ( !rsd1.is_bonded(rsd2) ) return;

	bool res1first = rsd2.seqpos() > rsd1.seqpos();
	if ( rsd1.has_upper_connect() && rsd1.connected_residue_at_resconn( rsd1.type().upper_connect_id() ) == rsd2.seqpos() ) {
		res1first = true;
	} else if ( rsd2.has_upper_connect() && rsd2.connected_residue_at_resconn( rsd2.type().upper_connect_id() ) == rsd1.seqpos() ) {
		res1first = false;
	}

	TR.Debug << "Score residue pair " << rsd1.seqpos() << "." << rsd1.name()
		<< " " << rsd2.seqpos() << "." << rsd2.name() << std::endl;

	conformation::Residue const &resLower = res1first?rsd1:rsd2;
	conformation::Residue const &resUpper = res1first?rsd2:rsd1;
	//if ( resLower.has_variant_type(core::chemical::CUTPOINT_LOWER) ) return;
	//if ( resUpper.has_variant_type(core::chemical::CUTPOINT_UPPER) ) return;

	//auto const & excl_info
	// ( static_cast< GenBondedExclInfo const & >( pose.data().get( pose::datacache::CacheableDataType::GEN_BONDED_EXCL_INFO ) ) );

	// TODO: also add logic for excl_improper...
	//utility::vector1< ResidueExclParamsCOP > rsd_pair_excl_info; rsd_pair_excl_info.resize( 3 );
	//rsd_pair_excl_info[1] = excl_info.get_residue_data( resLower.seqpos() );
	//rsd_pair_excl_info[2] = excl_info.get_residue_data( resUpper.seqpos() );
	//rsd_pair_excl_info[3] = excl_info.get_residue_pair_data( resLower.seqpos(), resUpper.seqpos() );

	utility::vector1< DerivVectorPair > dummy1, dummy2;
	Real E_bond  = eval_residue_pair_energy_and_derivative_bond( resLower, resUpper, dummy1, dummy2 );
	Real E_angle   = eval_residue_pair_energy_and_derivative_angle( resLower, resUpper, dummy1, dummy2 );
	Real E_torsion = eval_residue_pair_energy_and_derivative_torsion( resLower, resUpper, dummy1, dummy2); //, rsd_pair_excl_info );
	Real E_improper= eval_residue_pair_energy_and_derivative_improper( resLower, resUpper, dummy1, dummy2 );

	TR.Debug << "SUMM: seqpos1/2/bond/ang/tors/improp: " << I(3,int(resLower.seqpos())) << " " << I(3,int(resUpper.seqpos()))
		<< " " << F(8,3,E_bond) << " " << F(8,3,E_angle)
		<< " " << F(8,3,E_torsion) << " " << F(8,3,E_improper) << std::endl;

	emap[ gen_bonded ]  += E_bond + E_angle + E_torsion + E_improper;
	emap[ gen_bonded_bond ]  += E_bond;
	emap[ gen_bonded_angle ]   += E_angle;
	emap[ gen_bonded_torsion ] += E_torsion;
	emap[ gen_bonded_improper ] += E_improper;

}

// derivatives (1b)
void
GenericBondedPotential::residue_derivatives(
	conformation::Residue const & rsd,
	pose::Pose const & /*pose*/,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {

	// TODO: also add logic for excl_improper...
	//auto const & excl_info
	// ( static_cast< GenBondedExclInfo const & >( pose.data().get( pose::datacache::CacheableDataType::GEN_BONDED_EXCL_INFO ) ) );

	//ResidueExclParamsCOP rsd_excl_info = excl_info.get_residue_data( rsd.seqpos() );

	eval_residue_energy_and_derivative_bond( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_bond ], true );
	eval_residue_energy_and_derivative_angle( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_angle ], true );
	//eval_residue_energy_and_derivative_torsion( rsd, atom_derivs, rsd_excl_info, weights[ gen_bonded ] + weights[ gen_bonded_torsion ], true );
	eval_residue_energy_and_derivative_torsion( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_torsion ], true );
	eval_residue_energy_and_derivative_improper( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_improper ], true );
}

// derivatives (2b)
void
GenericBondedPotential::residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs_r1,
	utility::vector1< DerivVectorPair > & atom_derivs_r2
) const {
	// sanity checks
	//if ( !rsd1.is_bonded(rsd2) ) return;

	bool res1first = rsd2.seqpos() > rsd1.seqpos();
	if ( rsd1.has_upper_connect() && rsd1.connected_residue_at_resconn( rsd1.type().upper_connect_id() ) == rsd2.seqpos() ) {
		res1first = true;
	} else if ( rsd2.has_upper_connect() && rsd2.connected_residue_at_resconn( rsd2.type().upper_connect_id() ) == rsd1.seqpos() ) {
		res1first = false;
	}

	conformation::Residue const &resLower = res1first?rsd1:rsd2;
	conformation::Residue const &resUpper = res1first?rsd2:rsd1;
	//if ( resLower.has_variant_type(core::chemical::CUTPOINT_LOWER) ) return;
	//if ( resUpper.has_variant_type(core::chemical::CUTPOINT_UPPER) ) return;

	utility::vector1< DerivVectorPair > & atom_derivs_rL  = res1first?atom_derivs_r1:atom_derivs_r2;
	utility::vector1< DerivVectorPair > & atom_derivs_rU  = res1first?atom_derivs_r2:atom_derivs_r1;

	//auto const & excl_info
	// ( static_cast< GenBondedExclInfo const & >( pose.data().get( pose::datacache::CacheableDataType::GEN_BONDED_EXCL_INFO ) ) );

	// TODO: also add logic for excl_improper...
	//utility::vector1< ResidueExclParamsCOP > rsd_pair_excl_info; rsd_pair_excl_info.resize( 3 );
	//rsd_pair_excl_info[1] = excl_info.get_residue_data( resLower.seqpos() );
	//rsd_pair_excl_info[2] = excl_info.get_residue_data( resUpper.seqpos() );
	//rsd_pair_excl_info[3] = excl_info.get_residue_pair_data( resLower.seqpos(), resUpper.seqpos() );

	eval_residue_pair_energy_and_derivative_bond(
		resLower, resUpper, atom_derivs_rL, atom_derivs_rU, weights[ gen_bonded ]+weights[ gen_bonded_bond ], true );
	eval_residue_pair_energy_and_derivative_angle(
		resLower, resUpper, atom_derivs_rL, atom_derivs_rU, weights[ gen_bonded ]+weights[ gen_bonded_angle ], true );
	//eval_residue_pair_energy_and_derivative_torsion(
	// resLower, resUpper, atom_derivs_rL, atom_derivs_rU, rsd_pair_excl_info, weights[ gen_bonded ]+weights[ gen_bonded_torsion ], true );
	eval_residue_pair_energy_and_derivative_torsion(
		resLower, resUpper, atom_derivs_rL, atom_derivs_rU, weights[ gen_bonded ]+weights[ gen_bonded_torsion ], true );
	eval_residue_pair_energy_and_derivative_improper(
		resLower, resUpper, atom_derivs_rL, atom_derivs_rU, weights[ gen_bonded ]+weights[ gen_bonded_improper ], true );
}

// bond-length (1b)
Real
GenericBondedPotential::eval_residue_energy_and_derivative_bond(
	conformation::Residue const & rsd,
	utility::vector1< DerivVectorPair > & atom_derivs,
	Real const weight,
	bool const calc_deriv
) const {
	Real totalscore = 0.0;

	// for each bond
	for ( Size atm_i=1; atm_i<=rsd.natoms(); ++atm_i ) {
		chemical::AtomIndices atm_nbrs = rsd.type().nbrs( atm_i );
		for ( Size j=1; j<=atm_nbrs.size(); ++j ) {
			Size atm_j = atm_nbrs[j];
			if ( atm_i >= atm_j ) continue; // only score each bond once -- use restype index to define ordering

			SpringParams const &param_ij = lookup_bond_params( rsd.atom_type_index(atm_i) , rsd.atom_type_index(atm_j) );
			if ( param_ij.is_null() ) continue;

			Vector const &xyz1( rsd.xyz( atm_i ) );
			Vector const &xyz2( rsd.xyz( atm_j ) );
			Real d = xyz1.distance( xyz2 );
			Real score = param_ij.energy( d );
			totalscore += score;

			if ( TR.Debug.visible() ) {
				std::string aname1 = rsd.atom_name( atm_i );
				std::string aname2 = rsd.atom_name( atm_j );
				TR.Debug << "BND1 seqpos/atmnames/k/d0/d/score: "
					<< " " << I(4,int(rsd.seqpos()))
					<< " " << A(9, aname1+"-"+aname2)
					<< " " << F(6,3,param_ij.k())
					<< " " << F(6,3,param_ij.delta())
					<< " " << F(6,3,d)
					<< " " << F(8,5,score)
					<< std::endl;
			}

			if ( calc_deriv ) {
				Real dE_dd = weight*param_ij.deriv( d );
				Vector f1( 0.0 ), f2( 0.0 );
				numeric::deriv::distance_f1_f2_deriv( xyz1, xyz2, d, f1, f2 );
				atom_derivs[ atm_i ].f1() += dE_dd * f1;
				atom_derivs[ atm_i ].f2() += dE_dd * f2;

				atom_derivs[ atm_j ].f1() -= dE_dd * f1;
				atom_derivs[ atm_j ].f2() -= dE_dd * f2;
			}
		}
	}

	return totalscore;
}

// bond-length (2b)
Real
GenericBondedPotential::eval_residue_pair_energy_and_derivative_bond(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< DerivVectorPair > & atom_derivs_r1,
	utility::vector1< DerivVectorPair > & atom_derivs_r2,
	Real const weight,
	bool const calc_deriv
) const {
	Real totalscore = 0.0;

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

	//TR.Debug << "nbonds: " << rsd1.seqpos() << " " << rsd2.seqpos()
	//<< " " << r1_resconn_ids.size() << std::endl;
	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );
		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		SpringParams const &param_ij = lookup_bond_params(
			rsd1.atom_type_index(resconn_atomno1) , rsd2.atom_type_index(resconn_atomno2) );

		/*
		TR.Debug << "connected atoms? "
		<< rsd1.atom_name( resconn_atomno1 ) << "-" << rsd2.atom_name( resconn_atomno2 )
		<< ": " << rsd1.atom_type(resconn_atomno1).atom_type_name()
		<< "-" << rsd2.atom_type(resconn_atomno2).atom_type_name()
		<< " " << !param_ij.is_null()
		<< std::endl;
		*/

		if ( param_ij.is_null() ) continue;

		Vector const &xyz1( rsd1.xyz( resconn_atomno1 ) );
		Vector const &xyz2( rsd2.xyz( resconn_atomno2 ) );
		Real d = xyz1.distance( xyz2 );
		Real score = param_ij.energy( d );
		totalscore += score;

		if ( TR.Debug.visible() ) {
			std::string aname1 = rsd1.atom_name( resconn_atomno1 );
			std::string aname2 = rsd2.atom_name( resconn_atomno2 );
			TR.Debug << "BND2 seqpos1,2/atmnames/k/d0/d/score: "
				<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
				<< " " << A(9, aname1+"-"+aname2)
				<< " " << F(6,3,param_ij.k())
				<< " " << F(6,3,param_ij.delta())
				<< " " << F(6,3,d)
				<< " " << F(8,5,score)
				<< std::endl;
		}

		if ( calc_deriv ) {
			Real dE_dd = weight*param_ij.deriv( d );
			Vector f1( 0.0 ), f2( 0.0 );
			numeric::deriv::distance_f1_f2_deriv( xyz1, xyz2, d, f1, f2 );
			atom_derivs_r1[ resconn_atomno1 ].f1() += dE_dd * f1;
			atom_derivs_r1[ resconn_atomno1 ].f2() += dE_dd * f2;

			atom_derivs_r2[ resconn_atomno2 ].f1() -= dE_dd * f1;
			atom_derivs_r2[ resconn_atomno2 ].f2() -= dE_dd * f2;
		}
	}

	return totalscore;
}

// bond-angle (1b)
Real
GenericBondedPotential::eval_residue_energy_and_derivative_angle(
	conformation::Residue const & rsd,
	utility::vector1< DerivVectorPair > & atom_derivs,
	Real const weight,
	bool const calc_deriv
) const {
	Real totalscore = 0.0;

	// for each angle in the residue
	for ( Size bondang = 1; bondang <= rsd.type().num_bondangles(); ++bondang ) {
		Size atm_i = ( rsd.type().bondangle( bondang ) ).key1();
		Size atm_j = ( rsd.type().bondangle( bondang ) ).key2();
		Size atm_k = ( rsd.type().bondangle( bondang ) ).key3();

		SpringParams const &param_ij = lookup_angle_params(
			rsd.atom_type_index(atm_i) , rsd.atom_type_index(atm_j) , rsd.atom_type_index(atm_k) );
		if ( param_ij.is_null() ) continue;

		Vector const &xyz1( rsd.xyz( atm_i ) );
		Vector const &xyz2( rsd.xyz( atm_j ) );
		Vector const &xyz3( rsd.xyz( atm_k ) );
		Real angle = numeric::angle_radians( xyz1, xyz2, xyz3 );
		Real score = param_ij.energy( angle );
		totalscore += score;

		if ( TR.Debug.visible() ) {
			std::string aname1 = rsd.atom_name( atm_i );
			std::string aname2 = rsd.atom_name( atm_j );
			std::string aname3 = rsd.atom_name( atm_k );
			TR.Debug << "ANG1 seqpos/atmnames/k/d0/d/score: "
				<< " " << I(4,int(rsd.seqpos()))
				<< " " << A(13, aname1+"-"+aname2+"-"+aname3)
				<< " " << F(6,3,param_ij.k())
				<< " " << F(6,3,param_ij.delta())
				<< " " << F(6,1,angle*180.0/3.14159216)
				<< " " << F(8,5,score)
				<< std::endl;
		}

		if ( calc_deriv ) {
			Real dE_dang = weight*param_ij.deriv( angle );

			Vector f1( 0.0 ), f2( 0.0 );
			numeric::deriv::angle_p1_deriv( xyz1, xyz2, xyz3, angle, f1, f2 );
			atom_derivs[ atm_i ].f1() += dE_dang * f1;
			atom_derivs[ atm_i ].f2() += dE_dang * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::angle_p2_deriv( xyz1, xyz2, xyz3, angle, f1, f2 );
			atom_derivs[ atm_j ].f1() += dE_dang * f1;
			atom_derivs[ atm_j ].f2() += dE_dang * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::angle_p1_deriv( xyz3, xyz2, xyz1, angle, f1, f2 );
			atom_derivs[ atm_k ].f1() += dE_dang * f1;
			atom_derivs[ atm_k ].f2() += dE_dang * f2;
		}
	}

	return totalscore;
}

// bond-angle (2b)
Real
GenericBondedPotential::eval_residue_pair_energy_and_derivative_angle(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< DerivVectorPair > & atom_derivs_r1,
	utility::vector1< DerivVectorPair > & atom_derivs_r2,
	Real const weight,
	bool const calc_deriv
) const {
	Real totalscore = 0.0;

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );
		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		// (a) compute bond-angle energies with two atoms on rsd1
		utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
			rsd1.type().atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
		for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
			Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();

			// angle is r1.res1_lower_atomno / r1.resconn_atomno1 / r2.resconn_atomno2
			SpringParams const &param_ij = lookup_angle_params(
				rsd1.atom_type_index(res1_lower_atomno) ,
				rsd1.atom_type_index(resconn_atomno1) ,
				rsd2.atom_type_index(resconn_atomno2) );

			if ( param_ij.is_null() ) continue;

			Vector const &xyz1( rsd1.xyz( res1_lower_atomno ) );
			Vector const &xyz2( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz3( rsd2.xyz( resconn_atomno2 ) );
			Real angle = numeric::angle_radians( xyz1, xyz2, xyz3 );
			Real score = param_ij.energy( angle );
			totalscore += score;

			if ( TR.Debug.visible() ) {
				std::string aname1 = rsd1.atom_name( res1_lower_atomno );
				std::string aname2 = rsd1.atom_name( resconn_atomno1 );
				std::string aname3 = rsd2.atom_name( resconn_atomno2 );
				TR.Debug << "ANG2 seqpos1,2/atmnames/k/d0/d/score: "
					<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
					<< " " << A(14, aname1+"-"+aname2+"-"+aname3)
					<< " " << F(6,3,param_ij.k())
					<< " " << F(6,3,param_ij.delta())
					<< " " << F(6,1,angle*180.0/3.14159216)
					<< " " << F(8,5,score)
					<< " ?" << param_ij.is_null()
					<< std::endl;
			}

			if ( calc_deriv ) {
				Real dE_dang = weight*param_ij.deriv( angle );

				Vector f1( 0.0 ), f2( 0.0 );
				numeric::deriv::angle_p1_deriv( xyz1, xyz2, xyz3, angle, f1, f2 );
				atom_derivs_r1[ res1_lower_atomno ].f1() += dE_dang * f1;
				atom_derivs_r1[ res1_lower_atomno ].f2() += dE_dang * f2;

				f1 = f2 = Vector(0.0);
				numeric::deriv::angle_p2_deriv( xyz1, xyz2, xyz3, angle, f1, f2 );
				atom_derivs_r1[ resconn_atomno1 ].f1() += dE_dang * f1;
				atom_derivs_r1[ resconn_atomno1 ].f2() += dE_dang * f2;

				f1 = f2 = Vector(0.0);
				numeric::deriv::angle_p1_deriv( xyz3, xyz2, xyz1, angle, f1, f2 );
				atom_derivs_r2[ resconn_atomno2 ].f1() += dE_dang * f1;
				atom_derivs_r2[ resconn_atomno2 ].f2() += dE_dang * f2;
			}
		}

		// (b) compute bond-angle energies with two atoms on rsd2
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2.type().atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
		for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
			Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

			// angle is r1.resconn_atomno1 / r2.resconn_atomno2 / r2.res2_lower_atomno
			SpringParams const &param_ij = lookup_angle_params(
				rsd1.atom_type_index(resconn_atomno1) ,
				rsd2.atom_type_index(resconn_atomno2) ,
				rsd2.atom_type_index(res2_lower_atomno) );

			if ( param_ij.is_null() ) continue;

			Vector const &xyz1( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz2( rsd2.xyz( resconn_atomno2 ) );
			Vector const &xyz3( rsd2.xyz( res2_lower_atomno ) );
			Real angle = numeric::angle_radians( xyz1, xyz2, xyz3 );
			Real score = param_ij.energy( angle );
			totalscore += score;

			if ( TR.Debug.visible() ) {
				std::string aname1 = rsd1.atom_name( resconn_atomno1 );
				std::string aname2 = rsd2.atom_name( resconn_atomno2 );
				std::string aname3 = rsd2.atom_name( res2_lower_atomno );
				TR.Debug << "ANG2 seqpos1,2/atmnames/k/d0/d/score: "
					<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
					<< " " << A(14, aname1+"-"+aname2+"-"+aname3)
					<< " " << F(6,3,param_ij.k())
					<< " " << F(6,3,param_ij.delta())
					<< " " << F(6,1,angle*180.0/3.14159216)
					<< " " << F(8,5,score)
					<< " ?" << param_ij.is_null()
					<< std::endl;
			}

			if ( calc_deriv ) {
				Real dE_dang = weight*param_ij.deriv( angle );

				Vector f1( 0.0 ), f2( 0.0 );
				numeric::deriv::angle_p1_deriv( xyz1, xyz2, xyz3, angle, f1, f2 );
				atom_derivs_r1[ resconn_atomno1 ].f1() += dE_dang * f1;
				atom_derivs_r1[ resconn_atomno1 ].f2() += dE_dang * f2;

				f1 = f2 = Vector(0.0);
				numeric::deriv::angle_p2_deriv( xyz1, xyz2, xyz3, angle, f1, f2 );
				atom_derivs_r2[ resconn_atomno2 ].f1() += dE_dang * f1;
				atom_derivs_r2[ resconn_atomno2 ].f2() += dE_dang * f2;

				f1 = f2 = Vector(0.0);
				numeric::deriv::angle_p1_deriv( xyz3, xyz2, xyz1, angle, f1, f2 );
				atom_derivs_r2[ res2_lower_atomno ].f1() += dE_dang * f1;
				atom_derivs_r2[ res2_lower_atomno ].f2() += dE_dang * f2;
			}
		}
	}

	return totalscore;
}

// bond-torsion (1b)
Real
GenericBondedPotential::eval_residue_energy_and_derivative_torsion(
	conformation::Residue const & rsd,
	utility::vector1< DerivVectorPair > & atom_derivs,
	//ResidueExclParamsCOP excl_params,
	Real const weight,
	bool const calc_deriv
) const {

	Real totalscore = 0.0;

	// for each dihedral angle
	for ( Size dihe = 1; dihe <= rsd.type().ndihe(); ++dihe ) {
		Size atm_i = ( rsd.type().dihedral( dihe ) ).key1();
		Size atm_j = ( rsd.type().dihedral( dihe ) ).key2();
		Size atm_k = ( rsd.type().dihedral( dihe ) ).key3();
		Size atm_l = ( rsd.type().dihedral( dihe ) ).key4();

		core::chemical::VD vd_j = rsd.type().atom_vertex(atm_j);
		core::chemical::VD vd_k = rsd.type().atom_vertex(atm_k);
		GenTorsionParams const &t_param_ij = lookup_tors_params(
			rsd.type().bond( vd_j, vd_k ),
			rsd.atom_type_index(atm_i), rsd.atom_type_index(atm_j),
			rsd.atom_type_index(atm_k), rsd.atom_type_index(atm_l) );

		if ( t_param_ij.is_null() ) continue;

		Vector const &xyz1( rsd.xyz( atm_i ) );
		Vector const &xyz2( rsd.xyz( atm_j ) );
		Vector const &xyz3( rsd.xyz( atm_k ) );
		Vector const &xyz4( rsd.xyz( atm_l ) );
		Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
		Real score = t_param_ij.energy( angle );
		totalscore += score;

		if ( TR.Debug.visible() ) {
			std::string aname1, aname2, aname3, aname4;
			aname1 = rsd.atom_name( atm_i );
			aname2 = rsd.atom_name( atm_j );
			aname3 = rsd.atom_name( atm_k );
			aname4 = rsd.atom_name( atm_l );

			TR.Debug << "TOR1 torsion_type/seqpos/atmnames/k123/angle/score: "
				<< " " << A(19, t_param_ij.torsion_type())
				<< " " << I(4,int(rsd.seqpos()))
				<< " " << A(19, aname1+"-"+aname2+"-"+aname3+"-"+aname4)
				<< " " << F(6,3,t_param_ij.k(1)) << "," << F(6,3,t_param_ij.k(2))
				<< "," << F(6,3,t_param_ij.k(3))<< "," << F(6,3,t_param_ij.k(4));
			if ( std::abs(t_param_ij.k(6)) > 1.0e-3 ) TR.Debug << "," << t_param_ij.k(6);
			TR.Debug << " " << F(6,1,angle*180.0/3.14159216)
				<< " " << F(6,2,score)
				<< std::endl;
		}

		if ( calc_deriv ) {
			Real dE_dtors = weight*t_param_ij.deriv( angle );

			Vector f1( 0.0 ), f2( 0.0 );
			numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
			atom_derivs[ atm_i ].f1() += dE_dtors * f1;
			atom_derivs[ atm_i ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
			atom_derivs[ atm_j ].f1() += dE_dtors * f1;
			atom_derivs[ atm_j ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
			atom_derivs[ atm_k ].f1() += dE_dtors * f1;
			atom_derivs[ atm_k ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
			atom_derivs[ atm_l ].f1() += dE_dtors * f1;
			atom_derivs[ atm_l ].f2() += dE_dtors * f2;
		}
	}

	return totalscore;
}

// bond-torsion (2b)
Real
GenericBondedPotential::eval_residue_pair_energy_and_derivative_torsion(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< DerivVectorPair > & atom_derivs_r1,
	utility::vector1< DerivVectorPair > & atom_derivs_r2,
	Real const weight,
	bool const calc_deriv
) const {
	Real totalscore = 0.0;

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		/// 1. iterate over all triples on residue 1 within 2 bonds of resconn_atomno1.
		utility::vector1< chemical::three_atom_set > const & rsd1_atoms_wi2_bonds_of_ii(
			rsd1.type().atoms_within_two_bonds_of_a_residue_connection( resconn_id1 ));

		for ( Size jj = 1; jj <= rsd1_atoms_wi2_bonds_of_ii.size(); ++jj ) {
			debug_assert( rsd1_atoms_wi2_bonds_of_ii[ jj ].key1() == resconn_atomno1 );

			Size const jj_atom2 = rsd1_atoms_wi2_bonds_of_ii[ jj ].key2();
			Size const jj_atom1 = rsd1_atoms_wi2_bonds_of_ii[ jj ].key3();

			// torsion is r1.jj_atom1 / r1.jj_atom2 / r1.resconn_atomno1 / r2.resconn_atomno2
			core::chemical::VD vd_j = rsd1.type().atom_vertex(jj_atom2);
			core::chemical::VD vd_k = rsd1.type().atom_vertex(resconn_atomno1);
			GenTorsionParams const &t_param_ij = lookup_tors_params(
				rsd1.type().bond( vd_j, vd_k ),
				rsd1.atom_type_index(jj_atom1) ,
				rsd1.atom_type_index(jj_atom2) ,
				rsd1.atom_type_index(resconn_atomno1),
				rsd2.atom_type_index(resconn_atomno2)
			);
			if ( t_param_ij.is_null() ) continue;

			Vector const &xyz1( rsd1.xyz( jj_atom1 ) );
			Vector const &xyz2( rsd1.xyz( jj_atom2 ) );
			Vector const &xyz3( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz4( rsd2.xyz( resconn_atomno2 ) );
			Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
			Real score = t_param_ij.energy( angle );
			totalscore += score;

			if ( TR.Debug.visible() ) {
				std::string aname1, aname2, aname3, aname4;
				aname1 = rsd1.atom_name( jj_atom1 );
				aname2 = rsd1.atom_name( jj_atom2 );
				aname3 = rsd1.atom_name( resconn_atomno1 );
				aname4 = rsd2.atom_name( resconn_atomno2 );

				TR.Debug << "TOR2 seqpos1,2/atmnames/k123/angle/score: "
					<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
					<< " " << A(19, aname1+"-"+aname2+"-"+aname3+"-"+aname4)
					<< " " << F(6,3,t_param_ij.k(1)) << "," << F(6,3,t_param_ij.k(2))
					<< "," << F(6,3,t_param_ij.k(3))<< "," << F(6,3,t_param_ij.k(4))
					<< " " << F(6,1,angle*180.0/3.14159216)
					<< " " << F(6,2,score)
					<< std::endl;
			}

			if ( calc_deriv ) {
				Real dE_dtors = weight*t_param_ij.deriv( angle );

				Vector f1( 0.0 ), f2( 0.0 );
				numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				atom_derivs_r1[ jj_atom1 ].f1() += dE_dtors * f1;
				atom_derivs_r1[ jj_atom1 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				atom_derivs_r1[ jj_atom2 ].f1() += dE_dtors * f1;
				atom_derivs_r1[ jj_atom2 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				atom_derivs_r1[ resconn_atomno1 ].f1() += dE_dtors * f1;
				atom_derivs_r1[ resconn_atomno1 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				atom_derivs_r2[ resconn_atomno2 ].f1() += dE_dtors * f1;
				atom_derivs_r2[ resconn_atomno2 ].f2() += dE_dtors * f2;
			}
		}

		/// 2. iterate over all triples on residue 2 within 2 bonds of resconn_atomno2.
		utility::vector1< chemical::three_atom_set > const & rsd2_atoms_wi2_bonds_of_ii(
			rsd2.type().atoms_within_two_bonds_of_a_residue_connection( resconn_id2 ));

		for ( Size jj = 1; jj <= rsd2_atoms_wi2_bonds_of_ii.size(); ++jj ) {
			debug_assert( rsd2_atoms_wi2_bonds_of_ii[ jj ].key1() == resconn_atomno2 );
			Size const jj_atom3 = rsd2_atoms_wi2_bonds_of_ii[ jj ].key2();
			Size const jj_atom4 = rsd2_atoms_wi2_bonds_of_ii[ jj ].key3();

			core::chemical::VD vd_j = rsd2.type().atom_vertex(resconn_atomno2);
			core::chemical::VD vd_k = rsd2.type().atom_vertex(jj_atom3);
			// torsion is r1.resconn_atomno1 / r2.resconn_atomno2 / r1.jj_atom3 / r2.jj_atom4
			GenTorsionParams const &t_param_ij = lookup_tors_params(
				rsd2.type().bond( vd_j, vd_k ),
				rsd1.atom_type_index(resconn_atomno1) ,
				rsd2.atom_type_index(resconn_atomno2) ,
				rsd2.atom_type_index(jj_atom3),
				rsd2.atom_type_index(jj_atom4)
			);
			if ( t_param_ij.is_null() ) continue;

			Vector const &xyz1( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz2( rsd2.xyz( resconn_atomno2 ) );
			Vector const &xyz3( rsd2.xyz( jj_atom3 ) );
			Vector const &xyz4( rsd2.xyz( jj_atom4 ) );
			Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
			Real score = t_param_ij.energy( angle );
			totalscore += score;

			if ( TR.Debug.visible() ) {
				std::string aname1, aname2, aname3, aname4;
				aname1 = rsd1.atom_name( resconn_atomno1 );
				aname2 = rsd2.atom_name( resconn_atomno2 );
				aname3 = rsd2.atom_name( jj_atom3 );
				aname4 = rsd2.atom_name( jj_atom4 );

				TR.Debug << "TOR2 seqpos1,2/atmnames/k123/angle/score: "
					<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
					<< " " << A(19, aname1+"-"+aname2+"-"+aname3+"-"+aname4)
					<< " " << F(6,3,t_param_ij.k(1)) << "," << F(6,3,t_param_ij.k(2))
					<< "," << F(6,3,t_param_ij.k(3))<< "," << F(6,3,t_param_ij.k(4))
					<< " " << F(6,1,angle*180.0/3.14159216)
					<< " " << F(6,2,score)
					<< std::endl;
			}

			if ( calc_deriv ) {
				Real dE_dtors = weight*t_param_ij.deriv( angle );

				Vector f1( 0.0 ), f2( 0.0 );
				numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				atom_derivs_r1[ resconn_atomno1 ].f1() += dE_dtors * f1;
				atom_derivs_r1[ resconn_atomno1 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				atom_derivs_r2[ resconn_atomno2 ].f1() += dE_dtors * f1;
				atom_derivs_r2[ resconn_atomno2 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				atom_derivs_r2[ jj_atom3 ].f1() += dE_dtors * f1;
				atom_derivs_r2[ jj_atom3 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				atom_derivs_r2[ jj_atom4 ].f1() += dE_dtors * f1;
				atom_derivs_r2[ jj_atom4 ].f2() += dE_dtors * f2;
			}
		}

		/// 3. iterate over all pairs of pairs within 1 bond of either residue connection atom.
		utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
			rsd1.type().atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2.type().atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));

		for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
			Size const jj_term_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();
			for ( Size kk = 1; kk <= rsd2_atoms_wi1_bond_of_ii.size(); ++kk ) {
				Size const kk_term_atomno = rsd2_atoms_wi1_bond_of_ii[ kk ].key2();

				// no "Bond" associated with inter-residue connections
				//   assume double bond non-ring for protein, single bond non-ring for all others
				//   if we ever want the concept of "multi-residue ligands" this will need to be addressed
				core::chemical::Bond tempbond(
					1.0, rsd1.is_protein() ? core::chemical::DoubleBond : core::chemical::SingleBond);
				tempbond.ringness(core::chemical::BondNotInRing);

				GenTorsionParams const &t_param_ij = lookup_tors_params(
					tempbond,
					rsd1.atom_type_index(jj_term_atomno),
					rsd1.atom_type_index(resconn_atomno1),
					rsd2.atom_type_index(resconn_atomno2),
					rsd2.atom_type_index(kk_term_atomno)
				);
				if ( t_param_ij.is_null() ) continue;

				Vector const &xyz1( rsd1.xyz( jj_term_atomno ) );
				Vector const &xyz2( rsd1.xyz( resconn_atomno1 ) );
				Vector const &xyz3( rsd2.xyz( resconn_atomno2 ) );
				Vector const &xyz4( rsd2.xyz( kk_term_atomno ) );
				Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
				Real score = t_param_ij.energy( angle );
				totalscore += score;

				if ( TR.Debug.visible() ) {
					std::string aname1, aname2, aname3, aname4;
					aname1 = rsd1.atom_name( jj_term_atomno );
					aname2 = rsd1.atom_name( resconn_atomno1 );
					aname3 = rsd2.atom_name( resconn_atomno2 );
					aname4 = rsd2.atom_name( kk_term_atomno );

					TR.Debug << "TOR2 seqpos1,2/atmnames/k123/angle/score: "
						<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
						<< " " << A(19, aname1+"-"+aname2+"-"+aname3+"-"+aname4)
						<< " " << F(6,3,t_param_ij.k(1)) << "," << F(6,3,t_param_ij.k(2))
						<< "," << F(6,3,t_param_ij.k(3))<< "," << F(6,3,t_param_ij.k(4))
						<< " " << F(6,1,angle*180.0/3.14159216)
						<< " " << F(6,2,score)
						<< std::endl;
				}

				if ( calc_deriv ) {
					Real dE_dtors = weight*t_param_ij.deriv( angle );

					Vector f1( 0.0 ), f2( 0.0 );
					numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
					atom_derivs_r1[ jj_term_atomno ].f1() += dE_dtors * f1;
					atom_derivs_r1[ jj_term_atomno ].f2() += dE_dtors * f2;
					f1 = f2 = Vector(0.0);
					numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
					atom_derivs_r1[ resconn_atomno1 ].f1() += dE_dtors * f1;
					atom_derivs_r1[ resconn_atomno1 ].f2() += dE_dtors * f2;
					f1 = f2 = Vector(0.0);
					numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
					atom_derivs_r2[ resconn_atomno2 ].f1() += dE_dtors * f1;
					atom_derivs_r2[ resconn_atomno2 ].f2() += dE_dtors * f2;
					f1 = f2 = Vector(0.0);
					numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
					atom_derivs_r2[ kk_term_atomno ].f1() += dE_dtors * f1;
					atom_derivs_r2[ kk_term_atomno ].f2() += dE_dtors * f2;
				}
			}
		}
	}

	return totalscore;
}


// inproper-torsion (1b)
Real
GenericBondedPotential::eval_residue_energy_and_derivative_improper(
	conformation::Residue const & rsd,
	utility::vector1< DerivVectorPair > & atom_derivs,
	Real const weight,
	bool const calc_deriv
) const {
	Real totalscore = 0.0;

	// for each atom
	for ( Size atm_i=1; atm_i<=rsd.type().natoms(); ++atm_i ) {
		chemical::AtomIndices atm_nbrs = rsd.type().nbrs( atm_i );

		if ( atm_nbrs.size() != 3 ) continue;

		int atm_j = (int)atm_nbrs[1];
		int atm_k = (int)atm_nbrs[2];
		int atm_l = (int)atm_nbrs[3];

		// note: may shuffle atom j/k/l
		SpringParams const &param_ij = lookup_improper_params(
			rsd.atom_type_index(atm_i) , rsd.atom_type_index(atm_j) , rsd.atom_type_index(atm_k), rsd.atom_type_index(atm_l),
			atm_j, atm_k, atm_l
		);
		if ( param_ij.is_null() ) continue;

		Vector const &xyz1( rsd.xyz( atm_i ) );
		Vector const &xyz2( rsd.xyz( atm_j ) );
		Vector const &xyz3( rsd.xyz( atm_k ) );
		Vector const &xyz4( rsd.xyz( atm_l ) );

		Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
		Real score = param_ij.energy( angle );
		totalscore += score;

		if ( TR.Debug.visible() ) {
			Real k = param_ij.k();
			std::string aname1, aname2, aname3, aname4;
			aname1 = rsd.atom_name( atm_i );
			aname2 = rsd.atom_name( atm_j );
			aname3 = rsd.atom_name( atm_k );
			aname4 = rsd.atom_name( atm_l );

			TR.Debug << "IMP1 seqpos/atmnames/k123/angle/score: "
				<< " " << I(4,int(rsd.seqpos()))
				<< " " << A(19, aname1+"-"+aname2+"-"+aname3+"-"+aname4)
				<< " " << F(6,3,k)
				<< " " << F(6,1,angle*180.0/3.14159216)
				<< " " << F(6,2,score)
				<< std::endl;
		}

		if ( calc_deriv ) {
			Real dE_dtors = weight*param_ij.deriv( angle );

			Vector f1( 0.0 ), f2( 0.0 );
			numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
			atom_derivs[ atm_i ].f1() += dE_dtors * f1;
			atom_derivs[ atm_i ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
			atom_derivs[ atm_j ].f1() += dE_dtors * f1;
			atom_derivs[ atm_j ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
			atom_derivs[ atm_k ].f1() += dE_dtors * f1;
			atom_derivs[ atm_k ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
			atom_derivs[ atm_l ].f1() += dE_dtors * f1;
			atom_derivs[ atm_l ].f2() += dE_dtors * f2;
		}
	}

	return totalscore;
}

// inproper-torsion (2b)
Real
GenericBondedPotential::eval_residue_pair_energy_and_derivative_improper(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< DerivVectorPair > & atom_derivs_r1,
	utility::vector1< DerivVectorPair > & atom_derivs_r2,
	Real const weight,
	bool const calc_deriv
) const {
	Real totalscore = 0.0;

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );
		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		// (a) compute impropers centered on 'resconn_atomno1'
		utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
			rsd1.type().atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
		if ( rsd1_atoms_wi1_bond_of_ii.size() == 2 ) {
			Size const res1_lower_atomno1 = rsd1_atoms_wi1_bond_of_ii[ 1 ].key2();
			Size const res1_lower_atomno2 = rsd1_atoms_wi1_bond_of_ii[ 2 ].key2();

			// improper is r1.resconn_atomno1 / r1.res1_lower_atomno1 / r1.res1_lower_atomno2 / r2.resconn_atomno2
			int atm_j=res1_lower_atomno1, atm_k=res1_lower_atomno2, atm_l=-resconn_atomno2; // negative indicates "alternate" residue

			// note: may shuffle atom j/k/l
			SpringParams const &param_ij = lookup_improper_params(
				rsd1.atom_type_index(resconn_atomno1),
				rsd1.atom_type_index(res1_lower_atomno1),
				rsd1.atom_type_index(res1_lower_atomno2),
				rsd2.atom_type_index(resconn_atomno2),
				atm_j, atm_k, atm_l
			);
			if ( param_ij.is_null() ) continue;

			conformation::Residue const &rsd_j = (atm_j<0) ? rsd2 : rsd1;
			conformation::Residue const &rsd_k = (atm_k<0) ? rsd2 : rsd1;
			conformation::Residue const &rsd_l = (atm_l<0) ? rsd2 : rsd1;

			Vector const &xyz1( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz2( rsd_j.xyz( atm_j<0 ? -atm_j : atm_j ) );
			Vector const &xyz3( rsd_k.xyz( atm_k<0 ? -atm_k : atm_k ) );
			Vector const &xyz4( rsd_l.xyz( atm_l<0 ? -atm_l : atm_l ) );

			Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
			Real score = param_ij.energy( angle );
			totalscore += score;

			if ( TR.Debug.visible() ) {
				Real k = param_ij.k();
				std::string aname1, aname2, aname3, aname4;
				aname1 = rsd1.atom_name( resconn_atomno1 );
				aname2 = rsd_j.atom_name( atm_j<0 ? -atm_j : atm_j );
				aname3 = rsd_k.atom_name( atm_k<0 ? -atm_k : atm_k );
				aname4 = rsd_l.atom_name( atm_l<0 ? -atm_l : atm_l );

				TR.Debug << "IMP2 seqpos1,2/atmnames/k123/angle/score: "
					<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
					<< " " << A(19, aname1+"-"+aname2+"-"+aname3+"-"+aname4)
					<< " " << F(6,3,k)
					<< " " << F(6,1,angle*180.0/3.14159216)
					<< " " << F(6,2,score)
					<< std::endl;
			}

			if ( calc_deriv ) {
				Real dE_dtors = weight*param_ij.deriv( angle );
				utility::vector1< DerivVectorPair > &deriv_j = (atm_j<0) ? atom_derivs_r2 : atom_derivs_r1;
				utility::vector1< DerivVectorPair > &deriv_k = (atm_k<0) ? atom_derivs_r2 : atom_derivs_r1;
				utility::vector1< DerivVectorPair > &deriv_l = (atm_l<0) ? atom_derivs_r2 : atom_derivs_r1;

				Vector f1( 0.0 ), f2( 0.0 );
				numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				atom_derivs_r1[ resconn_atomno1 ].f1() += dE_dtors * f1;
				atom_derivs_r1[ resconn_atomno1 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				deriv_j[ atm_j<0 ? -atm_j : atm_j ].f1() += dE_dtors * f1;
				deriv_j[ atm_j<0 ? -atm_j : atm_j ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				deriv_k[ atm_k<0 ? -atm_k : atm_k ].f1() += dE_dtors * f1;
				deriv_k[ atm_k<0 ? -atm_k : atm_k ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				deriv_l[ atm_l<0 ? -atm_l : atm_l ].f1() += dE_dtors * f1;
				deriv_l[ atm_l<0 ? -atm_l : atm_l ].f2() += dE_dtors * f2;
			}
		}

		// (b) compute impropers centered on 'resconn_atomno2'
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2.type().atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
		if ( rsd2_atoms_wi1_bond_of_ii.size() == 2 ) {
			Size const res2_lower_atomno1 = rsd2_atoms_wi1_bond_of_ii[ 1 ].key2();
			Size const res2_lower_atomno2 = rsd2_atoms_wi1_bond_of_ii[ 2 ].key2();

			// improper is r2.resconn_atomno2 / r1.resconn_atomno1 / r2.res2_lower_atomno1 / r2.res2_lower_atomno2
			int atm_j=-resconn_atomno1, atm_k=res2_lower_atomno1, atm_l=res2_lower_atomno2; // negative indicates "alternate" residue

			// note: may shuffle atom j/k/l
			SpringParams const &param_ij = lookup_improper_params(
				rsd2.atom_type_index(resconn_atomno2),
				rsd1.atom_type_index(resconn_atomno1),
				rsd2.atom_type_index(res2_lower_atomno1),
				rsd2.atom_type_index(res2_lower_atomno2),
				atm_j, atm_k, atm_l
			);
			if ( param_ij.is_null() ) continue;

			conformation::Residue const &rsd_j = (atm_j<0) ? rsd1 : rsd2;
			conformation::Residue const &rsd_k = (atm_k<0) ? rsd1 : rsd2;
			conformation::Residue const &rsd_l = (atm_l<0) ? rsd1 : rsd2;

			Vector const &xyz1( rsd2.xyz( resconn_atomno2 ) );
			Vector const &xyz2( rsd_j.xyz( atm_j<0 ? -atm_j : atm_j ) );
			Vector const &xyz3( rsd_k.xyz( atm_k<0 ? -atm_k : atm_k ) );
			Vector const &xyz4( rsd_l.xyz( atm_l<0 ? -atm_l : atm_l ) );

			Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
			Real score = param_ij.energy( angle );
			totalscore += score;

			if ( TR.Debug.visible() ) {
				Real k = param_ij.k();
				std::string aname1, aname2, aname3, aname4;
				aname1 = rsd1.atom_name( resconn_atomno2 );
				aname2 = rsd_j.atom_name( atm_j<0 ? -atm_j : atm_j );
				aname3 = rsd_k.atom_name( atm_k<0 ? -atm_k : atm_k );
				aname4 = rsd_l.atom_name( atm_l<0 ? -atm_l : atm_l );

				TR.Debug << "IMP2 seqpos1,2/atmnames/k123/angle/score: "
					<< " " << I(4,int(rsd1.seqpos())) << "-" << I(4,int(rsd2.seqpos()))
					<< " " << A(19, aname1+"-"+aname2+"-"+aname3+"-"+aname4)
					<< " " << F(6,3,k)
					<< " " << F(6,1,angle*180.0/3.14159216)
					<< " " << F(6,2,score)
					<< std::endl;
			}

			if ( calc_deriv ) {
				Real dE_dtors = weight*param_ij.deriv( angle );
				utility::vector1< DerivVectorPair > &deriv_j = (atm_j<0) ? atom_derivs_r1 : atom_derivs_r2;
				utility::vector1< DerivVectorPair > &deriv_k = (atm_k<0) ? atom_derivs_r1 : atom_derivs_r2;
				utility::vector1< DerivVectorPair > &deriv_l = (atm_l<0) ? atom_derivs_r1 : atom_derivs_r2;

				Vector f1( 0.0 ), f2( 0.0 );
				numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				atom_derivs_r2[ resconn_atomno2 ].f1() += dE_dtors * f1;
				atom_derivs_r2[ resconn_atomno2 ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, angle, f1, f2 );
				deriv_j[ atm_j<0 ? -atm_j : atm_j ].f1() += dE_dtors * f1;
				deriv_j[ atm_j<0 ? -atm_j : atm_j ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				deriv_k[ atm_k<0 ? -atm_k : atm_k ].f1() += dE_dtors * f1;
				deriv_k[ atm_k<0 ? -atm_k : atm_k ].f2() += dE_dtors * f2;
				f1 = f2 = Vector(0.0);
				numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, angle, f1, f2 );
				deriv_l[ atm_l<0 ? -atm_l : atm_l ].f1() += dE_dtors * f1;
				deriv_l[ atm_l<0 ? -atm_l : atm_l ].f2() += dE_dtors * f2;
			}

		}
	}


	return totalscore;
}

void
GenericBondedPotential::modify_tors_params(
	std::string atm1, std::string atm2,
	std::string bondtype,
	std::string atm3, std::string atm4,
	core::Real k1, core::Real k2, core::Real k3, core::Real k4,
	core::Real f1, core::Real f2, core::Real f3, core::Real f4
) {

	int n_tor_changed(0);
	utility::vector1< core::Size > wildcard_vec(1, 0);

	utility::vector1< core::Size > indicesBT = bondorders_map(bondtype);
	utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
	utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
	utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
	utility::vector1< core::Size > const & indices4 = name_index_map[atm4];

	utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
	utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
	utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;
	utility::vector1< core::Size > const &indices4_loop = (indices4.size() == natomtypes) ? wildcard_vec : indices4;

	utility::vector1< core::Size > const & indicesX = name_index_map["X"];
	core::Size const multi_max( indicesX.size()/2 );
	core::Size multBT = indicesBT.size()*multi_max*multi_max*multi_max*multi_max;

	core::Size multiplicity = multBT + indices1.size() * indices2.size() * indices3.size() * indices4.size();

	std::string torsion_type_in = atm1+","+atm2+","+bondtype+","+atm3+","+atm4;
	// this logic is tricky
	//   if the torsion specified is _more specific_ than any in the input, it still should apply
	// we will simply add the new torsion to the database
	//   the overwritten torsion still exists in tors_pot_ but might not be pointed at by anything
	tors_pot_.push_back( GenTorsionParams( k1, k2, k3, k4, f1, f2, f3, f4, multiplicity, torsion_type_in ) );
	core::Size tgt_pot_idx1 = tors_pot_.size();
	core::Size tgt_pot_idx2 = tgt_pot_idx1;
	if ( f1 != 0 || f2 != 0 || f3 != 0 || f4 != 0 ) {
		tors_pot_.push_back( GenTorsionParams( k1, k2, k3, k4, -f1, -f2, -f3, -f4, multiplicity, torsion_type_in ) );
		tgt_pot_idx2 = tors_pot_.size();
	}

	for ( auto i0 : indicesBT ) {
		for ( auto i1 : indices1_loop ) {
			for ( auto i2 : indices2_loop ) {
				for ( auto i3 : indices3_loop ) {
					for ( auto i4 : indices4_loop ) {
						int64_t hashval = get_parameter_hash( 0, i1, i2, i3, i4 );
						auto it = tors_lookup_.find( hashval );

						if ( it == tors_lookup_.end() ) {
							// in theory everything should be covered in the input DB
							// still we can allow it...
							TR << "WARNING: inserting torsion type that doesn't exist in database!" << std::endl;
							tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx1) );
							hashval = get_parameter_hash( i0, i4, i3, i2, i1 );
							tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx2) );
						} else if ( multiplicity <= tors_pot_[ it->second ].multiplicity() ) {
							TR<< "Overwrite params of " << tors_pot_[ it->second ].torsion_type()
								<< " with params of " << torsion_type_in << std::endl;
							it->second = tgt_pot_idx1;
							hashval = get_parameter_hash( i0, i4, i3, i2, i1 );
							it = tors_lookup_.find( hashval );
							it->second = tgt_pot_idx2;
							n_tor_changed++;
						}
					}
				}
			}
		}
	}

	if ( n_tor_changed==0 ) {
		TR << "Warning: No params changed for torsion: " << torsion_type_in << std::endl;
	} else {
		TR << n_tor_changed << " torsion types changed for " << torsion_type_in << std::endl;
	}
}


} // scoring
} // core

