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
	rmIMPROPER,
	rmSPECIAL_TORSION
};


// pre allocated "null" params
const SpringParams GenericBondedPotential::null_sp_param;
const GenTorsionParams GenericBondedPotential::null_tors_param;
const SpecialGenTorsionParams GenericBondedPotential::null_sp_tors_param;


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
	if ( k1_ < 0 ) arg += -2.0*k1_;
	if ( k2_ < 0 ) arg += -2.0*k2_;
	if ( k3_ < 0 ) arg += -2.0*k3_;

	//TR << "E("<<value<<") = " << arg << " ("<<k1_<<","<<k2_<<","<<k3_<<","<<f1_<<","<<f2_<<","<<f3_<<")"<<std::endl;

	return arg;
}

Real
GenTorsionParams::deriv ( core::Real value ) const {
	Real arg = -k1_ * (sin( 1*value - f1_ ));
	arg +=  -2*k2_ * (sin( 2*value - f2_ ));
	arg +=  -3*k3_ * (sin( 3*value - f3_ ));

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
	} else if ( keyword == "f1" ) {
		return f1_;
	} else if ( keyword == "f2" ) {
		return f2_;
	} else if ( keyword == "f3" ) {
		return f3_;
	} else {
		return 0;//TR << keyword << " is invalid!" << std::endl;
	}
}

Real
SpecialGenTorsionParams::energy ( core::Real value ) const {
	Real arg = k1_ * (1+cos( 1*value - f1_ ));
	arg +=  k2_ * (1+cos( 2*value - f2_ ));
	arg +=  k3_ * (1+cos( 3*value - f3_ ));
	arg +=  k4_ * (1+cos( 4*value - f4_ ));
	arg +=  k8_ * (1+cos( 8*value - f8_ ));
	if ( k1_ < 0 ) arg += -2.0*k1_;
	if ( k2_ < 0 ) arg += -2.0*k2_;
	if ( k3_ < 0 ) arg += -2.0*k3_;
	if ( k4_ < 0 ) arg += -2.0*k4_;
	if ( k8_ < 0 ) arg += -2.0*k8_;

	return arg;
}

Real
SpecialGenTorsionParams::deriv ( core::Real value ) const {
	Real arg = -k1_ * (sin( 1*value - f1_ ));
	arg +=  -2*k2_ * (sin( 2*value - f2_ ));
	arg +=  -3*k3_ * (sin( 3*value - f3_ ));
	arg +=  -4*k4_ * (sin( 4*value - f4_ ));
	arg +=  -8*k8_ * (sin( 8*value - f8_ ));

	return arg;
}

Real
SpecialGenTorsionParams::get_params(std::string keyword) const {
	if ( keyword == "k1" ) {
		return k1_;
	} else if ( keyword == "k2" ) {
		return k2_;
	} else if ( keyword == "k3" ) {
		return k3_;
	} else if ( keyword == "k4" ) {
		return k4_;
	} else if ( keyword == "k8" ) {
		return k8_;
	} else if ( keyword == "f1" ) {
		return f1_;
	} else if ( keyword == "f2" ) {
		return f2_;
	} else if ( keyword == "f3" ) {
		return f3_;
	} else if ( keyword == "f4" ) {
		return f4_;
	} else if ( keyword == "f8" ) {
		return f8_;
	} else {
		return 0; //TR << keyword << " is invalid!" << std::endl;
	}
}


SpringParams const &
GenericBondedPotential::lookup_bond_params( Size type1, Size type2 ) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] ) return null_sp_param;

	auto it = bond_lookup_.find( get_parameter_hash(type1, type2) );
	if ( it == bond_lookup_.end() ) it = bond_lookup_.find( get_parameter_hash(type1, 0) );
	if ( it == bond_lookup_.end() ) it = bond_lookup_.find( get_parameter_hash(0, type2) );
	if ( it == bond_lookup_.end() ) it = bond_lookup_.find( get_parameter_hash(0, 0) );

	return (it == bond_lookup_.end()) ? null_sp_param : bond_pot_[it->second];
}

SpringParams const &
GenericBondedPotential::lookup_angle_params( Size type1, Size type2, Size type3 ) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3] ) return null_sp_param;

	auto it = angle_lookup_.find( get_parameter_hash(type1, type2, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, type2, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(type1, type2, 0) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(type1, 0, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, type2, 0) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(type1, 0, 0) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, 0, type3) );
	if ( it == angle_lookup_.end() ) it = angle_lookup_.find( get_parameter_hash(0, 0, 0) );

	return (it == angle_lookup_.end()) ? null_sp_param : angle_pot_[it->second];
}

GenTorsionParams const &
GenericBondedPotential::lookup_tors_params( Size type1, Size type2, Size type3, Size type4 ) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3]|| !defined_atom_types_[type4] ) return null_tors_param;

	auto it = tors_lookup_.find( get_parameter_hash(type1, type2, type3, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, type2, type3, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(type1, 0, type3, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(type1, type2, 0, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(type1, type2, type3, 0) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, 0, type3, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, type2, 0, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, type2, type3, 0) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(type1, 0, 0, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(type1, 0, type3, 0) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(type1, type2, 0, 0) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(type1, 0, 0, 0) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, type2, 0, 0) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, 0, type3, 0) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, 0, 0, type4) );
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, 0, 0, 0) );

	return tors_pot_[it->second];
}

SpecialGenTorsionParams const &
GenericBondedPotential::lookup_special_tors_params( Size type1, Size type2, Size type3, Size type4 ) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3]|| !defined_atom_types_[type4] ) return null_sp_tors_param;

	auto it = special_tors_lookup_.find( get_parameter_hash(type1, type2, type3, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, type2, type3, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(type1, 0, type3, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(type1, type2, 0, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(type1, type2, type3, 0) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, 0, type3, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, type2, 0, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, type2, type3, 0) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(type1, 0, 0, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(type1, 0, type3, 0) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(type1, type2, 0, 0) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(type1, 0, 0, 0) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, type2, 0, 0) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, 0, type3, 0) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, 0, 0, type4) );
	if ( it == special_tors_lookup_.end() ) it = special_tors_lookup_.find( get_parameter_hash(0, 0, 0, 0) );

	return (it == special_tors_lookup_.end()) ? null_sp_tors_param : special_tors_pot_[it->second];
}

// note!  improper torsions are order-specific, me might need to shuffle
SpringParams const &
GenericBondedPotential::lookup_improper_params(
	Size type1, Size type2, Size type3, Size type4, int &idx2, int &idx3, int &idx4
) const {
	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3]|| !defined_atom_types_[type4] ) return null_sp_param;

	Size i2in=idx2, i3in=idx3, i4in=idx4;

	auto it = improper_lookup_.find( get_parameter_hash(type1, type2, type3, type4) );
	if ( it != improper_lookup_.end() ) {
		return improper_pot_[it->second];
	}
	it = improper_lookup_.find( get_parameter_hash(type1, type2, type4, type3) );
	if ( it != improper_lookup_.end() ) {
		idx3 = i4in; idx4 = i3in;
		return improper_pot_[it->second];
	}


	it = improper_lookup_.find( get_parameter_hash(type1, type3, type2, type4) );
	if ( it != improper_lookup_.end() ) {
		idx2 = i3in; idx3 = i2in;
		return improper_pot_[it->second];
	}
	it = improper_lookup_.find( get_parameter_hash(type1, type3, type4, type2) );
	if ( it != improper_lookup_.end() ) {
		idx2 = i3in; idx3 = i4in; idx4 = i2in;
		return improper_pot_[it->second];
	}

	it = improper_lookup_.find( get_parameter_hash(type1, type4, type2, type3) );
	if ( it != improper_lookup_.end() ) {
		idx2 = i4in; idx3 = i2in; idx4 = i3in;
		return improper_pot_[it->second];
	}
	it = improper_lookup_.find( get_parameter_hash(type1, type4, type3, type2) );
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
	if ( option[ corrections::genpotential::set_special_torsion_params ].user() ) {
		modify_special_torsion_params_from_cmd_line();
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
	runtime_assert( ntypes < 65536 ); // we pack 4 atom indices into 64 bits

	// map atom names to ATS indices - handle wildcards
	ReadMode read_mode = rmNONE;
	Size linenum( 0 );
	natomtypes = 0;
	std::string fileline, tag;

	// special treatment of "fallbacks", e.g. 'X' in the params file
	//
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
		} else if ( tag.compare("SPECIAL_TORSION") == 0 ) {
			read_mode = rmSPECIAL_TORSION;
			continue;
		} else if ( tag.substr(0,1) == "#" ) {
			continue; // comment or empty line
		}

		if ( read_mode == rmATOM ) {       // in "ATOM" block
			core::Size atm_idx = ats->atom_type_index(tag);
			defined_atom_types_[ atm_idx ] = true;
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
					int64_t hashval = get_parameter_hash( i1, i2 );
					auto it = bond_lookup_.find( hashval );
					// use the MOST SPECIFIC potential
					if ( it == bond_lookup_.end() ) {
						bond_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
						hashval = get_parameter_hash( i2, i1 );
						bond_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
					} else if ( multiplicity < bond_pot_[ it->second ].multiplicity() ) {
						it->second = tgt_pot_idx;
						hashval = get_parameter_hash( i2, i1 );
						it = bond_lookup_.find( hashval );
						it->second = tgt_pot_idx;
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
						int64_t hashval = get_parameter_hash( i1, i2, i3 );
						auto it = angle_lookup_.find( hashval );
						// use the MOST SPECIFIC potential
						if ( it == angle_lookup_.end() ) {
							angle_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
							hashval = get_parameter_hash( i3, i2, i1 );
							angle_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
						} else if ( multiplicity < angle_pot_[ it->second ].multiplicity() ) {
							it->second = tgt_pot_idx;
							hashval = get_parameter_hash( i3, i2, i1 );
							it = angle_lookup_.find( hashval );
							it->second = tgt_pot_idx;
						}
					}
				}
			}

		} else if ( read_mode == rmTORSION ) { // in "TORSION" block
			std::string atm1, atm2, atm3, atm4;
			core::Real k1, k2, k3, f1=0.0, f2=0.0, f3=0.0;
			std::string dummy, extra;

			linestream >> atm1 >> atm2 >> atm3 >> atm4 >> dummy >> k1 >> k2 >> k3 >> extra;
			if ( extra == "phase" ) {
				linestream >> f1 >> f2 >> f3;
			} else {
				f1 = 0.0; f2 = 0.0; f3 = 0.0;
			}

			utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
			utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
			utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
			utility::vector1< core::Size > const & indices4 = name_index_map[atm4];
			core::Size multiplicity = indices1.size() * indices2.size() * indices3.size() * indices4.size();

			utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
			utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
			utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;
			utility::vector1< core::Size > const &indices4_loop = (indices4.size() == natomtypes) ? wildcard_vec : indices4;

			tors_pot_.push_back( GenTorsionParams( k1, k2, k3, f1, f2, f3, multiplicity ) );
			core::Size tgt_pot_idx1 = tors_pot_.size();
			core::Size tgt_pot_idx2 = tgt_pot_idx1;

			// flipping atom order reverses phase -- generate a new constraint if necessary
			if ( f1 != 0 || f2 != 0 || f3 != 0 ) {
				tors_pot_.push_back( GenTorsionParams( k1, k2, k3, -f1, -f2, -f3, multiplicity ) );
				tgt_pot_idx2 = tors_pot_.size();
			}

			for ( auto i1 : indices1_loop ) {
				for ( auto i2 : indices2_loop ) {
					for ( auto i3 : indices3_loop ) {
						for ( auto i4 : indices4_loop ) {
							int64_t hashval = get_parameter_hash( i1, i2, i3, i4 );
							auto it = tors_lookup_.find( hashval );

							// use the MOST SPECIFIC potential
							if ( it == tors_lookup_.end()  ) {
								tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx1) );
								hashval = get_parameter_hash( i4, i3, i2, i1 );
								tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx2) );
							} else if ( multiplicity < tors_pot_[ it->second ].multiplicity() ) {
								it->second = tgt_pot_idx1;
								hashval = get_parameter_hash( i4, i3, i2, i1 );
								it = tors_lookup_.find( hashval );
								it->second = tgt_pot_idx2;
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
							int64_t hashval = get_parameter_hash( i1, i2, i3, i4 );
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

		} else if ( read_mode == rmSPECIAL_TORSION ) { // in "SPECIAL_TORSION " block
			std::string atm1=tag, atm2, atm3, atm4;
			core::Real k1, k2, k3, k4, k8;
			std::string dummy;

			linestream >> atm2 >> atm3 >> atm4 >> dummy >> k1 >> k2 >> k3 >> k4 >> k8;

			utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
			utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
			utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
			utility::vector1< core::Size > const & indices4 = name_index_map[atm4];
			core::Size multiplicity = indices1.size() * indices2.size() * indices3.size() * indices4.size();

			utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
			utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
			utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;
			utility::vector1< core::Size > const &indices4_loop = (indices4.size() == natomtypes) ? wildcard_vec : indices4;

			// since phase = 0 no need to generate double params
			special_tors_pot_.push_back( SpecialGenTorsionParams( k1, k2, k3, k4, k8, 0.0, 0.0, 0.0, 0.0, 0.0, multiplicity ) );
			core::Size tgt_pot_idx = special_tors_pot_.size();

			for ( auto i1 : indices1_loop ) {
				for ( auto i2 : indices2_loop ) {
					for ( auto i3 : indices3_loop ) {
						for ( auto i4 : indices4_loop ) {
							int64_t hashval = get_parameter_hash( i1, i2, i3, i4 );
							auto it = special_tors_lookup_.find( hashval );
							// use the MOST SPECIFIC potential
							if ( it == special_tors_lookup_.end() ) {
								special_tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
								hashval = get_parameter_hash( i4, i3, i2, i1 );
								special_tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx) );
							} else if ( multiplicity < special_tors_pot_[ it->second ].multiplicity() ) {
								it->second = tgt_pot_idx;
								hashval = get_parameter_hash( i4, i3, i2, i1 );
								it = tors_lookup_.find( hashval );
								it->second = tgt_pot_idx;
							}
						}
					}
				}
			}

		} //end if read_mode
	} //end while

	TR << "Added total " << bond_pot_.size() << " bond parameters corresponding to "
		<< bond_lookup_.size() << " unique bonds." << std::endl;
	TR << "Added total " << angle_pot_.size() << " angle parameters corresponding to "
		<< angle_lookup_.size() << " unique angles." << std::endl;
	TR << "Added total " << tors_pot_.size() << " torsion parameters corresponding to "
		<< tors_lookup_.size() << " unique torsions." << std::endl;
	TR << "Added total " << improper_pot_.size() << " improper parameters corresponding to "
		<< improper_lookup_.size() << " unique impropers." << std::endl;
	TR << "Added total " << special_tors_pot_.size() << " special-torsion parameters corresponding to "
		<< special_tors_lookup_.size() << " unique torsions." << std::endl;

}

void
GenericBondedPotential::modify_torsion_params_from_cmd_line() {
	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;

	utility::vector1< std::string > const & mods( option[ corrections::genpotential::set_torsion_params ] );

	string const errmsg( "-corrections:genpotential:set_torsion_params format should be:: -corrections:genpotential:set_torsion_params <set1>:<atom1>:<param1>:"
		"<setting1> <set2>:<atom2>:<param2>:<setting2> ...; for example: '-corrections:genpotential:set_torsion_params fa_standard:C*:CS:CS:C*:k1:0.0:k2:0.0:k3:0.077 fa_standard:CD:CS:CS:CD:k1:0.435:k2:0.039:k3:0.070' ");

	for ( string const & mod : mods ) {
		// mod should look like (for example):  "fa_standard:OOC:LK_RADIUS:4.5"

		core::uint const pos1( mod.find( ":" ) );
		if ( pos1 == string::npos ) utility_exit_with_message( errmsg );
		string const atomset_tag( mod.substr( 0, pos1 ) );
		if ( atomset_tag != "fa_standard" ) continue;

		core::uint pos_starting( pos1 + 1 );
		core::uint const pos2( mod.substr( pos_starting ).find( ":" ) );
		if ( pos2 == string::npos ) utility_exit_with_message( errmsg );
		string const atom_name1( mod.substr( pos_starting, pos2 ) );
		if ( name_index_map.find( atom_name1 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name1 );
		}

		pos_starting = pos_starting + pos2 + 1;
		core::uint const pos3( mod.substr( pos_starting ).find( ":" ) );
		if ( pos3 == string::npos ) utility_exit_with_message( errmsg );
		std::string const atom_name2( mod.substr( pos_starting, pos3 ) );
		if ( name_index_map.find( atom_name2 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name2 );
		}

		pos_starting = pos_starting + pos3 +1;
		core::uint const pos4( mod.substr( pos_starting ).find( ":" ) );
		if ( pos4 == string::npos ) utility_exit_with_message( errmsg );
		string const atom_name3( mod.substr( pos_starting, pos4 ) );
		if ( name_index_map.find( atom_name3 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name3 );
		}

		pos_starting = pos_starting + pos4 + 1;
		core::uint const pos5( mod.substr( pos_starting ).find( ":" ) );
		if ( pos5 == string::npos ) utility_exit_with_message( errmsg );
		string const atom_name4( mod.substr( pos_starting, pos5 ) );
		if ( name_index_map.find( atom_name4 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name1 );
		}

		pos_starting = pos_starting + pos5 + 1;
		core::uint const pos6( mod.substr( pos_starting ).find( ":" ) );
		if ( pos6 == string::npos ) utility_exit_with_message( errmsg );
		string const param1( mod.substr( pos_starting, pos6 ) );
		if ( param1 != "k1" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k1, got " + param1 + "instead." );
		}

		pos_starting = pos_starting + pos6 + 1;
		core::uint const pos7( mod.substr( pos_starting ).find( ":" ) );
		if ( pos7 == string::npos ) utility_exit_with_message( errmsg );
		string stringsetting( mod.substr( pos_starting, pos7 ) );
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting1( double_of( stringsetting ) );

		pos_starting = pos_starting + pos7 + 1;
		core::uint const pos8( mod.substr( pos_starting ).find( ":" ) );
		if ( pos8 == string::npos ) utility_exit_with_message( errmsg );
		string const param2( mod.substr( pos_starting, pos8 ) );
		if ( param2 != "k2" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k2, got " + param2 + "instead." );
		}

		pos_starting = pos_starting + pos8 + 1;
		core::uint const pos9( mod.substr( pos_starting ).find( ":" ) );
		if ( pos9 == string::npos ) utility_exit_with_message( errmsg );
		stringsetting = mod.substr( pos_starting, pos9 );
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting2( double_of( stringsetting ) );

		pos_starting = pos_starting + pos9 + 1;
		core::uint const pos10( mod.substr( pos_starting ).find( ":" ) );
		if ( pos10 == string::npos ) utility_exit_with_message( errmsg );
		string const param3( mod.substr( pos_starting, pos10 ) );
		if ( param3 != "k3" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k3, got " + param3 + "instead." );
		}

		pos_starting = pos_starting + pos10 + 1;
		stringsetting = mod.substr( pos_starting);
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting3( double_of( stringsetting ) );

		modify_tors_params( atom_name1, atom_name2,
			atom_name3, atom_name4, setting1, setting2, setting3
		);
		TR << "modify_torsion_params_from_cmd_line: setting " << "fa_standard" << ' ' << atom_name1 << "," << atom_name2 << "," << atom_name3 << "," << atom_name4 << ':' << ' ' << param1 << ':' << setting1 << ' ' << param2 << ':' << setting2 << ' ' << param3 << ':' << setting3 << endl;
	}// end for
}

void
GenericBondedPotential::modify_special_torsion_params_from_cmd_line() {
	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;

	utility::vector1< std::string > const & mods( option[ corrections::genpotential::set_special_torsion_params ] );

	string const errmsg( "-corrections:genpotential:set_special torsion_params format should be:: -corrections:genpotential:set_torsion_params <set1>:<atom1>:<param1>:"
		"<setting1> <set2>:<atom2>:<param2>:<setting2> ...; for example: '-corrections:genpotential:set_torsion_params fa_standard:X:CRb:CRb:X:k1:0.000:k2:-0.226:k3:0.000:k4:0.093:k8:0.000 ");

	for ( string const & mod : mods ) {
		// mod should look like (for example):  "fa_standard:OOC:LK_RADIUS:4.5"

		core::uint const pos1( mod.find( ":" ) );
		if ( pos1 == string::npos ) utility_exit_with_message( errmsg );
		string const atomset_tag( mod.substr( 0, pos1 ) );
		if ( atomset_tag != "fa_standard" ) continue;

		core::uint pos_starting( pos1 + 1 );
		core::uint const pos2( mod.substr( pos_starting ).find( ":" ) );
		if ( pos2 == string::npos ) utility_exit_with_message( errmsg );
		string const atom_name1( mod.substr( pos_starting, pos2 ) );
		if ( name_index_map.find( atom_name1 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name1 );
		}

		pos_starting = pos_starting + pos2 + 1;
		core::uint const pos3( mod.substr( pos_starting ).find( ":" ) );
		if ( pos3 == string::npos ) utility_exit_with_message( errmsg );
		string const atom_name2( mod.substr( pos_starting, pos3 ) );
		if ( name_index_map.find( atom_name2 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name2 );
		}

		pos_starting = pos_starting + pos3 +1;
		core::uint const pos4( mod.substr( pos_starting ).find( ":" ) );
		if ( pos4 == string::npos ) utility_exit_with_message( errmsg );
		string const atom_name3( mod.substr( pos_starting, pos4 ) );
		if ( name_index_map.find( atom_name3 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name3 );
		}

		pos_starting = pos_starting + pos4 + 1;
		core::uint const pos5( mod.substr( pos_starting ).find( ":" ) );
		if ( pos5 == string::npos ) utility_exit_with_message( errmsg );
		string const atom_name4( mod.substr( pos_starting, pos5 ) );
		if ( name_index_map.find( atom_name1 ) == name_index_map.end() ) {
			utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name1 );
		}

		pos_starting = pos_starting + pos5 + 1;
		core::uint const pos6( mod.substr( pos_starting ).find( ":" ) );
		if ( pos6 == string::npos ) utility_exit_with_message( errmsg );
		string const param1( mod.substr( pos_starting, pos6 ) );
		if ( param1 != "k1" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k1, got " + param1 + "instead." );
		}
		pos_starting = pos_starting + pos6 + 1;
		core::uint const pos7( mod.substr( pos_starting ).find( ":" ) );
		if ( pos7 == string::npos ) utility_exit_with_message( errmsg );
		string stringsetting( mod.substr( pos_starting, pos7 ) );
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting1( double_of( stringsetting ) );

		pos_starting = pos_starting + pos7 + 1;
		core::uint const pos8( mod.substr( pos_starting ).find( ":" ) );
		if ( pos8 == string::npos ) utility_exit_with_message( errmsg );
		string const param2( mod.substr( pos_starting, pos8 ) );
		if ( param2 != "k2" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k2, got " + param2 + "instead." );
		}
		pos_starting = pos_starting + pos8 + 1;
		core::uint const pos9( mod.substr( pos_starting ).find( ":" ) );
		if ( pos9 == string::npos ) utility_exit_with_message( errmsg );
		stringsetting = mod.substr( pos_starting, pos9 );
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting2( double_of( stringsetting ) );

		pos_starting = pos_starting + pos9 + 1;
		core::uint const pos10( mod.substr( pos_starting ).find( ":" ) );
		if ( pos10 == string::npos ) utility_exit_with_message( errmsg );
		string const param3( mod.substr( pos_starting, pos10 ) );
		if ( param2 != "k3" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k3, got " + param3 + "instead." );
		}
		pos_starting = pos_starting + pos10 + 1;
		core::uint const pos11( mod.substr( pos_starting ).find( ":" ) );
		if ( pos11 == string::npos ) utility_exit_with_message( errmsg );
		stringsetting = mod.substr( pos_starting, pos11 );
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting3( double_of( stringsetting ) );

		pos_starting = pos_starting + pos11 + 1;
		core::uint const pos12( mod.substr( pos_starting ).find( ":" ) );
		if ( pos12 == string::npos ) utility_exit_with_message( errmsg );
		string const param4( mod.substr( pos_starting, pos12 ) );
		if ( param4 != "k4" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k4, got " + param4 + "instead." );
		}
		pos_starting = pos_starting + pos12 + 1;
		core::uint const pos13( mod.substr( pos_starting ).find( ":" ) );
		if ( pos13 == string::npos ) utility_exit_with_message( errmsg );
		stringsetting = mod.substr( pos_starting, pos13 );
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting4( double_of( stringsetting ) );


		pos_starting = pos_starting + pos13 + 1;
		core::uint const pos14( mod.substr( pos_starting ).find( ":" ) );
		if ( pos14 == string::npos ) utility_exit_with_message( errmsg );
		string const param5( mod.substr( pos_starting, pos14 ) );
		if ( param5 != "k8" ) {
			utility_exit_with_message( errmsg + ". Expecting parameter name: k8, got " + param5 + "instead." );
		}
		pos_starting = pos_starting + pos14 + 1;
		stringsetting = mod.substr( pos_starting);
		if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
		Real const setting5( double_of( stringsetting ) );

		modify_tors_params( atom_name1, atom_name2,
			atom_name3, atom_name4, setting1, setting2, setting3, setting4, setting5
		);

		TR << "modify_special_torsion_params_from_cmd_line: setting " << "fa_standard" << ' ' << atom_name1 << "," << atom_name2 << "," << atom_name3 << "," << atom_name4 << ':' << ' ' << param1 << ':' << setting1 << ' ' << param2 << ':' << setting2 << ' ' << param3 << ':' << setting3 << ' ' << param4 << ':' << setting4 << ' ' << param5 << ':' << setting5 << endl;

	} // end for
}

void
GenericBondedPotential::setup_for_scoring( pose::Pose & /*pose*/ ) const
{ }

// energies (1b)
void
GenericBondedPotential::residue_energy(
	conformation::Residue const & rsd,
	EnergyMap & emap
) const {
	utility::vector1< DerivVectorPair > dummy; // I wanna get rid of this...

	Real E_bonded  = eval_residue_energy_and_derivative_bond( rsd, dummy );
	Real E_angle   = eval_residue_energy_and_derivative_angle( rsd, dummy );
	Real E_torsion = eval_residue_energy_and_derivative_torsion( rsd, dummy );
	Real E_improper= eval_residue_energy_and_derivative_improper( rsd, dummy );

	emap[ gen_bonded ]  += E_bonded + E_angle + E_torsion;
	emap[ gen_bonded_bond ]  += E_bonded;
	emap[ gen_bonded_angle ]   += E_angle;
	emap[ gen_bonded_torsion ] += E_torsion;
	emap[ gen_bonded_improper ] += E_improper;
}

// energies (2b)
void
GenericBondedPotential::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	EnergyMap & emap
) const {
	// sanity checks
	if ( !rsd1.is_bonded(rsd2) ) return;
	bool res1first = rsd2.seqpos() > rsd1.seqpos();
	if ( rsd1.has_upper_connect() && rsd1.connected_residue_at_resconn( rsd1.type().upper_connect_id() ) == rsd2.seqpos() ) {
		res1first = true;
	} else if ( rsd2.has_upper_connect() && rsd2.connected_residue_at_resconn( rsd2.type().upper_connect_id() ) == rsd1.seqpos() ) {
		res1first = false;
	}

	conformation::Residue const &resLower = res1first?rsd1:rsd2;
	conformation::Residue const &resUpper = res1first?rsd2:rsd1;

	if ( resLower.has_variant_type(core::chemical::CUTPOINT_LOWER) ) return;
	if ( resUpper.has_variant_type(core::chemical::CUTPOINT_UPPER) ) return;

	utility::vector1< DerivVectorPair > dummy1, dummy2;
	Real E_bonded  = eval_residue_pair_energy_and_derivative_bond( resLower, resUpper, dummy1, dummy2 );
	Real E_angle   = eval_residue_pair_energy_and_derivative_angle( resLower, resUpper, dummy1, dummy2 );
	Real E_torsion = eval_residue_pair_energy_and_derivative_torsion( resLower, resUpper, dummy1, dummy2 );
	Real E_improper= eval_residue_pair_energy_and_derivative_improper( resLower, resUpper, dummy1, dummy2 );

	emap[ gen_bonded ]  += E_bonded + E_angle + E_torsion;
	emap[ gen_bonded_bond ]  += E_bonded;
	emap[ gen_bonded_angle ]   += E_angle;
	emap[ gen_bonded_torsion ] += E_torsion;
	emap[ gen_bonded_improper ] += E_improper;
}

// derivatives (1b)
void
GenericBondedPotential::residue_derivatives(
	conformation::Residue const & rsd,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {
	eval_residue_energy_and_derivative_bond( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_bond ], true );
	eval_residue_energy_and_derivative_angle( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_angle ], true );
	eval_residue_energy_and_derivative_torsion( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_torsion ], true );
	eval_residue_energy_and_derivative_improper( rsd, atom_derivs, weights[ gen_bonded ] + weights[ gen_bonded_improper ], true );
}

// derivatives (2b)
void
GenericBondedPotential::residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs_r1,
	utility::vector1< DerivVectorPair > & atom_derivs_r2
) const {
	eval_residue_pair_energy_and_derivative_bond(
		rsd1, rsd2, atom_derivs_r1, atom_derivs_r2, weights[ gen_bonded ]+weights[ gen_bonded_bond ], true );
	eval_residue_pair_energy_and_derivative_angle(
		rsd1, rsd2, atom_derivs_r1, atom_derivs_r2, weights[ gen_bonded ]+weights[ gen_bonded_angle ], true );
	eval_residue_pair_energy_and_derivative_torsion(
		rsd1, rsd2, atom_derivs_r1, atom_derivs_r2, weights[ gen_bonded ]+weights[ gen_bonded_torsion ], true );
	eval_residue_pair_energy_and_derivative_improper(
		rsd1, rsd2, atom_derivs_r1, atom_derivs_r2, weights[ gen_bonded ]+weights[ gen_bonded_improper ], true );
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
	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );
		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		SpringParams const &param_ij = lookup_bond_params(
			rsd1.atom_type_index(resconn_atomno1) , rsd2.atom_type_index(resconn_atomno2) );
		if ( param_ij.is_null() ) continue;

		Vector const &xyz1( rsd1.xyz( resconn_atomno1 ) );
		Vector const &xyz2( rsd2.xyz( resconn_atomno2 ) );
		Real d = xyz1.distance( xyz2 );
		Real score = param_ij.energy( d );
		totalscore += score;

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

			Vector const &xyz1( rsd1.xyz( res1_lower_atomno ) );
			Vector const &xyz2( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz3( rsd2.xyz( resconn_atomno2 ) );
			Real angle = numeric::angle_radians( xyz1, xyz2, xyz3 );
			Real score = param_ij.energy( angle );
			totalscore += score;

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

			Vector const &xyz1( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz2( rsd2.xyz( resconn_atomno2 ) );
			Vector const &xyz3( rsd2.xyz( res2_lower_atomno ) );
			Real angle = numeric::angle_radians( xyz1, xyz2, xyz3 );
			Real score = param_ij.energy( angle );
			totalscore += score;

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

		// lookup in special torsions first
		SpecialGenTorsionParams const &s_param_ij = lookup_special_tors_params(
			rsd.atom_type_index(atm_i) , rsd.atom_type_index(atm_j) , rsd.atom_type_index(atm_k), rsd.atom_type_index(atm_l) );
		// then in std torsions
		GenTorsionParams const &t_param_ij = lookup_tors_params(
			rsd.atom_type_index(atm_i) , rsd.atom_type_index(atm_j) , rsd.atom_type_index(atm_k), rsd.atom_type_index(atm_l) );

		if ( s_param_ij.is_null() && t_param_ij.is_null() ) continue;

		Vector const &xyz1( rsd.xyz( atm_i ) );
		Vector const &xyz2( rsd.xyz( atm_j ) );
		Vector const &xyz3( rsd.xyz( atm_k ) );
		Vector const &xyz4( rsd.xyz( atm_l ) );
		Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
		Real score;
		if ( !s_param_ij.is_null() ) {
			score = s_param_ij.energy( angle );
		} else {
			score = t_param_ij.energy( angle );
		}
		totalscore += score;

		if ( calc_deriv ) {
			Real dE_dtors;
			if ( !s_param_ij.is_null() ) {
				dE_dtors = weight*s_param_ij.deriv( angle );
			} else {
				dE_dtors = weight*t_param_ij.deriv( angle );
			}

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

		/// 1. iterate over all pairs of pairs within 1 bond of either residue connection atom.
		utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
			rsd1.type().atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2.type().atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));

		for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
			Size const jj_term_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();
			for ( Size kk = 1; kk <= rsd2_atoms_wi1_bond_of_ii.size(); ++kk ) {
				Size const kk_term_atomno = rsd2_atoms_wi1_bond_of_ii[ kk ].key2();

				// torsion is r1.jj_term_atomno / r1.resconn_atomno1 / r2.resconn_atomno2 / r2.kk_term_atomno
				SpecialGenTorsionParams const &s_param_ij = lookup_special_tors_params(
					rsd1.atom_type_index(jj_term_atomno) ,
					rsd1.atom_type_index(resconn_atomno1) ,
					rsd2.atom_type_index(resconn_atomno2),
					rsd2.atom_type_index(kk_term_atomno) );
				GenTorsionParams const &t_param_ij = lookup_tors_params(
					rsd1.atom_type_index(jj_term_atomno) ,
					rsd1.atom_type_index(resconn_atomno1) ,
					rsd2.atom_type_index(resconn_atomno2),
					rsd2.atom_type_index(kk_term_atomno) );
				if ( s_param_ij.is_null() && t_param_ij.is_null() ) continue;

				Vector const &xyz1( rsd1.xyz( jj_term_atomno ) );
				Vector const &xyz2( rsd1.xyz( resconn_atomno1 ) );
				Vector const &xyz3( rsd2.xyz( resconn_atomno2 ) );
				Vector const &xyz4( rsd2.xyz( kk_term_atomno ) );
				Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
				Real score;
				if ( !s_param_ij.is_null() ) {
					score = s_param_ij.energy( angle );
				} else {
					score = t_param_ij.energy( angle );
				}
				totalscore += score;

				if ( calc_deriv ) {
					Real dE_dtors;
					if ( !s_param_ij.is_null() ) {
						dE_dtors = weight*s_param_ij.deriv( angle );
					} else {
						dE_dtors = weight*t_param_ij.deriv( angle );
					}

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

		/// 2. iterate over all triples on residue 1 within 2 bonds of resconn_atomno1.
		utility::vector1< chemical::three_atom_set > const & rsd1_atoms_wi2_bonds_of_ii(
			rsd1.type().atoms_within_two_bonds_of_a_residue_connection( resconn_id1 ));

		for ( Size jj = 1; jj <= rsd1_atoms_wi2_bonds_of_ii.size(); ++jj ) {
			debug_assert( rsd1_atoms_wi2_bonds_of_ii[ jj ].key1() == resconn_atomno1 );

			Size const jj_atom2 = rsd1_atoms_wi2_bonds_of_ii[ jj ].key2();
			Size const jj_atom1 = rsd1_atoms_wi2_bonds_of_ii[ jj ].key3();

			// torsion is r1.jj_atom1 / r1.jj_atom2 / r1.resconn_atomno1 / r2.resconn_atomno2
			SpecialGenTorsionParams const &s_param_ij = lookup_special_tors_params(
				rsd1.atom_type_index(jj_atom1) ,
				rsd1.atom_type_index(jj_atom2) ,
				rsd1.atom_type_index(resconn_atomno1),
				rsd2.atom_type_index(resconn_atomno2));
			GenTorsionParams const &t_param_ij = lookup_tors_params(
				rsd1.atom_type_index(jj_atom1) ,
				rsd1.atom_type_index(jj_atom2) ,
				rsd1.atom_type_index(resconn_atomno1),
				rsd2.atom_type_index(resconn_atomno2) );
			if ( s_param_ij.is_null() && t_param_ij.is_null() ) continue;

			Vector const &xyz1( rsd1.xyz( jj_atom1 ) );
			Vector const &xyz2( rsd1.xyz( jj_atom2 ) );
			Vector const &xyz3( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz4( rsd2.xyz( resconn_atomno2 ) );
			Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
			Real score;
			if ( !s_param_ij.is_null() ) {
				score = s_param_ij.energy( angle );
			} else {
				score = t_param_ij.energy( angle );
			}
			totalscore += score;

			if ( calc_deriv ) {
				Real dE_dtors;
				if ( !s_param_ij.is_null() ) {
					dE_dtors = weight*s_param_ij.deriv( angle );
				} else {
					dE_dtors = weight*t_param_ij.deriv( angle );
				}

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

		/// 3. iterate over all triples on residue 2 within 2 bonds of resconn_atomno2.
		utility::vector1< chemical::three_atom_set > const & rsd2_atoms_wi2_bonds_of_ii(
			rsd2.type().atoms_within_two_bonds_of_a_residue_connection( resconn_id2 ));

		for ( Size jj = 1; jj <= rsd2_atoms_wi2_bonds_of_ii.size(); ++jj ) {
			debug_assert( rsd2_atoms_wi2_bonds_of_ii[ jj ].key1() == resconn_atomno2 );
			Size const jj_atom3 = rsd2_atoms_wi2_bonds_of_ii[ jj ].key2();
			Size const jj_atom4 = rsd2_atoms_wi2_bonds_of_ii[ jj ].key3();

			// torsion is r1.resconn_atomno1 / r2.resconn_atomno2 / r1.jj_atom3 / r2.jj_atom4
			SpecialGenTorsionParams const &s_param_ij = lookup_special_tors_params(
				rsd1.atom_type_index(resconn_atomno1) ,
				rsd2.atom_type_index(resconn_atomno2) ,
				rsd2.atom_type_index(jj_atom3),
				rsd2.atom_type_index(jj_atom4) );
			GenTorsionParams const &t_param_ij = lookup_tors_params(
				rsd1.atom_type_index(resconn_atomno1) ,
				rsd2.atom_type_index(resconn_atomno2) ,
				rsd2.atom_type_index(jj_atom3),
				rsd2.atom_type_index(jj_atom4) );
			if ( s_param_ij.is_null() && t_param_ij.is_null() ) continue;

			Vector const &xyz1( rsd1.xyz( resconn_atomno1 ) );
			Vector const &xyz2( rsd2.xyz( resconn_atomno2 ) );
			Vector const &xyz3( rsd2.xyz( jj_atom3 ) );
			Vector const &xyz4( rsd2.xyz( jj_atom4 ) );
			Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 );
			Real score;
			if ( !s_param_ij.is_null() ) {
				score = s_param_ij.energy( angle );
			} else {
				score = t_param_ij.energy( angle );
			}
			totalscore += score;

			if ( calc_deriv ) {
				Real dE_dtors;
				if ( !s_param_ij.is_null() ) {
					dE_dtors = weight*s_param_ij.deriv( angle );
				} else {
					dE_dtors = weight*t_param_ij.deriv( angle );
				}

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
	std::string atm3, std::string atm4,
	core::Real k1_in, core::Real k2_in, core::Real k3_in,
	core::Real f1_in, core::Real f2_in, core::Real f3_in
) {

	int n_tor_changed(0);
	utility::vector1< core::Size > wildcard_vec(1, 0);

	utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
	utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
	utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
	utility::vector1< core::Size > const & indices4 = name_index_map[atm4];
	core::Size multiplicity = indices1.size() * indices2.size() * indices3.size() * indices4.size();

	utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
	utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
	utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;
	utility::vector1< core::Size > const &indices4_loop = (indices4.size() == natomtypes) ? wildcard_vec : indices4;

	for ( auto i1 : indices1_loop ) {
		for ( auto i2 : indices2_loop ) {
			for ( auto i3 : indices3_loop ) {
				for ( auto i4 : indices4_loop ) {
					int64_t hashval = get_parameter_hash( i1, i2, i3, i4 );
					auto it = tors_lookup_.find( hashval );

					// use the MOST SPECIFIC potential
					if ( it == tors_lookup_.end() ) {
						TR << "Error: torsion type doesn't exist in database!" << std::endl;
					} else if ( multiplicity == tors_pot_[ it->second ].multiplicity() ) {

						if ( f1_in == 0 && f2_in == 0 && f3_in == 0 ) {
							f1_in = tors_pot_[ it->second ].get_params("f1");
							f2_in = tors_pot_[ it->second ].get_params("f2");
							f3_in = tors_pot_[ it->second ].get_params("f3");
						}
						tors_pot_[ it->second ] = GenTorsionParams( k1_in, k2_in, k3_in, f1_in, f2_in, f3_in, multiplicity );
						it = tors_lookup_.find( get_parameter_hash( i4, i3, i2, i1 ) );
						tors_pot_[ it->second ] = GenTorsionParams( k1_in, k2_in, k3_in, -f1_in, -f2_in, -f3_in, multiplicity );
						++n_tor_changed;
					}
				}
			}
		}
	}

	if ( n_tor_changed==0 ) {
		TR << "Warning: No params changed for torsion:" << atm1 << "-" << atm2 << "-" << atm3 << "-" << atm4 << std::endl;
	} else {
		TR << n_tor_changed << " torsion types changed for " << atm1 << "-" << atm2 << "-" << atm3 << "-" << atm4 << std::endl;
	}
}

void
GenericBondedPotential::modify_special_tors_params(
	std::string atm1, std::string atm2,
	std::string atm3, std::string atm4,
	core::Real k1_in, core::Real k2_in, core::Real k3_in,
	core::Real k4_in, core::Real k8_in,
	core::Real f1_in, core::Real f2_in, core::Real f3_in,
	core::Real f4_in, core::Real f8_in
) {
	int n_tor_changed(0);
	utility::vector1< core::Size > wildcard_vec(1, 0);

	utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
	utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
	utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
	utility::vector1< core::Size > const & indices4 = name_index_map[atm4];
	core::Size multiplicity = indices1.size() * indices2.size() * indices3.size() * indices4.size();

	utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
	utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
	utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;
	utility::vector1< core::Size > const &indices4_loop = (indices4.size() == natomtypes) ? wildcard_vec : indices4;

	for ( auto i1 : indices1_loop ) {
		for ( auto i2 : indices2_loop ) {
			for ( auto i3 : indices3_loop ) {
				for ( auto i4 : indices4_loop ) {
					int64_t hashval = get_parameter_hash( i1, i2, i3, i4 );
					auto it = special_tors_lookup_.find( hashval );

					// use the MOST SPECIFIC potential
					if ( it == special_tors_lookup_.end() ) {
						TR << "Error: torsion type doesn't exist in database!" << std::endl;
					} else if ( multiplicity == special_tors_pot_[ it->second ].multiplicity() ) {

						if ( f1_in == 0 && f2_in == 0 && f3_in == 0 && f4_in == 0 && f8_in == 0 ) {
							f1_in = special_tors_pot_[ it->second ].get_params("f1");
							f2_in = special_tors_pot_[ it->second ].get_params("f2");
							f3_in = special_tors_pot_[ it->second ].get_params("f3");
							f4_in = special_tors_pot_[ it->second ].get_params("f4");
							f8_in = special_tors_pot_[ it->second ].get_params("f8");
						}

						special_tors_pot_[ it->second ] = SpecialGenTorsionParams( k1_in, k2_in, k3_in, k4_in, k8_in, f1_in, f2_in, f3_in, f4_in, f8_in, multiplicity );
						it = special_tors_lookup_.find( get_parameter_hash( i4, i3, i2, i1 ) );
						special_tors_pot_[ it->second ] = SpecialGenTorsionParams( k1_in, k2_in, k3_in, k4_in, k8_in, -f1_in, -f2_in, -f3_in, -f4_in, -f8_in, multiplicity );
						++n_tor_changed;
					}
				}
			}
		}
	}
	if ( n_tor_changed==0 ) {
		TR << "Warning: No params changed for special torsion:" << atm1 << "-" << atm2 << "-" << atm3 << "-" << atm4 << std::endl;
	} else {
		TR << n_tor_changed << " special torsion types changed for " << atm1 << "-" << atm2 << "-" << atm3 << "-" << atm4 << std::endl;
	}
}


} // scoring
} // core

