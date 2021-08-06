// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/sap/SapDatabase.cc
/// @brief Hdf5 version of rosetta database
/// @details
/// @author Brian Coventry (bcov@uw.edu)



// unit headers
#include <core/pack/guidance_scoreterms/sap/SapDatabase.hh>
#include <core/pack/guidance_scoreterms/sap/util.hh>

// package headers

// project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>


// utility headers
#include <basic/database/open.hh>
#include <utility/pointer/memory.hh>
#include <utility/io/izstream.hh>

// basic headers
#include <basic/Tracer.hh>

// boost

// C++ headers

// C headers

static basic::Tracer TR( "core.sap.SapDatabase" );


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {


SapDatabase::SapDatabase() {

	load_hydrophobic_data();
	load_block_data();
	generate_max_sasa();
}

utility::pointer::shared_ptr< std::unordered_map< std::string, BlockParam > const >
SapDatabase::atomtype_to_block_param() const {
	return atomtype_to_block_param_;
}

std::pair< char, std::string >
SapDatabase::get_name1_name3( core::conformation::Residue const & res, bool warn ) const {
	std::pair< char, std::string > name1_name3( res.name1(), res.name3() );

	if ( ! res.type().is_l_aa() ) {
		// It's not strictly correct to interchange L and D like this but it should be close enough
		if ( res.type().is_d_aa() ) {
			core::chemical::AA l_aa = core::chemical::get_L_equivalent( res.aa() );
			name1_name3.first = core::chemical::oneletter_code_from_aa( l_aa );
			name1_name3.second = core::chemical::name_from_aa( l_aa );
		}
	}

	if ( name1_to_hydrophobic_.count( name1_name3.first ) ) {
		return name1_name3;
	}

	name1_name3.first = 0;
	name1_name3.second = "";
	if ( warn ) {
		TR.Warning << "Residue: " << res.seqpos() << " " << res.name3()
			<< " can't be used for SapScore but will still be used for SASA." << std::endl;
	}
	return name1_name3;
}

Real
SapDatabase::hydrophobic_weight( char aa ) const {
	return name1_to_hydrophobic_.at( aa );
}

Real
SapDatabase::max_sasa( char aa ) const {
	return name1_to_max_sasa_.at( aa );
}


// Development of hydrophobicity parameters to analyze proteins which bear post- or cotranslational modifications
// then you subtract 0.5 from scaled
void
SapDatabase::load_hydrophobic_data( ) {
	name1_to_hydrophobic_['A'] = 0.116;
	name1_to_hydrophobic_['C'] = 0.18;
	name1_to_hydrophobic_['D'] = -0.472;
	name1_to_hydrophobic_['E'] = -0.457;
	name1_to_hydrophobic_['F'] = 0.5;
	name1_to_hydrophobic_['G'] = 0.001;
	name1_to_hydrophobic_['H'] = -0.335;
	name1_to_hydrophobic_['I'] = 0.443;
	name1_to_hydrophobic_['K'] = -0.217;
	name1_to_hydrophobic_['L'] = 0.443;
	name1_to_hydrophobic_['M'] = 0.238;
	name1_to_hydrophobic_['N'] = -0.264;
	name1_to_hydrophobic_['P'] = 0.211;
	name1_to_hydrophobic_['Q'] = -0.249;
	name1_to_hydrophobic_['R'] = -0.5;
	name1_to_hydrophobic_['S'] = -0.141;
	name1_to_hydrophobic_['T'] = -0.05;
	name1_to_hydrophobic_['V'] = 0.325;
	name1_to_hydrophobic_['W'] = 0.378;
	name1_to_hydrophobic_['Y'] = 0.38;
}


void
SapDatabase::load_block_data( ) {

	atomtype_to_block_param_ = utility::pointer::make_shared< std::unordered_map< std::string, BlockParam > >();

	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/sap_sasa_calib.dat" );

	std::string tmp;
	stream >> tmp;
	runtime_assert( tmp == "atom_type" );
	stream >> tmp;
	runtime_assert( tmp == "full_block" );
	stream >> tmp;
	runtime_assert( tmp == "max_sasa" );
	stream >> tmp;
	runtime_assert( tmp == "no_block" );


	std::string atom_type;
	Real full_block;
	Real max_sasa;
	Real no_block;

	while ( !stream.eof() ) {

		stream >> atom_type;
		stream >> full_block;
		stream >> max_sasa;
		stream >> no_block;

		full_block /= SAP_BLOCK_STORE_SCALE;
		no_block /= SAP_BLOCK_STORE_SCALE;

		Real max_sasa_score = max_sasa;

		uint8_t byte_full_block = std::min<uint8_t>( 255, full_block );
		uint8_t byte_no_block = std::min<uint8_t>( 255, no_block );

		atomtype_to_block_param_->emplace( std::make_pair( atom_type, BlockParam( max_sasa_score, byte_no_block, byte_full_block ) ));

	}
}


void
SapDatabase::generate_max_sasa() {

	std::string name1s = "ACDEFGHIKLMNPQRSTVWY";

	utility::vector0<utility::vector1<Real>> const set_chis {
		utility::vector1<Real> {  },
		utility::vector1<Real> { 63.600,-60.000 },
		utility::vector1<Real> { -155.800,-6.900 },
		utility::vector1<Real> { 63.300,-177.200,-0.800 },
		utility::vector1<Real> { 60.700,91.100 },
		utility::vector1<Real> {  },
		utility::vector1<Real> { 62.300,74.600 },
		utility::vector1<Real> { -172.600,166.500 },
		utility::vector1<Real> { 62.600,-178.400,-179.600,-179.900 },
		utility::vector1<Real> { 70.700,165.700 },
		utility::vector1<Real> { 63.900,-172.400,72.000 },
		utility::vector1<Real> { -151.300,-28.100 },
		utility::vector1<Real> { 30.000,-33.900,24.824 },
		utility::vector1<Real> { 58.400,84.200,32.200 },
		utility::vector1<Real> { 62.500,176.900,176.600,85.600 },
		utility::vector1<Real> { 68.000,-180.000 },
		utility::vector1<Real> { -170.700,160.000 },
		utility::vector1<Real> { 63.300 },
		utility::vector1<Real> { 60.900,89.500 },
		utility::vector1<Real> { 61.500,90.400,180.000 },
		};


	std::stringstream ss;

	for ( Size iname1 = 0; iname1 < name1s.length(); iname1++ ) {
		char letter = name1s.at(iname1);

		ss << "AA" << letter << "AA";
	}

	pose::Pose pose;
	core::pose::make_pose_from_sequence( pose, ss.str(), "fa_standard" );

	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		pose.set_phi( seqpos, 180 );
		pose.set_psi( seqpos, 180 );
		pose.set_omega( seqpos, 180 );
	}

	for ( Size iname1 = 0; iname1 < name1s.length(); iname1++ ) {
		Size seqpos = iname1 * 5 + 3;

		for ( Size ichi = 1; ichi <= set_chis[iname1].size(); ichi++ ) {
			pose.set_chi( ichi, seqpos, set_chis[iname1][ichi] );
		}
	}


	select::residue_selector::ResidueSubset true_sub( pose.size(), true );
	core::id::AtomID_Map<Real> atom_sasa = sap_atom_sasa( pose, true_sub );


	for ( Size iname1 = 0; iname1 < name1s.length(); iname1++ ) {
		Size seqpos = iname1 * 5 + 3;

		conformation::Residue const & res = pose.residue(seqpos);

		Real sasa = 0;
		for ( Size iatom = 1; iatom <= res.natoms(); iatom++ ) {
			if ( res.atom_is_backbone( iatom ) ) continue;

			sasa += atom_sasa( seqpos, iatom );
		}
		name1_to_max_sasa_[res.name1()] = sasa;

	}

}

bool
SapDatabase::symm_debug() const {
	return symm_debug_;
}

bool
SapDatabase::symm_debug_force_map() const {
	return symm_debug_force_map_;
}


//namespaces
} //sap
} //guidance_scoreterms
} //pack
} //core


