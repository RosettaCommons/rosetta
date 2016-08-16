// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/raw_data/DecoyStruct.cc
///
/// @brief protein-specific "silent" file reader and writer for mini.
/// @author James Thompson, Monica Berrondo

// C++ Headers
#include <vector>
#include <string>
#include <map>

// mini headers


#include <core/io/raw_data/DecoyStruct.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <numeric/model_quality/rms.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


namespace core {
namespace io {
namespace raw_data {

using namespace ObjexxFCL::format;

DecoyStruct::DecoyStruct(
	core::pose::Pose const & pose,
	std::string tag, // = "empty_tag",
	bool fa // = false
) : fullatom_( fa ) {
	// tag information
	decoy_tag_ = tag;

	// conformation information
	sequence_ = pose.sequence();
	resize( pose.total_residue() );
	static const std::string important_atom = "CA";
	for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
		core::conformation::Residue resi = pose.residue(i);
		phi_  [i]     = resi.mainchain_torsion( 1 );
		psi_  [i]     = resi.mainchain_torsion( 2 );
		omega_[i]     = resi.mainchain_torsion( 3 );
		coords_[i]    = resi.xyz( important_atom );
		secstruct_[i] = pose.secstruct(i);
		if ( fullatom() ) {
			chi_[i] = resi.chi();
		} // if ( fullatom )
	} // for ( unsigned int i = 1; i <= pose.total_residue(); ++i )
} // DecoyStruct


/// @brief Resize this silent-struct to the appropriate number of residues.
void
DecoyStruct::resize(
	Size const nres_in
)
{
	nres_ = nres_in;
	secstruct_.resize( nres_ );
	phi_      .resize( nres_ );
	psi_      .resize( nres_ );
	omega_    .resize( nres_ );
	coords_   .resize( nres_ );
	chi_      .resize( nres_ );
}

// @brief Fill a Pose with the data in this DecoyStruct.
void DecoyStruct::fill_pose(
	core::pose::Pose & pose
) {
	using namespace core::chemical;
	ResidueTypeSetCOP residue_set;
	if ( fullatom() ) {
		residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	} else {
		residue_set = ChemicalManager::get_instance()->residue_type_set( CENTROID );
	}
	fill_pose( pose, *residue_set );
} // fill_pose

void DecoyStruct::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const& residue_set
) {

	// make_pose_match_sequence_( pose, sequence_, residue_set );
	core::pose::make_pose_from_sequence( pose, sequence_, residue_set );

	for ( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ) {
		pose.set_phi   ( seqpos, phi_   [seqpos] );
		pose.set_psi   ( seqpos, psi_   [seqpos] );
		pose.set_omega ( seqpos, omega_ [seqpos] );
		if ( fullatom_ ) {
			for ( Size j = 1; j <= pose.residue(seqpos).nchi(); ++j ) {
				pose.set_chi( j, seqpos, chi_[seqpos][j] );
			}
		}
		pose.set_secstruct( seqpos, secstruct_[seqpos] );
	}

} // fill_pose

/// @brief Print the conformation information contained in this object to the given ozstream.
void DecoyStruct::print_conformation( std::ostream& output ) const {
	using namespace ObjexxFCL;

	for ( Size i = 1; i <= nres_; ++i ) {
		output << I( 4, i ) << ' '
			<< secstruct_[i] << ' '
			<< F( 9, 3, phi_[i] )
			<< F( 9, 3, psi_[i] )
			<< F( 9, 3, omega_[i] )
			<< F( 9, 3, coords_[i].x() )
			<< F( 9, 3, coords_[i].y() )
			<< F( 9, 3, coords_[i].z() );
		if ( fullatom_ ) {
			for ( unsigned int chino = 1; chino <= 4; ++chino ) {
				Real chi_to_print = 0.0f;
				if ( chino <= chi_[i].size() ) {
					chi_to_print = chi_[i][chino];
				}
				output << F( 9, 3, chi_to_print );
			}
		}

		output << ' ' << decoy_tag_;
		output << std::endl;
	} // for ( Size i = 1; i <= nres; ++i )
} // print_conformation

Real DecoyStruct::get_debug_rmsd() {
	using namespace ObjexxFCL;

	pose::Pose temp_pose;
	FArray2D< Real > rebuilt_coords (3, coords_.size() ), original_coords( 3, coords_.size() );
	static std::string atom_name = "CA";

	// build temp_pose from coordinates
	fill_pose( temp_pose );

	for ( Size i = 1; i <= temp_pose.total_residue(); ++i ) {
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			rebuilt_coords (k,i) = temp_pose.residue(i).xyz( atom_name )[k-1];
			original_coords(k,i) = coords_[i][k-1];
		}
	}

	Real rmsd = numeric::model_quality::rms_wrapper( temp_pose.total_residue(), rebuilt_coords, original_coords );
	return rmsd;
}

} // namespace silent
} // namespace io
} // namespace core
