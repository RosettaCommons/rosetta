// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Symmetry Claimer
/// @brief  Basic Claimer for making Symmetry Claims
/// @author Justin Porter, Tatjana Braun

// Unit Headers
#include <protocols/topology_broker/SymmetryClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/SymmetryClaim.hh>
#include <protocols/topology_broker/claims/SequenceClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>
#include <protocols/topology_broker/TopologyBroker.hh>

// Project Headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/util/SwitchResidueTypeSet.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

//// C++ headers
#include <sstream>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

SymmetryClaimer::SymmetryClaimer() : symm_data_ ( /* NULL */ ){}

SymmetryClaimer::SymmetryClaimer( SymmetryClaimer const& src ) :
	TopologyClaimer( src ),
	symm_data_ ( src.symm_data_ )
{}

bool SymmetryClaimer::read_tag( std::string tag, std::istream& in ){

	if ( tag == "symmetry_definition" ) {
		std::string symdef_filename;
		in >> symdef_filename;
		tr.Debug << "Symmetry Definition filename read: " << symdef_filename << std::endl;
		symm_data_ = core::conformation::symmetry::SymmDataOP( new core::conformation::symmetry::SymmData() );
		symm_data_->read_symmetry_data_from_file(symdef_filename);
		return true;
	}

	return TopologyClaimer::read_tag( tag, in );
}

void SymmetryClaimer::generate_symmetry_claims( claims::SymmetryClaims& symm_claims ){
	symm_claims.push_back( claims::SymmetryClaimOP( new claims::SymmetryClaim( get_self_weak_ptr(), symm_data_,
		symm_data_->get_symmetry_name(),
		claims::DofClaim::INIT ) ) );
}

void SymmetryClaimer::symmetry_duplicate( claims::DofClaims& pre_accepted,
	core::pose::Pose& pose ){
	asymmetric_res_ = pose.total_residue();

	utility::vector1< claims::SequenceClaimOP > original_claims;
	for ( claims::DofClaims::iterator claim = pre_accepted.begin();
			claim != pre_accepted.end(); ++claim ) {
		claims::SequenceClaimOP seq_ptr( utility::pointer::dynamic_pointer_cast< claims::SequenceClaim >( *claim ) );
		runtime_assert( seq_ptr != 0 );
		original_claims.push_back( seq_ptr );
	}

	//Copy sequence claims (j-1) times to get j subunits in total.
	for ( Size j=1; j<symm_data_->get_subunits(); ++j ) {
		//Copy each pre accepted sequence claim (first n_asymm_claims in pre_accepted) and add them to pre_accepted
		for ( Size i=1; i<= original_claims.size(); ++i ) {
			claims::SequenceClaimOP old_claim = original_claims.at(i);
			std::ostringstream new_label_stream;
			new_label_stream << old_claim->label() << ":Symm" << j;

			pre_accepted.push_back( claims::DofClaimOP( new claims::SequenceClaim( get_self_weak_ptr(),
				old_claim->annotated_sequence(),
				new_label_stream.str() ) ) );
		}
	}

	//Use annotated sequence to duplicate
	std::string sequence = "";
	for ( Size i=1; i<=symm_data_->get_subunits(); ++i ) {
		sequence += pose.annotated_sequence();
	}

	core::pose::make_pose_from_sequence( pose, sequence,
		*( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )));

	//Build symmetry-related virtual residues
	for ( Size i=1; i<=symm_data_->get_num_virtual(); ++i ) {
		// create the new residue
		core::chemical::ResidueTypeSetCOP rsd_set( pose.conformation().residue(1).residue_type_set() );
		core::conformation::ResidueOP rsd( core::conformation::ResidueFactory::create_residue( rsd_set->name_map( "VRT" ) ) );

		std::string tag = symm_data_->get_virtual_num_to_id().find(i)->second;

		core::conformation::symmetry::VirtualCoordinate virt_coord( symm_data_->get_virtual_coordinates().find(tag)->second );
		rsd->set_xyz( "ORIG", virt_coord.get_origin() );
		rsd->set_xyz( "X", virt_coord.get_x().normalized() + virt_coord.get_origin() );
		rsd->set_xyz( "Y", virt_coord.get_y().normalized() + virt_coord.get_origin() );
		//append it to the end of the monomer i
		pose.append_residue_by_jump( *rsd, pose.total_residue() );

		pre_accepted.push_back( claims::DofClaimOP( new claims::SequenceClaim( get_self_weak_ptr(), "X", tag ) ));
	}

}

}
}
