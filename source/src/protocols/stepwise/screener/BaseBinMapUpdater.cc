// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/BaseBinMapUpdater.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/BaseBinMapUpdater.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.BaseBinMapUpdater" );

using namespace protocols::rotamer_sampler;
using namespace protocols::rotamer_sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	BaseBinMapUpdater::BaseBinMapUpdater( BaseBinMap & base_bin_map ):
		base_bin_map_( base_bin_map )
	{}

	//Destructor
	BaseBinMapUpdater::~BaseBinMapUpdater()
	{}

	// diagnostics.
	void
	BaseBinMapUpdater::get_update( rotamer_sampler::RotamerBaseOP sampler ){

		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ){
			RigidBodyRotamerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyRotamerWithResidueList * >( sampler.get() ) );
			update_base_bin_map( rigid_body_rotamer_with_copy_dofs.get_rigid_body_values() );
			return;
		}
		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyRotamerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyRotamerWithResidueAlternatives * >( sampler.get() ) );
			update_base_bin_map( rigid_body_rotamer_with_residue_alternatives.get_rigid_body_values() );
			return;
		}

		runtime_assert( sampler->type() == RIGID_BODY );
		RigidBodyRotamer & rigid_body_rotamer = *( static_cast< RigidBodyRotamer * >( sampler.get() ) );
		update_base_bin_map( rigid_body_rotamer.get_values() );
	}


	/////////////////////////////////////////////////////////////////////////////////////
	// diagnostics
	void
	BaseBinMapUpdater::update_base_bin_map( BaseBin const & base_bin ){
		std::map< BaseBin, int, compare_base_bin > ::const_iterator it = base_bin_map_.find( base_bin );
		if ( it == base_bin_map_.end() )	base_bin_map_[base_bin] = 0;
		base_bin_map_[base_bin] ++;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	// diagnostics
	void
	BaseBinMapUpdater::update_base_bin_map( utility::vector1< Real > const & rigid_body_values ){
		BaseBin base_bin;
		base_bin.centroid_x  = rigid_body_values[6];
		base_bin.centroid_y  = rigid_body_values[5];
		base_bin.centroid_z  = rigid_body_values[4];
		base_bin.euler_alpha = rigid_body_values[3];
		base_bin.euler_z     = rigid_body_values[2];
		base_bin.euler_gamma = rigid_body_values[1];
		update_base_bin_map( base_bin );
	}



} //screener
} //stepwise
} //protocols
