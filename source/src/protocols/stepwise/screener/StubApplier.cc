// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/StubApplier.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/StubApplier.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.StubApplier" );

using namespace core;
using namespace protocols::rotamer_sampler;
using namespace protocols::rotamer_sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	StubApplier::StubApplier( kinematics::Stub & stub ):
		stub_( stub )
	{}

	//Destructor
	StubApplier::~StubApplier()
	{}

	///////////////////////////////////////////////////////////////////
	void
	StubApplier::get_update( rotamer_sampler::RotamerBaseOP sampler ){

		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ){
			RigidBodyRotamerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyRotamerWithResidueList * >( sampler.get() ) );
			stub_ = rigid_body_rotamer_with_copy_dofs.get_stub();
			return;
		} else if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyRotamerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyRotamerWithResidueAlternatives * >( sampler.get() ) );
			stub_ = rigid_body_rotamer_with_residue_alternatives.get_stub();
			return;
		}

   runtime_assert( sampler->type() == RIGID_BODY );
	 RigidBodyRotamer & rigid_body_rotamer = *( static_cast< RigidBodyRotamer * >( sampler.get() ) );
	 stub_ = rigid_body_rotamer.get_stub();

	}

} //screener
} //stepwise
} //protocols
