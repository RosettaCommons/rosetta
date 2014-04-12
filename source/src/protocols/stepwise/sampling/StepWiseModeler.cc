// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/StepWiseModeler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/StepWiseModeler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Modeler.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModeler.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.StepWiseModeler" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {

	//Constructor
	StepWiseModeler::StepWiseModeler( rna::StepWiseRNA_ModelerOP stepwise_rna_modeler,
									 protein::StepWiseProteinModelerOP stepwise_protein_modeler ):
		stepwise_rna_modeler_( stepwise_rna_modeler ),
		stepwise_protein_modeler_( stepwise_protein_modeler ),
		moving_res_( 0 )
	{
	}

	//Destructor
	StepWiseModeler::~StepWiseModeler()
	{}

	StepWiseModeler::StepWiseModeler( StepWiseModeler const & src ):
		Mover( src )
	{
		*this = src;
	}

	StepWiseModelerOP
	StepWiseModeler::clone_modeler() const {
		return new StepWiseModeler( *this );
	}

	StepWiseModeler &
	StepWiseModeler::operator=( StepWiseModeler const & src ){
		stepwise_rna_modeler_ = src.stepwise_rna_modeler_->clone_modeler();
		stepwise_protein_modeler_ = src.stepwise_protein_modeler_->clone_modeler();
		return *this;
	}


	//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseModeler::apply( pose::Pose & pose ){
		Size const example_res = ( moving_res_ > 0 ) ? moving_res_ : 1;
		if ( pose.residue_type( example_res ).is_RNA() ){
			stepwise_rna_modeler_->apply( pose );
		} else {
			runtime_assert( pose.residue_type( example_res ).is_protein() );
			stepwise_protein_modeler_->apply( pose );
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseModeler::set_moving_res_and_reset( Size const moving_res ){
		moving_res_ = moving_res;
		stepwise_rna_modeler_->set_moving_res_and_reset( moving_res );
		stepwise_protein_modeler_->set_moving_res_and_reset( moving_res );
	}

	void
	StepWiseModeler::set_native_pose( pose::PoseCOP native_pose ){
		stepwise_rna_modeler_->set_native_pose( native_pose );
		stepwise_protein_modeler_->set_native_pose( native_pose );
	}

	void
	StepWiseModeler::set_working_minimize_res( utility::vector1< Size > const & working_minimize_res ){
		stepwise_rna_modeler_->set_working_minimize_res( working_minimize_res );
		stepwise_protein_modeler_->set_working_minimize_res( working_minimize_res );
	}

	void
	StepWiseModeler::set_skip_sampling( bool const & skip_sampling ){
		stepwise_rna_modeler_->set_skip_sampling( skip_sampling );
		stepwise_protein_modeler_->set_skip_sampling( skip_sampling );
	}

} //sampling
} //stepwise
} //protocols
