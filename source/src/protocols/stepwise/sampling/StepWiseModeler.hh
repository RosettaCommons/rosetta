// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/StepWiseModeler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_StepWiseModeler_HH
#define INCLUDED_protocols_stepwise_sampling_StepWiseModeler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/sampling/StepWiseModeler.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Modeler.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModeler.fwd.hh>
#include <core/pose/Pose.fwd.hh>

////////////////////////////////////////////////////////////////////////////////////
// A wrapper around protein or RNA modeler, as appropriate.
//
// Note that later we should be able to deprecate separate protein and RNA modeling,
//  and unify into a single function -- and that will be important for
//  modeling general RNA/DNA/protein/small-molecule problems.
//
//        -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {

	class StepWiseModeler: public protocols::moves::Mover  {

	public:

		//constructor
		StepWiseModeler( rna::StepWiseRNA_ModelerOP stepwise_rna_modeler,
										 protein::StepWiseProteinModelerOP stepwise_protein_modeler );

		StepWiseModeler( StepWiseModeler const & src );

		//destructor
		~StepWiseModeler();

	public:

		StepWiseModelerOP clone_modeler() const;

		StepWiseModeler & operator=( StepWiseModeler const & src );

		virtual void apply( pose::Pose & pose );

		virtual std::string get_name() const{ return "StepWiseModeler"; }

		void set_moving_res_and_reset( Size const moving_res );

		void set_native_pose( pose::PoseCOP );

		void set_working_minimize_res( utility::vector1< Size > const & setting );

		void set_skip_sampling( bool const & setting );

	private:

		rna::StepWiseRNA_ModelerOP stepwise_rna_modeler_;
		protein::StepWiseProteinModelerOP stepwise_protein_modeler_;
		Size moving_res_;
	};

} //sampling
} //stepwise
} //protocols

#endif
