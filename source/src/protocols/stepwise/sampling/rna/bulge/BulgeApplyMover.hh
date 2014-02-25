// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/bulge/BulgeApplyMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_rna_bulge_BulgeApplyMover_HH
#define INCLUDED_protocols_stepwise_sampling_rna_bulge_BulgeApplyMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/sampling/rna/bulge/BulgeApplyMover.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace bulge {

	class BulgeApplyMover: public protocols::moves::Mover {

	public:

		//constructor
		BulgeApplyMover( Size const moving_res );

		//destructor
		~BulgeApplyMover();

	public:

		virtual void apply( Pose & );

		virtual std::string get_name() const { return "BulgeApplyMover"; }

	private:

		Size const moving_res_;

	};

} //bulge
} //rna
} //sampling
} //stepwise
} //protocols

#endif
