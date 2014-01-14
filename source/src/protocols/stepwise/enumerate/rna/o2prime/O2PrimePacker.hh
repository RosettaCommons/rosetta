// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/o2prime/O2PrimePacker.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#ifndef INCLUDED_protocols_stepwise_rna_O2PrimePacker_HH
#define INCLUDED_protocols_stepwise_rna_O2PrimePacker_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/enumerate/rna/o2prime/O2PrimePacker.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

using namespace core;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace o2prime {

	class O2PrimePacker: public utility::pointer::ReferenceCount {

	public:

		//constructor
		O2PrimePacker( pose::Pose const & pose,
									 core::scoring::ScoreFunctionCOP const & scorefxn,
									 utility::vector1< core::Size > moving_res );

		//destructor
		~O2PrimePacker();

		void
		sample_o2prime_hydrogen();

		pose::Pose & pose();

		void
		copy_all_o2prime_torsions( core::pose::Pose & mod_pose ) const;

		void
		set_use_green_packer( bool const & setting ){ use_green_packer_ = setting ; }

		void set_partition_definition( ObjexxFCL::FArray1D < bool > const & setting ){ partition_definition_ = setting; }

	private:

		void
		initialize_o2prime_pack_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn );

		void
		initialize_o2prime_packer_task();

		void
		initialize_o2prime_green_packer();

	private:

		Pose const pose_with_original_HO2prime_torsion_;
		utility::vector1< Size > const moving_res_;
		Pose o2prime_pack_pose_;
		bool use_green_packer_;

		core::pack::task::PackerTaskOP o2prime_pack_task_;
		protocols::simple_moves::GreenPackerOP o2prime_green_packer_;

		core::scoring::ScoreFunctionOP o2prime_pack_scorefxn_;

		ObjexxFCL::FArray1D < bool > partition_definition_; // needed by green packer

	};

} //o2prime
} //rna
} //enumerate
} //stepwise
} //protocols

#endif
