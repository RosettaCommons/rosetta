// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/output/RNA_FragmentMonteCarloOutputter.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rna_denovo_output_RNA_FragmentMonteCarloOutputter_HH
#define INCLUDED_protocols_rna_denovo_output_RNA_FragmentMonteCarloOutputter_HH

#include <protocols/moves/Mover.hh>
#include <protocols/rna/denovo/output/RNA_FragmentMonteCarloOutputter.fwd.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <core/pose/rna/StubStubType.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <numeric/MathNTensor.fwd.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace rna {
namespace denovo {
namespace output {

class RNA_FragmentMonteCarloOutputter: public protocols::moves::Mover {

public:

	//constructor
	RNA_FragmentMonteCarloOutputter(  options::RNA_FragmentMonteCarloOptionsCOP options,
		core::pose::PoseCOP align_pose );

	//destructor
	~RNA_FragmentMonteCarloOutputter();

public:

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const{ return "RNA_FragmentMonteCarloOutputter"; }

	void
	output_running_info(
		core::Size const & r,
		core::Size const & i,
		core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP working_denovo_scorefxn );

	void
	finalize( core::scoring::ScoreFunctionCOP denovo_scorefxn );

	numeric::MathNTensorOP< core::Size, 6 > jump_histogram() const{ return jump_histogram_; }
	void set_jump_histogram( numeric::MathNTensorOP< core::Size, 6 > setting ) { jump_histogram_ = setting; }

private:

	void initialize( core::pose::PoseCOP align_pose );

	core::kinematics::RT
	get_output_jump_RT( core::pose::PoseCOP pose ) const;

	void
	get_output_jump_stub_stub( core::pose::Pose const & pose,
		core::kinematics::Stub & stub1,
		core::kinematics::Stub & stub2 ) const;

	void
	output_jump_information( core::pose::Pose const & pose );

private:

	options::RNA_FragmentMonteCarloOptionsCOP options_;

	utility::io::ozstream running_score_output_;
	numeric::MathNTensorOP< core::Size, 6 > jump_histogram_;
	utility::vector1< core::Real > jump_histogram_min_, jump_histogram_max_, jump_histogram_bin_width_;
	core::pose::rna::StubStubType output_stub_stub_type_;
	core::kinematics::RTOP reference_RT_;


};

} //output
} //denovo
} //rna
} //protocols

#endif
