// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_AddDeleteMonteCarlo.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_monte_carlo_RNA_AddDeleteMonteCarlo_hh
#define INCLUDED_protocols_stepwise_monte_carlo_RNA_AddDeleteMonteCarlo_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_AddDeleteMonteCarlo.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_TorsionMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_O2PrimeMover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>


namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_AddDeleteMonteCarlo: public protocols::moves::Mover {
public:


	RNA_AddDeleteMonteCarlo(  mover::AddOrDeleteMoverOP rna_add_or_delete_mover,
		RNA_TorsionMoverOP     rna_torsion_mover,
		RNA_O2PrimeMoverOP      rna_o2prime_mover,
		core::scoring::ScoreFunctionOP scorefxn );

	//destructor -- necessary? -- YES destructors are necessary.
	~RNA_AddDeleteMonteCarlo();

	// Undefinded, commenting out to fix PyRosetta build  void apply( core::pose::Pose & pose, Size const res_to_delete, protocols::stepwise::monte_carlo::MovingResidueCase const moving_residue_case  );

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void set_kT( core::Real const & setting ){ kT_ = setting; }

	void set_sample_range_small( core::Real const setting ){ sample_range_small_ = setting; }

	void set_sample_range_large( core::Real const setting ){ sample_range_large_ = setting; }

	void set_num_cycles( Size const setting ){ num_cycles_ = setting; }

	void set_output_period( Size const setting ){ output_period_ = setting; }

	void set_do_add_delete( bool const setting ){ do_add_delete_ = setting; }

	void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

private:

	void
	initialize_next_suite_atoms();

	void
	output_silent_file( core::pose::Pose & pose, Size const count );

private:

	mover::AddOrDeleteMoverOP rna_add_or_delete_mover_;
	RNA_TorsionMoverOP rna_torsion_mover_;
	RNA_O2PrimeMoverOP rna_o2prime_mover_;

	core::scoring::ScoreFunctionOP scorefxn_;

	Size num_cycles_, output_period_;
	core::Real sample_range_small_, sample_range_large_, kT_;

	bool do_add_delete_;
	std::string silent_file_;
	core::io::silent::SilentFileDataOP silent_file_data_;

	utility::vector1< std::string > next_suite_atoms_;

};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
