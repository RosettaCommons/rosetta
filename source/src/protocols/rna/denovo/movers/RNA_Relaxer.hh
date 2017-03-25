// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file loopRNA_relaxer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_Relaxer_hh
#define INCLUDED_protocols_rna_RNA_Relaxer_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

//Oooh.

//// C++ headers
#include <string>
#include <vector>

#include <protocols/rna/denovo/movers/RNA_FragmentMover.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rna {
namespace denovo {
namespace movers {


typedef utility::pointer::shared_ptr< RNA_FragmentMover > RNA_FragmentMoverOP;
typedef utility::pointer::shared_ptr< RNA_Minimizer > RNA_MinimizerOP;

/// @brief The RNA de novo structure modeling protocol
class RNA_Relaxer: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_Relaxer( RNA_FragmentMoverOP rna_fragment_mover, RNA_MinimizerOP rna_minimizer );

	~RNA_Relaxer();

	/// @brief Clone this object
	// virtual RNA_Relaxer* clone() const {
	//  return new RNA_Relaxer(*this);
	// }

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void
	simple_rmsd_cutoff_relax( bool const setting ){ simple_rmsd_cutoff_relax_ = setting; }

private:

	void
	make_fragment_moves( core::pose::Pose & pose );

	void
	find_fragment_by_simple_rmsd_cutoff( core::pose::Pose & pose );

	void
	lores_monte_carlo( core::pose::Pose & pose );

	//data
	RNA_FragmentMoverOP rna_fragment_mover_;
	RNA_MinimizerOP rna_minimizer_;

	Size const relax_cycles_;
	Size const max_frag_size_;
	Size const num_find_fragment_tries_;
	core::Real const rmsd_min_cutoff_;
	core::Real const rmsd_max_cutoff_;

	bool simple_rmsd_cutoff_relax_;

}; // class RNA_Relaxer


} //movers
} //denovo
} //rna
} //protocols

#endif
