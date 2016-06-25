// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_Minimizer_HH
#define INCLUDED_protocols_rna_RNA_Minimizer_HH

#include <protocols/moves/Mover.hh>
#include <protocols/farna/options/RNA_MinimizerOptions.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace protocols::farna::options;

namespace protocols {
namespace farna {
namespace movers {

/// @brief The RNA de novo structure modeling protocol
class RNA_Minimizer: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_Minimizer( RNA_MinimizerOptionsCOP options = 0 );

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new RNA_Minimizer(*this) );
	}

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	virtual void show(std::ostream & output=std::cout) const;

	void set_perform_minimizer_run( bool const setting ){ perform_minimizer_run_ = setting; }

	void set_options( RNA_MinimizerOptionsCOP setting ){ options_ = setting; }

	void set_include_default_linear_chainbreak( bool const setting){ include_default_linear_chainbreak_ = setting; }

	void set_min_type( std::string const & setting){ min_type_ = setting; }

	void
	set_atom_level_domain_map(toolbox::AtomLevelDomainMapCOP atom_level_domain_map  ) {
		atom_level_domain_map_input_ = atom_level_domain_map;
	}

	void
	set_score_function( core::scoring::ScoreFunctionCOP scorefxn );

	core::scoring::ScoreFunctionOP const &
	score_function() const{ return scorefxn_; }

	core::scoring::ScoreFunctionOP
	clone_scorefxn() const{ return scorefxn_->clone(); }

private:

	// Make this a Mover?
	void
	o2prime_trials( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP const & scorefxn ) const;

	void
	setup_movemap( core::kinematics::MoveMap & mm, core::pose::Pose & pose );

	void
	update_atom_level_domain_map_with_extra_minimize_res( core::pose::Pose const & pose );

	void
	update_atom_level_domain_map_to_move_rosetta_library_chunks();

private:

	RNA_MinimizerOptionsCOP options_;
	core::Real const coord_sdev_;
	core::Real const coord_cst_weight_;
	bool perform_minimizer_run_;
	bool include_default_linear_chainbreak_;
	bool close_loops_;
	std::string min_type_;

	toolbox::AtomLevelDomainMapCOP atom_level_domain_map_input_;
	toolbox::AtomLevelDomainMapOP atom_level_domain_map_;

	core::scoring::ScoreFunctionOP scorefxn_;
}; // class RNA_Minimizer

std::ostream &operator<< ( std::ostream &os, RNA_Minimizer const &mover );

} //movers
} //farna
} //protocols

#endif
