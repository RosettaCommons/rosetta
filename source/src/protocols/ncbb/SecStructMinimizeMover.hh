// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ncbb/SecStructMinimizeMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_ncbb_SecStructMinimizeMover_HH
#define INCLUDED_protocols_ncbb_SecStructMinimizeMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/ncbb/SecStructMinimizeMover.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace ncbb {

class SecStructMinimizeMover: public protocols::moves::Mover {

public:

	//constructor
	SecStructMinimizeMover():
		constrain_( false ),
		score_fxn_( new core::scoring::ScoreFunction )
	{
		Mover::type("SecStructMinimizeMover");
	}

	SecStructMinimizeMover(
		core::scoring::ScoreFunctionOP const & score_fxn,
		std::string const & dihedral_pattern,
		std::string const & alpha_beta_pattern
	):
		constrain_( false ),
		score_fxn_( score_fxn ),
		dihedral_pattern_( dihedral_pattern ),
		alpha_beta_pattern_( alpha_beta_pattern )
	{}

	//destructor
	~SecStructMinimizeMover();

	protocols::moves::MoverOP fresh_instance() const { return SecStructMinimizeMoverOP( new SecStructMinimizeMover ); }
	protocols::moves::MoverOP clone() const { return protocols::moves::MoverOP( new SecStructMinimizeMover(
		score_fxn_,
		dihedral_pattern_,
		alpha_beta_pattern_ ) ); }

	void add_dihedral_constraints_to_pose( Pose & pose, Size number_dihedral_sets, utility::vector1< char > uniqs );

	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

public:

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const { return "SecStructMinimizeMover"; }

	void set_scorefunction( core::scoring::ScoreFunctionOP const & score_fxn ) {
		score_fxn_ = score_fxn;
	}

	core::scoring::ScoreFunctionOP
	get_scorefunction() {
		return score_fxn_;
	}

	void set_dihedral_pattern( std::string const & dihedral_pattern ) {
		dihedral_pattern_ = dihedral_pattern;
	}

	std::string
	get_dihedral_pattern() {
		return dihedral_pattern_;
	}

	void set_alpha_beta_pattern( std::string const & alpha_beta_pattern ) {
		alpha_beta_pattern_ = alpha_beta_pattern;
	}

	std::string
	get_alpha_beta_pattern() {
		return alpha_beta_pattern_;
	}

	void set_dihedrals( utility::vector1< core::Real > const & dihedrals ) {
		dihedrals_ = dihedrals;
	}

	utility::vector1< core::Real >
	get_dihedrals() {
		return dihedrals_;
	}

	void set_constrain( bool const constrain ) {
		constrain_ = constrain;
	}

	bool
	get_constrain() {
		return constrain_;
	}

private:

	bool constrain_;
	core::scoring::ScoreFunctionOP score_fxn_;
	std::string dihedral_pattern_;
	std::string alpha_beta_pattern_;
	utility::vector1< core::Real > dihedrals_;
};

} //ncbb
} //protocols

#endif
