// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MiniRelax.hh
/// @author James Thompson

#ifndef INCLUDED_protocols_relax_MiniRelax_hh
#define INCLUDED_protocols_relax_MiniRelax_hh

#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/MiniRelax.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace relax {

class MiniRelax : public RelaxProtocolBase {

public:
	typedef RelaxProtocolBase parent;

public:
	MiniRelax( core::scoring::ScoreFunctionOP scorefxn_in );
	MiniRelax( MiniRelax const & other );

	~MiniRelax() override;
	protocols::moves::MoverOP clone() const override;
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	core::scoring::ScoreFunctionOP scorefxn_;
};

}
} // protocols

#endif
