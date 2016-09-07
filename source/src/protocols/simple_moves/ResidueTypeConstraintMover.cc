// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueTypeConstraintMover.cc
/// @brief Assigns a ResidueTypeConstraint to a pose.
/// @author Doo Nam Kim (doonam.kim@gmail.com) (All I did is just making a simple mover using Sarel's ResidueTypeConstraint)

#include <protocols/simple_moves/ResidueTypeConstraintMover.hh>
#include <protocols/simple_moves/ResidueTypeConstraintMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.ResidueTypeConstraintMover" );

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace scoring;
using namespace constraints;

using namespace utility::tag;

std::string
ResidueTypeConstraintMoverCreator::keyname() const
{
	return ResidueTypeConstraintMoverCreator::mover_name();
}

protocols::moves::MoverOP
ResidueTypeConstraintMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ResidueTypeConstraintMover );
}

std::string
ResidueTypeConstraintMoverCreator::mover_name()
{
	return "ResidueTypeConstraintMover";
}

// constructors
ResidueTypeConstraintMover::ResidueTypeConstraintMover()
: protocols::moves::Mover( ResidueTypeConstraintMoverCreator::mover_name() )
{
}

ResidueTypeConstraintMover::ResidueTypeConstraintMover( std::string const & type )
: protocols::moves::Mover(type)
{
}

// destructor
ResidueTypeConstraintMover::~ResidueTypeConstraintMover()= default;


void
ResidueTypeConstraintMover::apply( Pose & pose )
{
	int AA_name3_length_1 = std::count(AA_name3_.begin(), AA_name3_.end(), ',');
	std::string delimiter = ",";
	for ( Size resnum=1; resnum<=pose.total_residue(); ++resnum ) {
		int substr_index = 0;
		for ( int i =0; i<=(AA_name3_length_1); i++ ) {
			size_t pos = AA_name3_.find(delimiter);

			std::string favored_res = AA_name3_.substr(substr_index, pos);
			substr_index = substr_index + pos + 1;
			// essentially just same as "substr_index = substr_index + 4";

			ResidueTypeConstraintOP res_favor( new core::scoring::constraints::ResidueTypeConstraint(pose, resnum, favored_res, favor_bonus_) );

			pose.add_constraint(res_favor);
		}
	}
}

std::string
ResidueTypeConstraintMover::get_name() const {
	return ResidueTypeConstraintMoverCreator::mover_name();
}

protocols::moves::MoverOP ResidueTypeConstraintMover::clone() const { return protocols::moves::MoverOP( new protocols::simple_moves::ResidueTypeConstraintMover( *this ) ); }
protocols::moves::MoverOP ResidueTypeConstraintMover::fresh_instance() const { return protocols::moves::MoverOP( new ResidueTypeConstraintMover ); }

void
ResidueTypeConstraintMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	if ( tag->hasOption("AA_name3") ) {
		AA_name3_ = tag->getOption<std::string>("AA_name3"); // for example: ASP,GLU
	} else {
		TR << "please specifiy AA_name3 like SER,THR" << std::endl;
	}
	favor_bonus_ = tag->getOption<Real>("favor_bonus", 0.5);
	// since this mover does not necessarily favor native sequence only, I named it "favor_bonus" instead of "native_bonus" to make it more general
	// positively higher bonus gives more favorable selection to (a) specified residue(s)

}
} // moves
} // protocols
