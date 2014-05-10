// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMoverBase.cc
/// @brief implementation of PerturbChiSidechainMover class and functions
/// @author Oliver Lange ( oliver.lange@tum.de )


#include <protocols/simple_moves/sidechain_moves/PerturbChiSidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbChiSidechainMoverCreator.hh>

// Procols Headers
#include <basic/datacache/DataMap.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/prof.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Utility
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// C++ Headers
#include <sstream>
#include <fstream>
#include <utility/fixedsizearray1.hh>
using namespace core;
using namespace core::pose;

static numeric::random::RandomGenerator RG(18615125);
static basic::Tracer TR("protocols.simple_moves.sidechain_moves.PerturbChiSidechainMover");

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

using namespace chemical;
using namespace conformation;


std::string
PerturbChiSidechainMoverCreator::keyname() const {
	return PerturbChiSidechainMoverCreator::mover_name();
}

protocols::moves::MoverOP
PerturbChiSidechainMoverCreator::create_mover() const {
	return new PerturbChiSidechainMover;
}

std::string
PerturbChiSidechainMoverCreator::mover_name() {
	return "PerturbChiSidechain";
}



PerturbChiSidechainMover::PerturbChiSidechainMover() {
	protocols::moves::Mover::type( "PerturbChiSidechain" );
	set_defaults();
}

PerturbChiSidechainMover::PerturbChiSidechainMover(
	pack::dunbrack::RotamerLibrary const & rotamer_library
) :	Parent( rotamer_library ) {
	set_defaults();
}

PerturbChiSidechainMover::PerturbChiSidechainMover(
  PerturbChiSidechainMover const & mover
) : Parent ( mover ) {
	set_defaults();
}

protocols::moves::MoverOP
PerturbChiSidechainMover::clone() const {
	return new protocols::simple_moves::sidechain_moves::PerturbChiSidechainMover(*this);
}

void
PerturbChiSidechainMover::set_defaults() {
	magnitude_ = 10;
	gaussian_ = false;
}

void
PerturbChiSidechainMover::parse_my_tag(
  utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	pose::Pose const & /*pose*/
) {
	magnitude_ = tag->getOption<Real>( "magnitude", magnitude_ );
	gaussian_ = tag->getOption<bool>( "gaussian", gaussian_ );
}

std::string
PerturbChiSidechainMover::get_name() const {
	return "PerturbChiSidechainMover";
}

void
PerturbChiSidechainMover::make_chi_move(
	conformation::Residue const&,
	ChiVector const& old_chi,
	ChiVector& new_chi
) {
	new_chi.resize( old_chi.size() );
	for ( Size i = 1; i <= old_chi.size(); i++) {
		if ( !gaussian_ ) {
			Real rand = RG.uniform();
			new_chi[ i ] = basic::periodic_range( (( 2.0*rand-1.0 )*magnitude_ + old_chi[ i ]) , 360.0 );
		} else {
			new_chi[ i ] = basic::periodic_range( old_chi[i] + RG.gaussian()*magnitude_, 360.0 );
		}
	}
}

} // sidechain_moves
} // simple_moves
} // protocols
