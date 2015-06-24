// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMoverBase.cc
/// @brief implementation of PerturbRotamerSidechainMover class and functions
/// @author Oliver Lange ( oliver.lange@tum.de )


#include <protocols/simple_moves/sidechain_moves/PerturbRotamerSidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbRotamerSidechainMoverCreator.hh>
#include <protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMover.hh>
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
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
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

static thread_local basic::Tracer TR( "protocols.simple_moves.sidechain_moves.PerturbRotamerSidechainMover" );

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

using namespace chemical;
using namespace conformation;


std::string
PerturbRotamerSidechainMoverCreator::keyname() const {
	return PerturbRotamerSidechainMoverCreator::mover_name();
}

protocols::moves::MoverOP
PerturbRotamerSidechainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PerturbRotamerSidechainMover );
}

std::string
PerturbRotamerSidechainMoverCreator::mover_name() {
	return "PerturbRotamerSidechain";
}


PerturbRotamerSidechainMover::PerturbRotamerSidechainMover() {
	protocols::moves::Mover::type( "PerturbRotamerSidechain" );
	set_defaults();
}

PerturbRotamerSidechainMover::PerturbRotamerSidechainMover(
	pack::dunbrack::RotamerLibrary const & rotamer_library
) :	Parent( rotamer_library ) {
	protocols::moves::Mover::type( "PerturbRotamerSidechain" );
	set_defaults();
}

PerturbRotamerSidechainMover::PerturbRotamerSidechainMover(
  PerturbRotamerSidechainMover const & mover
) : Parent ( mover ) {
	set_defaults();
}

protocols::moves::MoverOP
PerturbRotamerSidechainMover::clone() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::sidechain_moves::PerturbRotamerSidechainMover(*this) );
}

void
PerturbRotamerSidechainMover::set_defaults() {
	temperature_ = 1;
}

void
PerturbRotamerSidechainMover::parse_my_tag(
  utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	pose::Pose const & /*pose*/
) {

}

std::string
PerturbRotamerSidechainMover::get_name() const {
	return "PerturbRotamerSidechainMover";
}

void
PerturbRotamerSidechainMover::make_chi_move(
	Residue const& residue,
	ChiVector const& old_chi,
	ChiVector& new_chi
) {

	RotamerList rotamers;
	build_rotamer_list( residue, false, rotamers );

	// find the rotamer that has the highest probability of proposing the previous chi angles
	Real max_rot_prob = 0;
	Size max_rot_num = 0;
	//Real rot_prob_normalize (0);

	for (Size ii = 1; ii <= rotamers.size() ; ++ii) {
		//TR << "rotamer.size is " << rotamers.size() << "  rotamer number is " << ii << std::endl;
		Real rot_prob( rotamers[ii].chi_probability( old_chi, temperature_ )*rotamers[ ii ].probability());
		//Real rot_prob( rotamers[ii].chi_probability( old_chi, temperature_ ));
		if ( rot_prob > max_rot_prob ) {
			max_rot_num = ii;
			max_rot_prob = rot_prob;
		}
		//TR << "number is " << ii << "  rot_prob is " << rot_prob << std::endl;
	}

	new_chi=old_chi;
	//TR << " max_rot_num is " << max_rot_num << std::endl;
	rotamers[ max_rot_num ].assign_random_chi( new_chi, numeric::random::rg(), temperature_ );
	//TR << " new chi is " << new_chi << std::endl;
}

///all angles in degree
Real PerturbRotamerSidechainMover::compute_proposal_density(
  Residue const & new_residue,
	Size const,
	chemical::ResidueType const &,
	ChiVector const& old_chi /* in degree */
) const {

	utility::vector1<Real> const & new_chi( new_residue.chi() );

	RotamerList rotamers;
	build_rotamer_list( new_residue, false /*no filtering*/, rotamers );

	//Real rot_density;
	Real within_rot_density;
	//compute_rotdensities( rotamers, old_chi, new_chi, within_rot_density );
	compute_rotdensities( rotamers, old_chi, new_chi, within_rot_density );
	return within_rot_density;
}

///all angles in degree
void
PerturbRotamerSidechainMover::build_rotamer_list(
	Residue const& residue,
	bool filter_low_probabilities,
	RotamerList& rotamers
) const {
	using namespace pack::dunbrack;
	using namespace pack::rotamers;

	rotamers.clear();
	rotamers.reserve( 200 ); //no idea but that should be plenty.

	SingleResidueRotamerLibraryCOP residue_rotamer_library(
		core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( residue.type() )
	);

	SingleResidueDunbrackLibraryCOP residue_dunbrack_library(
		utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const >( residue_rotamer_library )
	);

	if ( !residue_dunbrack_library ) return;

	//amw
	//get phi-psi dependent rotamer library
	Real const phi( residue_dunbrack_library->get_phi_from_rsd( residue ) );
	Real const psi( residue_dunbrack_library->get_psi_from_rsd( residue ) );
	utility::fixedsizearray1< Real, 5 > bbs;
	//for ( Size i = 1; i <= residue.type().mainchain_torsions().size() - 1; ++i )
	bbs[ 1 ] = phi;
	bbs[2] = psi; 
	//residue_dunbrack_library->get_bb_from_rsd( i, residue );
	//RotamerList rotamers_raw( residue_dunbrack_library->get_all_rotamer_samples( phi, psi ) );
	RotamerList rotamers_raw( residue_dunbrack_library->get_all_rotamer_samples( bbs ) );

	//make short list of most probable rotamers
	if ( !filter_low_probabilities ) {
		rotamers=rotamers_raw;
	}	else {
		utility::vector1< DunbrackRotamerSampleData > most_probable_rotamers;
		most_probable_rotamers.reserve( rotamers.size() );
		Real probability_threshold = 0.01;
		Size i = 1;
		while ( i <= rotamers.size() &&
			rotamers[ i ].probability() > probability_threshold ) {
			rotamers.push_back( rotamers[ i ] );
			i++;
		}
	}
}

void
PerturbRotamerSidechainMover::compute_rotdensities(
	RotamerList const& rotamers,
  ChiVector const& old_chi,
	ChiVector const& new_chi,
	Real& within_rot_density
) const {

	Real max_new_rot_prob(0);
	Real max_old_rot_prob(0);
	Real norm_prob(0);
	//Real const inv_nrot( 1.0 / rotamers.size() );

	for (Size jj=1; jj <= rotamers.size(); ++jj) {
		norm_prob+=rotamers[jj].probability();
	}
	//if ( std::abs( norm_prob - 1.0 ) > 0.00001 ) {
	//TR.Warning << "ALARM: probs are not normalized correctly: " << norm_prob << std::endl;
		//}

	for ( Size ii = 1; ii <= rotamers.size(); ++ii ) {
		//for each rotamer evaluate the density at our new chi angles
		Real well_prob( rotamers[ii].probability() );
		Real const new_rot_prob( rotamers[ii].chi_probability( new_chi, temperature() )*well_prob );

		Real const old_rot_prob( rotamers[ii].chi_probability( old_chi, temperature() )*well_prob );
		if (old_rot_prob > max_old_rot_prob) {
			max_old_rot_prob = old_rot_prob;
			max_new_rot_prob = new_rot_prob;
		}
	}
	within_rot_density = max_new_rot_prob;
}


} // sidechain_moves
} // simple_moves
} // protocols
