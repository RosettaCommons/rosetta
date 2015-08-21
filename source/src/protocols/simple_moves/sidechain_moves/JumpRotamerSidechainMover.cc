// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMoverBase.cc
/// @brief implementation of JumpRotamerSidechainMover class and functions
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#include <protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMoverCreator.hh>

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

static thread_local basic::Tracer tr( "protocols.simple_moves.sidechain_moves.JumpRotamerSidechainMover" );

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

using namespace chemical;
using namespace conformation;


std::string
JumpRotamerSidechainMoverCreator::keyname() const {
	return JumpRotamerSidechainMoverCreator::mover_name();
}

protocols::moves::MoverOP
JumpRotamerSidechainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new JumpRotamerSidechainMover );
}

std::string
JumpRotamerSidechainMoverCreator::mover_name() {
	return "JumpRotamerSidechain";
}


JumpRotamerSidechainMover::JumpRotamerSidechainMover() {
	protocols::moves::Mover::type( "JumpRotamerSidechain" );
	set_defaults();
}

JumpRotamerSidechainMover::JumpRotamerSidechainMover(
	pack::dunbrack::RotamerLibrary const & rotamer_library
) : Parent( rotamer_library ) {
	set_defaults();
}

JumpRotamerSidechainMover::JumpRotamerSidechainMover(
	JumpRotamerSidechainMover const & mover
) : Parent ( mover ) {
	set_defaults();
}

protocols::moves::MoverOP
JumpRotamerSidechainMover::clone() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::sidechain_moves::JumpRotamerSidechainMover(*this) );
}

void
JumpRotamerSidechainMover::set_defaults() {
	sample_rotwells_unif_=false;
}

void
JumpRotamerSidechainMover::parse_my_tag(
	utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	pose::Pose const & /*pose*/
) {

}

std::string
JumpRotamerSidechainMover::get_name() const {
	return "JumpRotamerSidechainMover";
}


void
JumpRotamerSidechainMover::make_chi_move(
	conformation::Residue const& residue,
	ChiVector const&,
	ChiVector& new_chi
) {

	RotamerList rotamers;
	build_rotamer_list( residue, sample_rotwells_unif_, rotamers );

	/// temper the well probabilities
	utility::vector1< Real > rot_probs;
	Real normalize;
	bool is_tempered( std::abs( temperature() - 1.0 ) > 0.01 );
	if ( is_tempered ) {
		compute_tempered_rotamer_probabilities( rotamers, temperature(), rot_probs, normalize );
	}

	/// select a random rotamer
	Size rotnum;
	Real rand = numeric::random::rg().uniform();
	//if ( rand <=0 ){
	//tr.Debug << "RG.uniform is " << rand << std::endl;}
	Real const inv_nrot( 1.0/rotamers.size() );
	Real rot_prob_normalize (0);
	runtime_assert( !is_tempered );
	for ( rotnum=1; rotnum <=rotamers.size() &&rand > 0; ++rotnum ) {
		rot_prob_normalize+=rotamers [rotnum].probability();
	}
	//tr.Debug << "rot_prob_normalized is" << rot_prob_normalize << std::endl;
	for ( rotnum=1; rotnum <= rotamers.size() && rand > 0; ++rotnum ) {
		if ( sample_rotwells_unif_ ) {
			rand -= inv_nrot; //pick any rotamer with uniform probability
		} else {
			rand -= ( is_tempered ? rot_probs[ rotnum ]*normalize : ( rotamers[ rotnum ].probability()/rot_prob_normalize ));
		}
	}
	if ( rand <= 0 ) {
		rotnum -= 1;
	}
	runtime_assert( rotnum >= 1 && rotnum <= rotamers.size() );
	//rotamer_sample_data[rotnum].assign_random_chi(last_chi_angles_,numeric::random::rg());
	rotamers[rotnum].assign_random_chi( new_chi, numeric::random::rg(), temperature() );
}

///all angles in degree
Real JumpRotamerSidechainMover::compute_proposal_density(
	Residue const & new_residue,
	Size const,
	chemical::ResidueType const &,
	utility::vector1<Real> const &
) const {
	utility::vector1<Real> const & new_chi( new_residue.chi() );
	using namespace pack::dunbrack;
	RotamerList rotamers;
	build_rotamer_list( new_residue, sample_rotwells_unif_, rotamers );

	Real rot_density;
	compute_rotdensities( rotamers, new_chi, rot_density );

	return rot_density;
}

void
JumpRotamerSidechainMover::compute_rotdensities(
	RotamerList const& rotamers,
	ChiVector const& new_chi,
	Real& rot_density
) const {
	rot_density=0;
	//Real max_new_rot_prob(0);
	//Real max_old_rot_prob(0);

	utility::vector1< Real > rot_probs;
	Real normalize;
	bool is_tempered( std::abs( temperature() - 1.0 ) > 0.01 );
	if ( is_tempered ) {
		compute_tempered_rotamer_probabilities( rotamers, temperature(), rot_probs, normalize );
	}

	Real norm_prob=0;
	for ( RotamerList::const_iterator it = rotamers.begin(); it!=rotamers.end(); ++it ) {
		norm_prob += it->probability();
	}
	if ( std::abs( norm_prob - 1.0 ) > 0.00001 ) {
		tr.Warning << "ALARM: probs are not normalized correctly: " << norm_prob << std::endl;
	}
	Real const inv_nrot( 1.0 / rotamers.size() );
	for ( Size ii = 1; ii <= rotamers.size(); ++ii ) {
		//for each rotamer evaluate the density at our new chi angles
		Real const within_well_prob( rotamers[ii].chi_probability( new_chi, temperature() ) );

		if ( sample_rotwells_unif_ ) {
			rot_density += inv_nrot * within_well_prob;
		} else {
			Real const well_prot( is_tempered ? rot_probs[ ii ]*normalize : rotamers[ii].probability()/norm_prob );
			rot_density += exp( log( well_prot ) + log( within_well_prob ) );
		}
	}
}

void
JumpRotamerSidechainMover::compute_tempered_rotamer_probabilities(
	RotamerList const& rotamers,
	core::Real,
	utility::vector1< Real >& rot_probs,
	core::Real& normalize
) const {
	rot_probs.clear();
	rot_probs.resize( rotamers.size() );
	Real const inv_T( 1.0/Parent::temperature() );
	Real sum_prob( 0.0 );
	for ( Size rotnum=1; rotnum <= rotamers.size(); ++rotnum ) {
		sum_prob+=std::pow( rotamers[ rotnum ].probability(), inv_T );
	}
	normalize = 1.0/sum_prob;
}


} // sidechain_moves
} // simple_moves
} // protocols
