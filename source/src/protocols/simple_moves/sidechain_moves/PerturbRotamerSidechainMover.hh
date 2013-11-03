// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/PerturbRotamerSidechainMover.hh
/// @brief definition of PerturbRotamerSidechainMover class and functions
/// @author Oliver Lange (oliver.lange@tum.de) adapted from Colin A. Smith's code


#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_PerturbRotamerSidechainMover_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_PerturbRotamerSidechainMover_hh

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/PerturbRotamerSidechainMover.fwd.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMoverBase.hh>

// Protocols Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>

// Core Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/id/DOF_ID_Range.fwd.hh>

#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

/// @brief class for non-discrete side chain sampling using Dunbrack rotamer probabilities/distributions
class PerturbRotamerSidechainMover : public protocols::simple_moves::sidechain_moves::SidechainMoverBase {

public:
	typedef SidechainMoverBase Parent;
	/// @brief default constructor
	PerturbRotamerSidechainMover();

	/// @brief constructor with user provided rotamer library
	PerturbRotamerSidechainMover(
		core::pack::dunbrack::RotamerLibrary const & rotamer_library
	);

	PerturbRotamerSidechainMover(
		PerturbRotamerSidechainMover const & mover
	);

	virtual
	protocols::moves::MoverOP
	clone() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual	void
	make_chi_move(
	  core::conformation::Residue const& residue,
		ChiVector const& old_chi,
		ChiVector&  new_chi
	);

	virtual std::string get_name() const;

	virtual core::Real
	compute_proposal_density(
		core::conformation::Residue const & new_residue,
		core::Size const resnum,
		core::chemical::ResidueType const & old_res_type,
		ChiVector const & old_chi
	) const;

	void set_temperature( core::Real setting ) {
		temperature_ = setting;
	}

	core::Real temperature() const {
		return temperature_;//temperature_;
	}

protected:
	void set_defaults();
public:
	typedef utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > RotamerList;
	void build_rotamer_list(
    core::conformation::Residue const&,
		bool filter, /* filter out low probabilities p<0.01 */
		RotamerList&
	) const;

	//helper function
	void compute_rotdensities(
		RotamerList const&,
		ChiVector const& old_chi,
		ChiVector const& new_chi,
		core::Real& within_rot_density
	) const;


private:
	//the dunbrack rotamers std-deviation can be scaled ... this is similar to a temperature;
	// 1.0 means the standard dunbrack probabilities are used --- does that mean cryo-temperature ?
	core::Real temperature_;
}; //PerturbRotamerSidechainMover


} // sidechain_moves
} // simple_moves
} // protocols

#endif
