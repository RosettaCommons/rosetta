// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/PerturbChiSidechainMover.hh
/// @brief definition of PerturbChiSidechainMover class and functions
/// @author Oliver Lange (oliver.lange@tum.de) adapted from Colin A. Smith's code


#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_PerturbChiSidechainMover_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_PerturbChiSidechainMover_hh

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/PerturbChiSidechainMover.fwd.hh>
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
class PerturbChiSidechainMover : public protocols::simple_moves::sidechain_moves::SidechainMoverBase {

public:
	typedef SidechainMoverBase Parent;
	/// @brief default constructor
	PerturbChiSidechainMover();

	/// @brief constructor with user provided rotamer library
	PerturbChiSidechainMover(
		core::pack::dunbrack::RotamerLibrary const & rotamer_library
	);

	PerturbChiSidechainMover(
		PerturbChiSidechainMover const & mover
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

	virtual
	void
	make_chi_move(
		core::conformation::Residue const& residue,
		ChiVector const& old_chi,
		ChiVector&  new_chi
	);

	virtual std::string get_name() const;

	virtual core::Real
	compute_proposal_density(
		core::conformation::Residue const &,
		core::Size const,
		core::chemical::ResidueType const &,
		ChiVector const &
	) const { return 1.0; };

	void set_magnitude( core::Real setting ) {
		magnitude_ = setting;
	}
	void set_gaussian( bool setting ) {
		gaussian_ = setting;
	}

protected:
	void set_defaults();
private:
	core::Real magnitude_;
	bool gaussian_;
}; //PerturbChiSidechainMover


} // sidechain_moves
} // simple_moves
} // protocols

#endif
