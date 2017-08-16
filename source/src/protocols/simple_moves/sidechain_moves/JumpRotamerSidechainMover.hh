// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMover.hh
/// @brief definition of JumpRotamerSidechainMover class and functions
/// @author Oliver Lange (oliver.lange@tum.de) adapted from Colin A. Smith's code


#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_JumpRotamerSidechainMover_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_JumpRotamerSidechainMover_hh

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMover.fwd.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbRotamerSidechainMover.hh>

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
class JumpRotamerSidechainMover : public protocols::simple_moves::sidechain_moves::PerturbRotamerSidechainMover {

public:
	typedef PerturbRotamerSidechainMover Parent;
	/// @brief default constructor
	JumpRotamerSidechainMover();

	/// @brief constructor with user provided rotamer library
	JumpRotamerSidechainMover(
		core::pack::dunbrack::RotamerLibrary const & rotamer_library
	);

	JumpRotamerSidechainMover(
		JumpRotamerSidechainMover const & mover
	);

	protocols::moves::MoverOP
	clone() const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	core::Real
	compute_proposal_density(
		core::conformation::Residue const & new_residue,
		core::Size const resnum,
		core::chemical::ResidueType const & old_res_type,
		ChiVector const & old_chi
	) const override;

	void
	make_chi_move(
		core::conformation::Residue const& residue,
		ChiVector const& old_chi,
		ChiVector&  new_chi
	) override;


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	void set_defaults();

	void compute_tempered_rotamer_probabilities(
		RotamerList const&,
		core::Real temperature,
		utility::vector1< core::Real >& prob,
		core::Real& normalize
	) const;

	//helper function
	void compute_rotdensities(
		RotamerList const&,
		ChiVector const& new_chi,
		core::Real& rot_density
	) const;

private:
	bool sample_rotwells_unif_;
}; //JumpRotamerSidechainMover


} // simple_moves
} // sidechain_moves
} // protocols

#endif
