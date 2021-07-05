// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/relax/LocalRelax.hh
/// @brief A relax protocol that iteratively cart relaxes clustered subsets of residues
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_relax_LocalRelax_hh
#define INCLUDED_protocols_relax_LocalRelax_hh

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>




namespace protocols {
namespace relax {


class LocalRelax : public protocols::moves::Mover {
private:
	core::Size NCYC_, NEXP_;
	core::Real K_, max_iter_;
	bool ramp_cart_;
	bool verbose_;

	utility::vector1< core::Real > ramp_schedule_;

	core::scoring::ScoreFunctionOP pack_sfxn_, min_sfxn_;

public:
	LocalRelax();


	/// @brief one cycle of local optimization
	void
	optimization_loop(
		Pose & pose,
		core::pack::task::PackerTaskOP ptask,
		core::kinematics::MoveMapOP mm,
		core::Real, core::Real);

	/// @brief get matrix of interacting residues
	utility::vector1< utility::vector1<bool> >
	get_neighbor_graph(Pose const & pose);

	/// @brief RS integration
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief return a fresh instance of this class in an owning pointer
	protocols::moves::MoverOP clone() const override {
		return utility::pointer::make_shared< LocalRelax >(*this);
	}

	void apply( core::pose::Pose & pose) override;

	void set_pack_sfxn(core::scoring::ScoreFunctionOP const & sf) { pack_sfxn_ = sf; }
	void set_min_sfxn(core::scoring::ScoreFunctionOP const & sf) { min_sfxn_ = sf; }
	void set_sfxn(core::scoring::ScoreFunctionOP const & sf) { min_sfxn_ = sf; pack_sfxn_ = sf; }
	void set_K(core::Real const K) { K_ = K; }
	void set_max_iter(core::Real const max_iter) { max_iter_ = max_iter; }
	void set_ncyc(core::Size const ncyc) { NCYC_ = ncyc; }
	void set_nexp(core::Size const nexp) { NEXP_ = nexp; }

	core::scoring::ScoreFunctionOP get_pack_sfxn() const { return pack_sfxn_; }
	core::scoring::ScoreFunctionOP get_min_sfxn() const { return min_sfxn_; }
	core::Real get_K() const { return K_; }
	core::Real get_max_iter() const { return max_iter_; }
	core::Size get_ncyc() const { return NCYC_; }
	core::Size get_nexp() const { return NEXP_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};


}
} // protocols

#endif
