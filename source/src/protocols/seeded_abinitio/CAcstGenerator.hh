// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/seeded_abinitio/CAcstGenerator.hh
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_CAcstGenerator_hh
#define INCLUDED_protocols_seeded_abinitio_CAcstGenerator_hh

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace seeded_abinitio {

class CAcstGenerator : public protocols::moves::Mover {
public:
	typedef core::pose::Pose Pose;

	CAcstGenerator();

	void apply( core::pose::Pose & pose ) override;

	// XRW TEMP  std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	~CAcstGenerator() override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::scoring::constraints::ConstraintSetOP ca_cst_;

	/// determines whether constraints for the areas should be which will be "replaced" by the seeds
	bool add_cst_seed_;

	/// container for the cutpoints, since there shouldnt be constraints to cutpoints
	utility::vector1< core::Size > cut_points_;

	/// container with residues from seeds that should have constraints
	utility::vector1< core::Size > seed_exceptions_;

	/// stddeviation for the harmonic CA constraints
	core::Real stddev_;

	/// container that has the seed information
	protocols::loops::Loops all_seeds_;

	/// residues for which no constraints should be derrived
	protocols::loops::Loops clear_seeds_;

	/// user specified which chain to gather the constraints from
	core::Size from_chain_;

	/// user specified to which chain of the input chain is applied to
	core::Size to_chain_;

	///user specified a template
	bool template_presence_;

	/// template pdb
	core::pose::PoseOP template_pdb_;

	/// the chain/pose that the user actually wants to read the constraints from
	core::pose::PoseOP curr_pose_;

	/// replace constraints or add onto them
	bool replace_;

	/// sequence separation after which the pair constraints are added
	core::Size seq_separation_;

	/// distance separation
	core::Real distance_cutoff_;
};

}
}

#endif
