// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MotifDnaPacker.hh
/// @author sthyme
/// @brief
/// @details

#ifndef INCLUDED_protocols_motifs_MotifDnaPacker_hh
#define INCLUDED_protocols_motifs_MotifDnaPacker_hh

// Unit Headers
#include <protocols/motifs/MotifDnaPacker.fwd.hh>

// Project Headers
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/dna/DnaInterfacePacker.fwd.hh>
#include <protocols/dna/DnaDesignDef.fwd.hh>
#include <protocols/dna/PDBOutput.fwd.hh>
#include <protocols/motifs/MotifSearch.fwd.hh>

#include <list>
#include <map>
#include <set>

#include <utility/vector1.hh>

namespace protocols {
namespace motifs {

class MotifDnaPacker : public protocols::moves::Mover {

public:
	MotifDnaPacker();

	MotifDnaPacker(
		core::scoring::ScoreFunctionOP,
		bool minimize = false,
		std::string filename_root = "dnapacker"
	);

	~MotifDnaPacker() override; // important for properly "releasing" owning pointer data

	protocols::moves::MoverOP fresh_instance() const override;
	protocols::moves::MoverOP clone() const override;

	void apply( core::pose::Pose & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	//void minimize_dna( core::pose::Pose & );

private:
	void init_standard( core::pose::Pose & );
	void init_options();

	void
	run_motifs(
		core::pose::Pose & pose,
		utility::vector1< core::Size > & design_positions,
		std::set< core::Size > & src_pos,
		std::map< core::Size, core::pack::rotamer_set::Rotamers > & rotamer_map,
		std::map< core::Size, std::set< std::string > > & types_map,
		std::list< std::string > & info_lines,
		core::pack::task::TaskFactoryOP taskfactory
	);

	// aromatic_motifs and expand_motifs could likely be combined at some point
	// if an input type list was added that defaults to all 20 amino acids
	void
	expand_motifs(
		core::pose::Pose & pose,
		utility::vector1< core::Size > & design_positions,
		std::set< core::Size > & src_pos,
		std::map< core::Size, core::pack::rotamer_set::Rotamers > & rotamer_map,
		std::map< core::Size, std::set< std::string > > & types_map,
		std::list< std::string > & info_lines,
		core::pack::task::TaskFactoryOP taskfactory
	);

	void
	aromatic_motifs(
		core::pose::Pose & pose,
		utility::vector1< core::Size > & design_positions,
		std::set< core::Size > & src_pos,
		std::map< core::Size, core::pack::rotamer_set::Rotamers > & rotamer_map,
		std::map< core::Size, std::set< std::string > > & types_map,
		std::list< std::string > & info_lines,
		core::pack::task::TaskFactoryOP taskfactory
	);

	void
	motif_expansion_inner_loop(
		core::pose::Pose & pose,
		std::set< core::Size > current_pos,
		std::map< core::Size, core::pack::rotamer_set::Rotamers >::const_iterator it,
		std::set< std::string >::const_iterator it2,
		std::map< core::Size, core::pack::rotamer_set::Rotamers > & rotamer_map,
		std::list< std::string > & info_lines,
		core::pack::task::TaskFactoryOP taskfactory
	);

private:
	protocols::dna::DnaInterfacePackerOP dna_packer_;
	protocols::motifs::MotifSearchOP motif_search_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::dna::PDBOutputOP pdboutput_;
	core::pose::PoseCOP starting_pose_;
	protocols::dna::DnaDesignDefOPs targeted_dna_;

	// user defined options
	bool run_motifs_;
	bool expand_motifs_;
	bool aromatic_motifs_;
	bool minimize_dna_;
	core::Real special_rotweight_;
	core::Real num_repacks_;
	bool flex_dna_sugar_;
	std::string filename_root_;
	bool dna_design_;

};

} // namespace motifs
} // namespace protocols

#endif
