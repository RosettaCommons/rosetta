// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/DrugDesignMover.hh
/// @brief do MonteCarlo drug design in a protein context
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Yidan Tang (yidan.tang@vanderbilt.edu)

#ifndef INCLUDED_protocols_drug_design_DrugDesignMover_hh
#define INCLUDED_protocols_drug_design_DrugDesignMover_hh

// Unit header
#include <protocols/drug_design/DrugDesignMover.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/chemistries/Chemistry.fwd.hh>
#include <core/chemical/PoseResidueTypeSet.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>

// Project Headers
#include <protocols/filters/Filter.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// External headers

// C/C++ headers
#include <string>

namespace protocols {
namespace drug_design {

class DrugDesignMover : public protocols::moves::Mover {
	// TODO: Checkpointing
public:
	/// @brief default constructor
	DrugDesignMover();

	/// @brief destructor
	~DrugDesignMover() override;

	/// @brief create copy constructor
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of objectt
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief apply DrugDesignMover
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

public: // accessors

	/// @brief Return the redocker being used
	protocols::moves::MoverCOP redocker() const;

	/// @brief Return the scoring filter
	protocols::filters::FilterCOP scorer() const;

	/// @brief Return the scoring temperature
	core::Real temperature() const {
		return temperature_;
	}

	/// @brief Returns maximum number of trials
	core::Size maxtrials() const { return maxtrials_; }

	/// @brief What chain is the ligand supposed to be?
	char chain() const { return chain_; }

	/// @brief What's the prefix of the filename to which intermediate debug structures are dumped?
	std::string debug_prefix() const { return debug_prefix_; }

public: // setters

	/// @brief Add a chemistry to use
	void add_chemistry( protocols::chemistries::ChemistryOP chemistry, core::Real weight = 1.0 );

	/// @brief Add a chemistry to use
	void add_before_chemistry( protocols::chemistries::ChemistryOP chemistry );

	/// @brief Add a chemistry to use
	void add_after_chemistry( protocols::chemistries::ChemistryOP chemistry );

	/// @brief set redocking mover
	void redocker( protocols::moves::MoverOP const & mover );

	/// @brief Set the scoring filter in use
	void scorer( protocols::filters::FilterOP const & scorer );

	/// @brief Set the filter to use before redocking
	void prefilter( protocols::filters::FilterOP const & setting );

	/// @brief Set the filter to use after redocking
	void postfilter( protocols::filters::FilterOP const & setting );

	/// @brief set temperature
	void temperature( core::Real const temp ) {
		temperature_ = temp;
	}

	/// @brief set max trials of MC trials
	void maxtrials( core::Size const ntrial ) {
		maxtrials_ = ntrial;
	}

	/// @brief set the ligand chain to use for design
	void chain( char chain ) {
		chain_ = chain;
	}

	/// @brief What's the prefix of the filename to which intermediate debug structures are dumped?
	void debug_prefix( std::string const & setting ) {
		debug_prefix_ = setting;
	}

public: // Other functions

	/// @brief Is position n a valid designable residue according to the settings
	bool
	check_design_position( core::pose::Pose const & pose, core::Size n );

	/// @brief Determine which residue on this pose should be designed
	core::Size
	find_design_position( core::pose::Pose const & pose );

	/// @brief parse xml file
	void
	parse_my_tag( TagCOP tag, basic::datacache::DataMap & data ) override;

	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

protected:

	/// @brief Return a unique string representing which job this is
	std::string
	get_jobname() const;

	/// @brief Subclass access to the chemistries
	utility::vector1< protocols::chemistries::ChemistryOP > & chemistries() { return chemistries_; }

	/// @brief Pre-process the ResidueType
	/// * Apply the always chemistries
	bool // Returns true if we can't use this residuetype
	pre_process_restype(
		core::chemical::MutableResidueTypeOP restype, // Can/will be modified!
		core::chemical::IndexVDMapping & index_vd_mapping, // Can/will be modified!
		core::pose::Pose const & pose // for context
	);

	/// @brief Post-process the residuetype
	/// * Apply the always chemistries
	/// * Fix up some after issues
	bool // Returns true if we can't use this residuetype
	post_process_restype(
		core::chemical::MutableResidueTypeOP restype, // Can/will be modified!
		core::chemical::IndexVDMapping & index_vd_mapping, // Can/will be modified!
		std::string const & new_name,
		core::pose::Pose const & pose // for context
	);

	/// @brief Take a new residuetype, place it on the pose,
	/// and see if it passes filters.
	/// @details Differs from place_residue_type() in that it also filters
	bool // Returns true if we can't use this residuetype
	emplace_residue_type(
		core::pose::Pose & pose, // Can/will be modified!
		core::chemical::MutableResidueTypeOP restype, // Can/will be modified!
		core::chemical::IndexVDMapping const & index_VD_mapping
	);

	/// @brief Find a new residue name, if the name already exists in the residue type set.
	std::string
	find_new_res_name( std::string original_name, core::Size iteration, core::Size subiteration = 0 ) const;

	/// @brief Dump the ligand to the given file (appending)
	void
	dump_molecule( core::chemical::MutableResidueType const & restype, std::string const & stage) const;

	/// @brief Dump the ligand to the given file (appending)
	void
	dump_molecule( core::conformation::Residue const & residue, std::string const & stage) const;

	protocols::chemistries::ChemistryOP
	chemistry_from_subtag( utility::tag::TagCOP const subtag, basic::datacache::DataMap & data  ) const;

//	/// @brief Return updated temperature if using simulated annealing
//	core::Real
//	update_temperature( core::Real currT, core::Real curr_accept_ratio );

private: // Data

	/// @brief What chemistries can be done in the design
	utility::vector1< protocols::chemistries::ChemistryOP > chemistries_;

	/// @brief The weight to use each chemistry at
	utility::vector1< core::Real > weights_;

	/// @brief What chemistries are always applied before the randomly chosen one.
	utility::vector1< protocols::chemistries::ChemistryOP > before_chemistries_;

	/// @brief What chemistries are always applied after the randomly chosen one.
	utility::vector1< protocols::chemistries::ChemistryOP > after_chemistries_;

	/// @brief Mover which is used for redocking
	protocols::moves::MoverOP redocker_;

	/// @brief The filter used to do scoring
	protocols::filters::FilterOP scorer_;
	
	/// @brief if set true, a ligand efficiency interface score will be used for MC acceptance
	bool lig_efficiency_;

	/// @brief acceptance criterion temperature
	core::Real temperature_;
	
	/// @brief if set true, temperature varies throughout MC trials
	bool simulated_annealing_;
	/// @brief parameters for simulated annealing, in the order: min, max, step, OFF_after_n_trials
	utility::vector0< core::Real > temperature_params;
	
//	/// @brief if use simulated annealing, the desired acceptance ratio
//	core::Real target_accept_ratio_;
	
//	/// @brief if use simulated annealing, the update interval for temperature
//	core::Real temp_update_interval_;
	
//	/// @brief if use simulated annealing, the update formula for temperature
//	std::string cooling_scheme_;
	
//	/// @brief if use simulated annealing, use current accept ratio or cumulative ratio
//	bool useCurrRatio_;

	/// @brief the number of MC trials to run in the complete run.
	core::Real maxtrials_;

	/// @brief What chain is the ligand?
	char chain_;

	/// @brief The pass/fail filter to use prior to redocking
	protocols::filters::FilterOP prefilter_;

	/// @brief The pass/fail filter to use after redocking
	protocols::filters::FilterOP postfilter_;

	/// @brief If non-empty, dump intermediate structures to files based on this prefix.
	std::string debug_prefix_;

	/// @brief TODO: HACK - needed to make sure that restypes don't go out of scope. Need to handle this better.
	utility::vector1< core::chemical::ResidueTypeCOP > restypes_;
	core::chemical::PoseResidueTypeSetOP restypeset_;
};

} // namespace drug_design
} // namespace protocols

#endif
