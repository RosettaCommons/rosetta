// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/denovo_design/task_operations/DesignBySecondaryStructure.hh
/// @brief Design residues with secondary structures that don't match the desired secondary structure
/// specified in the blueprint.
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_devel_denovo_design_task_operations_designbysecondarystructure_hh
#define INCLUDED_devel_denovo_design_task_operations_designbysecondarystructure_hh

// unit headers
#include <devel/denovo_design/task_operations/DesignBySecondaryStructure.fwd.hh>

// protocol headers
#include <protocols/denovo_design/filters/PsiPredInterface.fwd.hh>
#include <devel/denovo_design/task_operations/HighestEnergyRegion.hh>


// project headers
#include <protocols/ss_prediction/SS_predictor.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

namespace devel {
namespace denovo_design {
namespace task_operations {

class DesignBySecondaryStructureOperation : public HighestEnergyRegionOperation {
public:
	/// @brief default constructor
	DesignBySecondaryStructureOperation();

	/// @brief value constructor
	DesignBySecondaryStructureOperation( std::string const bp_file,
		std::string const cmd,
		bool const prevent_native,
		bool const prevent_bad_point_mutations );

	/// @brief copy constructor
	DesignBySecondaryStructureOperation( DesignBySecondaryStructureOperation const & rval );

	/// @brief destructor
	virtual ~DesignBySecondaryStructureOperation();

	/// @brief make clone
	virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief apply
	virtual void apply( Pose const & pose, core::pack::task::PackerTask & task ) const;

	/// @brief Runs the calculation and caches residues to design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "SSPrediction"; }

public:
	void parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & );

	/// @brief opens the passed blueprint file and determines the desired secondary structure
	void initialize_blueprint_ss( std::string const blueprint_file );

private:
	std::string blueprint_ss_; // the blueprint secondary structure
	std::string pred_ss_; // cache of the predicted secondary structure
	utility::vector1< core::Real > psipred_prob_; // cache of the predicted probability of secondary structure

	/// @brief If set, the amino acid residue already in the pose will be disallowed default=false
	bool prevent_native_aa_;
	/// @brief If set, all mutations at all positions will be scanned in one pass, and those that cause worse psipred secondary structure agreement will be disallowed (default=false)
	bool prevent_bad_point_mutants_;
	/// @brief the object which directly communicates with psipred and parses psipred output
	protocols::denovo_design::filters::PsiPredInterfaceOP psipred_interface_;
	/// @brief the svm secondary structure predictor
	protocols::ss_prediction::SS_predictorOP ss_predictor_;
};


} // task_operations
} // denovo_design
} // devel

#endif
