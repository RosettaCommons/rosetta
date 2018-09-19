// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/ConservativeDesignOperation.hh
/// @brief TaskOperation to allow only conservative mutations at designable residues
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_task_operations_ConservativeDesignOperation_hh
#define INCLUDED_protocols_task_operations_ConservativeDesignOperation_hh

#include <protocols/task_operations/ConservativeDesignOperationCreator.hh>
#include <protocols/task_operations/ConservativeDesignOperation.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace task_operations {

/// @brief A TaskOperation that sets the allowed amino acids of designable residues to the native amino acid's conservative mutations.
///
/// @details Default is to act on all designable residues.  Use limit_to_positions to limit this.
/// Default is to replace the allowed_aas with these conservative mutations.
/// Data is loaded from database/sequence/resinfo.db.
///
class ConservativeDesignOperation : public core::pack::task::operation::TaskOperation {

public:

	/// @brief Default constructor.  Will use native aa from apply (Changes each pack if tf is passed).
	ConservativeDesignOperation();

	/// @brief Constructor with setting of data source.
	ConservativeDesignOperation(std::string data_source);

	~ConservativeDesignOperation() override;

	ConservativeDesignOperation(ConservativeDesignOperation const & src);

	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;

	/// @brief Used to parse an xml-like tag to load parameters and properties.
	void parse_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "ConservativeDesignOperation"; }

public:

	////////////////////////////////////////////////////////////////////////////
	// Class Options
	//
	//

	/// @brief Limit to a subset of residue positions, already set to designable.
	void
	limit_to_positions( utility::vector1< Size > const & positions );

	void
	include_residue( core::Size resnum );

	/// @brief Clear any set positions.
	void
	clear_positions();

	/// @brief Add to the allowed amino acids list instead of replacing it.  Default false.
	void
	add_to_allowed_aas( bool setting ){
		add_to_allowed_aas_ = setting;
	}

	/// @brief Include native amino acid in the allowed_aas list.  Default true.
	void
	include_native_aa( bool setting ){
		include_native_aa_ = setting;
	}

	/// @brief Will use native residues from the this pose to determine conserved aa.
	void
	use_pose_sequence_as_native( core::pose::Pose const & pose );

	void
	set_native_sequence( std::string const & seq ){
		pose_sequence_ = seq;
	}

	/// @brief Set the source of the data used to define what is conservative.
	/// Options are: chothia_76 and the Blosum matrices from 30 to 100; designated as blosum30, 62, etc.
	/// Default is blosum62.  The higher the number, the more conservative the set of mutations (numbers are sequence identity cutoffs)
	void
	set_data_source( std::string const & data_source );

public:
	ConservativeDesignOperation & operator =( ConservativeDesignOperation const & rhs);

	core::pack::task::operation::TaskOperationOP clone() const override;

protected:
	bool
	skip_resid( core::Size resid, core::pack::task::PackerTask const & task, std::string const & seq ) const;

private:

	void load_data_from_db();

	//void load_data_from_file();

	void init_for_equal_operator_and_copy_constructor( ConservativeDesignOperation & lhs, ConservativeDesignOperation const & rhs);

	void set_defaults();

	//void read_cmd_line_options();


private:

	utility::vector1< utility::vector1<bool> > conserved_mutations_; //AA index in alphabetical order -> conservative mutations.

	//An AND operation will be performed between these two, unless positions_ is empty.
	utility::vector1< core::Size > positions_;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

	bool add_to_allowed_aas_;
	bool include_native_aa_;

	std::string pose_sequence_;
	std::string data_source_;

};

} //task_operations
} //protocols


#endif //INCLUDED_protocols_antibody_design_ConservativeDesignOperation_hh

