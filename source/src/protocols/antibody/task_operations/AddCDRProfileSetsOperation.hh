// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/task_operations/AddCDRProfileSetsOperation.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_task_operations_AddCDRProfileSetsOperation_hh
#define INCLUDED_protocols_antibody_task_operations_AddCDRProfileSetsOperation_hh

#include <protocols/antibody/task_operations/AddCDRProfileSetsOperation.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>

#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace antibody {
namespace task_operations {
			
///@brief Add Cluster-based sets of mutations as a TaskOperation.
/// Essentially samples full sequences of CDRs within a particular CDR cluster randomly each time the packer is called.
/// Does this for each CDR.
///
/// See protocols/toolbox/task_operations/MutationSetDesignOperation for more info.
/// If a CDR has an unknown cluster or there are no data for that particular CDR, will skip that CDR.
///
/// CDR definitions used are North/Dunbrack as the clusters are defined using it.
///
///
///@details
/// Note that by default, a data cutoff of 10 is set.  If the cluster has less than 10 sequences it will be skipped.
/// Use the set_cutoff function to change this.
///
class AddCDRProfileSetsOperation : public core::pack::task::operation::TaskOperation {
public:

	AddCDRProfileSetsOperation();
	
	AddCDRProfileSetsOperation(AntibodyInfoCOP ab_info);
	
	AddCDRProfileSetsOperation(AntibodyInfoCOP ab_info, utility::vector1<bool> const & cdrs);
	
	AddCDRProfileSetsOperation(AntibodyInfoCOP ab_info, utility::vector1< bool> const & cdrs, bool limit_only_to_length);
	
	
	AddCDRProfileSetsOperation(AddCDRProfileSetsOperation const & src);

	virtual ~AddCDRProfileSetsOperation();

	core::pack::task::operation::TaskOperationOP
	clone() const;
	
	/// @brief Configure from a RosettaScripts XML tag.
	virtual void
	parse_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & );
	
	
	///////////////////////
	
	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;
	
	///@brief Pre load the data according to options instead of at each apply.  
	utility::vector1<bool>
	pre_load_data(core::pose::Pose const & pose);
	
	///@brief Set a single CDR with which to use this task op on.
	/// Default is all of them.
	void
	set_cdr_only(CDRNameEnum cdr);
	
	///@brief Set the CDRs with which to use this task op on.
	/// Default is all of them.
	void
	set_cdrs( utility::vector1< bool > const & cdrs);
	
	///@brief Set the class to sample all CDR sequences within a particular length, 
	/// instead of by cluster, which is the default.
	/// Default False.
	void
	set_limit_only_to_length( bool limit_only_to_length);
	
	
	void
	set_defaults();
	
	///@brief Force the use of the antibody database that houses only the North data.
	/// This is the db distributed with Rosetta.  If a current one is present in the database, it will use that instead.
	/// Default False.
	void
	set_force_north_paper_db(bool force_north_db);
	
	///@brief Use the light chain type if held by AntibodyInfo.  
	/// Default is true.
	void
	set_use_light_chain_type(bool use_light_chain_type);
	
	///@brief Use cluster outliers as defined using DihedralDistance and RMSD.
	/// Default false. 
	void
	set_use_outliers( bool use_outliers);
	
	///@brief Set the cutoff.  Will not add the profile set if the total is less than or equal to this number.
	/// Default is 10.
	void
	set_cutoff( core::Size cutoff);
	
	///@brief Set the number of times a sequence each chosen.  Increase this number to increase variability of design.
	/// Default 1 round
	void
	set_picking_rounds(core::Size rounds);
	
	///@brief Add to the current set of amino acids in the task or replace them? 
	/// Default False
	void
	set_add_to_current(bool add_to_current);
	
	///@brief Include the poses current residue type in the allowed amino acids.
	/// Default True.  
	void
	set_include_native_type(bool use_native);
	
private:
	AntibodyInfoCOP ab_info_;
	utility::vector1<bool> cdrs_;
	
	bool limit_only_to_length_;
	bool force_north_paper_db_;
	
	bool use_light_chain_type_;
	bool use_outliers_;
	
	core::Size cutoff_;
	
	//Profile Options
	core::Size picking_rounds_;
	bool keep_task_allowed_aas_;
	bool include_native_restype_;
	
	CDRDBSequenceSet sequences_;
	bool pre_loaded_data_;
		
	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;

};

} //task_operations
} //antibody
} //protocols

#endif	//INCLUDED_ _AddCDRProfileSetsOperation_hh



