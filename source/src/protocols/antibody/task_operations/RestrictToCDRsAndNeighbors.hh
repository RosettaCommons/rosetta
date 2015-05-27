// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file /RestrictToCDRsAndNeighbors.hh
/// @brief TaskOperation to restrict to Specific CDRs and neighbors.
/// @author Jared Adolf-Bryfogle


#ifndef INCLUDED_protocols_antibody_task_operations_RestrictToCDRsAndNeighbors_hh
#define INCLUDED_protocols_antibody_task_operations_RestrictToCDRsAndNeighbors_hh

#include <protocols/antibody/task_operations/RestrictToCDRsAndNeighbors.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace antibody {
namespace task_operations {


///@brief Task Operation to restrict packing/design to specific CDRs and neighbors. 
/// See DisableAntibodyRegionOperation and DisableCDRsOperation further restrict TaskFactory.
///
///@detals See options for control of the design of CDRs, antigen, and framework.
/// By default, restricts to all packing of all CDRs and neighbors.
///  See options for control of which CDRs, including whether to only restrict design and control of whether we design 
///  neighbor antigen and/or framework residues.
///
class RestrictToCDRsAndNeighbors : public core::pack::task::operation::TaskOperation {
public:

	///Default Constructor. 
	RestrictToCDRsAndNeighbors();
	
	///@brief Constructor with AntibodyInfo
	RestrictToCDRsAndNeighbors(AntibodyInfoCOP ab_info);
	
	///@brief Constructor specifying CDRs to restrict to.
	RestrictToCDRsAndNeighbors(AntibodyInfoCOP ab_info, utility::vector1< bool > cdrs);
	
	///@brief Constructor with more options
	RestrictToCDRsAndNeighbors(
			AntibodyInfoCOP ab_info,
			utility::vector1<bool> cdrs,
			bool allow_cdr_design);
	
	///@brief Constructor with most options
	RestrictToCDRsAndNeighbors(
			AntibodyInfoCOP ab_info,
			utility::vector1<bool> cdrs,
			bool allow_cdr_design,
			bool allow_neighbor_framework_design,
			bool allow_neighbor_antigen_design);
	
	RestrictToCDRsAndNeighbors(RestrictToCDRsAndNeighbors const & src);

	virtual ~RestrictToCDRsAndNeighbors();
	core::pack::task::operation::TaskOperationOP
	clone() const;
	
	/// @brief Configure from a RosettaScripts XML tag.
	virtual void
	parse_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & );
	
	
	//////////////////////////
	
	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;
		
	void
	set_cdrs(utility::vector1<bool> const & cdrs);
	
	void
	set_cdr_only(CDRNameEnum cdr);
	
	///@brief Set the distance for the detection of neighbor residues. 6A default
	void
	set_neighbor_distance(core::Real neighbor_dis);
	
	
	///@brief Allow design of CDRs?
	void
	set_allow_design_cdr(bool allow_cdr_design);
	
	///@brief Allow design of neighbor framework residues?
	void
	set_allow_design_neighbor_framework(bool allow_framework_design);
	
	///@brief Allow design of neighbor antigen residues?
	void
	set_allow_design_neighbor_antigen(bool allow_antigen_design);
	
	///@brief Set the size of the stem - the number of residues going into the framework.  This will be included as part of the loop when determining neighbors.
	/// However, the residues of the stem itself count as part of the framework
	void
	set_stem_size(core::Size stem_size);
	

	void
	set_defaults();

private:
	
	AntibodyInfoCOP ab_info_;
	utility::vector1<bool> cdrs_;
	core::Real neighbor_dis_;
	
	bool design_cdrs_;
	bool design_antigen_;
	bool design_framework_;
	
	core::Size stem_size_;
	
	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;
	CDRDefinitionEnum cdr_definition_;
};

			
} //task_operations
} //antibody
} //protocols

#endif	//INCLUDED_ _RestrictToCDRsAndNeighbors_hh



