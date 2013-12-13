// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/CrystalContactsOperation.cc 
/// @brief  Exclude crystal contacts from design 
/// @author Patrick Conway (ptconway@uw.edu) 

// Unit Headers
#include <protocols/toolbox/task_operations/CrystalContactsOperation.hh>
#include <protocols/toolbox/task_operations/CrystalContactsOperationCreator.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/sasa.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ Headers

static basic::Tracer TR("protocols.toolbox.task_operations.CrystalContactsOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

core::pack::task::operation::TaskOperationOP
CrystalContactsOperationCreator::create_task_operation() const
{
	return new CrystalContactsOperation;
}


	CrystalContactsOperation::CrystalContactsOperation( core::Real all_gap, core::Real polar_gap, core::Real max_buried_sasa, bool invert ):
	all_gap_(all_gap),					// add this to all calculated distances
	polar_gap_(polar_gap),				// if either residue is polar - add this to calculated distances
	max_buried_sasa_(max_buried_sasa),  // ignore buried residues as defined by maximum allowed sasa
	invert_(invert)						// design residues in contact
{}

CrystalContactsOperation::~CrystalContactsOperation() {}

core::pack::task::operation::TaskOperationOP CrystalContactsOperation::clone() const
{
	return new CrystalContactsOperation( *this );
}

void
CrystalContactsOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core;
	using namespace basic;
	using namespace pose;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace utility;
	
	if( !is_symmetric(pose) )
		utility_exit_with_message( "Cannot evaluate crystal contacts on an asymmetric pose" );
	
	// calc sc sasa (mc can be exposed as long as sc is buried)
	Pose asymm_pose;
	extract_asymmetric_unit( pose, asymm_pose, true );
	utility::vector1< core::Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa_sc( asymm_pose, rsd_sasa, true /*normalize*/);
	
	// find crystal contacts - iterate through all residue pairs between subunit 1 and the other subunits
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	std::set<Size> contacts;
	for(Size ir=1; ir<=sym_info->num_independent_residues(); ir++) {
		TR.Debug << "sasa: " << ir << " " << rsd_sasa[ir] << std::endl;
		if (rsd_sasa[ir] < max_buried_sasa_) continue;						// if residue has no sc SASA, ignore
		if(!pose.residue(ir).is_protein()) continue;
		std::string atom_i = (pose.residue(ir).name3() == "GLY") ? "CA" : "CB";
		
		for(Size jr=sym_info->num_independent_residues()+1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
			if(!pose.residue(jr).is_protein()) continue;
			std::string atom_j = (pose.residue(jr).name3() == "GLY") ? "CA" : "CB";
			// decide if residue pair are contacts
			Real contact_distance = pose.residue(ir).nbr_radius() + pose.residue(jr).nbr_radius() + all_gap_ + polar_gap_;
			if(pose.residue(ir).xyz(atom_i).distance_squared(pose.residue(jr).xyz(atom_j)) <= contact_distance*contact_distance) {
				contacts.insert(ir);
				TR.Debug << "contact: " << ir << " " << jr << " " << contact_distance << std::endl;   
				break;
			}
		}
	}
	
	std::string design_select = "";
	for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		// if find residue in contacts or greater than asymm.nres, exclude
		// use != invert_ to flip Boolean, exclude residues not in contact or greater than asymm nres
		if (( ((contacts.find(ir) != contacts.end()) != invert_ )) || ir > sym_info->num_independent_residues() ) {	
			task.nonconst_residue_task(ir).prevent_repacking();
		}
		else {
			design_select += ObjexxFCL::string_of(ir)+"+";
		}
	}
	TR << "sele resi " << design_select << std::endl;
}

void
CrystalContactsOperation::parse_tag( TagCOP tag, DataMap & )
{
	all_gap_ = tag->getOption<core::Real>("all_gap", 2);
	polar_gap_ = tag->getOption<core::Real>("polar_gap", 1);
	max_buried_sasa_ = tag->getOption<core::Real>("max_buried_sasa", 0.01);
	invert_ = tag->getOption< bool >("invert",0);
}

void
CrystalContactsOperation::parse_def( utility::lua::LuaObject const & def)
{
	all_gap_ = def["all_gap"] ? def["all_gap"].to<core::Real>() : 2;
	polar_gap_ = def["polar_gap"] ? def["polar_gap"].to<core::Real>() : 1;
	max_buried_sasa_ = def["max_buried_sasa"] ? def["max_buried_sasa"].to<core::Real>() : 0.01;
}

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
