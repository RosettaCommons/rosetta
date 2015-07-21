// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/LinkResidues.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/LinkResidues.hh>
#include <protocols/toolbox/task_operations/LinkResiduesCreator.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pose/selection.hh>

// Utility Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>


using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.LinkResidues" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace std;

LinkResidues::LinkResidues() {
}

LinkResidues::~LinkResidues() {}

core::pack::task::operation::TaskOperationOP
LinkResiduesCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new LinkResidues );
}

core::pack::task::operation::TaskOperationOP LinkResidues::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new LinkResidues( *this ) );
}

void LinkResidues::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
		using namespace std;
		using namespace utility;
		using namespace core;
		core::pack::rotamer_set::RotamerLinksOP links( new core::pack::rotamer_set::RotamerLinks );
		Size nres = pose.total_residue();
		links->resize(nres);
		utility::vector1< utility::vector1< Size > > equiv_pos;
		//all positions are equivalent to themselves
		for (Size ii = 1; ii<= nres ; ++ii){
				utility::vector1<Size> list;
				list.push_back(ii);
				equiv_pos.push_back(list);
		}
		//add the groups in. Each group contains multiple sets of residues.
		for (Size ii=1; ii<=allGroups_.size(); ++ii) {
				vector1< string> const grp_s( utility::string_split( allGroups_[ii] , ',' ) );
				vector1< vector1 <Size> > grp_i;
				//get residues in number format into grp_res
				for(Size kk=1; kk<=grp_s.size(); ++kk){
						vector1<Size> set_i = core::pose::get_resnum_list_ordered( grp_s[kk], pose );
						grp_i.push_back(set_i);
				}
				//check that all sets in the group have the same number of residues
				Size numResInSet = grp_i[1].size();
				for(Size kk=2; kk<=grp_i.size(); ++kk){
						runtime_assert(grp_i[kk].size() == numResInSet);
				}
				//go through the sets and set them to be equal
				for(Size kk=1; kk<=grp_i.size(); ++kk){
						for(Size ll=1; ll<=grp_i.size(); ++ll){
								Size numResInSet = grp_i[kk].size();
								for(Size mm=1; mm<=numResInSet; ++mm){
										Size source = grp_i[kk][mm];
										Size target = grp_i[ll][mm];
										equiv_pos[source].push_back(target);
								}
						}
				}
		}
		//print out equivalent residues
		for(Size ii=1; ii<=equiv_pos.size(); ++ii){
				TR.Debug << ii <<" ";
				for(Size kk=1; kk<=equiv_pos[ii].size(); ++kk){
						TR.Debug << equiv_pos[ii][kk] <<",";
				}
				TR.Debug << std::endl;
		}
		//set up links
		for(Size ii=1; ii<=equiv_pos.size(); ++ii){
				if(equiv_pos[ii].size()!= 0)
						links->set_equiv(ii,equiv_pos[ii]);
		}
		task.rotamer_links( links );
}

void
LinkResidues::add_group( std::string group  ){
	allGroups_.push_back(group);
}

void
LinkResidues::parse_tag( TagCOP tag , DataMap & )
{
  // now parse ncs groups <<< subtags
  utility::vector1< TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< TagCOP >::const_iterator tag_it;
	for( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ){
		if( (*tag_it)->getName() == "LinkGroup" || (*tag_it)->getName() == "linkgroup" ){
				std::string grp_i = (*tag_it)->getOption<std::string>( "group" );
				TR.Debug << "Adding grp " << grp_i << std::endl;
				allGroups_.push_back( grp_i );
		}
	}
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
