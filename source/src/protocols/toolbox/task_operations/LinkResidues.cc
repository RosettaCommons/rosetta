// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/LinkResidues.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/LinkResidues.hh>
#include <protocols/toolbox/task_operations/LinkResiduesCreator.hh>

#include <core/id/SequenceMapping.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pose/selection.hh>

// Utility Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <set>


using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.LinkResidues" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace std;

core::pack::task::operation::TaskOperationOP
LinkResiduesCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new LinkResidues );
}

void LinkResiduesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LinkResidues::provide_xml_schema( xsd );
}

std::string LinkResiduesCreator::keyname() const
{
	return LinkResidues::keyname();
}

LinkResidues::LinkResidues() {
}

LinkResidues::~LinkResidues() {}

core::pack::task::operation::TaskOperationOP LinkResidues::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new LinkResidues( *this ) );
}

void LinkResidues::remap_allowed_residues_to_template(core::pose::Pose const & pose, core::pack::rotamer_set::RotamerLinksOP & links) const{
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace core::id;
	SequenceMappingOP seqmap( new core::id::SequenceMapping() );
	for ( Size ii=1; ii<=pose.total_residue(); ++ii ) {
		links->set_template(ii,ii);
	}
	if ( template_group_!= "" ) {
		vector1< Size > template_set = core::pose::get_resnum_list_ordered( template_group_, pose );
		//get positions
		for ( Size ii=1; ii<=allGroups_.size(); ++ii ) {
			vector1< string > const grp_s( utility::string_split( allGroups_[ii] , ',' ) );
			vector1< vector1< Size > > grp_i;
			//get residues in number format into grp_res
			for ( Size kk=1; kk<=grp_s.size(); ++kk ) {
				vector1< Size > set_i = core::pose::get_resnum_list_ordered( grp_s[kk], pose );
				grp_i.push_back(set_i);
			}
			//check that all sets in the group have the same number of residues
			Size numResInSet = grp_i[1].size();
			for ( Size kk=2; kk<=grp_i.size(); ++kk ) {
				runtime_assert(grp_i[kk].size() == numResInSet);
			}
			//go through the sets and set them to be equal
			for ( Size kk=1; kk<=grp_i.size(); ++kk ) {
				Size numResInSet = grp_i[kk].size();
				if ( grp_i[kk].size() != template_set.size() ) {
					utility_exit_with_message("template and groups must all be the same size. Can not continue");
				}
				for ( Size mm=1; mm<=numResInSet; ++mm ) {
					Size source = template_set[mm];
					Size target = grp_i[kk][mm];
					links->set_template(target,source);
				}
			}
		}
	}
}


void LinkResidues::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace std;
	using namespace utility;
	using namespace core;
	core::pack::rotamer_set::RotamerLinksOP links( new core::pack::rotamer_set::RotamerLinks );
	Size nres = pose.size();
	links->resize(nres);
	remap_allowed_residues_to_template(pose,links);
	utility::vector1< set < Size > > equiv_pos;
	//all positions are equivalent to themselves
	for ( Size ii = 1; ii<= nres ; ++ii ) {
		set<Size> tmp_set;
		tmp_set.insert(ii);
		equiv_pos.push_back(tmp_set);
	}
	//add the groups in. Each group contains multiple sets of residues.
	for ( Size ii=1; ii<=allGroups_.size(); ++ii ) {
		vector1< string > const grp_s( utility::string_split( allGroups_[ii] , ',' ) );
		vector1< vector1< Size > > grp_i;
		//get residues in number format into grp_res
		for ( Size kk=1; kk<=grp_s.size(); ++kk ) {
			vector1< Size > set_i = core::pose::get_resnum_list_ordered( grp_s[kk], pose );
			grp_i.push_back(set_i);
		}
		//check that all sets in the group have the same number of residues
		Size numResInSet = grp_i[1].size();
		for ( Size kk=2; kk<=grp_i.size(); ++kk ) {
			runtime_assert(grp_i[kk].size() == numResInSet);
		}
		//go through the sets and set them to be equal
		for ( Size kk=1; kk<=grp_i.size(); ++kk ) {
			for ( Size ll=1; ll<=grp_i.size(); ++ll ) {
				Size numResInSet = grp_i[kk].size();
				for ( Size mm=1; mm<=numResInSet; ++mm ) {
					Size source = grp_i[kk][mm];
					Size target = grp_i[ll][mm];
					equiv_pos[source].insert(target);
				}
			}
		}
	}
	//print out equivalent residues
	for ( Size ii=1; ii<=equiv_pos.size(); ++ii ) {
		TR.Debug << ii <<" ";
		for ( set<Size>::iterator equiv_pos_itr=equiv_pos[ii].begin(); equiv_pos_itr!=equiv_pos[ii].end(); ++equiv_pos_itr ) {
			TR.Debug << *equiv_pos_itr <<",";
		}
		TR.Debug << std::endl;
	}
	//set up links
	for ( Size ii=1; ii<=equiv_pos.size(); ++ii ) {
		if ( equiv_pos[ii].size()!= 0 ) {
			utility::vector1< Size> equiv_pos_vector;
			equiv_pos_vector.assign(equiv_pos[ii].begin(),equiv_pos[ii].end());
			links->set_equiv(ii,equiv_pos_vector);
		}
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
	for ( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "LinkGroup" || (*tag_it)->getName() == "linkgroup" ) {
			std::string grp_i = (*tag_it)->getOption<std::string>( "group" );
			TR.Debug << "Adding grp " << grp_i << std::endl;
			allGroups_.push_back( grp_i );
		}
	}
	if ( tag->hasOption("template_group") ) {
		template_group_ = tag->getOption<std::string>("template_group");
	} else {
		template_group_ = "";
	}

}

std::string linkres_subelement_ctname( std::string const & name ) { return "linkres_" + name + "Type"; }
std::string linkres_subelement_group() { return "linkres_link_group"; }

void LinkResidues::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & linkres_subelement_ctname );
	AttributeList subattributes; subattributes.push_back( XMLSchemaAttribute::required_attribute( "group", xs_string ));
	subelements.add_simple_subelement( "LinkGroup", subattributes );
	subelements.add_simple_subelement( "linkgroup", subattributes );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen
		.element_name( keyname() )
		.complex_type_naming_func( & complex_type_name_for_task_op )
		.add_attribute( optional_name_attribute() )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
