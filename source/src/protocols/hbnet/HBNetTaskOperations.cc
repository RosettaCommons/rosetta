// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hbnet/HBNetTaskOperations.cc
/// @brief task operations for HBNet
/// @author Scott Boyken (sboyken@gmail.com)

#include <protocols/hbnet/HBNetTaskOperations.hh>
#include <protocols/hbnet/ConstrainHBondNetworkCreator.hh>

#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>

#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/pointer/access_ptr.hh>

// utility includes
#include <utility/string_util.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace hbnet {

static basic::Tracer TR("protocols.hbnet.HBNetTaskOperations");

/// @brief sets allowed residue types to constrain HBNet residues in downstream design;
///         add this Taskop to any design movers downstream of HBNet
ConstrainHBondNetwork::ConstrainHBondNetwork(){}

//destructor
ConstrainHBondNetwork::~ConstrainHBondNetwork(){}

core::pack::task::operation::TaskOperationOP
ConstrainHBondNetwork::clone() const
{
	ConstrainHBondNetworkOP new_clone( new ConstrainHBondNetwork( *this ) ); // why this?
	return new_clone;
}

/// @brief sets allowed residue types to constrain HBNet residues in downstream design;
///         add this Taskop to any design movers downstream of HBNet
void
ConstrainHBondNetwork::apply(
	Pose const & pose,
	PackerTask & task ) const
{
	core::scoring::constraints::ConstraintSetCOP cst_op(pose.constraint_set()); //needs to be COP

	protocols::toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	protocols::toolbox::match_enzdes_util::EnzdesCstCacheCOP cst_cache(0);
	protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP cstio;
	if ( enz_obs ) {
		cst_cache = enz_obs->cst_cache();
	}
	if ( cst_cache ) {
		cstio = cst_cache->enzcst_io();
	}

	if ( !cst_cache && !(cst_op->has_residue_pair_constraints()) ) {
		TR << "No constraints found; returning" << std::endl;
		return;
	}
	std::string hb_info("ConstrainHBondNetwork task operation will constrain the following hbnet residues: ");

	for ( core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i ) {
		//NEED TO CHECK FOR LIGAND CASE OR UNNATURALS
		if ( cst_cache && cst_cache->contains_position( i ) ) {
			if ( pose.residue( i ).is_protein() ) {
				hb_info = hb_info + utility::to_string( i ) + "(";

				utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
				utility::vector1< std::string > allowed_name3 = cstio->allowed_res_name3_at_position( pose, i );
				for ( core::Size j =1; j <= allowed_name3.size(); ++j ) {
					allowed_aas[ core::chemical::aa_from_name( allowed_name3[ j ] ) ] = true;
					hb_info = hb_info + allowed_name3[ j ] + ",";
				}
				hb_info = hb_info + "), ";
				task.nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task.nonconst_residue_task( i ).prevent_repacking();
			}
		} else if ( cst_op->residue_pair_constraints_exists(i) ) {
			if ( pose.residue( i ).is_protein() ) {
				hb_info = hb_info + utility::to_string( i ) + "(" + pose.residue(i).name1() + "), ";
				utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
				allowed_aas[ core::chemical::aa_from_name( pose.residue(i).name3() ) ] = true;
				task.nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task.nonconst_residue_task( i ).prevent_repacking();
			}
		}
	}
	TR << hb_info << std::endl;
}

void
ConstrainHBondNetwork::parse_tag( TagCOP , DataMap & )
{

}

std::string
ConstrainHBondNetwork::keyname() {
	return "ConstrainHBondNetwork";
}

void
ConstrainHBondNetwork::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	//TODO: provide a proper XML schema for this TaskOperation to replace the following line.
	core::pack::task::operation::task_op_schema_empty( xsd, keyname() );
}

core::pack::task::operation::TaskOperationOP
ConstrainHBondNetworkCreator::create_task_operation() const
{
	//return new ConstrainHBondNetwork;
	return core::pack::task::operation::TaskOperationOP( new ConstrainHBondNetwork );
}

void
ConstrainHBondNetworkCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const
{
	ConstrainHBondNetwork::provide_xml_schema( xsd );
}

}//namespace hbnet
}//namespace protocols
