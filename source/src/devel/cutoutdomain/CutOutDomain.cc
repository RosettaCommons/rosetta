// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/cutoutdomain/CutOutdomain.cc
/// @brief
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit headers
#include <devel/cutoutdomain/CutOutDomain.hh>
#include <devel/cutoutdomain/CutOutDomainCreator.hh>

// Package headers
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/moves/Mover.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/util.hh>
#include <core/pose/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace devel {
namespace cutoutdomain {

static basic::Tracer TR( "devel.cutoutdomain.CutOutDomain" );

CutOutDomain::CutOutDomain()
: protocols:: moves::Mover("CutOutDomain")
{
}

CutOutDomain::~CutOutDomain() = default;

void
CutOutDomain::apply( core::pose::Pose & pose )
{
	core::pose::Pose Temp_pose;
	core::import_pose::pose_from_file( Temp_pose, source_pdb_name_ , core::import_pose::PDB_file);
	core::Size from = find_nearest_res(Temp_pose,pose,start_res_, 1/*chain*/ );
	TR<<from<<std::endl;
	core::Size to  = find_nearest_res(Temp_pose,pose,end_res_, 1/*chain*/ );
	//TR<<to<<std::endl;
	TR<<"Start resdiue on the source pose "<<source_pdb_name_<<" : "<<Temp_pose.residue(start_res_).name1()<<start_res_<<std::endl;
	TR<<"Start resdiue on the target pose "<<pose.pdb_info()->name()<<" : "<<pose.residue(from).name1()<<from<<std::endl;
	TR<<"End resdiue on the source pose "<<source_pdb_name_<<" : "<<Temp_pose.residue(end_res_).name1()<<end_res_<<std::endl;
	TR<<"End resdiue on the target pose "<<pose.pdb_info()->name()<<" : "<<pose.residue(to).name1()<<to<<std::endl;
	pose.conformation().delete_residue_range_slow( to+delta_c_ter_+1,pose.size() );
	pose.conformation().delete_residue_range_slow( 1,from+1-delta_n_ter_ );
	//hack so Rosetta doesn't call scoring function at the end
	pose.dump_pdb(suffix_);
}

void
CutOutDomain::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	start_res_ = tag->getOption<core::Size>( "start_res", 1 );
	end_res_ = tag->getOption<core::Size>( "end_res", 1 );
	source_pdb_name( tag->getOption< std::string >( "source_pdb" ));
	suffix( tag->getOption< std::string >( "suffix" ,""));
	delta_n_ter_ = tag->getOption<core::Size>( "delta_n_ter", 0 );
	delta_c_ter_ = tag->getOption<core::Size>( "delta_c_ter", 0 );
}

core::Size
CutOutDomain::find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain/*=0*/ ){
	core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
	core::Size i;
	for ( i = 1; i < target.size(); ++i ) {
		if ( target.residue( i ).is_ligand() ) continue;
		if ( chain && target.residue( i ).chain() != chain ) continue;
		// TR<<"the residue examnied is:"<<i<<target.residue(i).name1()<<std::endl;
		core::Real const dist( source.residue( res ).xyz( "CA" ).distance( target.residue( i ).xyz( "CA" ) ) );
		TR<<"distance between source res:"<<res<<source.residue( res ).name1()<<" and target res: "<<i<<target.residue(i).name1();
		TR<<" is: " <<dist<<std::endl;
		if ( dist <= min_dist ) {
			min_dist = dist;
			nearest_res = i;
		}
	}
	static basic::Tracer TR( "This is nearest_res" );

	TR<<nearest_res<<std::endl;
	if ( min_dist <= 10.0 ) return nearest_res;
	else return 0;
}

std::string CutOutDomain::get_name() const {
	return mover_name();
}

std::string CutOutDomain::mover_name() {
	return "CutOutDomain";
}

void CutOutDomain::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"start_res", xsct_non_negative_integer,
		"begin residue on the template pdb",
		"1");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"end_res", xsct_non_negative_integer,
		"end residue on the template pdb",
		"1");

	attlist + XMLSchemaAttribute::required_attribute(
		"source_pdb", xs_string,
		"name of the pdb to be cut");

	attlist + XMLSchemaAttribute(
		"suffix", xs_string,
		"suffix of the outputted structure");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"delta_n_ter", xsct_non_negative_integer,
		"do not delete n N-terminal residues",
		"0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"delta_c_ter", xsct_non_negative_integer,
		"do not delete n C-terminal residues",
		"0");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"This mover takes a template pdb and cuts the active pose accroding to "
		"start and end position relative to the template pose",
		attlist );
}

std::string CutOutDomainCreator::keyname() const {
	return CutOutDomain::mover_name();
}

protocols::moves::MoverOP
CutOutDomainCreator::create_mover() const {
	return protocols::moves::MoverOP( new CutOutDomain );
}

void CutOutDomainCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CutOutDomain::provide_xml_schema( xsd );
}

} //cutoutdomain
} //devel
