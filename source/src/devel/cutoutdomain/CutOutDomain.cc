// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/io/pdb/pose_io.hh>
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

namespace devel {
namespace cutoutdomain {

static basic::Tracer TR( "devel.cutoutdomain.CutOutDomain" );

std::string
CutOutDomainCreator::keyname() const
{
	return CutOutDomainCreator::mover_name();
}

protocols::moves::MoverOP
CutOutDomainCreator::create_mover() const {
	return new CutOutDomain;
}

std::string
CutOutDomainCreator::mover_name()
{
	return "CutOutDomain";
}

CutOutDomain::CutOutDomain()
	: protocols:: moves::Mover("CutOutDomain")
{
}

CutOutDomain::~CutOutDomain() {}

void
CutOutDomain::apply( core::pose::Pose & pose )
{
	core::pose::Pose Temp_pose;
	core::import_pose::pose_from_pdb( Temp_pose, source_pdb_name_ );
	core::Size from = find_nearest_res(Temp_pose,pose,start_res_, 1/*chain*/ );
		TR<<from<<std::endl;
		core::Size to  = find_nearest_res(Temp_pose,pose,end_res_, 1/*chain*/ );
		//TR<<to<<std::endl;
		TR<<"Start resdiue on the source pose "<<source_pdb_name_<<" : "<<Temp_pose.residue(start_res_).name1()<<start_res_<<std::endl;
		TR<<"Start resdiue on the target pose "<<pose.pdb_info()->name()<<" : "<<pose.residue(from).name1()<<from<<std::endl;
		TR<<"End resdiue on the source pose "<<source_pdb_name_<<" : "<<Temp_pose.residue(end_res_).name1()<<end_res_<<std::endl;
		TR<<"End resdiue on the target pose "<<pose.pdb_info()->name()<<" : "<<pose.residue(to).name1()<<to<<std::endl;
		pose.conformation().delete_residue_range_slow( to+delta_c_ter_+1,pose.total_residue() );
		pose.conformation().delete_residue_range_slow( 1,from+1-delta_n_ter_ );
		//hack so Rosetta doesn't call scoring function at the end
		pose.dump_pdb(suffix_);
}


std::string
CutOutDomain::get_name() const {
	return CutOutDomainCreator::mover_name();
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
  for( i = 1; i < target.total_residue(); ++i ){
		if( target.residue( i ).is_ligand() ) continue;
		if( chain && target.residue( i ).chain() != chain ) continue;
		// TR<<"the residue examnied is:"<<i<<target.residue(i).name1()<<std::endl;
    core::Real const dist( source.residue( res ).xyz( "CA" ).distance( target.residue( i ).xyz( "CA" ) ) );
    TR<<"distance between source res:"<<res<<source.residue( res ).name1()<<" and target res: "<<i<<target.residue(i).name1();
    TR<<" is: " <<dist<<std::endl;
    if( dist <= min_dist ){
      min_dist = dist;
      nearest_res = i;
    }
  }
  static basic::Tracer TR("This is nearest_res");

      TR<<nearest_res<<std::endl;
  if( min_dist <= 10.0 ) return nearest_res;
  else return 0;
}
} //cutoutdomain
} //devel
