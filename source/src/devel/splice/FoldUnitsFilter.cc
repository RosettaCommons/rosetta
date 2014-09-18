// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/FoldUnitsFilter.cc
/// @brief
/// @author Sarel Fleishman


//Unit Headers
#include <devel/splice/FoldUnitsFilter.hh>
#include <devel/splice/FoldUnitsFilterCreator.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <boost/foreach.hpp>
#include <stdio.h>
#include <fstream>
#include <protocols/simple_filters/StemFinderFilter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/PDBInfo.hh>

namespace devel{
namespace splice {

using namespace core;
using namespace std;

static thread_local basic::Tracer TR( "devel.splice.FoldUnitsFilter" );

protocols::filters::FilterOP
FoldUnitsFilterCreator::create_filter() const { return new FoldUnitsFilter; }

std::string
FoldUnitsFilterCreator::keyname() const { return "FoldUnits"; }

//default ctor
FoldUnitsFilter::FoldUnitsFilter() :
protocols::filters::Filter( "FoldUnits" ),
minimal_length_( 20 ),
maximal_length_( 40 ),
ends_distance_( 4.5 ),
filename_( "FoldUnits.out" )
{
}

FoldUnitsFilter::~FoldUnitsFilter() {}

void
FoldUnitsFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	minimal_length( tag->getOption< Size >( "minimal_length", 20 ) );
	maximal_length( tag->getOption< Size >( "maximal_length", 40 ) );
	ends_distance( tag->getOption< Real >( "ends_distance", 4.5 ) );
	filename_ = tag->getOption< string >("filename", "FoldUnits.out" );

	TR<<"options minimal_length: "<<minimal_length()<<" maximal_length: "<<maximal_length()<<" ends_distance: "<<ends_distance()<<" filename: "<<filename()<<std::endl;
}

bool
FoldUnitsFilter::apply( core::pose::Pose const & pose ) const {
	compute( pose );
	return true;
}

void
FoldUnitsFilter::report( std::ostream & , core::pose::Pose const & ) const {
}

core::Real
FoldUnitsFilter::report_sm( core::pose::Pose const & ) const {
	return 0.0;
}

std::string
aa_sequence( pose::Pose const & pose ){
	std::string seq( "" );

	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( !pose.conformation().residue( i ).is_protein() )
			continue;
		seq += pose.conformation().residue( i ).name1();
	}
	return seq;
}

std::string
input_file_name(){
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string const input_file_name( job->input_tag() );
  core::Size const wheres_period( input_file_name.find_last_of( "." ) );
	std::string const pdb_name( input_file_name.substr(0, wheres_period ) );
	return pdb_name;
}

void
FoldUnitsFilter::write_to_file( std::string const pdb, core::Size const from_res, core::Size const to_res, std::string const dssp, std::string const sequence, core::pose::Pose const & pose ) const{
	std::ofstream data;
	data.open( filename_.c_str(), std::ios::app );
	runtime_assert( data );
	data<<pdb<<' '<<from_res<<' '<<to_res<<' '<<dssp<<' '<<sequence<<' ';
	for( Size resi = from_res; resi <= to_res; ++resi ){
		Real const phi = pose.phi( resi );
		Real const psi = pose.psi( resi );
		Real const omega  = pose.omega( resi );
		data<<phi<<' '<<psi<<' '<<omega<<' ';
	}
	data<<std::endl;

	data.close();
}

core::Real
atom_distance( core::pose::Pose const & p1, core::Size const r1, std::string const a1,
               core::pose::Pose const & p2, core::Size const r2, std::string const a2 ){
	  return( p1.conformation().residue( r1 ).xyz( a1 ).distance( p2.conformation().residue( r2 ).xyz( a2 ) ) );
}


bool
find_cut_in_segment( core::pose::Pose const & pose, core::Size const from_res, core::Size const to_res ){
	for( Size resi = from_res; resi < to_res; ++resi ){
		if( atom_distance( pose, resi, "C", pose, resi + 1, "N" ) >= 1.5 )
			return true;
	}
	return false;
}

/// @brief find residue pairs within contact and within the same chain that are at the beginning and end of 2ary structural elements. Output to a database file pdb name, residue boundaries, dssp, aa sequence and phi/psi/omega
core::Real
FoldUnitsFilter::compute(
	core::pose::Pose const & pose
) const {
	scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	string const sec_struct( dssp.get_dssp_secstruct() );
	string const aa_seq( aa_sequence( pose ) );
	string const pdb( input_file_name() );
	pose::PDBInfoCOP pdb_info( pose.pdb_info() );
	TR<<"secstruct: "<<sec_struct<<'\n'<<"sequence: "<<aa_seq<<std::endl;

	for( Size resi = 2; resi <= pose.total_residue(); ++resi ){
		TR<<"resi: "<<resi<<' ';
		if( !pose.conformation().residue( resi ).is_protein() )
			continue;
		if( !(sec_struct[ resi - 1 ] != 'L' && sec_struct[ resi - 2 ] == 'L' ) )
			continue;
		for( Size resj = resi + minimal_length(); resj <= resi + maximal_length() && resj < pose.total_residue(); ++resj ){
			if( !pose.conformation().residue( resj ).is_protein() )
				continue;
			if( pdb_info->chain( resi ) != pdb_info->chain( resj ) )
				continue;
			if( !(sec_struct[ resj - 1 ] != 'L' && sec_struct[ resj ] == 'L' ) )
				continue;

			Real const dist( protocols::simple_filters::res_res_min_distance( pose, resi, pose, resj ) );
			if( dist >= ends_distance() )
				continue;

			if( find_cut_in_segment( pose, resi, resj ) )
				continue;
			write_to_file( pdb, resi, resj, sec_struct.substr( resi - 1, resj - resi + 1 ), aa_seq.substr( resi - 1, resj - resi + 1 ), pose );
		}
	}
	return 0.0;
}

}
}
