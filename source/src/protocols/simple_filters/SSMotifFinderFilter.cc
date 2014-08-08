// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SSMotifFinderFilter.cc
/// @brief
/// @author Ravit Netzer & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/SSMotifFinderFilter.hh>
#include <protocols/simple_filters/SSMotifFinderFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <fstream>
#include <utility/io/izstream.hh>
#include <sstream>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/Jump.hh>
#include <core/import_pose/import_pose.hh>
//#include <boost/regex.hpp>
//#include <boost/regex/v4/regex_search.hpp>
#include <boost/foreach.hpp>
//#include <boost/xpressive/xpressive.hpp>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/Edge.hh>
#include <core/pose/PDBInfo.hh>

namespace protocols{
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.SSMotifFinder" );

protocols::filters::FilterOP
SSMotifFinderFilterCreator::create_filter() const { return new SSMotifFinder; }

std::string
SSMotifFinderFilterCreator::keyname() const { return "SSMotifFinder"; }

//default ctor
SSMotifFinder::SSMotifFinder() :
protocols::filters::Filter( "SSMotifFinder" ),
from_res_( 0 ),
to_res_( 0 ),
template_pose_( NULL ),
template_stem1_( 0 ),
template_stem2_ ( 0 ),
rmsd_( 0.0 ),
filename_( "" ),
pdbname_( "" )
{
}

SSMotifFinder::~SSMotifFinder() {}

void
SSMotifFinder::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	from_res( tag->getOption< core::Size >( "from_res" ) );
	to_res( tag->getOption< core::Size >( "to_res" ) );
	rmsd( tag->getOption< core::Real >( "rmsd" ) );
	filename( tag->getOption< std::string >( "filename" ) );
	pdbname( tag->getOption< std::string >( "pdbname" ) );

	std::string template_pose_filename( tag->getOption< std::string >( "template_pose" ) );
	template_pose_ = new core::pose::Pose;
	core::import_pose::pose_from_pdb( *template_pose_, template_pose_filename );
	template_stem1(core::pose::parse_resnum(tag->getOption<std::string>("template_stem1", "0"), *template_pose_));
	template_stem2(core::pose::parse_resnum(tag->getOption<std::string>("template_stem2", "0"), *template_pose_));

	jump_ = compute_jump( *template_pose_, template_stem1(), template_stem2() );

	TR<<"SSMotifFinder with options: from_res "<<from_res()<<", to_res "<<to_res()<<", template_stem1 "<<template_stem1()<<", template_stem2 "<<template_stem2()<<", rmsd "<<rmsd()<<", filename "<<filename()<<", template_pose "<<template_pose_filename<<", pdbname "<<pdbname()<<std::endl;
}

core::Real
atom_distance( core::conformation::Residue const & r1, std::string const a1,
               core::conformation::Residue const & r2, std::string const a2 ){
	  return( r1.xyz( a1 ).distance( r2.xyz( a2 ) ) );
}

core::Real
two_res_rmsd( core::conformation::Residue const & r1, core::conformation::Residue const & r2 ){
	core::Real rmsd( 0.0 );

	rmsd += pow( atom_distance( r1, "N", r2, "N"   ), 2.0 );
	rmsd += pow( atom_distance( r1, "O", r2, "O"   ), 2.0 );
	rmsd += pow( atom_distance( r1, "CA", r2, "CA" ), 2.0 );
	rmsd += pow( atom_distance( r1, "C", r2, "C"   ), 2.0 );
	if ( (r1.name3()!="GLY") && (r2.name3()!="GLY") ){
		rmsd += pow( atom_distance( r1, "CB", r2, "CB" ), 2.0 );
		rmsd /= 5.0;
	}
	else
		rmsd /= 4.0;
	rmsd = pow( rmsd, 0.5 );
	return rmsd;
}

void
write_to_file( std::string const filename, core::Size const stem1, core::Size const stem2, core::Real const rmsd, std::string const pdbname, core::pose::Pose const & pose ){
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	char const chain1( pdb_info->chain( stem1 ) ), chain2( pdb_info->chain( stem2 ) );
	int const resnum1( pdb_info->number( stem1 )), resnum2( pdb_info->number( stem2 ) );

  std::ofstream stem_file;
	stem_file.open( filename.c_str(), std::ios::app );
	runtime_assert( stem_file );
	stem_file<<pdbname<<' '<<resnum1<<chain1<<' '<<resnum2<<chain2<<' '<<rmsd<<'\n';
}

bool
SSMotifFinder::apply( core::pose::Pose const & pose ) const {
	using namespace std;
	using core::Size;
	using core::Real;
	using namespace boost;
	using utility::vector1;

  core::scoring::dssp::Dssp pose_dssp( pose );
	pose_dssp.dssp_reduced(); // switch to simplified H E L notation
	string const sec_struct( pose_dssp.get_dssp_secstruct() );

	core::scoring::dssp::Dssp template_dssp( *template_pose_ );
	template_dssp.dssp_reduced(); // switch to simplified H E L notation
	string const template_sec_struct( template_dssp.get_dssp_secstruct() );


//	xpressive::match_results< std::string::const_iterator > regex_results;
//	vector1< pair< Size/*start*/, Size/*length*/ > > submatches;
//	submatches.clear();
//	TR<<"Found following submatches: ";
//	while( regex_search( sec_struct.begin(), sec_struct.end(), regex_results, motif_ ) ){
//		Size const start( regex_results.position() );
//		Size const length( regex_results.length() );
//		submatches.push_back( pair< Size, Size >( start, length ) );

//		TR<<start<<": "<<sec_struct.substr( start, length )<<'\n';
//	}
//	TR<<std::endl;

	vector1< pair< Size/*stem1*/, Size/*stem2*/ > > stems;
	stems.clear();
	for( Size stem1 = 1; stem1 <= pose.total_residue() - from_res(); ++stem1 ){
		if( sec_struct[ stem1 - 1 ] != template_sec_struct[ template_stem1() - 1 ] ||
				!pose.conformation().residue( stem1 ).is_protein() )
			continue;
		for( Size stem2 = stem1 + from_res(); stem2 <= stem1 + to_res() && stem2 <= pose.total_residue(); ++stem2 ){
			if( sec_struct[ stem2 - 1 ] != template_sec_struct[ template_stem2() - 1 ] ||
			    !pose.conformation().residue( stem2 ).is_protein() )
				continue;
			stems.push_back( pair< Size, Size >( stem1, stem2 ) );
		}
	}
//	sort( stems.begin(), stems.end() );
//	vector1< pair< Size, Size > >::iterator it = unique( stems.begin(), stems.end() );
//	stems.resize( std::distance( stems.begin(), it ) );
	TR<<"List of potential pairs: ";
	for( vector1< pair< Size, Size > >::const_iterator vit = stems.begin(); vit != stems.end(); ++vit )
		TR<<vit->first<<' '<<vit->second<<'\n';
	TR<<std::endl;

	for( vector1< pair< Size, Size > >::const_iterator vit = stems.begin(); vit != stems.end(); ++vit ){
		Size const stem1( vit->first ), stem2( vit->second );
		TR<<"Testing stems: "<<stem1<<' '<<stem2<<std::endl;
		core::conformation::Residue const res_stem2( pose.conformation().residue( stem2 ) );
		core::pose::Pose copy_pose( pose );
		copy_pose.conformation().append_residue_by_jump( res_stem2, stem1 );
		TR<<"new foldtree: "<<copy_pose.fold_tree()<<std::endl;
		copy_pose.set_jump( copy_pose.num_jump(), jump() );
		Real const motif_rmsd = two_res_rmsd( res_stem2, copy_pose.conformation().residue( copy_pose.total_residue() ) );
		TR<<"rmsd="<<motif_rmsd<<std::endl;
		if( motif_rmsd <= rmsd() )
			write_to_file( filename_, stem1, stem2, motif_rmsd, pdbname_, pose );
	}

	return true;
}

void
SSMotifFinder::report( std::ostream &, core::pose::Pose const & ) const {
}

core::Real
SSMotifFinder::report_sm( core::pose::Pose const & ) const {
	return 0.0;
}

core::kinematics::Jump
SSMotifFinder::jump() const{ return jump_; }

core::kinematics::Jump
SSMotifFinder::compute_jump( core::pose::Pose const & pose, core::Size const start, core::Size const end ) const {
	using namespace core::kinematics;

	runtime_assert( pose.conformation().num_chains() == 1 );
	core::pose::Pose copy_pose( pose );
	core::conformation::Residue const res( pose.conformation().residue( end ) );
	copy_pose.conformation().append_residue_by_jump( res, start );
//	copy_pose.dump_pdb( "template_pose_with_stem2.pdb" );

/*	core::kinematics::FoldTree ft;
	ft.clear();
	runtime_assert( end <= pose.total_residue() );
	runtime_assert( end > start );
  ft.add_edge( 1, start, -1 );
	ft.add_edge( start, end, 1 );
	ft.add_edge( end, start + 1, -1 );
	ft.add_edge( end, pose.total_residue(), -1 );

	core::pose::Pose copy_pose( pose );
	copy_pose.fold_tree( ft );*/
	TR<<"Template ft: "<<copy_pose.fold_tree()<<std::endl;
	return( copy_pose.jump( 1 ) );
}

}
}
