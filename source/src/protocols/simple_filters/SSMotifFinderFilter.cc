// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <fstream>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/Jump.hh>
#include <core/import_pose/import_pose.hh>
//#include <boost/regex.hpp>
//#include <boost/regex/v4/regex_search.hpp>
//#include <boost/xpressive/xpressive.hpp>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/toolbox/superimpose.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.SSMotifFinder" );



//default ctor
SSMotifFinder::SSMotifFinder() :
	protocols::filters::Filter( "SSMotifFinder" ),
	from_res_( 0 ),
	to_res_( 0 ),
	template_pose_( /* NULL */ ),
	template_stem1_( 0 ),
	template_stem2_ ( 0 ),
	rmsd_( 0.0 ),
	filename_( "" ),
	pdbname_( "" )
{
}

SSMotifFinder::~SSMotifFinder() = default;

void
SSMotifFinder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) {
	from_res( tag->getOption< core::Size >( "from_res" ) );
	to_res( tag->getOption< core::Size >( "to_res" ) );
	rmsd( tag->getOption< core::Real >( "rmsd" ) );
	filename( tag->getOption< std::string >( "filename" ) );
	pdbname( tag->getOption< std::string >( "pdbname" ) );

	std::string template_pose_filename( tag->getOption< std::string >( "template_pose" ) );
	template_pose_ = utility::pointer::make_shared< core::pose::Pose >();
	core::import_pose::pose_from_file( *template_pose_, template_pose_filename , core::import_pose::PDB_file);
	template_stem1(core::pose::parse_resnum(tag->getOption<std::string>("template_stem1", "0"), *template_pose_));
	template_stem2(core::pose::parse_resnum(tag->getOption<std::string>("template_stem2", "0"), *template_pose_));

	jump_ = compute_jump( *template_pose_, template_stem1(), template_stem2() );

	TR<<"SSMotifFinder with options: from_res "<<from_res()<<", to_res "<<to_res()<<", template_stem1 "<<template_stem1()<<", template_stem2 "<<template_stem2()<<", rmsd "<<rmsd()<<", filename "<<filename()<<", template_pose "<<template_pose_filename<<", pdbname "<<pdbname()<<std::endl;
}

core::Real
atom_distance( core::conformation::Residue const & r1, std::string const & a1,
	core::conformation::Residue const & r2, std::string const & a2 ){
	return( r1.xyz( a1 ).distance( r2.xyz( a2 ) ) );
}

core::Real
res_rmsd( utility::vector1< core::Size > const pose_res_for_rmsd, utility::vector1< core::Size > const template_res_for_rmsd, core::pose::Pose const & copy_pose, core::pose::Pose const & template_pose ){
	//two_res_rmsd( core::conformation::Residue const & r1, core::conformation::Residue const & r2 ){
	core::Real rmsd( 0.0 );
	core::Size count = 0.0;
	runtime_assert( pose_res_for_rmsd.size() == template_res_for_rmsd.size() );

	for ( core::Size i = 1; i <= pose_res_for_rmsd.size(); i++ ) {
		core::conformation::Residue r1 = copy_pose.conformation().residue( pose_res_for_rmsd[i] );
		core::conformation::Residue r2 = template_pose.conformation().residue( template_res_for_rmsd[i] );
		rmsd += pow( atom_distance( r1, "N", r2, "N"   ), 2.0 );
		rmsd += pow( atom_distance( r1, "O", r2, "O"   ), 2.0 );
		rmsd += pow( atom_distance( r1, "CA", r2, "CA" ), 2.0 );
		rmsd += pow( atom_distance( r1, "C", r2, "C"   ), 2.0 );
		count = count + 4;
		//if ( (r1.name3()!="GLY") && (r2.name3()!="GLY") ){
		//rmsd += pow( atom_distance( r1, "CB", r2, "CB" ), 2.0 );
		//count++;
		//}
	}
	rmsd /= (core::Real) count;
	rmsd = pow( rmsd, 0.5 );
	return rmsd;
}

void
write_to_file( std::string const & filename, core::Size const stem1, core::Size const stem2, core::Real const rmsd, std::string const & pdbname, core::pose::Pose const & pose ){
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	std::string const chain1( pdb_info->chain( stem1 ) ), chain2( pdb_info->chain( stem2 ) );
	int const resnum1( pdb_info->number( stem1 )), resnum2( pdb_info->number( stem2 ) );

	std::ofstream stem_file;
	stem_file.open( filename.c_str(), std::ios::app );
	runtime_assert( stem_file.good() );
	stem_file<<pdbname<<' '<<resnum1<<chain1<<' '<<resnum2<<chain2<<' '<<rmsd<<'\n';
}

utility::vector1< numeric::xyzVector< core::Real > >
coords( core::pose::Pose const & pose, utility::vector1< core::Size > const positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	for ( core::Size const pos : positions ) {
		coords.push_back( pose.residue( pos ).xyz( "N" ) );
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
		coords.push_back( pose.residue( pos ).xyz( "C" ) );
		coords.push_back( pose.residue( pos ).xyz( "O" ) );
		//if (pose.residue( pos ).name3()!="GLY")
		//coords.push_back( pose.residue( pos ).xyz( "CB" ) );
	}
	return coords;
}

void
SSMotifFinder::superimpose_pose_on_template( core::pose::Pose const & template_pose, core::pose::Pose & copy_pose, core::Size const pose_stem1, core::Size const pose_stem2) const{
	using namespace protocols::toolbox;

	utility::vector1< core::Size > pose_positions, template_positions;
	pose_positions.clear(); template_positions.clear();
	template_positions.push_back( template_stem1() -1 );
	template_positions.push_back( template_stem1() );
	template_positions.push_back( template_stem1() +1 );
	template_positions.push_back( template_stem2() -1 );
	template_positions.push_back( template_stem2() );
	template_positions.push_back( template_stem2() +1 );
	pose_positions.push_back( pose_stem1 -1 );
	pose_positions.push_back( pose_stem1 );
	pose_positions.push_back( pose_stem1 +1 );
	pose_positions.push_back( pose_stem2 -1 );
	pose_positions.push_back( pose_stem2 );
	pose_positions.push_back( pose_stem2 +1 );

	utility::vector1< numeric::xyzVector< core::Real > > init_coords( coords( copy_pose, pose_positions ) ), ref_coords( coords( template_pose, template_positions )
	);

	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;

	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );

	apply_superposition_transform( copy_pose, rotation, to_init_center, to_fit_center );
}

bool
SSMotifFinder::apply( core::pose::Pose const & pose ) const {
	using namespace std;
	using core::Size;
	using core::Real;
	using utility::vector1;

	core::scoring::dssp::Dssp pose_dssp( pose );
	pose_dssp.dssp_reduced(); // switch to simplified H E L notation
	string const sec_struct( pose_dssp.get_dssp_secstruct() );

	core::scoring::dssp::Dssp template_dssp( *template_pose_ );
	template_dssp.dssp_reduced(); // switch to simplified H E L notation
	string const template_sec_struct( template_dssp.get_dssp_secstruct() );


	// xpressive::match_results< std::string::const_iterator > regex_results;
	// vector1< pair< core::Size/*start*/, core::Size/*length*/ > > submatches;
	// submatches.clear();
	// TR<<"Found following submatches: ";
	// while( regex_search( sec_struct.begin(), sec_struct.end(), regex_results, motif_ ) ){
	//  core::Size const start( regex_results.position() );
	//  core::Size const length( regex_results.length() );
	//  submatches.push_back( pair< core::Size, core::Size >( start, length ) );

	//  TR<<start<<": "<<sec_struct.substr( start, length )<<'\n';
	// }
	// TR<<std::endl;

	vector1< pair< core::Size/*stem1*/, core::Size/*stem2*/ > > stems;
	stems.clear();
	for ( core::Size stem1 = 2; stem1 <= pose.size() - from_res() - 1 ; ++stem1 ) { //changed this line from for(core::Size stem1 = 1; stem1 <= pose.size() - from_res() ; ++stem1 )
		if ( sec_struct[ stem1 - 1 ] != template_sec_struct[ template_stem1() - 1 ] ||
				!pose.conformation().residue( stem1 ).is_protein() ) {
			continue;
		}
		for ( core::Size stem2 = stem1 + from_res(); stem2 <= stem1 + to_res() && stem2 <= pose.size() - 1; ++stem2 ) { //changed this line from for( core::Size stem2 = stem1 + from_res(); stem2 <= stem1 + to_res() && stem2 <= pose.size() ; ++stem2 )
			if ( sec_struct[ stem2 - 1 ] != template_sec_struct[ template_stem2() - 1 ] ||
					!pose.conformation().residue( stem2 ).is_protein() ) {
				continue;
			}
			stems.push_back( pair< core::Size, core::Size >( stem1, stem2 ) );
		}
	}
	// sort( stems.begin(), stems.end() );
	// vector1< pair< core::Size, core::Size > >::iterator it = unique( stems.begin(), stems.end() );
	// stems.resize( std::distance( stems.begin(), it ) );
	TR<<"List of potential pairs: ";
	for ( auto const & stem : stems ) {
		TR<< stem.first<<' '<< stem.second<<'\n';
	}
	TR<<std::endl;

	for ( auto const & stem : stems ) {
		core::Size const stem1( stem.first ), stem2( stem.second );
		TR<<"Testing stems: "<<stem1<<' '<<stem2<<std::endl;
		//core::conformation::Residue const res_stem2( pose.conformation().residue( stem2 ) );
		core::pose::Pose copy_pose( pose );
		//copy_pose.conformation().append_residue_by_jump( res_stem2, stem1 );
		//TR<<"new foldtree: "<<copy_pose.fold_tree()<<std::endl;
		//copy_pose.set_jump( copy_pose.num_jump(), jump() );
		superimpose_pose_on_template( *template_pose_, copy_pose, stem1, stem2);

		vector1< core::Size > pose_res_for_rmsd, template_res_for_rmsd;
		pose_res_for_rmsd.clear(); template_res_for_rmsd.clear();
		pose_res_for_rmsd.push_back( stem1 - 1 );
		pose_res_for_rmsd.push_back( stem1 );
		pose_res_for_rmsd.push_back( stem1 + 1 );
		pose_res_for_rmsd.push_back( stem2 - 1 );
		pose_res_for_rmsd.push_back( stem2 );
		pose_res_for_rmsd.push_back( stem2 + 1 );

		template_res_for_rmsd.push_back( template_stem1() - 1 );
		template_res_for_rmsd.push_back( template_stem1() );
		template_res_for_rmsd.push_back( template_stem1() + 1 );
		template_res_for_rmsd.push_back( template_stem2() - 1 );
		template_res_for_rmsd.push_back( template_stem2() );
		template_res_for_rmsd.push_back( template_stem2() + 1 );

		Real const motif_rmsd = res_rmsd( pose_res_for_rmsd, template_res_for_rmsd, copy_pose, *template_pose_ );
		//Real const motif_rmsd = two_res_rmsd( res_stem2, copy_pose.conformation().residue( copy_pose.size() ) );

		TR<<"rmsd="<<motif_rmsd<<std::endl;
		if ( motif_rmsd <= rmsd() ) {
			write_to_file( filename_, stem1, stem2, motif_rmsd, pdbname_, pose );
		}
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
	// copy_pose.dump_pdb( "template_pose_with_stem2.pdb" );

	/* core::kinematics::FoldTree ft;
	ft.clear();
	runtime_assert( end <= pose.size() );
	runtime_assert( end > start );
	ft.add_edge( 1, start, -1 );
	ft.add_edge( start, end, 1 );
	ft.add_edge( end, start + 1, -1 );
	ft.add_edge( end, pose.size(), -1 );

	core::pose::Pose copy_pose( pose );
	copy_pose.fold_tree( ft );*/
	TR<<"Template ft: "<<copy_pose.fold_tree()<<std::endl;
	return( copy_pose.jump( 1 ) );
}

std::string SSMotifFinder::name() const {
	return class_name();
}

std::string SSMotifFinder::class_name() {
	return "SSMotifFinder";
}

void SSMotifFinder::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("from_res", xsct_non_negative_integer, "Starting residue number; default to 0")
		+ XMLSchemaAttribute::required_attribute("to_res", xsct_non_negative_integer, "Ending residue number; default to 0")
		+ XMLSchemaAttribute::required_attribute("rmsd", xsct_real, "RMSD to base filter on")
		+ XMLSchemaAttribute::required_attribute("filename", xs_string, "file name to write to")
		+ XMLSchemaAttribute::required_attribute("pdb_name", xs_string, "name of input pdb")
		+ XMLSchemaAttribute::required_attribute("template_pose", xs_string, "filename of input template")
		+ XMLSchemaAttribute::attribute_w_default("template_stem1", xsct_refpose_enabled_residue_number, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default("template_stem2", xsct_refpose_enabled_residue_number, "XRW TO DO", "0");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SSMotifFinderFilterCreator::keyname() const {
	return SSMotifFinder::class_name();
}

protocols::filters::FilterOP
SSMotifFinderFilterCreator::create_filter() const {
	return utility::pointer::make_shared< SSMotifFinder >();
}

void SSMotifFinderFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SSMotifFinder::provide_xml_schema( xsd );
}


}
}
