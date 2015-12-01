// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/setup/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <core/pose/PDBInfo.hh>
#include <core/init/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/conversions.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

//RNA stuff.
#include <protocols/farna/util.hh>
#include <protocols/farna/BasePairStepLibrary.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

#include <utility/excn/Exceptions.hh>
//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

using namespace core;
using namespace core::pose::rna;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;
using protocols::farna::MAX_BULGE_LENGTH;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( IntegerVector, exclude_res )
OPT_KEY( Boolean,       general_bps )
OPT_KEY( Boolean,       use_lores_base_pair_classification )

///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions_test( ){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	utility::vector1 < std::string >  infiles  = option[ in::file::s ]();
	std::string outfile  = option[ basic::options::OptionKeys::rna::farna::vall_torsions ]();
	utility::vector1< core::Size > const exclude_res_list = option[exclude_res]();

	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ](); // old syntax -- should deprecate.

	utility::io::ozstream torsions_out( outfile );

	for ( Size n = 1; n <= infiles.size(); n++ ) {
		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, *rsd_set, infiles[n] );
		/////////////////////////////////////////
		protocols::farna::make_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		protocols::farna::create_rna_vall_torsions( pose, torsions_out, exclude_res_list );

		std::cout << "***********************************************************" << std::endl;
		std::cout << "Put torsions from PDB file " <<  infiles[n] << " into " << outfile << std::endl;
		std::cout << "***********************************************************" << std::endl;

	}


}

////////////////////////////////////////////////////////////////////////////////////
bool
check_for_contacts( pose::Pose & pose, Size const i,
	Vector const & atom_j, Vector const dir_vector,
	char & edge_i, char & orientation )
{
	static Real const CONTACT_CUTOFF2( 3.0 * 3.0 );

	scoring::rna::RNA_ScoringInfo  & rna_scoring_info( scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	scoring::rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	//Figure out jump.
	conformation::Residue const & rsd_i( pose.residue( i ) );

	bool found_contact( false );
	for ( Size ii=rsd_i.first_sidechain_atom()+1 ; ii<= rsd_i.nheavyatoms(); ++ii ) {
		//  if ( rsd_i.atom_name( ii ) == " O2'" ) continue;
		Real const dist2( (rsd_i.xyz( ii ) - atom_j ).length_squared() ) ;
		if ( dist2 < CONTACT_CUTOFF2 ) {
			//   std::cout << dist2 << " " <<  i << " " << rsd_i.atom_name( ii ) << std::endl;
			found_contact = true;
			break;
		}
	}

	if ( !found_contact ) return false;

	Vector const & centroid_i( base_centroids[i] );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	Vector d_ij = atom_j - centroid_i;
	Real const dist_x = dot_product( d_ij, x_i );
	Real const dist_y = dot_product( d_ij, y_i );
	// Real const dist_z = dot_product( d_ij, z_i );
	// Real const rho2 = dist_x*dist_x + dist_y*dist_y;

	Real const zeta = numeric::conversions::degrees( std::atan2( dist_y, dist_x) );
	if ( std::abs(zeta) < 60.0 ) edge_i = 'W';  //Watson-Crick edge
	else if ( zeta > +60.0 )   edge_i = 'H'; // Hoogsteen edge
	else                       edge_i = 'S'; // Sugar edge

	orientation = dot( dir_vector, z_i ) > 0 ?  'A' : 'P';

	return true;
}

/////////////////////////////////////////////////////////////////////////////////
void
check_for_contacts_and_output_jump_o2prime( pose::Pose & pose, Size const i, Size const j,
	utility::io::ozstream & dataout ){

	char edge_i, orientation;

	conformation::Residue const & rsd_i( pose.residue( i ) );
	conformation::Residue const & rsd_j( pose.residue( j ) );

	std::string const atom_name = " O2'";
	Vector const atom_vector = rsd_j.xyz( atom_name );
	Vector const dir_vector  = rsd_j.xyz( " O2'" ) - rsd_j.xyz( " C2'" );
	if ( !check_for_contacts( pose, i, atom_vector, dir_vector, edge_i, orientation) ) return;

	char const edge_j = '2';

	kinematics::Stub const stub1( rsd_i.xyz( rsd_i.chi_atoms(1)[4] ),
		rsd_i.xyz( rsd_i.chi_atoms(1)[3] ),
		rsd_i.xyz( rsd_i.chi_atoms(1)[2] ) );

	kinematics::Stub const stub2( rsd_j.xyz( atom_name ),
		rsd_j.xyz( " C2'" ),
		rsd_j.xyz( " C3'" ) );

	dataout << "PAIR " <<
		I(5, i) << ' ' << edge_i << ' ' <<
		I(5, j) << ' ' << edge_j << "   " <<
		orientation << "   " <<
		pose.residue(i).name1() << ' ' << pose.residue(j).name1() << " " <<
		rsd_i.atom_name( rsd_i.chi_atoms(1)[4] ) <<  " " <<
		atom_name <<  " " <<
		kinematics::Jump( stub1, stub2) <<
		std::endl;

}


/////////////////////////////////////////////////////////////////////////////////
void
check_for_contacts_and_output_jump_phosphate( pose::Pose & pose, Size const i, Size const j,
	utility::io::ozstream & dataout ){

	char edge_i, orientation;

	conformation::Residue const & rsd_i( pose.residue( i ) );
	conformation::Residue const & rsd_j( pose.residue( j ) );
	conformation::Residue const & prev_rsd( pose.residue( j-1 ) );


	Vector dir_vector  = cross( rsd_j.xyz( " OP1" ) - rsd_j.xyz( " P  " ),
		rsd_j.xyz( " OP2" ) - rsd_j.xyz( " P  " ) );

	std::string atom_name = " OP2";
	Vector atom_vector = rsd_j.xyz( atom_name );
	if ( !check_for_contacts( pose, i, atom_vector, dir_vector, edge_i, orientation) ) {
		atom_name = " OP1";
		atom_vector = rsd_j.xyz( atom_name );
		if ( !check_for_contacts( pose, i, atom_vector, dir_vector, edge_i, orientation) ) return;
	}


	char const edge_j = 'P';

	kinematics::Stub const stub1( rsd_i.xyz( rsd_i.chi_atoms(1)[4] ),
		rsd_i.xyz( rsd_i.chi_atoms(1)[3] ),
		rsd_i.xyz( rsd_i.chi_atoms(1)[2] ) );

	kinematics::Stub const stub2_fwd( rsd_j.xyz( atom_name ),
		rsd_j.xyz( " P  " ),
		prev_rsd.xyz( " O3'" ) );

	kinematics::Stub const stub2_back( rsd_j.xyz( atom_name ),
		rsd_j.xyz( " P  " ),
		rsd_j.xyz( " O5'" ) );

	dataout << "PAIR " <<
		I(5, i) << ' ' << edge_i << ' ' <<
		I(5, j) << ' ' << edge_j << "   " <<
		orientation << "   " <<
		pose.residue(i).name1() << ' ' << pose.residue(j).name1() << " " <<
		rsd_i.atom_name( rsd_i.chi_atoms(1)[4] ) <<  " " <<
		atom_name <<  " " <<
		kinematics::Jump( stub1, stub2_fwd ) <<
		kinematics::Jump( stub1, stub2_back ) <<
		std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////
// JUMP extractor.
void
create_bp_jump_database_test( ){

	using namespace chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::farna;

	utility::vector1< core::Size > const exclude_res_list = option[exclude_res]();

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	std::string infile  = option[ in::file::s ][1];
	std::string outfile  = option[ out::file::o ];

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, *rsd_set, infile );
	make_phosphate_nomenclature_matches_mini( pose );

	// core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );

	// Fill base pairing information... these are
	// all functions used in scoring... see RNA_BaseBaseEnergy.cc
	ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	(*lores_scorefxn)( pose );
	lores_scorefxn->show( std::cout, pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	EnergyBasePairList scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	utility::io::ozstream dataout( outfile );

	for ( EnergyBasePairList::const_iterator it = scored_base_pair_list.begin();
			it != scored_base_pair_list.end(); ++it ) {

		BasePair const base_pair = it->second;

		int const i = base_pair.res1();
		int const k = base_pair.edge1();

		int const j = base_pair.res2();
		int const m = base_pair.edge2();

		if ( is_num_in_list(i, exclude_res_list) || is_num_in_list(j, exclude_res_list) ) {
			continue;
		}

		char const orientation = ( base_pair.orientation() == 1) ? 'A' : 'P';

		char const edge_i = core::chemical::rna::get_edge_from_num( k );
		char const edge_j = core::chemical::rna::get_edge_from_num( m );

		//Figure out jump.
		conformation::Residue const & rsd_i( pose.residue( i ) );
		conformation::Residue const & rsd_j( pose.residue( j ) );
		kinematics::Stub const stub_i( rsd_i.xyz( rsd_i.chi_atoms(1)[4] ),
			rsd_i.xyz( rsd_i.chi_atoms(1)[3] ),
			rsd_i.xyz( rsd_i.chi_atoms(1)[2] ) );
		kinematics::Stub const stub_j( rsd_j.xyz( rsd_j.chi_atoms(1)[4] ),
			rsd_j.xyz( rsd_j.chi_atoms(1)[3] ),
			rsd_j.xyz( rsd_j.chi_atoms(1)[2] ) );

		dataout << "PAIR " <<
			I(5, i) << ' ' << edge_i << ' ' <<
			I(5, j) << ' ' << edge_j << "   " <<
			orientation << "   " <<
			pose.residue(i).name1() << ' ' << pose.residue(j).name1() << " " <<
			rsd_i.atom_name( rsd_i.chi_atoms(1)[4] ) <<  " " <<
			rsd_j.atom_name( rsd_j.chi_atoms(1)[4] ) <<  " " <<
			kinematics::Jump( stub_i, stub_j) <<
			std::endl;

	}

	//How about 2' and Phosphate jumps?
	// Look at each base.
	core::scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		// Neighboring residues making base-phosphate or base-2'OH contacts?
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> ( *iru ) );
			Size const j( edge->get_other_ind(i) );

			if ( is_num_in_list(i, exclude_res_list) || is_num_in_list(j, exclude_res_list) ) {
				continue;
			}

			//   EnergyGraph const & energy_graph( pose.energies().energy_graph() );

			check_for_contacts_and_output_jump_o2prime( pose, i, j, dataout );

			if ( j > 1  ) check_for_contacts_and_output_jump_phosphate( pose, i, j, dataout );

		}
	}


	dataout.close();

	std::cout << "***********************************************************" << std::endl;
	std::cout << "Put jumps from PDB file " <<  infile << " into " << outfile << std::endl;
	std::cout << "***********************************************************" << std::endl;

}

///////////////////////////////////////////////////////////////
std::string
get_bps_seq( utility::vector1< Size > const & base_pair_res, std::string const & sequence  ){
	runtime_assert( base_pair_res.size() == 4 );

	std::string bps_seq;
	bps_seq += sequence[ base_pair_res[1]-1 ];
	bps_seq += sequence[ base_pair_res[2]-1 ];

	bps_seq += "_";

	bps_seq += sequence[ base_pair_res[3]-1 ];
	int seq_sep = int( base_pair_res[ 4 ] ) - int( base_pair_res[ 3 ] );
	if ( seq_sep > 0 &&  seq_sep <= (MAX_BULGE_LENGTH+1) ) {
		for ( Size i = base_pair_res[3]+1; i < base_pair_res[4]; i++ ) bps_seq += 'n';  // bulge
	} else {
		bps_seq += "x"; // greater than 3 nts, or possibly different chain
	}
	bps_seq += sequence[ base_pair_res[4]-1 ];

	return bps_seq;
}

///////////////////////////////////////////////////////////
std::string
get_bps_tag( utility::vector1< Size > const &  base_pair_res, std::string const & infile, pose::Pose const & pose ){

	using namespace ObjexxFCL::format;

	pose::PDBInfoCOP pdb_info = pose.pdb_info();

	std::stringstream pdb_tag_stream;
	pdb_tag_stream << infile;
	for ( Size n = 1; n <= base_pair_res.size(); n++ ) {
		pdb_tag_stream << "_";
		pdb_tag_stream << pdb_info->chain( base_pair_res[1] ) << pdb_info->number( base_pair_res[n] );
	}

	return pdb_tag_stream.str();
}



///////////////////////////////////////////////////////////
std::string
get_bps_sub_dir( int const seq_sep ){
	if ( seq_sep == 1 ) {
		return "standard/";
	} else if ( seq_sep >= 2 && seq_sep <= (MAX_BULGE_LENGTH+1) ) {
		return "bulge_"+ObjexxFCL::format::I( 1, seq_sep - 1 )+"nt/";
	} else {
		return "long_distance/";
	}
}


///////////////////////////////////////////////////////////
std::string
get_bps_sub_dir( utility::vector1< Size > const &  base_pair_res ){
	int seq_sep = int( base_pair_res[ 4 ] ) - int( base_pair_res[ 3 ] );
	return get_bps_sub_dir( seq_sep );
}

///////////////////////////////////////////////////////////////
void
write_base_pair_step_to_silent_struct( pose::Pose const & pose,
	utility::vector1< int > base_pair_res,
	std::string const & intag,
	bool const use_sub_directories = true )
{
	using namespace core::io::silent;
	static SilentFileData silent_file_data;

	pose::Pose bps_pose;
	core::pose::pdbslice( bps_pose, pose, base_pair_res );

	std::string bps_seq = get_bps_seq( base_pair_res, pose.sequence() );
	std::string bps_tag = get_bps_tag( base_pair_res, intag, pose );

	std::cout << "Found base pair step! " << bps_seq << " " << bps_tag << std::endl;

	std::string bps_sub_dir;
	if ( use_sub_directories ) {
		bps_sub_dir = get_bps_sub_dir( base_pair_res );
		if ( !utility::file::file_exists( bps_sub_dir ) ) utility::file::create_directory_recursive( bps_sub_dir );
	}
	std::string const silent_file = bps_sub_dir + bps_seq + ".out";
	BinarySilentStruct s( bps_pose, bps_tag );
	silent_file_data.write_silent_struct( s, silent_file,  false /* score_only */ );
}

///////////////////////////////////////////////////////////////
void
write_base_pair_step_to_silent_struct( pose::Pose const & pose, Size const & i, Size const & j,
	std::string const & intag )
{
	utility::vector1< int > const base_pair_res = utility::tools::make_vector1( i, i+1, j-1, j);
	write_base_pair_step_to_silent_struct( pose, base_pair_res, intag, false /* use_sub_directories */ );
}

///////////////////////////////////////////////////////////////
// Following just gets Watson-Crick and G*U pairs
void
output_canonical_base_pair_steps( pose::Pose & pose,
	std::string const & intag )
{
	std::map< Size, Size > partner;
	protocols::farna::figure_out_base_pair_partner( pose, partner, false );

	for ( std::map< Size, Size >::const_iterator it = partner.begin();
			it != partner.end(); ++it ) {

		Size const i = it->first;
		Size const j = it->second;

		if ( i < pose.total_residue() &&
				partner.find( i+1 ) != partner.end() && j > 1 && partner[ i+1 ] == j-1 &&
				!pose.fold_tree().is_cutpoint(i) && !pose.fold_tree().is_cutpoint(j-1) ) {

			write_base_pair_step_to_silent_struct( pose, i, j, intag );

		}
	}

}

///////////////////////////////////////////////////////////////
// Following outputs all base pair steps (including neighboring non-canonicals)
void
output_general_base_pair_steps( pose::Pose const & pose,
	std::string const & intag )
{

	utility::vector1< core::pose::rna::BasePair > base_pair_list;
	if ( basic::options::option[ use_lores_base_pair_classification ]() ) {
		base_pair_list= protocols::farna::classify_base_pairs_lores( pose );
	} else {
		base_pair_list= classify_base_pairs( pose );
	}

	std::map< std::pair< Size, Size >, bool > partnered;
	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		BasePair const base_pair = base_pair_list[ n ];

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		runtime_assert( pose.residue( i ).is_RNA() );
		runtime_assert( pose.residue( j ).is_RNA() );

		partnered[ std::make_pair( i, j ) ] = true;
		partnered[ std::make_pair( j, i ) ] = true;
	}

	for ( Size i = 1; i < pose.total_residue(); i++ ) {
		for ( Size j = 1; j <= pose.total_residue(); j++ ) {

			if ( partnered[ std::make_pair( i,   j )   ] ) {

				vector1< bool > outputted_bps( pose.total_residue(), false );

				// Looking for i+1 to also base pair to something.
				// base-pair-step, and then bulges: single, double, triple
				for ( int n = 0; n <= MAX_BULGE_LENGTH; n++ ) {

					if ( int(j) > (n+1) && std::abs( int(i) - int(j) ) > (int(n)+1) &&
							partnered[ std::make_pair( i+1, j-n-1 ) ] &&
							!pose.fold_tree().is_cutpoint( i ) ) {

						bool found_cut( false );
						for ( int q = 1; q <= n+1; q++ ) {
							if ( pose.fold_tree().is_cutpoint( j-q ) ) {
								found_cut = true;
							}
						}
						if ( found_cut ) continue;

						write_base_pair_step_to_silent_struct( pose, utility::tools::make_vector1( i, i+1, j-n-1, j), intag, true /*create subdir*/ );

						outputted_bps[ j-n-1 ] = true;

					}
				}


				// scan for adjacent base pair that has one partner distant in sequence; may even involve totally different chain!
				for ( Size k = 1; k <= pose.total_residue(); k++ ) {
					if ( outputted_bps[ k ] ) continue;
					if ( partnered[ std::make_pair( i+1, k ) ] ) {
						write_base_pair_step_to_silent_struct( pose, utility::tools::make_vector1( i, i+1, k, j), intag, true /*create subdir*/ );
					}
				}

			}

		}
	}

}


///////////////////////////////////////////////////////////////
void
create_base_pair_step_database_test( ){

	using namespace chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::farna;
	using namespace protocols::stepwise;

	utility::vector1< core::Size > const exclude_res_list = option[exclude_res]();

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	utility::vector1< std::string > const infiles  = setup::load_s_and_l();

	for ( Size n = 1; n <= infiles.size(); n++ ) {
		std::string const & infile = infiles[ n ];

		std::string intag = infile;
		Size pos( intag.find( ".pdb" ) );
		intag.replace( pos, 4, "" );

		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, *rsd_set, std::string( option[ in::path::pdb ](1) ) + infile );
		core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );
		std::string const sequence = pose.sequence();

		if ( option[ general_bps ]() ) {
			output_general_base_pair_steps( pose, intag );
		} else {
			output_canonical_base_pair_steps( pose, intag );
		}

		std::cout << "***********************************************************" << std::endl;
		std::cout << "Put base pair steps from " <<  infile << " into .out files"  << std::endl;
		std::cout << "***********************************************************" << std::endl;

	}


}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::rna::farna;

	if ( option[ jump_database ] ) {
		create_bp_jump_database_test();
	} else if ( option[ vall_torsions ].user() ) {
		create_rna_vall_torsions_test();
	} else if ( option[ bps_database ]() ) {
		create_base_pair_step_database_test();
	} else {
		std::cout << std::endl;
		std::cout << "Please specify: " << std::endl;
		std::cout << "  -vall_torsions   for generating torsions library " << std::endl;
		std::cout << "  -jump_database   for generating database of rigid-body orientations " << std::endl;
		std::cout << "  -bps_database    for generating database of Watson-Crick base pair steps " << std::endl;
		std::cout << std::endl;
	}

	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		utility::vector1< Size > blank_size_vector;


		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb1> [... <more pdbs>] -vall_torsions [name of output torsions file to create] " << std::endl;
		std::cout <<              "              " << argv[0] << "  -s <pdb1> [... <more pdbs>] -jump_database [name of rigid-body orientation database to create] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		NEW_OPT( exclude_res, "Residues excluded for database creation (works for one file only)", blank_size_vector );
		NEW_OPT( general_bps, "For bps_database, output base pair steps involving noncanonicals", false );
		NEW_OPT( use_lores_base_pair_classification, "Use loose base-pair classifier (same as used in FARNA)", false );
		option.add_relevant( basic::options::OptionKeys::rna::farna::vall_torsions );
		option.add_relevant( basic::options::OptionKeys::rna::farna::jump_database );
		option.add_relevant( basic::options::OptionKeys::rna::farna::bps_database );


		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);


		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
