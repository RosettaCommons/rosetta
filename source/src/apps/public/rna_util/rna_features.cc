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
/// @author Rhiju, rhiju@stanford.edu

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  RNA features -- produce a huge dump of features that can be regressed against 'deep chemical profiling' data
//  that is being collected in the Das lab.  -- initially craeted on May 7, 2013 by Rhiju.
//
//  Uses new RDAT 'feature' output format.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/sasa.hh>
#include <protocols/farna/util.hh>
#include <protocols/farna/RNA_BasePairClassifier.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>

#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/util.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef
// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <numeric/xyz.functions.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

using namespace core;
//using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace core::chemical::rna;

using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////////////////
std::string get_WC_atom( core::chemical::AA const & res_type ){
	using namespace core::chemical;
	std::string WC_atom( "" );
	if ( res_type == na_rad ) WC_atom = " N1 ";
	if ( res_type == na_rcy ) WC_atom = " N3 ";
	if ( res_type == na_rgu ) WC_atom = " N1 ";
	if ( res_type == na_ura ) WC_atom = " N3 ";
	return WC_atom;
}


///////////////////////////////////////////////////////////////////////////
void
save_feature( 	vector1< std::string > & feature_names,
								vector1< Real > & feature_vals,
								Size & feature_counter,
								std::string const feature_name,
								Real const feature_val ){

	feature_counter++;

	if ( feature_counter > feature_names.size() ) feature_names.push_back( feature_name );
	runtime_assert( feature_counter <= feature_names.size() );

	feature_vals.push_back( feature_val );
}

///////////////////////////////////////////////////////////////////////////////
void
update_edge_paired(  Size const i, Size const k,
										 vector1< bool > & wc_edge_paired,
										 vector1< bool > & hoogsteen_edge_paired,
										 vector1< bool > & sugar_edge_paired ) {

	if (k == WATSON_CRICK ) wc_edge_paired[i] = true;
	else if (k == HOOGSTEEN ) hoogsteen_edge_paired[i] = true;
	else {
		runtime_assert( k == SUGAR );
		sugar_edge_paired[i] = true;
	}

}

///////////////////////////////////////////////////////////////////////////////
void
get_rna_base_pairing_status( core::pose::Pose & pose,
														 vector1< bool > & wc_edge_paired,
														 vector1< bool > & hoogsteen_edge_paired,
														 vector1< bool > & sugar_edge_paired,
														 vector1< bool > & is_bulged ){

	using namespace core::scoring;

	//initialize
	is_bulged.clear();
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		wc_edge_paired.push_back( false );
		hoogsteen_edge_paired.push_back( false );
		sugar_edge_paired.push_back( false );
	}

	// score pose with low res scorefunction -- will determine base pairs
	vector1< core::scoring::rna::Base_pair > base_pair_list;

	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	scorefxn->set_weight( hbond_sc, 1.0 ); // used for HBondSet
	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	protocols::farna::classify_base_pairs( pose, base_pair_list, is_bulged ); // could also get 'energies' for each base pair...

	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {
		scoring::rna::Base_pair const base_pair = base_pair_list[ n ];
		update_edge_paired(  base_pair.res1, base_pair.edge1, wc_edge_paired, hoogsteen_edge_paired, sugar_edge_paired );
		update_edge_paired(  base_pair.res2, base_pair.edge2, wc_edge_paired, hoogsteen_edge_paired, sugar_edge_paired );
	}

}

///////////////////////////////////////////////////////////////////////
bool
check_hbonded( pose::Pose const & pose,
							 Size const & i /*residue number*/,
							 std::string const & atom_name,
							 scoring::hbonds::HBondSetOP hbond_set,
							 bool is_acceptor ){

	using namespace core::scoring::hbonds;
	// need to know what atoms could be in this nt.
	// go through list of donors and list of acceptors

	// check in HBond List
	for (Size n = 1; n <= hbond_set->nhbonds(); n++ ) {
		HBond const & hbond( hbond_set->hbond( n ) );
		if ( is_acceptor && hbond.acc_res() == i && pose.residue_type(i).atom_name( hbond.acc_atm() ) == atom_name )			return true;
		if ( !is_acceptor && hbond.don_res() == i && pose.residue_type(i).atom_name( hbond.don_hatm() ) == atom_name )			return true;
	}
	return false;

}

///////////////////////////////////////////////////////////////////////////////
Size
rna_features_from_pose( utility::io::ozstream & out, pose::Pose & pose )
{

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::id;
	using namespace protocols::farna;
	using namespace chemical::rna;

	vector1< std::string > feature_names;
	vector1< char > chains;
	vector1< char > seqchars;
	vector1< int > resnums;
	vector1< vector1< Real > > all_feature_vals;

	Size res_count( 0 ), num_features( 0 );
	Size const nres = pose.total_residue();
	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
	vector1<bool> wc_edge_paired, hoog_edge_paired, sugar_edge_paired, is_bulged;
	get_rna_base_pairing_status( pose, wc_edge_paired, hoog_edge_paired, sugar_edge_paired, is_bulged );

	//	(*scorefxn)(pose); // may need to create a dummy scorefxn with hbond_sc and score the pose...
	hbonds::HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
	hbond_options->use_hb_env_dep( false );
	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( hbond_options ) );
	hbonds::fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

	// sasa calculation
	AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;
	Real const probe_radius( 1.4 );
	scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, true );

	vector1< char  > nt_names = utility::tools::make_vector1( 'a', 'c', 'g', 'u' );

	Size const num_nt = nt_names.size();
	runtime_assert( num_nt == 4 );

	for (Size i = 1; i <= nres; i++) {

		Residue const & rsd = pose.residue( i );
		if ( !rsd.is_RNA() ) continue;
		res_count++;

		seqchars.push_back( rsd.name1() );
		resnums.push_back( pdb_info->number( i ) );
		chains.push_back( pdb_info->chain( i ) );

		vector1< Real > feature_vals;
		Size feature_counter( 0 );

		// is_a, is_c, etc.
		vector1< bool > is_nt;
		vector1< std::string > is_nt_tag;
		for ( Size m = 1; m <= num_nt; m++ ){
			is_nt.push_back( rsd.name1() == nt_names[m] );
			is_nt_tag.push_back( "is_" + std::string(&nt_names[m],1) /* man, converting char to string is a pain!*/ );
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m], is_nt[m] );
 		}

		// probably makes sense to explicitly take product with is_a, etc. above -- in case classifier is not
		// smart about leveraging products of features.
		for ( Size m = 1; m <= num_nt; m++ ){
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_wc_edge_paired"   , wc_edge_paired[i] * is_nt[m] );
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_hoog_edge_paired" , hoog_edge_paired[i] * is_nt[m]);
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_sugar_edge_paired", sugar_edge_paired[i] * is_nt[m]);
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_is_bulged", is_bulged[i] * is_nt[m] );
		}


		// What is hydrogen bonded?
		for ( Size m = 1; m <= num_nt; m++ ){

			// need an instance of this nucleotide to play around with.
			// This better have RNA residue types in it.
			ResidueTypeSet const & rsd_set = rsd.residue_type_set();
			ResidueTypeCOPs const & rsd_types = rsd_set.aa_map( chemical::aa_from_oneletter_code( nt_names[m] ) );
			ResidueType const & rsd_type = *rsd_types[1];

			bool ade_n1_hbonded = is_nt[1] && check_hbonded( pose, i, " N1 ", hbond_set, true /*acceptor*/ );

			// what is hydrogen bonded? Acceptors.
			AtomIndices accpt_pos = rsd_type.accpt_pos();
			for ( Size k = 1; k <= accpt_pos.size(); k++ ){
				std::string atom_name = rsd_type.atom_name( accpt_pos[ k ] );
				bool hbonded = is_nt[m] && check_hbonded( pose, i, atom_name, hbond_set, true /*acceptor*/ );
				ObjexxFCL::strip_whitespace(atom_name);
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded", hbonded);

				// special for DMS
				if (m==1) save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded"+"_and_N1_hbonded", hbonded && ade_n1_hbonded);
			}

			// what is hydrogen bonded? Donor hydrogens.
			AtomIndices Hpos_polar = rsd_type.Hpos_polar();
			for ( Size k = 1; k <= Hpos_polar.size(); k++ ){
				std::string atom_name = rsd_type.atom_name( Hpos_polar[ k ] );
				bool hbonded = is_nt[m] && check_hbonded( pose, i, atom_name, hbond_set, false /*acceptor --> check for polar hydrogen names*/ );
				ObjexxFCL::strip_whitespace(atom_name);
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded", hbonded);

				// special for DMS
				if (m==1) save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded"+"_and_N1_hbonded", hbonded && ade_n1_hbonded);
			}
		}


		//sasa
		for ( Size m = 1; m <= num_nt; m++ ){

			// need an instance of this nucleotide to play around with.
			// This better have RNA residue types in it.
			ResidueTypeSet const & rsd_set = rsd.residue_type_set();
			ResidueTypeCOPs const & rsd_types = rsd_set.aa_map( chemical::aa_from_oneletter_code( nt_names[m] ) );
			ResidueType const & rsd_type = *rsd_types[1];
			for ( Size k = 1; k <= rsd_type.nheavyatoms(); k++ ){
				std::string atom_name = rsd_type.atom_name( k );
				Real sasa_value = 0.0;
				if ( rsd.type().has( atom_name ) ){
					AtomID atom_id = pose::named_atom_id_to_atom_id( NamedAtomID( atom_name, i ), pose );
					sasa_value = atom_sasa[ atom_id ];
				}
				ObjexxFCL::strip_whitespace(atom_name);
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_SASA", sasa_value);
			}
		}

		if ( num_features == 0 ) num_features = feature_counter; // initialize num_features.
		runtime_assert( feature_vals.size() == num_features ); // sanity check.

		all_feature_vals.push_back( feature_vals );
	}

	runtime_assert( all_feature_vals.size() == res_count );

	/////////////////////////////////////////////////////////////////////////////////////////////
	// following should probably go in its own block of code -- perhaps used by an RDAT class.
	/////////////////////////////////////////////////////////////////////////////////////////////
	// output features to RDAT file.
	for ( Size n = 1; n <= num_features; n++ ){
		out << "ANNOTATION_DATA:" << n << '\t'<< "feature:" << feature_names[n] << std::endl;
	}
	out << std::endl;

	// SEQPOS line
	out << "SEQPOS";
	for ( Size n = 1; n <= res_count; n++ ){
		out << "\t" << seqchars[n] << resnums[n];
	}
	out << std::endl;

	// REACTIVITY line [or should this be called feature_value?]
	for ( Size k = 1; k <= num_features; k++ ){

		out << "REACTIVITY:"<< k;
		for ( Size n = 1; n <= res_count; n++ ){
			runtime_assert( all_feature_vals[ n ].size() == num_features );
			out << '\t' << all_feature_vals[ n ][ k ];
		}
		out << std::endl;

	}

	return res_count;

}

///////////////////////////////////////////////
// Perhaps should create an RDAT class inside
// Rosetta that can handle all this, including
// readin of data!
///////////////////////////////////////////////
void
create_rdat_header( utility::io::ozstream & out, pose::Pose & pose ){
	std::string const rdat_version_num_string = "0.33";
	out << "RDAT_VERSION\t" << rdat_version_num_string << std::endl;

	std::string pdb_name = pose.pdb_info()->name();
	out << "NAME\t" << pdb_name << std::endl;
	out << "SEQUENCE\t" << pose.sequence() << std::endl;

	// later can replace this with actual structure...
	std::string structure;
	for ( Size n = 1; n <= pose.sequence().size(); n++ ) structure += '.';
	out << "STRUCTURE\t" << structure << std::endl;
	out << std::endl;
	out << "COMMENT\tGenerated by rna_features in Rosetta." << std::endl;
	out << std::endl;
	out << "ANNOTATION\tsequenceSource:PDB:" << pdb_name << std::endl;
	out << std::endl;

	// offset -- could be an issue -- need to devise a better solution, like ALL_SEQPOS or something, which is an alternative.
	int offset = pose.pdb_info()->number(1) - 1;
	out << "OFFSET\t" << offset << std::endl;

}

///////////////////////////////////////////////
void
rhiju_pdbstats()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::import_pose;

	vector1 < std::string> pdb_files( option[ in::file::s ]() );
	std::string const file_path( option[ in::path::pdb ]( 1 ) );

	if ( option[ in::file::l ].user() ) {

		std::string const pdb_list(  option[ in::file::l ](1) );
		utility::io::izstream instream( pdb_list );
		if (!instream){
			std::cerr  << "Can't find list file " << pdb_list << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			return;
		}

		pdb_files.clear();
		std::string pdb_file, line;
		while ( 	getline( instream, line )  ) {
			std::istringstream line_stream( line );
			line_stream >> pdb_file;
			pdb_files.push_back( pdb_file );
		}

	}

	// later hope to replace this with fa_standard, which should soon include RNA, DNA, & protein.
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	Size count( 0 );
	std::string outfile;
	pose::Pose pose;
	Size total_residues( 0 );

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {

		std::string const pdb_file = pdb_files[n];

		std::string pdb_file_load = pdb_file;
		if ( file_path  != "./" ) pdb_file_load = file_path + '/' + pdb_file;
		pose_from_pdb( pose, *rsd_set, pdb_file_load );

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << std::endl;

		if ( option[out::file::o].user() ) outfile  = option[ out::file::o ];
		else outfile = pdb_file + ".rdat";

		utility::io::ozstream out( outfile );

		create_rdat_header( out, pose );

		total_residues += rna_features_from_pose( out, pose );

		std::cout << "Creating output RDAT file: " << outfile << std::endl;
		out.close();

	}

	std::cout << "FINISHED ==> TOTAL RESIDUES Processed: " << total_residues << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {
		using namespace basic::options;

		core::init::init(argc, argv);
		rhiju_pdbstats();

		exit( 0 );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
