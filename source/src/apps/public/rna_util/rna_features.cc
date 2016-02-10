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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.hh>
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.hh>
#include <core/io/rna/RDAT.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>

#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <core/import_pose/import_pose.hh>

#include <devel/init.hh>

#include <protocols/farna/util.hh>
#include <protocols/viewer/viewers.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/constants.hh>

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
#include <basic/options/keys/chemical.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <numeric/xyz.functions.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

using namespace core;
using namespace core::conformation;
//using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace core::chemical::rna;

using utility::vector1;
using utility::tools::make_vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;


typedef  numeric::xyzMatrix< core::Real > Matrix;


///////////////////////////////////////////////////////////////////////////
void
save_feature(  vector1< std::string > & feature_names,
	vector1< core::Real > & feature_vals,
	Size & feature_counter,
	std::string const & feature_name,
	core::Real const feature_val ){

	feature_counter++;

	if ( feature_counter > feature_names.size() ) feature_names.push_back( feature_name );
	runtime_assert( feature_counter <= feature_names.size() );

	feature_vals.push_back( feature_val );
}

///////////////////////////////////////////////////////////////////////////////
Size
rna_features_from_pose( core::io::rna::RDAT & rdat, pose::Pose & pose )
{
	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna::data;
	using namespace core::kinematics;
	using namespace core::id;
	using namespace protocols::farna;
	using namespace core::chemical::rna;

	vector1< std::string > feature_names;
	vector1< char > chains;
	vector1< char > seqchars;
	vector1< int > resnums;
	vector1< vector1< core::Real > > all_feature_vals;

	Size res_count( 0 ), num_features( 0 );
	Size const nres = pose.total_residue();
	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();

	// sasa calculation
	AtomID_Map< core::Real > atom_sasa;
	utility::vector1< core::Real > rsd_sasa;
	core::Real const probe_radius( 1.4 );
	scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, true );

	RNA_DMS_Potential & rna_dms_potential  = ScoringManager::get_instance()->get_RNA_DMS_Potential();
	pose.update_residue_neighbors();
	rna_dms_potential.initialize( pose );

	RNA_DMS_LowResolutionPotential & rna_dms_low_resolution_potential = ScoringManager::get_instance()->get_RNA_DMS_LowResolutionPotential();
	// rna_dms_low_resolution_potential.set_careful_base_pair_classifier( false );
	rna_dms_low_resolution_potential.initialize( pose );

	vector1< char  > nt_names = utility::tools::make_vector1( 'a', 'c', 'g', 'u' );
	Size const num_nt = nt_names.size();
	runtime_assert( num_nt == 4 );

	for ( Size i = 1; i <= nres; i++ ) {

		Residue const & rsd = pose.residue( i );
		if ( !rsd.is_RNA() ) continue;
		res_count++;

		seqchars.push_back( rsd.name1() );
		resnums.push_back( pdb_info->number( i ) );
		chains.push_back( pdb_info->chain( i ) );

		vector1< core::Real > feature_vals;
		Size feature_counter( 0 );

		// is_a, is_c, etc.
		vector1< bool > is_nt;
		vector1< std::string > is_nt_tag;
		for ( Size m = 1; m <= num_nt; m++ ) {
			is_nt.push_back( rsd.name1() == nt_names[m] );
			is_nt_tag.push_back( "is_" + std::string(&nt_names[m],1) /* man, converting char to string is a pain!*/ );
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m], is_nt[m] );
		}

		// look for O2' in WC sector -- probably explains rest of DMS protections. Note that this is general, so
		// perhaps should not be in a dms_low_resolution_potential, but instead in an rna_base_pair_info object?
		bool wc_near_o2prime = rna_dms_low_resolution_potential.get_wc_near_o2prime( pose, i );

		// probably makes sense to explicitly take product with is_a, etc. above -- in case classifier is not
		// smart about leveraging products of features.
		for ( Size m = 1; m <= num_nt; m++ ) {
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_wc_edge_paired"   ,
				rna_dms_low_resolution_potential.wc_edge_paired()[i]    * is_nt[m] );
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_hoog_edge_paired" ,
				rna_dms_low_resolution_potential.hoog_edge_paired()[i]  * is_nt[m]);
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_sugar_edge_paired",
				rna_dms_low_resolution_potential.sugar_edge_paired()[i] * is_nt[m]);
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_is_bulged",
				rna_dms_low_resolution_potential.is_bulged()[i] * is_nt[m] );
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_wc_near_o2prime"   , wc_near_o2prime * is_nt[m] );
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_wc_edge_paired_or_near_o2prime",
				( wc_near_o2prime || rna_dms_low_resolution_potential.wc_edge_paired()[i] ) * is_nt[m] );
		}


		// What is hydrogen bonded?
		bool ade_n1_hbonded = is_nt[1] && rna_dms_potential.check_hbonded( pose, i, " N1 ", true /*acceptor*/ );

		for ( Size m = 1; m <= num_nt; m++ ) {

			// need an instance of this nucleotide to play around with.
			// This better have RNA residue types in it.
			ResidueTypeSet const & rsd_set = *rsd.residue_type_set();
			ResidueType const & rsd_type = *( rsd_set.get_representative_type_aa(  core::chemical::aa_from_oneletter_code( nt_names[m] ) ) );

			// what is hydrogen bonded? Acceptors.
			AtomIndices accpt_pos = rsd_type.accpt_pos();
			for ( Size k = 1; k <= accpt_pos.size(); k++ ) {
				std::string atom_name = rsd_type.atom_name( accpt_pos[ k ] );
				bool hbonded = is_nt[m] && rna_dms_potential.check_hbonded( pose, i, atom_name, true /*acceptor*/ );
				ObjexxFCL::strip_whitespace(atom_name);
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded", hbonded);

				// special for DMS
				if ( m==1 ) save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded"+"_and_N1_hbonded", hbonded && ade_n1_hbonded);
			}

			if ( m == 1 ) {
				bool ade_n1_chbonded( false );
				if ( is_nt[1] ) ade_n1_chbonded = rna_dms_potential.check_chbonded( pose, i, " N1 " );
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+"_and_N1_chbonded", ade_n1_chbonded );
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+"_and_N1_hbonded_or_chbonded",
					ade_n1_chbonded || ade_n1_hbonded );
			}


			// what is hydrogen bonded? Donor hydrogens.
			AtomIndices Hpos_polar = rsd_type.Hpos_polar();
			for ( Size k = 1; k <= Hpos_polar.size(); k++ ) {
				std::string atom_name = rsd_type.atom_name( Hpos_polar[ k ] );
				bool hbonded = is_nt[m] && rna_dms_potential.check_hbonded( pose, i, atom_name, false /*acceptor --> check for polar hydrogen names*/ );
				ObjexxFCL::strip_whitespace(atom_name);
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded", hbonded);

				// special for DMS
				if ( m==1 ) save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_hbonded"+"_and_N1_hbonded", hbonded && ade_n1_hbonded);
			}
		}


		//sasa
		for ( Size m = 1; m <= num_nt; m++ ) {

			// need an instance of this nucleotide to play around with.
			// This better have RNA residue types in it.
			ResidueTypeSet const & rsd_set = *rsd.residue_type_set();
			ResidueType const & rsd_type = *( rsd_set.get_representative_type_aa( core::chemical::aa_from_oneletter_code( nt_names[m] ) ) );
			for ( Size k = 1; k <= rsd_type.nheavyatoms(); k++ ) {
				std::string atom_name = rsd_type.atom_name( k );
				core::Real sasa_value = 0.0;
				if ( rsd.type().has( atom_name ) ) {
					AtomID atom_id = pose::named_atom_id_to_atom_id( NamedAtomID( atom_name, i ), pose );
					sasa_value = atom_sasa[ atom_id ];
				}
				ObjexxFCL::strip_whitespace(atom_name);
				save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+ "_and_"+atom_name+"_SASA", sasa_value);
			}
		}

		// chi, delta
		for ( Size m = 1; m <= num_nt; m++ ) {
			core::Real chi = is_nt[m] ? pose.chi( i ) : 0.0;
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+"_and_chi", chi);
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+"_and_is_syn", (chi < 0.0) );
		}
		for ( Size m = 1; m <= num_nt; m++ ) {
			core::Real delta = is_nt[m] ? pose.delta( i ) : 0.0;
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+"_and_delta", delta);
			save_feature( feature_names, feature_vals, feature_counter, is_nt_tag[m]+"_and_is_south", ( delta > 120.0) );
		}

		// occlusion of pseudo-methyl at N1 -- towards fitting a potential.
		utility::vector1< Distance > probe_dists = make_vector1( 0.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.5 );
		utility::vector1< Distance > const shells = make_vector1( 0.0, 2.0, 4.0, 6.0, 8.0 );
		//utility::vector1< core::Real > const shells = make_vector1( 0.0, 3.0, 6.0, 9.0 );
		for ( Size m = 1; m <= probe_dists.size(); m++ ) {
			utility::vector1< core::Real > occupancy_densities( (shells.size()-1), 0.0 ); // atoms per A^3.
			core::Real const probe_dist = probe_dists[ m ];
			if ( is_nt[ 1 ] ) {
				core::Vector const probe_xyz = rna_dms_potential.get_probe_xyz( pose.residue( i ), probe_dist );
				rna_dms_potential.get_occupancy_densities( occupancy_densities, pose, i /*for exclusion*/, probe_xyz, shells );
			}
			for ( Size k = 1; k <= occupancy_densities.size(); k++ ) {
				std::string const tag =  is_nt_tag[1]+"_and_pseudomethyldist"+F(3,1,probe_dist)+"_shell_"+I( 1, shells[k] ) + "-" + I( 1, shells[k+1]);
				save_feature( feature_names, feature_vals, feature_counter, tag, occupancy_densities[k] );
			}
		}

		// actual computation of a potential around pseudo-methyl atom.
		vector1< ScoreFunctionOP > probe_scorefxns = make_vector1( rna_dms_potential.get_probe_scorefxn( false /* soft_rep*/, false /* just_fa_atr*/  ),
			rna_dms_potential.get_probe_scorefxn( true, false  ),
			rna_dms_potential.get_probe_scorefxn( true, true ) );
		vector1< std::string > scorefxntags = make_vector1( "hard", "soft", "justatrrep_soft" );
		for ( Size k = 1; k <= probe_scorefxns.size(); k++ ) {
			core::Real pseudo_methyl_plus_oxygen_energy( 0 );
			for ( Size m = 1; m <= probe_dists.size(); m++ ) {
				core::Real const probe_dist = probe_dists[ m ];
				core::Real binding_energy( 0.0 );
				if ( is_nt[ 1 ] ) {
					core::Vector const probe_xyz = rna_dms_potential.get_probe_xyz( pose.residue( i ), probe_dist );
					binding_energy = rna_dms_potential.get_binding_energy( i, probe_xyz, *probe_scorefxns[k] );
				}
				if ( probe_dist == 3.5 || probe_dist == 5.5 ) pseudo_methyl_plus_oxygen_energy += binding_energy;
				std::string const tag =  is_nt_tag[1]+"_and_pseudomethyldist"+F(3,1,probe_dist)+"_"+scorefxntags[k]+"_energy";
				if ( is_nt[ 1 ] ) std::cout << pose.pdb_info()->number( i ) << " " << tag << ": " << binding_energy << std::endl; // some output since this takes so long.
				save_feature( feature_names, feature_vals, feature_counter, tag, binding_energy );
			}

			if ( is_nt[ 1 ] ) std::cout << std::endl;
			std::string const tag =  is_nt_tag[1]+"_and_pseudomethyl_plus_oxygen_"+scorefxntags[k]+"_energy";
			save_feature( feature_names, feature_vals, feature_counter, tag, pseudo_methyl_plus_oxygen_energy );
		}

		if ( num_features == 0 ) num_features = feature_counter; // initialize num_features.
		runtime_assert( feature_vals.size() == num_features ); // sanity check.

		all_feature_vals.push_back( feature_vals );
	}

	runtime_assert( all_feature_vals.size() == res_count );

	rdat.fill_data_from_features( seqchars, resnums, feature_names, all_feature_vals );

	return res_count;
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
		if ( !instream ) {
			std::cerr  << "Can't find list file " << pdb_list << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			//return;
		}

		pdb_files.clear();
		std::string pdb_file, line;
		while (  getline( instream, line )  ) {
			std::istringstream line_stream( line );
			line_stream >> pdb_file;
			pdb_files.push_back( pdb_file );
		}

	}

	// later hope to replace this with fa_standard, which should soon include RNA, DNA, & protein.
	ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	Size count( 0 );
	std::string outfile;
	pose::Pose pose;
	Size total_residues( 0 );

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {

		std::string const pdb_file = pdb_files[n];

		std::string pdb_file_load = pdb_file;
		if ( file_path  != "./" ) pdb_file_load = file_path + '/' + pdb_file;
		pose_from_file( pose, *rsd_set, pdb_file_load , core::import_pose::PDB_file);

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << std::endl;
		std::cout << "Read in pose sequence: " << pose.annotated_sequence() << std::endl;
		pose.dump_pdb( "dump.pdb" );

		if ( option[out::file::o].user() ) outfile  = option[ out::file::o ];
		else outfile = pdb_file + ".rdat";

		io::rna::RDAT rdat;
		rdat.fill_header_information( pose );
		total_residues += rna_features_from_pose( rdat, pose );
		rdat.output_rdat_to_file( outfile );

	}

	std::cout << "FINISHED ==> TOTAL RESIDUES Processed: " << total_residues << std::endl;

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rhiju_pdbstats();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {
		using namespace basic::options;

		core::init::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
