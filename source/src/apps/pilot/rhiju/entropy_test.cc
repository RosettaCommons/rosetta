// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/scoring/dunbrack/DunbrackRotamer.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/viewer/viewers.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

//StepWiseProtein!
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/StepWiseProteinPoseSetup.hh>
#include <protocols/stepwise/StepWisePoseSampleGenerator.hh>
#include <protocols/stepwise/StepWiseDoNothingSampleGenerator.hh>
#include <protocols/stepwise/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/protein/StepWiseProteinMainChainSampleGenerator.hh>
#include <protocols/stepwise/protein/StepWiseProteinPoseSetup.hh>
#include <protocols/stepwise/protein/StepWiseProteinScreener.hh>
#include <protocols/stepwise/protein/util.hh>
#include <protocols/stepwise/protein/StepWiseProteinConnectionSampler.hh>
#include <protocols/stepwise/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/protein/MainChainTorsionSet.hh>

//clustering
#include <protocols/cluster/cluster.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/optimizeH.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>

#include <core/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/pose_stream/PoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>
#include <utility/tools/make_vector1.hh>


#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

//Job dsitributor
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <vector>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/cluster.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

typedef  numeric::xyzMatrix< Real > Matrix;
//typedef std::map< std::string, core::pose::PoseOP > PoseList;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, each_aa )
OPT_KEY( Boolean, connectivity )
OPT_KEY( Boolean, dist_matrix )
OPT_KEY( StringVector, silent_file_parent_build_Cterm )
OPT_KEY( StringVector, silent_file_parent_build_Nterm )


/////////////////////////////////////////////////////////////
Real
sidechain_sample( pose::Pose & pose,
									utility::vector1< Real > const & temperatures,
									utility::vector1< Real > & deltaG_diff_tot,
									Real & min_delta_score_tot,
									bool ignore_current_sidechain
									)
{

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::dunbrack;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace ObjexxFCL::format;

	static core::scoring::dunbrack::RotamerLibrary rotamer_library_(core::scoring::ScoringManager::get_instance()->get_RotamerLibrary());

	Pose pose_start = pose;
	static const ScoreFunctionOP scorefxn = get_score_function();

	Real const start_score = (*scorefxn)( pose );
	min_delta_score_tot = 0.0;

	deltaG_diff_tot.clear();
	for ( Size i = 1; i <= temperatures.size(); i++ ) deltaG_diff_tot.push_back( 0.0 );

	for ( Size n = 1; n <= pose.total_residue(); n++ ) {

		pose = pose_start;

		// Figure out how many alternative rotamers there might be.
		core::chemical::ResidueTypeCOP residue_type( pose.residue( n ).type() );

		SingleResidueRotamerLibraryCAP residue_rotamer_library	( rotamer_library_.get_rsd_library(*residue_type) );
		SingleResidueDunbrackLibraryCAP	residue_dunbrack_library( dynamic_cast< core::scoring::dunbrack::SingleResidueDunbrackLibrary const * >(residue_rotamer_library.get())	);

		if ( !residue_dunbrack_library){
			//std::cout << " no rotamers for position " << n << std::endl;
			continue;
		}

		utility::vector1< DunbrackRotamerSampleData> rotamer_samples = residue_dunbrack_library->get_all_rotamer_samples(pose.phi(n), pose.psi(n) );

		Size const nrots( rotamer_samples.size() );
		//		std::cout << "FOUND " << nrots << " rotamers for position  " << n << " with amino acid " << residue_type->name1() << std::endl;

		// Cycle through rotamers -- gobbledygook from protocols/SidechainMover.cc
		Size closest_rot( 0 );
		Real closest_chi_dist( 999999999999999999999.99999);
		utility::vector1< Real > delta_scores;

		Real min_delta_score = 0.0;

		for (Size rotnum = 1; rotnum <= nrots ; rotnum++) {

			DunbrackRotamerSampleData const & rotamer_sample_data = rotamer_samples[ rotnum ];

			Real4 const & chi_means( rotamer_sample_data.chi_mean() );

			// closest rotamer check?
			Real chi_dist( 0.0 );

			for ( Size chinum = 1; chinum <= rotamer_sample_data.nchi(); ++chinum) {
				Real const chi_val = chi_means[ chinum ];
				//std::cout << " " << F(8,3,chi_val);
				pose.set_chi(chinum, n, chi_val );

				chi_dist += numeric::principal_angle_degrees( chi_val - pose_start.chi( chinum, n ) );

			}

			if ( chi_dist < closest_chi_dist ) {
				closest_rot = rotnum;
				closest_chi_dist = chi_dist;
			}

			Real const score = (*scorefxn)( pose );
			//std::cout << " --> " << F(8,3,score) << std::endl;
			Real const delta_score = score - start_score;
			delta_scores.push_back( delta_score );

			if ( delta_score < min_delta_score ) min_delta_score = delta_score;
		}

		// Replace closest rotamer score with current score.
		if ( !ignore_current_sidechain ) delta_scores[ closest_rot ] = 0.0;

		for ( Size i = 1; i <= temperatures.size(); i++ ) {
			Real const temperature = temperatures[ i ];
			Real Z( 0.0 );
			for (Size rotnum = 1; rotnum <= nrots ; rotnum++) {
				Z += exp( -1.0 * ( delta_scores[ rotnum ] - min_delta_score ) / temperature );
			}
			Real const delG_diff = - temperature * log( Z ) + min_delta_score;
			deltaG_diff_tot[ i ] += delG_diff;
		}

		min_delta_score_tot += min_delta_score;

	}

	pose = pose_start;

	return ( start_score );

}

////////////////////////////////////////////////////////////////
void
output_stuff( Real const start_score, Real const min_delta_score_tot,
							utility::vector1< Real > const & temperatures,
							utility::vector1< Real > const & deltaG_diff_tot,
							std::ostream & out )
{
		out << start_score << "   " << min_delta_score_tot << "   ";
		for ( Size i = 1; i <= temperatures.size(); i++ ) out << ' ' << deltaG_diff_tot[ i ];
		out << std::endl;
}
////////////////////////////////////////////////////////////////
void
entropy_calculate_test() {

	using namespace core::chemical;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::io::pose_stream;
	using namespace protocols::swa;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	PoseInputStreamOP input;

	if ( option[ in::file::s].user() ) {
		input = new PDBPoseInputStream( option[ in::file::s ]() );
	} else if ( option[ in::file::silent].user() ) {
		input = new SilentFilePoseInputStream( option[ in::file::silent ](), true /*order by energy*/  );
	}

	std::string const outfile = option[ out::file::o ];
	utility::io::ozstream out( outfile );

	utility::vector1< Real > temperatures;
	for ( Size i = 1; i <= 40; i++ ) temperatures.push_back( i * 0.05 );

	Pose pose;
	while( input->has_another_pose() ) {

		///////////////////////////////////////////////////////////////
		input->fill_pose( pose, *rsd_set );

		utility::vector1< Real > deltaG_diff_tot;
		Real min_delta_score_tot( 0.0 );
		Real const start_score = sidechain_sample( pose, temperatures, deltaG_diff_tot, min_delta_score_tot, false /*ignore_current_sidechain*/ );

		output_stuff( start_score, min_delta_score_tot, temperatures, deltaG_diff_tot, std::cout );
		output_stuff( start_score, min_delta_score_tot, temperatures, deltaG_diff_tot, out );

	}
	out.close();

}

///////////////////////////////////////////////////////////////
void
each_aa_test(){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::swa;
	using namespace protocols::stepwise::protein;

	utility::vector1< Size > blank_size_vector;
	Size cutpoint_closed = 0;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	std::string aas = "ACDEFGHIKLMNPQRSTVWY";
	//	std::string aas = "R";

	utility::vector1< Real > temperatures;
	for ( Size i = 1; i <= 40; i++ ) temperatures.push_back( i * 0.05 );

	std::string const outfile = "each_aa.txt";
	utility::io::ozstream out( outfile );

	for ( Size n = 1; n <= aas.size(); n++ ) {

		///////////////////////////////////////
		// setup pose
		Pose pose;
		std::string const seq = aas.substr(n-1,1);

		// old style
    //ResidueOP rsd1( ResidueFactory::create_residue( *(rsd_set->aa_map( aa_from_oneletter_code( aas[n-1] ) )[1] ) ) );
    //pose.append_residue_by_bond( *rsd1, true /*build_ideal_geometry*/ );
		//chemical::add_variant_type_to_pose_residue( pose, "N_ACETYLATION", 1 );
		//chemical::add_variant_type_to_pose_residue( pose, "C_METHYLAMIDATION", 1 );

		// new style
		utility::vector1< InputStreamWithResidueInfoOP > input_streams; // empty
		utility::vector1< core::Size > const moving_res_list = utility::tools::make_vector1( 1 );
		StepWiseProteinPoseSetupOP stepwise_pose_setup = new StepWiseProteinPoseSetup( moving_res_list, seq,
																																		 input_streams,
																																		 blank_size_vector, blank_size_vector /*cutpoint_closed*/ );
		stepwise_pose_setup->set_add_peptide_plane_variants( true );
		stepwise_pose_setup->set_remove_nterminus_variant( true );
		stepwise_pose_setup->set_remove_cterminus_variant( true );

		stepwise_pose_setup->apply( pose );
		pose.dump_pdb( seq + ".pdb" );

		/////////////////////////////////////////////
		//get ready to sample backbone residues.
		working_parameters::StepWiseWorkingParametersOP & working_parameters = stepwise_pose_setup->working_parameters();
		StepWiseProteinScreener stepwise_screener( working_parameters );
		stepwise_screener.apply( pose );
		StepWisePoseSampleGeneratorOP sample_generator;
		sample_generator = new StepWiseProteinMainChainSampleGenerator( stepwise_screener.which_torsions(),
																																		stepwise_screener.main_chain_torsion_set_lists_real() );
		std::cout << "Using StepWiseProteinMainChainSampleGenerator. Num poses: " << stepwise_screener.main_chain_torsion_set_lists_real().size() << std::endl;


		////////////////////////////////////////////////////////////
		// Not necessary, but allows "packing" of one residue
		StepWiseProteinPacker stepwise_packer( moving_res_list, new StepWiseDoNothingSampleGenerator() );
		stepwise_packer.set_scorefxn( get_score_function() );

		// keep track of free energies vs. temperature.
		utility::vector1< utility::vector1<Real >  > all_deltaG_diff_tot;
		utility::vector1< Real > all_start_score,	all_min_delta_score_tot, all_min_score;


		// sample backbone residues, and sidechain.
		sample_generator->reset();
		while( sample_generator->has_another_sample() ){

			sample_generator->get_next_sample( pose );

			// might as well pack it.
			stepwise_packer.apply( pose );

			utility::vector1< Real > deltaG_diff_tot;
			Real min_delta_score_tot( 0.0 );
			Real const start_score = sidechain_sample( pose, temperatures, deltaG_diff_tot, min_delta_score_tot, true /*ignore_current_sidechain*/ );

			output_stuff( start_score, min_delta_score_tot, temperatures, deltaG_diff_tot, std::cout );

			all_deltaG_diff_tot.push_back( deltaG_diff_tot );
			all_start_score.push_back( start_score );
			all_min_delta_score_tot.push_back( min_delta_score_tot );
			all_min_score.push_back( min_delta_score_tot + start_score );

		}

		// First need to figure out min score.
		Real min_score_overall = all_min_score[1];
		for ( Size i = 1; i <= all_min_score.size(); i++ ){	 if ( all_min_score[i] < min_score_overall ) min_score_overall = all_min_score[i]; }

		utility::vector1< Real > deltaG_diff_overall;

		for ( Size j = 1; j <= temperatures.size(); j++ ) {
			Real Z_across_backbone_rotamers = 0.0;

			for ( Size i = 1; i <= all_min_score.size(); i++ ){
				Real const deltaG_diff = all_deltaG_diff_tot[i][j] + all_start_score[i] - min_score_overall;
				Z_across_backbone_rotamers += exp( - deltaG_diff / temperatures[j] );
			}

			deltaG_diff_overall.push_back(  -temperatures[j] * log( Z_across_backbone_rotamers ) );

		}

		std::cout << "OVERALL: " << std::endl;
		output_stuff( min_score_overall, 0.0, temperatures, deltaG_diff_overall, std::cout );
		output_stuff( min_score_overall, 0.0, temperatures, deltaG_diff_overall, out );

	}

	out.close();
	std::cout << std::endl;
	std::cout << "Output delG for each amino acid in " << outfile << std::endl;


}

///////////////////////////////////////////////////////////////
void
output_connect( std::ostream & out,
								Size const & count_child,
								Real const & best_rmsd,
								Size const & count_parent,
								Real const build_direction ){

	out << count_child << "      " << build_direction << " " << count_parent << "   " << best_rmsd << std::endl;

}

///////////////////////////////////////////////////////////////
void
setup_atom_id_map( std::map< id::AtomID, id::AtomID > & atom_id_map,
									 pose::Pose const & pose_child,
									 pose::Pose const & pose_parent,
									 int const build_direction ){

	using namespace core::id;

	utility::vector1< std::string > atom_names;
	atom_names.push_back( " N  " );
	atom_names.push_back( " C  " );
	atom_names.push_back( " CA " );
	atom_names.push_back( " O  " );

	Size parent_count = 0;
	for ( Size i = 1; i <= pose_child.total_residue(); i++ ){
		if ( build_direction == +1 && i == pose_child.total_residue() ) continue;
		if ( build_direction == -1 && i == 1                          ) continue;
		parent_count++;
		if ( pose_child.sequence()[i-1] != pose_parent.sequence()[parent_count-1] ) utility_exit_with_message( "Mismatch between parent and child." );

		for ( Size n = 1; n <= atom_names.size(); n++ ){
			atom_id_map[ AtomID( pose_child.residue_type(i).atom_index( atom_names[ n ] ), i ) ] =
									 AtomID( pose_parent.residue_type(parent_count).atom_index( atom_names[ n ] ), parent_count );
		}
	}

}

///////////////////////////////////////////////////////////////
void
find_close_parents( pose::Pose const & pose_child,
										Size const count_child,
										utility::vector1< pose::PoseOP > & pose_parent_list,
										int const build_direction,
										std::ostream & out,
										Size & total_connect,
										Real & best_rmsd, Size & best_count, int & best_parent )
{
	using namespace core::pose;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	Real rmsd_cutoff( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) rmsd_cutoff = option[ OptionKeys::cluster::radius ]();

	std::map< AtomID, AtomID > atom_id_map;
	setup_atom_id_map( atom_id_map, pose_child, *(pose_parent_list[1]), build_direction );

	for ( Size n = 1; n <= pose_parent_list.size(); n++ ){

		Real const rmsd = rms_at_corresponding_atoms( pose_child, *(pose_parent_list[n]), atom_id_map );
		if ( rmsd < rmsd_cutoff ) {
			output_connect( out, count_child, rmsd, n, build_direction );
			total_connect++;
		}

		if ( rmsd < best_rmsd ) {
			best_rmsd = rmsd;
			best_count = n;
			best_parent = build_direction;
		}

	}

}

///////////////////////////////////////////////////////////////
void
fill_pose_list( 	utility::vector1< pose::PoseOP > & pose_list,
									core::io::pose_stream::PoseInputStreamOP input_parent,
									core::chemical::ResidueTypeSetCAP rsd_set ){

	using namespace core::pose;

	input_parent->reset();

	while( input_parent->has_another_pose() ) {
		PoseOP pose_op = new Pose;
		input_parent->fill_pose( *pose_op, *rsd_set );
		pose_list.push_back( pose_op );
	}

}


///////////////////////////////////////////////////////////////
void
connectivity_test(){

	using namespace core::chemical;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::io::pose_stream;
	using namespace protocols::swa;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	utility::vector1< std::string > const & silent_files = option[ in::file::silent ]();

	PoseInputStreamOP input_child, input_parent_build_Cterm, input_parent_build_Nterm;
	input_child = new SilentFilePoseInputStream( silent_files, true );

	if ( option[ silent_file_parent_build_Cterm ].user() ){
		input_parent_build_Cterm = new SilentFilePoseInputStream( option[ silent_file_parent_build_Cterm ](), true );
	}
	if ( option[ silent_file_parent_build_Nterm ].user() ){
		input_parent_build_Nterm = new SilentFilePoseInputStream( option[ silent_file_parent_build_Nterm ](), true );
	}

	Pose pose_child, pose_parent;

	input_child->reset();

	std::string const outfile = option[ out::file::o ];
	utility::io::ozstream out( outfile );
	Size count_child( 0 );
	utility::vector1< PoseOP > pose_list_build_Cterm, pose_list_build_Nterm;
	fill_pose_list( pose_list_build_Cterm, input_parent_build_Cterm, rsd_set );
	fill_pose_list( pose_list_build_Nterm, input_parent_build_Nterm, rsd_set );

	while( input_child->has_another_pose() ){

		count_child++;
		input_child->fill_pose( pose_child, *rsd_set );

		Real best_rmsd( 99999.999999 );
		Size best_count( 0 ), total_connect( 0 );
		int best_parent( 0 );


		if ( pose_list_build_Cterm.size() > 0 ) {
			find_close_parents( pose_child, count_child,
													pose_list_build_Cterm,
													+1, out, total_connect,
													best_rmsd, best_count, best_parent );
		}

		if ( pose_list_build_Nterm.size() > 0 ) {
			find_close_parents( pose_child, count_child,
													pose_list_build_Nterm,
													-1, out, total_connect,
													best_rmsd, best_count, best_parent );
		}

		if ( total_connect == 0 )  output_connect( out, count_child, best_rmsd, best_count, best_parent );
		std::cout << "Finished child pose " << count_child << " --> best_rmsd " << best_rmsd << std::endl;

	}

	out.close();

}


///////////////////////////////////////////////////////////////
Vector
get_xyz( pose::Pose const & pose, Size const i ){

	if ( pose.residue_type( i ).has( " CB " ) )  return pose.residue( i ) .xyz( " CB " );
	return pose.residue( i ) .xyz( " CA " );

}

///////////////////////////////////////////////////////////////
void
output_distance_matrix( pose::Pose const & pose, std::ostream & out ){

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		Vector const v_i = get_xyz( pose, i );

		for ( Size j = 1; j <= pose.total_residue(); j++ ) {
			Vector const v_j = get_xyz( pose, j );
			out << ' ' << ( v_i - v_j ).length();
		}

	}

	out << std::endl;

}

///////////////////////////////////////////////////////////////
void
dist_matrix_test(){

	using namespace core::chemical;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::io::pose_stream;
	using namespace protocols::swa;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	utility::vector1< std::string > const & silent_files = option[ in::file::silent ]();

	PoseInputStreamOP input = new SilentFilePoseInputStream( silent_files, true );

	std::string const outfile = option[ out::file::o ];
	utility::io::ozstream out( outfile );

	Pose pose;

	while( input->has_another_pose() ) {
		input->fill_pose( pose, *rsd_set );
		output_distance_matrix( pose, out );
	}

	out.close();


}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ each_aa ]() ){
		each_aa_test();
	} else if ( option[ connectivity ]() ){
		connectivity_test();
	} else if ( option[ dist_matrix ]() ){
		dist_matrix_test();
	} else {
		entropy_calculate_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


	using namespace core::options;

	utility::vector1< std::string > blank_string_vector;

	NEW_OPT( each_aa, "Scan through each AA", false );
	NEW_OPT( connectivity, "Figure out connectivity between models", false );
	NEW_OPT( dist_matrix, "Output distance matrix for pose", false );
	NEW_OPT( silent_file_parent_build_Cterm, "silent file for parent", blank_string_vector );
	NEW_OPT( silent_file_parent_build_Nterm, "silent file for parent", blank_string_vector );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

	exit( 0 );

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
