// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Minimizer
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh> //Sept 26, 2011
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/screener/ChainClosableScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/simple_moves/ConstrainToIdealMover.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009
#include <core/conformation/Conformation.hh>
#include <protocols/rna/RNA_LoopCloser.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>

#include <time.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of RNA (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static numeric::random::RandomGenerator RG(3892128);
static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_Minimizer" ) ;

namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_Minimizer::StepWiseRNA_Minimizer(
																						 utility::vector1 < core::pose::PoseOP > const & pose_list,
																						 StepWiseRNA_JobParametersCOP & job_parameters ):
	pose_list_( pose_list ),
	job_parameters_( job_parameters ),
	silent_file_( "silent_file.txt" ),
	verbose_( false ),
	native_screen_( false ),
	native_screen_rmsd_cutoff_( 3.0 ),
	perform_electron_density_screen_( false ),
	rm_virt_phosphate_( false ),
	native_edensity_score_cutoff_(  - 1.0 ),
	centroid_screen_( true ),
	perform_o2prime_pack_( true ),
	output_before_o2prime_pack_( false ),
	perform_minimize_( true ),
	num_pose_minimize_( 999999 ), //Feb 02, 2012
	minimize_and_score_sugar_( true ),
	rename_tag_( true ),
	output_minimized_pose_list_( true ),
	base_centroid_screener_( 0 ), //Owning-pointer.
	vary_bond_geometry_frequency_( 0.0 ),
	allow_variable_bond_geometry_( false )
{
	set_native_pose( job_parameters_->working_native_pose() );
}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseRNA_Minimizer::~StepWiseRNA_Minimizer()
{}

/////////////////////
std::string
StepWiseRNA_Minimizer::get_name() const {
	return "StepWiseRNA_Minimizer";
}
			
//////////////////////////////////////////

void
StepWiseRNA_Minimizer::apply( core::pose::Pose & pose ) {

	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::optimization;
	using namespace protocols::rna;
	using namespace protocols::swa;

	//Real const bulge_weight = scorefxn_->get_weight( num_stacks );
	//scorefxn_->set_weight( num_stacks, 0.0);
	original_geometry_weight_ = scorefxn_->get_weight( rna_bond_geometry );
	output_title_text( "Enter StepWiseRNA_Minimizer::apply", TR.Debug );
	output_parameters();
	clock_t const time_start( clock() );

	Size const gap_size(  job_parameters_->gap_size() );
	bool const close_chainbreak = gap_size == 0 && ( job_parameters_->moving_res() > 0 );
	RNA_LoopCloser rna_loop_closer;

	SilentFileData silent_file_data;

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025 );
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	/////////////////////////////////////////////////////////
	if ( pose_list_.size() == 0 ){
		TR.Debug << "pose_list_.size() == 0, early exit from StepWiseRNA_Minimizer::apply" << std::endl;
		output_empty_minimizer_silent_file();
		return;
	}

	pose::Pose dummy_pose = ( *pose_list_[1] );
	Size const nres =  dummy_pose.total_residue();
	bool o2prime_pack_verbose( verbose_ );
	utility::vector1< Size > const working_moving_res = get_working_moving_res( nres ); // will be used for o2prime packing.

	//Check scorefxn
	TR.Debug << "check scorefxn" << std::endl;
	scorefxn_->show( TR.Debug, dummy_pose );
	
	bool vary_bond_geometry( false );
	if ( allow_variable_bond_geometry_ && ( RG.uniform() < vary_bond_geometry_frequency_ ) ) {
		vary_bond_geometry = true;
	}

	if ( move_map_list_.size() == 0 ) {
		move_map_list_ = get_default_movemap( dummy_pose );
	}

	if ( ! minimize_and_score_sugar_ ){
		TR.Debug << "WARNING: minimize_and_score_sugar_ is FALSE, freezing all DELTA, NU_1 and NU_2 torsions." << std::endl;
		for ( Size n = 1; n <= move_map_list_.size(); n++ ){
			Freeze_sugar_torsions( move_map_list_[n], nres );
		}
	}

	for ( Size n = 1; n <= move_map_list_.size(); n++ ){
		TR.Debug << "OUTPUT move_map_list[" << n << "]:" << std::endl;
		output_movemap( move_map_list_[n], dummy_pose, TR.Debug );
	}

	minimized_pose_list_.clear();

	for ( Size pose_ID = 1; pose_ID <= pose_list_.size(); pose_ID++ ){

		if ( minimized_pose_list_.size() >= num_pose_minimize_ ){
			TR.Debug << "WARNING MAX num_pose_minimize_( " << num_pose_minimize_ << " ) EXCEEDED, EARLY BREAK." << std::endl;
			break;
		}

		TR.Debug << "Minimizing pose_ID = " << pose_ID << std::endl;

		std::string tag = tag_from_pose( *pose_list_[pose_ID] );
		pose = ( *pose_list_[pose_ID] ); //This actually creates a hard copy.....need this to get the graphic working..

		if ( rm_virt_phosphate_ ){ //Fang's electron density code
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", ii );
			}
		}

		if ( verbose_ && output_before_o2prime_pack_ ) output_pose_wrapper( tag, 'B', pose, silent_file_data, silent_file_ + "_before_o2prime_pack" );

		remove_virtual_O2Prime_hydrogen( pose );

		if ( perform_o2prime_pack_ ) o2prime_trials( pose, scorefxn_, get_surrounding_O2prime_hydrogen( pose, working_moving_res, o2prime_pack_verbose ) );

		if ( verbose_ && !output_before_o2prime_pack_ ) output_pose_wrapper( tag, 'B', pose, silent_file_data, silent_file_ + "_before_minimize" );

		utility::vector1< Size > moving_chainbreaks = figure_out_moving_chain_break_res( pose, move_map_list_[ 1 ] );

		///////Minimization/////////////////////////////////////////////////////////////////////////////
		if ( close_chainbreak ){
			rna_loop_closer.apply( pose, moving_chainbreaks ); //This doesn't do anything if rna_loop_closer was already applied during sampling stage...May 10,2010
			if ( verbose_ ) output_pose_wrapper( tag, 'C', pose, silent_file_data, silent_file_ + "_after_loop_closure_before_minimize" );
		}

		for ( Size round = 1; round <= move_map_list_.size(); round++ ){
			core::kinematics::MoveMap mm = move_map_list_[round];

			moving_chainbreaks = figure_out_moving_chain_break_res( pose, mm );
			
			if ( vary_bond_geometry ) {
				TR << "Performing variable geometry minimization..." << std::endl;
				scorefxn_->set_weight( rna_bond_geometry, 8.0 );
				protocols::simple_moves::ConstrainToIdealMover CTIMover;
				core::kinematics::MoveMapOP mmop(mm.clone());
				CTIMover.set_movemap(mmop);
				CTIMover.set_AllowInsert(allow_insert_->clone());
				CTIMover.apply(pose);
				core::kinematics::MoveMap new_mm = (*CTIMover.get_movemap());
				//CartesianMinimizer cartesian_minimizer;
				if ( perform_minimize_ ) minimizer.run( pose, new_mm, *(scorefxn_ ), options );
				scorefxn_->set_weight( rna_bond_geometry, original_geometry_weight_ );
			}

			if ( perform_minimize_ ) minimizer.run( pose, mm, *( scorefxn_ ), options );
			if ( perform_o2prime_pack_ ) o2prime_trials( pose, scorefxn_, get_surrounding_O2prime_hydrogen( pose, working_moving_res, o2prime_pack_verbose ) );

			if ( close_chainbreak ){
				rna_loop_closer.apply( pose, moving_chainbreaks );
				if ( perform_o2prime_pack_ ) o2prime_trials( pose, scorefxn_, get_surrounding_O2prime_hydrogen( pose, working_moving_res, o2prime_pack_verbose ) );
				if ( perform_minimize_ ) minimizer.run( pose, mm, *( scorefxn_ ), options );

			}

		}

		////////////////////////////////Final screening //////////////////////////////////////////////
		if ( ! pass_all_pose_screens( pose, tag, silent_file_data ) ) continue;

		//////////////////////////////////////////////////////////////////////////////////////////////
		tag_into_pose( pose, tag );
		minimized_pose_list_.push_back( pose.clone() );

		TR.Debug << "SO FAR: pose_ID = " << pose_ID  << " | minimized_pose_list_.size() = " << minimized_pose_list_.size();
		TR.Debug << " | time taken = " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
		TR.Debug << "-------------------------------------------------------------------------------------------------" << std::endl;
		TR.Debug << "-------------------------------------------------------------------------------------------------" << std::endl;

		}

	TR.Debug << "FINAL minimized_pose_list_.size() = " << minimized_pose_list_.size() << std::endl;

	if ( minimized_pose_list_.size() == 0 ){
		TR.Debug << "After finish minimizing,  minimized_pose_list_.size() == 0!" << std::endl;
		output_empty_minimizer_silent_file();
	}

	std::sort( minimized_pose_list_.begin(), minimized_pose_list_.end(), sort_pose_by_score );

	if ( output_minimized_pose_list_ ) output_minimized_pose_list();

	// output best scoring pose.
	if ( minimized_pose_list_.size() > 0 ) pose = *( minimized_pose_list_[1] );

	TR.Debug << "--------------------------------------------------------------------" << std::endl;
	TR.Debug << "Total time in StepWiseRNA_Minimizer: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
	TR.Debug << "--------------------------------------------------------------------" << std::endl;

	output_title_text( "Exit StepWiseRNA_Minimizer::apply", TR.Debug );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// following is deprecated -- delete after testing figure_out_moving_chainbreak_res
core::Size
StepWiseRNA_Minimizer::figure_out_actual_five_prime_chain_break_res( pose::Pose const & pose ) const{

	using namespace ObjexxFCL;

	Size five_prime_chain_break_res = 0;

	if ( job_parameters_->is_simple_full_length_job_params() ){
		//March 03, 2012. If job_parameters_->is_simple_full_length_job_params()==true, then
		//actual five_prime_chain_break_res varies from pose to pose and might not correspond to the one given in the job_params!

		Size const gap_size(  job_parameters_->gap_size() );

		if ( gap_size != 0 ) utility_exit_with_message( "job_parameters_->is_simple_full_length_job_params() == true but gap_size != 0" );

		Size num_cutpoint_lower_found = 0;

		for ( Size lower_seq_num = 1; lower_seq_num < pose.total_residue(); lower_seq_num++ ){

			Size const upper_seq_num = lower_seq_num + 1;

			if ( pose.residue( lower_seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ){
				if ( pose.residue( upper_seq_num ).has_variant_type( chemical::CUTPOINT_UPPER ) == false ){
					std::string error_message = "";
					error_message += "seq_num " + string_of( lower_seq_num ) + " is a CUTPOINT_LOWER ";
					error_message += "but seq_num " + string_of( upper_seq_num ) + " is not a cutpoint CUTPOINT_UPPER??";
					utility_exit_with_message( error_message );
				}

				five_prime_chain_break_res = lower_seq_num;
				num_cutpoint_lower_found++;

			}
		}

		if ( num_cutpoint_lower_found != 1 ) utility_exit_with_message( "num_cutpoint_lower_found = ( " + string_of( num_cutpoint_lower_found ) + " ) != 1" );

	} else{ //Default
		five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
	}

	return five_prime_chain_break_res;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_pose_wrapper( std::string const & tag,
																								 pose::Pose & pose,
																								 std::string const & out_silent_file ) const {

	core::io::silent::SilentFileData silent_file_data;
	runtime_assert( tag.size() > 1 );
	std::string tag_local = tag;
	output_pose_wrapper( tag_local, tag[0], pose, silent_file_data, out_silent_file );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_pose_wrapper( std::string & tag,
																								 char tag_first_char,
																								 pose::Pose & pose,
																								 core::io::silent::SilentFileData & silent_file_data,
																								 std::string const out_silent_file ) const {

	using namespace core::io::silent;

	//tag_first_char=B for before minimize
	//tag_first_char=C for consistency_check
	//tag_first_char=M for minimize

	tag[0] = tag_first_char;
	TR.Debug << "tag = " << tag <<  std::endl;
	( *scorefxn_ )( pose ); //Score pose to ensure that that it has a score to be output
	output_data( silent_file_data, out_silent_file, tag, false, pose, get_native_pose(), job_parameters_ );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_empty_minimizer_silent_file() const {

	if ( utility::file::file_exists( silent_file_ ) ) utility_exit_with_message( silent_file_ + " already exist!" );

	std::ofstream outfile;
	outfile.open( silent_file_.c_str() ); //Opening the file with this command removes all prior content..

	outfile << "StepWiseRNA_Minimizer:: num_pose_outputted == 0, empty silent_file!\n"; //specific key signal to SWA_sampling_post_process.py

	outfile.flush();
	outfile.close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::Freeze_sugar_torsions( core::kinematics::MoveMap & mm, Size const nres ) const {

	 using namespace core::id;

	 TR.Debug << "Freeze pose sugar torsions, nres = " << nres << std::endl;

	 for ( Size i = 1; i <= nres; i++ ){

			mm.set( TorsionID( i , id::BB,  4 ), false ); //delta_i
			mm.set( TorsionID( i , id::CHI, 2 ), false ); //nu2_i
			mm.set( TorsionID( i , id::CHI, 3 ), false );	//nu1_i

	 }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////

//Cannot pass pose in as constant due to the setPoseExtraScores function
bool
StepWiseRNA_Minimizer::pass_all_pose_screens( core::pose::Pose & pose, std::string const in_tag, core::io::silent::SilentFileData & silent_file_data ) const{

	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace protocols::rna;
	using namespace core::optimization;


	Size const gap_size( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */

	bool pass_screen = true;

	//March 03, 2012. Don't need this for post-process analysis. For full_length_job_params, is_prepend, moving_res and five_prime_chain_break_res are not well defined.
	if ( job_parameters_->is_simple_full_length_job_params() == false ){

		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
		bool const is_prepend(  job_parameters_->is_prepend() );
		Size const moving_res(  job_parameters_->working_moving_res() );

		if ( gap_size == 0 && five_prime_chain_break_res > 0 ){
			//Sept 19, 2010 screen for weird delta/sugar conformation. Sometime at the chain closure step, minimization will severely screw up the pose sugar
			//if the chain is not a closable position. These pose will have very bad score and will be discarded anyway. The problem is that the clusterer check
			//for weird delta/sugar conformation and call utility_exit_with_message() if one is found (purpose being to detect bad silent file conversion)
			//So to prevent code exit, will just screen for out of range sugar here.
			//Oct 28, 2010 ...ok check for messed up nu1 and nu2 as well. Note that check_for_messed_up_structure() also check for messed up delta but the range is smaller than below.
			if ( check_for_messed_up_structure( pose, in_tag ) == true ){
				TR.Debug << "gap_size == 0, " << in_tag << " discarded: messed up structure " << std::endl;
				pass_screen = false;
			}

			conformation::Residue const & five_prime_rsd = pose.residue( five_prime_chain_break_res );
			Real const five_prime_delta = numeric::principal_angle_degrees( five_prime_rsd.mainchain_torsion( DELTA ) );



			if ( ( five_prime_rsd.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) == false ) && ( five_prime_rsd.has_variant_type( "VIRTUAL_RIBOSE" ) == false ) ){
				if ( ( five_prime_delta > 1.0 && five_prime_delta < 179.00 ) == false ){

					TR.Debug << "gap_size == 0, " << in_tag << " discarded: five_prime_chain_break_res = " << five_prime_chain_break_res << " five_prime_CB_delta = " << five_prime_delta << " is out of range " << std::endl;
					pass_screen = false;
				}
			}

			conformation::Residue const & three_prime_rsd = pose.residue( five_prime_chain_break_res + 1 );
			Real const three_prime_delta = numeric::principal_angle_degrees( three_prime_rsd.mainchain_torsion( DELTA ) );

			if ( ( three_prime_rsd.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) == false ) && ( three_prime_rsd.has_variant_type( "VIRTUAL_RIBOSE" ) == false ) ){
				if ( ( three_prime_delta > 1.0 && three_prime_delta < 179.00 ) == false ){

					TR.Debug << "gap_size == 0, " << in_tag << " discarded: three_prime_chain_break_res = " << ( five_prime_chain_break_res + 1 ) << " three_prime_CB_delta = " << three_prime_delta << " is out of range " << std::endl;
					pass_screen = false;
				}
			}
		}


		if ( centroid_screen_ ){
			if ( base_centroid_screener_ == 0 ) utility_exit_with_message( "base_centroid_screener_ == 0!" );

			if ( !base_centroid_screener_->update_base_stub_list_and_check_that_terminal_res_are_unstacked( pose, true /*reinitialize*/ ) ){
				TR.Debug << in_tag << " discarded: fail check_that_terminal_res_are_unstacked	" << std::endl;
				pass_screen = false;
			}
		}

		screener::ChainClosableScreener chain_closable_screener( five_prime_chain_break_res, gap_size );
		if ( gap_size == 1 &&  !chain_closable_screener.check_screen( pose ) ) pass_screen = false;

		if ( native_screen_ && get_native_pose() ){	//Before have the (&& !is_chain_break condition). Parin Dec 21, 2009

			Real const rmsd = suite_rmsd( *get_native_pose(), pose,  moving_res, is_prepend, true /*ignore_virtual_atom*/ );
			Real const loop_rmsd =	 rmsd_over_residue_list( *get_native_pose(), pose, job_parameters_, true /*ignore_virtual_atom*/ );

			if ( rmsd > native_screen_rmsd_cutoff_ || loop_rmsd > native_screen_rmsd_cutoff_ ){
				TR.Debug << in_tag << " discarded: fail native_rmsd_screen. rmsd = " << rmsd << " loop_rmsd = " << loop_rmsd << " native_screen_rmsd_cutoff_ = " << native_screen_rmsd_cutoff_ << std::endl;
				pass_screen = false;
			}
		}
	}

	if ( perform_electron_density_screen_ ){ //Fang's electron density code

		//Note to Fang: I moved this inside this funciton. While creating new pose is time-consuming..this does be effect code speed since the rate-limiting step here is the minimizer.run(). Parin S.
		//setup native score screening
		pose::Pose native_pose = *( job_parameters_->working_native_pose() ); //Hard copy!

		for ( Size i = 1; i <= native_pose.total_residue(); ++i ) {
			pose::remove_variant_type_from_pose_residue( native_pose, "VIRTUAL_PHOSPHATE", i );
		}
		/////////////////////////////////////////////////////////

		bool const pass_native = native_edensity_score_screener( pose, native_pose );
		if ( pass_native == false ){
			TR.Debug << in_tag << " discarded: fail native_edensity_score_screening" << std::endl;
			pass_screen = false;
		}
	}

	//March 20, 2011..This is neccessary though to be consistent with FARFAR. Cannot have false low energy state that are due to empty holes in the structure.
	//Feb 22, 2011. Warning this is inconsistent with the code in SAMPLERER:
	//In Both standard and floating base sampling, user_input_VDW_bin_screener_ is not used for gap_size==0 or internal case.
	//This code is buggy for gap_size==0 case IF VDW_screener_pose contains the residue right at 3' or 5' of the building_res
	//Internal case should be OK.
	if ( user_input_VDW_bin_screener_ && user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ){

		//Convert to using physical_pose instead of bin for screening (as default) on June 04, 2011
		utility::vector1 < core::Size > const & working_global_sample_res_list = job_parameters_->working_global_sample_res_list();
		bool const pass_physical_pose_VDW_rep_screen = user_input_VDW_bin_screener_->VDW_rep_screen_with_act_pose( pose, working_global_sample_res_list, true /*local verbose*/ );

		if ( pass_physical_pose_VDW_rep_screen == false ){
			TR.Debug << in_tag << " discarded: fail physical_pose_VDW_rep_screen" << std::endl;
			pass_screen = false;
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::string temp_tag = in_tag;

	if ( verbose_ ) output_pose_wrapper( temp_tag, 'S', pose, silent_file_data, silent_file_ + "_screen" );

	return pass_screen;

}
////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_Minimizer::native_edensity_score_screener( pose::Pose & pose, pose::Pose & native_pose ) const{

	using namespace core::scoring;

	//get the score of elec_dens_atomwise only
	static ScoreFunctionOP eden_scorefxn = new ScoreFunction;
	eden_scorefxn -> set_weight( elec_dens_atomwise, 1.0 );
	core::Real pose_score = ( ( *eden_scorefxn )( pose ) );
	core::Real native_score = ( ( *eden_scorefxn )( native_pose ) );
	core::Size nres = pose.total_residue();
	core::Real native_score_cutoff = native_score / ( static_cast < double > ( nres ) ) *
	( 1 - native_edensity_score_cutoff_ );
	TR.Debug << "pose_score = " << pose_score << std::endl;
	TR.Debug << "native_score = " << native_score << std::endl;
	if ( pose_score > native_score - native_score_cutoff ) {
		TR.Debug << "Fail native edensity score screening!" << std::endl;
		return false;
	} else {
		return true;
	}
}


////////////////////////////////////////////////////////////////////////////////////////
utility::vector1 < core::Size >
StepWiseRNA_Minimizer::get_working_moving_res( Size const & nres ) const {
	utility::vector1 < core::Size > working_moving_res;

	utility::vector1 < core::Size > const & working_fixed_res = job_parameters_->working_fixed_res();
	utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();
	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
		if ( !working_fixed_res.has_value( seq_num ) || working_moving_partition_pos.has_value( seq_num ) ){
			working_moving_res.push_back( seq_num );
		}
	}
	//	output_seq_num_list( "working_moving_res", working_moving_res );
	return working_moving_res;
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1 < core::kinematics::MoveMap >
StepWiseRNA_Minimizer::get_default_movemap( core::pose::Pose const & pose ) {

	utility::vector1 < core::kinematics::MoveMap > move_map_list;

	core::kinematics::MoveMap mm;
	figure_out_moving_residues( mm, pose );

	// Allow sugar torsions to move again (RD 01/31/2010), now
	// that rotamers have been pre-optimized for sugar closure, and
	// sugar_close is turned back on.
	//		Freeze_sugar_torsions(mm, nres); //Freeze the sugar_torsions!

	move_map_list.push_back( mm );
	return move_map_list;
}

////////////////////////////////////////////////////////////////////////////////////
// This is similar to code in RNA_Minimizer.cc
void
StepWiseRNA_Minimizer::figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose )
{
	using namespace core::id;

	utility::vector1< core::Size > const & fixed_res( job_parameters_->working_fixed_res() );
	allow_insert_ = new toolbox::AllowInsert( pose );
	for ( Size i = 1; i <= fixed_res.size(); i++ ) allow_insert_->set( fixed_res[ i ], false ) ;
	update_allow_insert_with_extra_minimize_res( pose, allow_insert_, extra_minimize_res_ );
	figure_out_swa_rna_movemap( mm, pose, allow_insert_ );
}

////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_minimized_pose_list(){

	io::silent::SilentFileData silent_file_data;

	for ( Size pose_ID = 1; pose_ID <= minimized_pose_list_.size(); pose_ID++ ){
		pose::Pose & pose = *(minimized_pose_list_[pose_ID]);
		std::string tag;
		if ( rename_tag_ ){
			tag = "S_" + ObjexxFCL::lead_zero_string_of( pose_ID /* start with zero */, 6 );
		} else{
			tag = tag_from_pose( pose );
			tag[0] = 'M';  //M for minimized, these are the poses that will be clustered by master node.
		}
		output_pose_wrapper( tag, tag[0], pose, silent_file_data, silent_file_ );
	}
}

////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_move_map_list( utility::vector1 < core::kinematics::MoveMap > const & move_map_list ){
	move_map_list_ = move_map_list;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_silent_file( std::string const & silent_file ){
	silent_file_ = silent_file;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
core::io::silent::SilentFileDataOP &
StepWiseRNA_Minimizer::silent_file_data(){
	return sfd_;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener ){
	base_centroid_screener_ = screener;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_native_edensity_score_cutoff( core::Real const & setting ){
	native_edensity_score_cutoff_ = setting;
	perform_electron_density_screen_ = ( native_edensity_score_cutoff_ > -0.99999 || native_edensity_score_cutoff_ < -1.00001 );
}
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_parameters(){
	output_boolean( " verbose_ = ", verbose_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( " native_screen_ = ", native_screen_, TR.Debug ); TR.Debug << std::endl;
	TR.Debug << " native_screen_rmsd_cutoff_ = " << native_screen_rmsd_cutoff_ << std::endl;
	output_boolean( " perform_electron_density_screen_ = ", perform_electron_density_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( " rm_virt_phosphate_ = ", rm_virt_phosphate_, TR.Debug ); TR.Debug << std::endl;
	TR.Debug << " native_edensity_score_cutoff_ = " << native_edensity_score_cutoff_ << std::endl;
	output_boolean( " centroid_screen_ = ", centroid_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( " perform_o2prime_pack_ = ", perform_o2prime_pack_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( " output_before_o2prime_pack_ = ", output_before_o2prime_pack_, TR.Debug ); TR.Debug << std::endl;
	TR.Debug << " ( Upper_limit ) num_pose_minimize_ = " << num_pose_minimize_ << " pose_list.size() = " << pose_list_.size() <<  std::endl;
	output_boolean( " minimize_and_score_sugar_ = ", minimize_and_score_sugar_, TR.Debug ); TR.Debug << std::endl;
	if (user_input_VDW_bin_screener_) output_boolean( " user_inputted_VDW_bin_screener_ = ", user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose(), TR.Debug ); TR.Debug << std::endl;
	output_seq_num_list( " working_global_sample_res_list = ", job_parameters_->working_global_sample_res_list(), TR.Debug );
	output_boolean( " perform_minimize_ = ", perform_minimize_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( " rename_tag_ = ", rename_tag_, TR.Debug ); TR.Debug << std::endl;
}


}
}
}
