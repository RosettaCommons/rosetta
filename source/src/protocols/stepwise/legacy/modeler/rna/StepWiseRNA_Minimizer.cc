// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseRNA_Minimizer
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @details
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


//////////////////////////////////
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_OutputData.hh> //Sept 26, 2011
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh> //Sept 26, 2011
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/movemap/util.hh>
#include <protocols/stepwise/modeler/util.hh>
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
#include <core/chemical/rna/util.hh> // for DELTA
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009
#include <core/conformation/Conformation.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>

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
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::rna;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise modeler of RNA (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.legacy.modeler.rna.StepWiseRNA_Minimizer" );

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_Minimizer::StepWiseRNA_Minimizer(
		utility::vector1 < core::pose::PoseOP > const & pose_list,
		stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP & working_parameters ) :
	pose_list_( pose_list ),
	working_parameters_( working_parameters ),
	silent_file_( "silent_file.txt" ),
	perform_electron_density_screen_( false ), // updated below
	base_centroid_checker_( 0 ), //Owning-pointer.
	vary_bond_geometry_frequency_( 0.0 ),
	allow_variable_bond_geometry_( false )
{
	set_native_pose( working_parameters_->working_native_pose() );
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
	using namespace protocols::farna;
	using namespace protocols::stepwise;

	//Real const bulge_weight = scorefxn_->get_weight( num_stacks );
	//scorefxn_->set_weight( num_stacks, 0.0);
	original_geometry_weight_ = scorefxn_->get_weight( rna_bond_geometry );
	stepwise::modeler::rna::output_title_text( "Enter StepWiseRNA_Minimizer::apply", TR.Debug );
	output_parameters();
	clock_t const time_start( clock() );

	Size const gap_size(  working_parameters_->gap_size() );
	bool const close_chainbreak = gap_size == 0 && ( working_parameters_->moving_res() > 0 );
	RNA_LoopCloser rna_loop_closer;
	SilentFileData silent_file_data;
	perform_electron_density_screen_ = ( options_->native_edensity_score_cutoff() > -0.99999 || options_->native_edensity_score_cutoff() < -1.00001 );

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025 );
	bool const use_nblist( true );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	/////////////////////////////////////////////////////////
	if ( pose_list_.size() == 0 ){
		TR.Debug << "pose_list_.size() == 0, early exit from StepWiseRNA_Minimizer::apply" << std::endl;
		output_empty_minimizer_silent_file();
		return;
	}

	pose::Pose dummy_pose = ( *pose_list_[1] );
	Size const nres =  dummy_pose.size();
	bool o2prime_pack_verbose( options_->verbose() );
	utility::vector1< Size > const working_moving_res = get_working_moving_res( nres ); // will be used for o2prime packing.

	//Check scorefxn
	TR.Debug << "check scorefxn" << std::endl;
	scorefxn_->show( TR.Debug, dummy_pose );

	bool vary_bond_geometry( false );
	if ( allow_variable_bond_geometry_ && ( numeric::random::rg().uniform() < vary_bond_geometry_frequency_ ) ) {
		vary_bond_geometry = true;
	}

	if ( !move_map_ ) move_map_ = get_default_movemap( dummy_pose );
	if ( ! options_->minimize_and_score_sugar() )	freeze_sugar_torsions( *move_map_, nres );
	//	output_movemap( *move_map_, dummy_pose, TR.Debug );

	minimized_pose_list_.clear();

	for ( Size pose_ID = 1; pose_ID <= pose_list_.size(); pose_ID++ ){

		if ( options_->num_pose_minimize() > 0 &&
				 minimized_pose_list_.size() >= options_->num_pose_minimize() ){
			TR.Debug << "WARNING MAX options_->num_pose_minimize()( " << options_->num_pose_minimize() << " ) EXCEEDED, EARLY BREAK." << std::endl;
			break;
		}

		 TR.Debug << "Minimizing pose_ID = " << pose_ID << std::endl;

		std::string tag = tag_from_pose( *pose_list_[pose_ID] );
		pose = ( *pose_list_[pose_ID] ); //This actually creates a hard copy.....need this to get the graphic working..

		if ( options_->rm_virt_phosphate() ){ //Fang's electron density code
			for ( Size ii = 1; ii <= pose.size(); ++ii ) {
				pose::remove_variant_type_from_pose_residue( pose, chemical::VIRTUAL_PHOSPHATE, ii );
			}
		}

		if ( options_->verbose() && options_->minimizer_output_before_o2prime_pack() ) output_pose_wrapper( tag, 'B', pose, silent_file_data, silent_file_ + "_before_o2prime_pack" );

		remove_virtual_O2Prime_hydrogen( pose );

		if ( options_->minimizer_perform_o2prime_pack() ) o2prime_trials( pose, scorefxn_, get_surrounding_O2prime_hydrogen( pose, working_moving_res, o2prime_pack_verbose ) );

		if ( options_->verbose() && !options_->minimizer_output_before_o2prime_pack() ) output_pose_wrapper( tag, 'B', pose, silent_file_data, silent_file_ + "_before_minimize" );

		utility::vector1< Size > moving_chainbreaks = figure_out_moving_chain_break_res( pose, *move_map_ );

		///////Minimization/////////////////////////////////////////////////////////////////////////////
		if ( close_chainbreak ){
			rna_loop_closer.apply( pose, moving_chainbreaks ); //This doesn't do anything if rna_loop_closer was already applied during modeler stage...May 10,2010
			if ( options_->verbose() ) output_pose_wrapper( tag, 'C', pose, silent_file_data, silent_file_ + "_after_loop_closure_before_minimize" );
		}

		// set up move_map on the fly if variants are different in different poses! E.g. packed phosphates.
		move_map_ = get_default_movemap( pose );
		core::kinematics::MoveMap & mm = *move_map_;
		//		output_movemap( *move_map_, dummy_pose, TR );

		moving_chainbreaks = figure_out_moving_chain_break_res( pose, mm );

 		if ( vary_bond_geometry ) {
			TR << "Performing variable geometry minimization..." << std::endl;
			ScoreFunctionOP scorefxn_vbg = scorefxn_->clone();
			scorefxn_vbg->set_weight( rna_bond_geometry, 8.0 );
			setup_vary_rna_bond_geometry( mm, pose, atom_level_domain_map_ );
			if ( !options_->skip_minimize() ) minimizer.run( pose, new_mm, *(scorefxn_vbg ), options );
		}
		if ( !options_->skip_minimize() ) minimizer.run( pose, mm, *( scorefxn_ ), options );

		if ( options_->minimizer_perform_o2prime_pack() ) o2prime_trials( pose, scorefxn_, get_surrounding_O2prime_hydrogen( pose, working_moving_res, o2prime_pack_verbose ) );

		if ( close_chainbreak ){
			rna_loop_closer.apply( pose, moving_chainbreaks );
			if ( options_->minimizer_perform_o2prime_pack() ) o2prime_trials( pose, scorefxn_, get_surrounding_O2prime_hydrogen( pose, working_moving_res, o2prime_pack_verbose ) );
			if ( !options_->skip_minimize() ) minimizer.run( pose, mm, *( scorefxn_ ), options );

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

	if ( options_->output_minimized_pose_list() ) output_minimized_pose_list();

	// output best scoring pose.
	if ( minimized_pose_list_.size() > 0 ) pose = *( minimized_pose_list_[1] );

	TR.Debug << "--------------------------------------------------------------------" << std::endl;
	TR.Debug << "Total time in StepWiseRNA_Minimizer: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
	TR.Debug << "--------------------------------------------------------------------" << std::endl;

	stepwise::modeler::rna::output_title_text( "Exit StepWiseRNA_Minimizer::apply", TR.Debug );

}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_pose_wrapper(
	std::string const & tag,
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
	output_data( silent_file_data, out_silent_file, tag, false, pose, get_native_pose(), working_parameters_ );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_empty_minimizer_silent_file() const {

	if ( utility::file::file_exists( silent_file_ ) ) utility_exit_with_message( silent_file_ + " already exist!" );

	std::ofstream outfile;
	outfile.open( silent_file_.c_str() ); //Opening the file with this command removes all prior content..

	outfile << "StepWiseRNA_Minimizer:: num_pose_outputted == 0, empty silent_file!\n"; //specific key signal to SWA_modeler_post_process.py

	outfile.flush();
	outfile.close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::freeze_sugar_torsions( core::kinematics::MoveMap & mm, Size const nres ) const {

	 using namespace core::id;

	 TR.Debug << "Freeze pose sugar torsions, nres = " << nres << std::endl;

	 for ( Size i = 1; i <= nres; i++ ){

			mm.set( TorsionID( i , id::BB,  4 ), false ); //delta_i
			mm.set( TorsionID( i , id::CHI, 2 ), false ); //nu2_i
			mm.set( TorsionID( i , id::CHI, 3 ), false );	//nu1_i

	 }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////

//Cannot pass pose in as constant due to the setPoseExtraScore function
bool
StepWiseRNA_Minimizer::pass_all_pose_screens( core::pose::Pose & pose, std::string const in_tag, core::io::silent::SilentFileData & silent_file_data ) const {

	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace protocols::farna;
	using namespace core::optimization;


	Size const gap_size( working_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */

	bool pass_screen = true;

	//March 03, 2012. Don't need this for post-process analysis. For full_length_job_params, is_prepend, moving_res and five_prime_chain_break_res are not well defined.
	if ( working_parameters_->is_simple_full_length_job_params() == false ){

		Size const five_prime_chain_break_res = working_parameters_->five_prime_chain_break_res();
		bool const is_prepend(  working_parameters_->is_prepend() );
		Size const moving_res(  working_parameters_->working_moving_res() );

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

			if ( ( five_prime_rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) == false ) && ( five_prime_rsd.has_variant_type( chemical::VIRTUAL_RIBOSE ) == false ) ){
				if ( ( five_prime_delta > 1.0 && five_prime_delta < 179.00 ) == false ){

					TR.Debug << "gap_size == 0, " << in_tag << " discarded: five_prime_chain_break_res = " << five_prime_chain_break_res << " five_prime_CB_delta = " << five_prime_delta << " is out of range " << std::endl;
					pass_screen = false;
				}
			}

			conformation::Residue const & three_prime_rsd = pose.residue( five_prime_chain_break_res + 1 );
			Real const three_prime_delta = numeric::principal_angle_degrees( three_prime_rsd.mainchain_torsion( DELTA ) );

			if ( ( three_prime_rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) == false ) && ( three_prime_rsd.has_variant_type( chemical::VIRTUAL_RIBOSE ) == false ) ){
				if ( ( three_prime_delta > 1.0 && three_prime_delta < 179.00 ) == false ){
					TR.Debug << "gap_size == 0, " << in_tag << " discarded: three_prime_chain_break_res = " << ( five_prime_chain_break_res + 1 ) << " three_prime_CB_delta = " << three_prime_delta << " is out of range " << std::endl;
					pass_screen = false;
				}
			}
		}


		if ( base_centroid_checker_ != 0 ){

			if ( !base_centroid_checker_->update_base_stub_list_and_check_that_terminal_res_are_unstacked( pose, true /*reinitialize*/ ) ){
				TR.Debug << in_tag << " discarded: fail check_that_terminal_res_are_unstacked	" << std::endl;
				pass_screen = false;
			}
		}

		checker::RNA_ChainClosableGeometryChecker chain_closable_geometry_checker( five_prime_chain_break_res, gap_size );
		if ( gap_size == 1 &&  !chain_closable_geometry_checker.check_screen( pose ) ) pass_screen = false;


		if ( options_->rmsd_screen() && get_native_pose() ){	//Before have the (&& !is_chain_break condition). Parin Dec 21, 2009
			Real const rmsd = suite_rmsd( *get_native_pose(), pose,  moving_res, is_prepend, true /*ignore_virtual_atom*/ );
			Real const loop_rmsd =	 rmsd_over_residue_list( *get_native_pose(), pose, working_parameters_, true /*ignore_virtual_atom*/ );

			if ( rmsd > options_->rmsd_screen() || loop_rmsd > options_->rmsd_screen() ){
				TR.Debug << in_tag << " discarded: fail rmsd_screen. rmsd = " << rmsd << " loop_rmsd = " << loop_rmsd << " options_->sampler_native_screen_rmsd_cutoff() = " << options_->rmsd_screen() << std::endl;
				pass_screen = false;
			}
		}
	}

	if ( perform_electron_density_screen_ ){ //Fang's electron density code

		//Note to Fang: I moved this inside this funciton. While creating new pose is time-consuming..this does be effect code speed since the rate-limiting step here is the minimizer.run(). Parin S.
		//setup native score screening
		pose::Pose native_pose = *( working_parameters_->working_native_pose() ); //Hard copy!

		for ( Size i = 1; i <= native_pose.size(); ++i ) {
			pose::remove_variant_type_from_pose_residue( native_pose, chemical::VIRTUAL_PHOSPHATE, i );
		}
		/////////////////////////////////////////////////////////

		bool const pass_native = native_edensity_score_checker( pose, native_pose );
		if ( pass_native == false ){
			TR.Debug << in_tag << " discarded: fail native_edensity_score_screening" << std::endl;
			pass_screen = false;
		}
	}

	//March 20, 2011..This is neccessary though to be consistent with FARFAR. Cannot have false low energy state that are due to empty holes in the structure.
	//Feb 22, 2011. Warning this is inconsistent with the code in SAMPLERER:
	//In Both standard and floating base modeler, user_input_VDW_bin_checker_ is not used for gap_size==0 or internal case.
	//This code is buggy for gap_size==0 case IF VDW_checker_pose contains the residue right at 3' or 5' of the building_res
	//Internal case should be OK.
	if ( user_input_VDW_bin_checker_ && user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ){

		//Convert to using physical_pose instead of bin for screening (as default) on June 04, 2011
		utility::vector1 < core::Size > const & working_global_sample_res_list = working_parameters_->working_global_sample_res_list();
		bool const pass_physical_pose_VDW_rep_screen = user_input_VDW_bin_checker_->VDW_rep_screen_with_act_pose( pose, working_global_sample_res_list, true /*local verbose*/ );

		if ( pass_physical_pose_VDW_rep_screen == false ){
			TR.Debug << in_tag << " discarded: fail physical_pose_VDW_rep_screen" << std::endl;
			pass_screen = false;
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::string temp_tag = in_tag;

	if ( options_->verbose() ) output_pose_wrapper( temp_tag, 'S', pose, silent_file_data, silent_file_ + "_screen" );

	return pass_screen;

}
////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_Minimizer::native_edensity_score_checker( pose::Pose & pose, pose::Pose & native_pose ) const{

	using namespace core::scoring;

	//get the score of elec_dens_atomwise only
	static ScoreFunctionOP eden_scorefxn = new ScoreFunction;
	eden_scorefxn -> set_weight( elec_dens_atomwise, 1.0 );
	core::Real pose_score = ( ( *eden_scorefxn )( pose ) );
	core::Real native_score = ( ( *eden_scorefxn )( native_pose ) );
	core::Size nres = pose.size();
	core::Real native_score_cutoff = native_score / ( static_cast < double > ( nres ) ) *
	( 1 - options_->native_edensity_score_cutoff() );
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

	utility::vector1 < core::Size > const & working_fixed_res = working_parameters_->working_fixed_res();
	utility::vector1 < core::Size > const & working_moving_partition_res = working_parameters_->working_moving_partition_res();
	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
		if ( !working_fixed_res.has_value( seq_num ) || working_moving_partition_res.has_value( seq_num ) ){
			working_moving_res.push_back( seq_num );
		}
	}
	//	stepwise::modeler::rna::output_seq_num_list( "working_moving_res", working_moving_res );
	return working_moving_res;
}

////////////////////////////////////////////////////////////////////////////////////////
core::kinematics::MoveMapOP
StepWiseRNA_Minimizer::get_default_movemap( core::pose::Pose const & pose ) {
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	figure_out_stepwise_movemap( *mm, atom_level_domain_map_, pose, working_parameters_->working_fixed_res(), working_extra_minimize_res_ );
	return mm;
}

////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_minimized_pose_list(){

	io::silent::SilentFileData silent_file_data;

	for ( Size pose_ID = 1; pose_ID <= minimized_pose_list_.size(); pose_ID++ ){
		pose::Pose & pose = *(minimized_pose_list_[pose_ID]);
		std::string tag;
		if ( options_->minimizer_rename_tag() ){
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
StepWiseRNA_Minimizer::set_move_map( core::kinematics::MoveMapOP move_map ) {
	move_map_ = move_map;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_silent_file( std::string const & silent_file ){
	silent_file_ = silent_file;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
core::io::silent::SilentFileDataOP &
StepWiseRNA_Minimizer::silent_file_data(){
	return sfd_;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker ){
	base_centroid_checker_ = checker;
}

//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::output_parameters(){
	output_boolean( " verbose = ", options_->verbose(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( " sampler_native_rmsd_screen = ", options_->rmsd_screen(), TR.Debug ); TR.Debug << std::endl;
	TR.Debug << " sampler_native_screen_rmsd_cutoff = " << options_->rmsd_screen() << std::endl;
	output_boolean( " perform_electron_density_screen_ = ", perform_electron_density_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( " rm_virt_phosphate = ", options_->rm_virt_phosphate(), TR.Debug ); TR.Debug << std::endl;
	TR.Debug << " native_edensity_score_cutoff = " << options_->native_edensity_score_cutoff() << std::endl;
	output_boolean( " centroid_screen = ", (base_centroid_checker_ != 0), TR.Debug ); TR.Debug << std::endl;
	output_boolean( " minimizer_perform_o2prime_pack = ", options_->minimizer_perform_o2prime_pack(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( " minimizer_output_before_o2prime_pack = ", options_->minimizer_output_before_o2prime_pack(), TR.Debug ); TR.Debug << std::endl;
	TR.Debug << " ( Upper_limit ) num_pose_minimize = " << options_->num_pose_minimize() << " pose_list.size() = " << pose_list_.size() <<  std::endl;
	output_boolean( " minimize_and_score_sugar = ", options_->minimize_and_score_sugar(), TR.Debug ); TR.Debug << std::endl;
	if (user_input_VDW_bin_checker_) output_boolean( " user_inputted_VDW_bin_checker_ = ", user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose(), TR.Debug ); TR.Debug << std::endl;
	stepwise::modeler::rna::output_seq_num_list( " working_global_sample_res_list = ", working_parameters_->working_global_sample_res_list(), TR.Debug );
	output_boolean( " perform_minimize = ", !options_->skip_minimize(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( " minimizer_rename_tag = ", options_->minimizer_rename_tag(), TR.Debug ); TR.Debug << std::endl;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Minimizer::set_options( protocols::stepwise::modeler::options::StepWiseModelerOptionsOP options ){
	options_ = options;
}


} //rna
} //modeler
} //legacy
} //stepwise
} //protocols
