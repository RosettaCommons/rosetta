// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/mpi_refinement/MPI_Refinement.cc
/// @brief
/// @author Mike Tyka
#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/mpi_refinement/MPI_Refinement.hh>
#include <protocols/mpi_refinement/MultiObjective.hh>
#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
#include <protocols/mpi_refinement/util.hh>

#include <protocols/relax/WorkUnit_BatchRelax.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
//#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// Utility headers
#include <numeric/random/random.hh>

// C/C++ headers
#ifndef _WIN32 // REQUIRED FOR WINDOWS
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <utility/vector1.hh>

using namespace ObjexxFCL::format;

namespace protocols {
namespace mpi_refinement {

using namespace protocols::wum;

static THREAD_LOCAL basic::Tracer TR("MPI.LHR.A");

MPI_Refinement::MPI_Refinement( char machine_letter ):
	MPI_WorkUnitManager( machine_letter ),
	max_lib_size_(50),
	max_ref_lib_size_(50),
	save_state_interval_(100000), // don't save... original value is 300
	last_save_state_(0),
	totaltime_loophash_(0),
	n_loophash_(0),
	totaltime_batchrelax_(0),
	n_batchrelax_(0),
	total_structures_(0),
	total_structures_relax_(0),
	total_metropolis_(1),
	total_metropolis_accepts_(0),
	ident_string_("ident")
{
	set_defaults();  // constructors must must must call this!
}

void
MPI_Refinement::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	max_lib_size_               = option[ OptionKeys::lh::max_lib_size ]();
	max_ref_lib_size_           = option[ OptionKeys::lh::max_ref_lib_size ]();
	save_state_interval_        = option[ OptionKeys::lh::mpi_save_state_interval ]();
	mpi_feedback_               = option[ OptionKeys::lh::mpi_feedback ]();
	mpi_metropolis_temp_        = option[ OptionKeys::lh::mpi_metropolis_temp ]();
	rms_limit_                  = option[ OptionKeys::lh::rms_limit ]();
	//objective_function_         = option[ OptionKeys::lh::objective_function ]();

	core::chemical::ResidueTypeSetCOP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	native_given_ = false;
	if ( option[ in::file::native ].user() ) {
		native_given_ = true;
		core::import_pose::pose_from_file( native_pose_, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file);
	}

	// initial "reference" structure should be provided as the first at -in:file:s
	core::import_pose::pose_from_file( pose0_, *rsd_set, option[ in::file::s ](1) , core::import_pose::PDB_file);

	fobj_ = MultiObjectiveOP( new MultiObjective() );

	mpi_resume_ = "";
	if ( option[ OptionKeys::lh::mpi_resume ].user() ) {
		mpi_resume_ = option[ OptionKeys::lh::mpi_resume ]();
	}
	jobname_ = option[ OptionKeys::lh::jobname ]();
	sim_replace_obj_ = option[ OptionKeys::lh::sim_replace_obj ]();
	// Make ident string:
	// Make a medley of a string and the PID
	ident_string_ =  option[ lh::jobname ];
	TR << "IDENT: " << ident_string_ << std::endl;

	// make sure the state saves are randomly staggered - otherwise all the masters dump several hundred megs at once!
	last_save_state_ = time(nullptr)  + core::Size( numeric::random::rg().uniform()  * (core::Real) save_state_interval_) ;
	TR.Debug << "Interlace dumps: " << last_save_state_ << "  " << time(nullptr)  << "  " << last_save_state_ - time(nullptr) << "  " << save_state_interval_ << "  " << std::endl;

	starttime_ = time(nullptr);
}

void
MPI_Refinement::load_structures_from_cmdline_into_library(
	protocols::wum::SilentStructStore &library )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;

	start_timer( TIMING_IO_READ  );
	TR << "Reading in structures..." << std::endl;

	core::chemical::ResidueTypeSetCOP rsd_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

	core::Size count = 0;
	SilentStructStore temp_lib;

	// Rescore input structures for fair comparison with whatever comes back from relax calls
	//scorefxn = core::scoring::getScoreFunction();
	core::scoring::ScoreFunctionOP scorefxn = fobj_->get_scorefxn( 1 ); // this should be rosetta standard...
	core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn );
	scorefxn->set_weight( core::scoring::coordinate_constraint, 10.0 ); //10.0 will allow about 0.1 ang diff...

	while ( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		core::pose::Pose pose_save( pose );

		// Minimize pose prior to store
		if ( option[ OptionKeys::lh::mpi_packmin_init ]() ) {
			protocols::wum::SilentStructStore library_tmp;
			core::pose::set_ss_from_phipsi( pose );
			core::io::silent::SilentFileOptions opts;
			core::io::silent::SilentStructOP ss = option[ OptionKeys::lh::bss]() ?
				core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary",opts) :
				core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out(opts);

			ss->fill_struct( pose );
			TR << "STARTLIB, before min" << std::endl;
			fobj_->add_objective_function_info( library_tmp );
			print_library( library_tmp, "MasterLIB Init");

			if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				ConstraintSetOP cstset =
					ConstraintIO::get_instance()->read_constraints( get_cst_fa_file_option(),
					ConstraintSetOP( new ConstraintSet ), pose );
				TR << "Read constraints: " << ( cstset->has_constraints() ? "YES" : "NONE" ) << std::endl;
				pose.constraint_set( cstset );
			} else { // strong coordinate constraint
				protocols::relax::AtomCoordinateCstMover coord_cst_mover;
				coord_cst_mover.cst_sd( 1.0 );
				coord_cst_mover.apply( pose );
			}
			ramp_minpack_pose( pose, scorefxn, true, false ); // this is cartmin!
			// For debug
			core::Real rmsd = core::scoring::CA_rmsd( pose, pose_save );
			TR << "rmsd after minimizing: " << rmsd << std::endl;
		}

		(*scorefxn)(pose);

		core::pose::set_ss_from_phipsi( pose );
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentStructOP ss = option[ OptionKeys::lh::bss]() ?
			core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary",opts) :
			core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out(opts);

		ss->fill_struct( pose );
		ss->add_energy( "parent_score", 0 );  // fa score of parent
		ss->add_energy( "nuse", 0 );   // lh count is the number of times this structure has been used to generate new structures
		ss->add_energy( "round", 0 );     // round is the number of consecutive loophashes this structure has gone through
		ss->add_energy( "iter", 0 );     // round is the number of consecutive loophashes this structure has gone through
		ss->add_energy( "expire", 0 );    // how many times has this structure been retired to the temporor ? >1 means it's been sent out again
		//if(  ss->get_string_value("usid") == "" ){
		//  ss->add_string_value( "usid", "_" + ObjexxFCL::string_of(count)); // A unique structure id (usid) that will identify the structure until it is modified.
		//}
		//ss->add_string_value( "husid", ss->get_string_value("usid") ); // history of usids

		ss->add_energy( "state", 0 );     // state: 0 init, 1 loophashed, 2 relaxed

		if ( native_given_ ) add_poseinfo_to_ss( *ss, native_pose_, "" );
		add_poseinfo_to_ss( *ss, pose0_, "_i" );

		library.add( ss );
		count ++;
	}

	if ( library.size() == 0 ) {
		TR.Error << "Error reading starting structures: 0 valid structures read. Check your input files and parameters" << std::endl;
		utility_exit_with_message( "Error reading starting structures: 0 valid structures read. Check your input files and parameters" );
	}

	fobj_->add_objective_function_info( library );
	TR << "Added " << library.size() << " structures to library " << std::endl;

	runtime_assert(  library_central_.size() <= max_lib_size_ );
	start_timer( TIMING_CPU );
}

void
MPI_Refinement::save_state(std::string prefix ){
	start_timer( TIMING_IO_WRITE );
	long starttime = time(nullptr);
	write_queues_to_file( prefix + "." + ObjexxFCL::string_of(mpi_rank()) );
	library_central_.serialize_to_file( prefix + "." + ObjexxFCL::string_of(mpi_rank())+ ".lib.library_central" );
	long endtime = time(nullptr);
	TR << "Saved state: " << endtime - starttime << "s " << inbound().size() << " + " << outbound().size() << " + " <<  library_central_.size() << " + " << ( inbound().size() + outbound().size() + library_central_.size() ) << std::endl;
	start_timer( TIMING_CPU );
}

void
MPI_Refinement::save_state_auto(){
	if ( (core::Size)(last_save_state_ + save_state_interval_ ) < (core::Size)time(nullptr) ) {
		TR << "Saving state.. " << std::endl;
		last_save_state_ = time(nullptr);
		save_state( ident_string_ );
	}
}

void
MPI_Refinement::load_state(std::string prefix ){
	start_timer( TIMING_IO_READ );
	inbound().clear();
	outbound().clear();
	read_queues_from_file( prefix + "." + ObjexxFCL::string_of(mpi_rank()) );
	library_central_.clear();
	library_central_.read_from_file( prefix +  "." + ObjexxFCL::string_of(mpi_rank()) + ".lib.library_central" );
	start_timer( TIMING_CPU );
}

void
MPI_Refinement::print_stats(){
	static int lasttime = 0;
	if ( (time(nullptr) - lasttime) < 300 ) return;
	lasttime = time(nullptr);

	TR.Debug << "STATL: "
		<< wall_time() << "s  "
		<< total_structures_ << "  "
		<< total_structures_relax_
		<< " Acc: " << F(5,3, core::Real(total_metropolis_accepts_)/ core::Real(total_metropolis_) )
		<< " CPU: " << int((totaltime_batchrelax_+totaltime_loophash_)/3600) << " hrs  "
		<< " r/l: " << F(5,2, float(totaltime_batchrelax_)/(float(totaltime_loophash_)+0.1))
		<< " LHav: " << int(totaltime_loophash_/(n_loophash_+1)) << "s "
		<< " BRav: " << int(totaltime_batchrelax_/(n_batchrelax_+1)) << "s "
		<< " MEM: " << int(library_central_.mem_footprint()/1024) << " kB "  << std::endl;
}

bool
MPI_Refinement::add_structure_to_library( core::io::silent::SilentStructOP ss, std::string add_algorithm ){
	// reset the nuse to 0
	ss->add_energy( "nuse", 0 );
	//ss->add_energy( "ltime", time(NULL) );

	bool result = false;

	// Fill in missing scores for input structure
	fobj_->add_objective_function_info( ss, library_central() );

	// Add native quality info before adding
	//if( native_given_ ) add_nativeinfo_to_ss( ss, native_pose_ );

	// default algorithm is set by class wide parameter
	if ( add_algorithm == "" ) add_algorithm = mpi_feedback_;

	if     ( add_algorithm == "no" )             result = false;
	else if ( add_algorithm == "add_n_limit" )    result = add_structure_to_library_direct( *ss );
	else if ( add_algorithm == "add_n_replace" )  result = add_structure_to_library_add_n_replace( *ss );
	else if ( add_algorithm == "single_replace" ) result = add_structure_to_library_single_replace( *ss );

	else {
		utility_exit_with_message( "FATAL ERROR:  Unknown adding algorithm: '" + add_algorithm + "'" );
	}

	// Update mean/stdev after every replacement
	// this will make writing very slow!
	//if( result ) fobj_->add_objective_function_info( library_central_.store() );

	return result;
}

bool
MPI_Refinement::add_structure_to_library_direct( core::io::silent::SilentStruct &pss )
{
	TR.Debug << "Add: Direct addition called." << std::endl;
	library_central_.add( pss );
	return true;
}

bool
MPI_Refinement::add_structure_to_library_add_n_replace( core::io::silent::SilentStruct &pss )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// just hard coded for now
	std::string const dist_measure( "Sscore" );

	core::Real new_struct_score = fobj_->get_fobj( pss, sim_replace_obj_ );
	//pss.add_energy( "nuse", 0 );
	start_timer( TIMING_CPU );
	TRDEBUG << "Checking for zero size .. " << std::endl;

	if ( library_central_.size() == 0 ) {
		add_structure_to_library_direct( pss );
		return true;
	}

	library_central_.sort_by( sim_replace_obj_ );
	// if energy is wrse then worst member - ignore structure
	if ( new_struct_score > fobj_->get_fobj( *library_central_.store().back(), 1 ) ) {
		TR << "Add: NReplace, Ignoring struc (E vs max(Elib)): " << new_struct_score << " > ";
		TR << fobj_->formatted_objs_values( *library_central_.store().back() );
		TR << std::endl;
		return false;
	}

	// now find the closest DIST
	core::Real closest_dist = 1000000;
	//core::io::silent::SilentStructOP closest_struct;
	SilentStructStore::iterator closest_struct;

	//for( core::Size i = 0; i < library_central_.size(); ++i ){
	//core::io::silent::SilentStructOP ss = library_central_.get_struct( i );
	for ( auto jt =  library_central_.begin(),
			end = library_central_.end(); jt != end; ++jt ) {

		core::Real the_dist( 0.0 );
		core::Real dumm;
		// slow down?
		core::io::silent::SilentStructOP ss( pss.clone() );
		if ( dist_measure.compare("rmsd") == 0 ) {
			dumm = 1.0 - CA_Sscore( ss, *jt, the_dist );
		} else if ( dist_measure.compare("Sscore") == 0 ) {
			the_dist = 1.0 - CA_Sscore( ss, *jt, dumm );
		} else if ( dist_measure.compare("looprmsd") == 0 ) {}

		TRDEBUG << "The dist: " << the_dist << std::endl;
		if ( the_dist < closest_dist ) {
			closest_dist = the_dist;
			closest_struct = jt;
		}
	}

	// we need to add one more option for this...
	bool close_to_existing( false );
	if ( dist_measure.compare("Sscore") == 0 ) {
		if ( closest_dist < 0.05 ) close_to_existing = true;

	} else if ( dist_measure.compare("rmsd") == 0 ||
			dist_measure.compare("looprmsd") == 0 ) {
		if ( closest_dist < 0.5 ) close_to_existing = true;

	}

	if ( close_to_existing ) {
		// replace if lower in energy
		core::Real energy_old = fobj_->get_fobj( *(*closest_struct), sim_replace_obj_ );

		if ( new_struct_score < energy_old ) {
			TR << "Add: NReplace ACC, byE: " << format_silent_struct( pss ) << "|";
			TR << " " << F(8,3, new_struct_score );
			TR << " " << F(8,3, energy_old );
			TR << " " << F(8,3, new_struct_score - energy_old);
			TR << std::endl;
			library_central_.erase( closest_struct );
			library_central_.add( pss );
			return true;

		} else {
			TR << "Add: NReplace REJ, byE: " << format_silent_struct( pss ) << "|";
			TR << " " << F(8,3, new_struct_score );
			TR << " " << F(8,3, energy_old );
			TR << " " << F(8,3, new_struct_score - energy_old);
			TR << std::endl;
			return false;
		}

	} else {
		// add
		library_central_.add( pss );
		TR << "Add: NReplace Acc, new: " << format_silent_struct( pss ) << std::endl;
		return true;
	}

	// Nothing can reach here, but let's just return something :)
	return false;
}

bool
MPI_Refinement::add_structure_to_library_single_replace( core::io::silent::SilentStruct &pss )
{

	core::Real ssid = pss.get_energy("ssid");

	// now find the library structure with the same ssid
	find_SilentStructOPs predic("ssid", ssid);
	auto ssid_match = std::find_if( library_central_.begin(), library_central_.end(), predic );

	if ( ssid_match == library_central_.end() ) {
		TR << "Add: SingleReplace, ssid expired: " + ObjexxFCL::string_of( ssid ) + " - rejected structure" << std::endl;
		return false;
	}

	bool replace_it = false;
	// to get library replacement you must either: be of a more advanced round or have lower energy
	core::Size dround( pss.get_energy("round") - (*ssid_match)->get_energy("round") );
	core::Real dobj( fobj_->get_fobj( pss, 1 ) - fobj_->get_fobj( *ssid_match, 1 ) );
	if ( dround > 1 || dobj < 0.0 ) {
		replace_it = true;
		//TR << "SingleReplace Filter Pass on " << ssid << ", dround/dobj: " << dround << " " << dobj << std::endl;
	} else {
		TR.Debug << "Add: SingleReplace Filter Fail on " << ssid << ", dround/dobj(vsMax): " << dround;
		TR.Debug << " " << F(8,3,dobj);
		TR.Debug << std::endl;
	}

	if ( replace_it ) {
		core::Real new_energy = fobj_->get_fobj( pss, 1 );
		core::Real old_energy = fobj_->get_fobj( *ssid_match, 1 );

		bool metropolis_replace = false;

		core::Real energy_diff_T = 0;
		if ( mpi_metropolis_temp_ > 0.0 ) energy_diff_T = old_energy - new_energy;

		if ( ( energy_diff_T >= 0.0 ) ) metropolis_replace = true; // energy of new is simply lower
		else if ( energy_diff_T/mpi_metropolis_temp_ > (-10.0) ) { // exp(-10) ~ 0, so check this before trying to evaluate exp(-100) and getting NaN
			core::Real random_float = numeric::random::rg().uniform();
			if ( random_float < exp( energy_diff_T/mpi_metropolis_temp_ ) )  metropolis_replace = true;
		}
		total_metropolis_++;
		if ( metropolis_replace ) {
			total_metropolis_accepts_++;
			TR << "Add: SingleReplace ACC: " << format_silent_struct( *ssid_match) << " with " << format_silent_struct( pss ) << "  " << -energy_diff_T << std::endl;
			*ssid_match = pss.get_self_ptr();
			return true;
		} else {
			TR << "Add: SingleReplace REJ: " << format_silent_struct( *ssid_match) << " with " << format_silent_struct( pss ) << "  " << -energy_diff_T << std::endl;
		}
	}

	return false;
}

void
MPI_Refinement::print_summary( std::string const prefix )
{
	//Should be sorted by fobj in order to show min/max status
	core::io::silent::SilentStructOP minss, maxss;

	protocols::wum::SilentStructStore lib_loc( library_central_);
	TR << prefix;
	for ( core::Size iobj = 1; iobj <= fobj_->nobjs(); ++iobj ) {
		TR << "/min_" << fobj_->fobjnames( iobj );
	}
	TR << "/Time(min): ";

	for ( core::Size iobj = 1; iobj <= fobj_->nobjs(); ++iobj ) {
		std::string objname = fobj_->fobjnames( iobj );
		lib_loc.sort_by( objname );
		TR << " " << F(8,2, lib_loc.get_struct( 0 )->get_energy( objname ) );
	}

	core::Real endtime = time(nullptr);
	core::Real dt_in_min = (endtime - starttime_)/60.0;
	TR << " " << F(6,1,dt_in_min);
	TR << std::endl;
}

void
MPI_Refinement::print_library( protocols::wum::SilentStructStore &library,
	std::string const prefix )
{
	TR << prefix << " N:" << library.store().size() << "        tag";
	TR << fobj_->formatted_objs_names();

	TR << "   Front Round Nuse ";
	if ( native_given_ ) TR << " GDTTM    GDTHA";
	TR << "  IGDTHA" << std::endl;

	for ( SilentStructStore::const_iterator it = library.begin(); it != library.end(); ++ it ) {
		TR << prefix << ": " << format_silent_struct( *it ) << std::endl;
	}
}

// iterate through the structure store and add strucutres to the central library according to the algorithm specified or the default algorithm
bool
MPI_Refinement::add_structures_to_library( SilentStructStore &new_structs, std::string add_algorithm )
{
	bool result = false;

	if ( add_algorithm.compare("NSGAII") == 0 ) {
		result = fobj_->update_library_NSGAII( library_central_, new_structs, max_lib_size() );

	} else {
		// First, clear marks on old members
		//for( SilentStructStore::iterator jt = library_central_.begin();
		//   jt != library_central_.end(); jt ++ )
		//(*jt)->add_string_value( "added", "" );

		for ( SilentStructStore::const_iterator it = new_structs.begin();
				it != new_structs.end(); ++it ) {
			runtime_assert( *it );
			//(*it)->add_string_value( "added", "o" );
			bool local_result = add_structure_to_library( *it, add_algorithm );
			result = result || local_result;
		}

		if ( result ) limit_library();
	}

	return result;
}

void
MPI_Refinement::limit_library()
{
	// now shave off the worst structres
	// By default this is done by "score" - could be changed to other Objective functions
	// so as to keep diversity
	library_central_.sort_by( sim_replace_obj_ );
	core::Size npop( 0 );
	while ( library_central_.size() > max_lib_size_ ) {
		npop++;
		core::io::silent::SilentStructOP ss = library_central_.get_struct( library_central_.size() - 1 );
		TR << "Limit: " << format_silent_struct( ss ) << std::endl;
		library_central_.store().pop_back();
	}
	if ( npop > 0 ) {
		TR << "Limit library: " << npop << " structures removed from library" << std::endl;
	}
}


void
MPI_Refinement::shave_library( SilentStructStore &new_structs,
	std::string const scorename,
	core::Real const frac ) const
{
	new_structs.sort_by( scorename );
	core::Size ntot = new_structs.size();
	core::Size n_to_pop = core::Size( new_structs.size()*frac );
	core::Size npop( 0 );

	while ( npop < n_to_pop ) {
		new_structs.store().pop_back();
		npop++;
	}
	TR << "Shave store: " << npop << " out of " << ntot << ", returning " << new_structs.size() << std::endl;
}

void
MPI_Refinement::dump_structures( const SilentStructStore &new_structs,
	bool score_only,
	std::string prefix ) const
{
	start_timer( TIMING_IO_WRITE  );
	core::io::silent::SilentFileOptions silent_options;
	core::io::silent::SilentFileData sfd( silent_options );
	std::string filename = jobname_ + "." + prefix + ObjexxFCL::string_of( mpi_rank() ) + ".out";

	core::Size istr( 0 );
	for ( auto const & new_struct : new_structs ) {
		istr++;
		sfd.write_silent_struct( *new_struct, filename, score_only );
	}
	core::Real write_time = start_timer( TIMING_CPU );
	TRDEBUG << "Write time: " << write_time << std::endl;
}

void
MPI_Refinement::send_sortedpick_library_structs( core::Size dest_rank,
	core::Size nsend,
	std::string const scorename,
	bool const weighted,
	bool const ,
	core::Real const kT )
{
	if ( library_central_.size() == 0 ) {
		TR.Error << "Have no structure to send" << std::endl;
		return;
	}

	// Get copy of library_central
	protocols::wum::SilentStructStore structs_loc( library_central_ );
	structs_loc.sort_by( scorename );
	core::Size const n = structs_loc.size();

	// get weighted probability based on scores ( exp(-(E - Emin)/(Emax-Emin) ))
	core::Real const Emax = structs_loc.get_struct( n - 1 )->get_energy( scorename );
	core::Real const Emin = structs_loc.get_struct( 0 )->get_energy( scorename );
	core::Real denom = kT*(Emax - Emin);
	if ( denom == 0.0 ) denom = 1.0;

	utility::vector0< core::Real > accprob( n, 0.0 ); //accummulative probability
	core::Real psum( 0.0 );
	core::Size i( 0 );

	for ( protocols::wum::SilentStructStore::const_iterator it = structs_loc.begin();
			it != structs_loc.end(); ++it, ++i ) {
		core::Real score = (*it)->get_energy( scorename );
		core::Real prob;
		if ( weighted ) {
			prob = std::exp( -score/denom );
			if ( prob <= 0.0000001 ) prob = 0.0000001;
		} else {
			if ( i <= nsend ) {
				prob = 1.0;
			} else {
				prob = 0.0;
			}
		}
		psum += prob;
		accprob[i] = psum;
	}

	// normalize
	bool assign_fail( false );
	for ( core::Size j = 0; j < n; ++j ) {
		accprob[j] /= psum;
		if ( !(accprob[j] > 0.0) || !(std::abs(accprob[j] ) <= 1.0 ) ||
				accprob[j] != accprob[j] ) {
			assign_fail = true;
		}
	}

	if ( assign_fail ) {
		psum = 0.0;
		for ( core::Size j = 0; j < n; ++j ) {
			accprob[j] = (core::Real)(j+1);
			psum += (core::Real)(j+1);
		}

		for ( core::Size j = 0; j < n; ++j ) accprob[j] /= psum;
	}

	// report
	TR << "Accum. prob. assigned: ";
	for ( core::Size j = 0; j < n; ++j ) TR << " " << F(6,3,accprob[j]);
	TR << std::endl;

	WorkUnit_SilentStructStoreOP resultfeedback( new WorkUnit_SilentStructStore( ) );
	resultfeedback->set_wu_type( "resultfeedback" );

	// not to allow duplication
	utility::vector1< core::Size > picked;
	TR << "Random number picked: ";
	for ( core::Size ipick = 1; ipick <= nsend; ++ipick ) {
		core::Size k( 0 );
		core::Size iter( 0 );
		while ( true ) {
			core::Real rnum = numeric::random::rg().uniform();
			iter++;

			for ( core::Size j = 0; j < n; ++j ) {
				if ( !(accprob[j] <= rnum) ) {
					k = j;
					break;
				}
			}

			if ( !picked.contains( k ) ) {
				TR << " " << F(6,4,rnum);
				picked.push_back( k );
				break;
			}
		}

		core::io::silent::SilentStructOP new_struct = structs_loc.get_struct( k )->clone();
		core::Size nuse = new_struct->get_energy("nuse")+1;
		new_struct->add_energy( "nuse", nuse );
		resultfeedback->decoys().add( new_struct );
	}
	TR << std::endl;

	send_MPI_workunit( resultfeedback, dest_rank );
}

void
MPI_Refinement::send_random_library_structs( core::Size dest_rank, core::Size nsend ){
	if ( library_central_.size() == 0 ) {
		TR.Error << "Have no structure to send" << std::endl;
		return;
	}

	WorkUnit_SilentStructStoreOP resultfeedback( new WorkUnit_SilentStructStore( ) );
	resultfeedback->set_wu_type( "resultfeedback" );

	// not to allow duplication
	utility::vector1< core::Size > picked;
	for ( core::Size i = 1; i <= nsend; ++i ) {

		core::Size irand( 0 );
		while ( true ) {
			irand = (core::Size)(numeric::random::rg().uniform() * library_central_.size());
			if ( irand >= library_central_.size() ) irand = library_central_.size() - 1;

			if ( nsend >= library_central_.size() ) {
				break; // just send whatever if pool is smaller
			} else if ( !picked.contains( irand ) ) {
				picked.push_back( irand );
				break;
			}
		}

		core::io::silent::SilentStructOP new_struct = library_central_.get_struct( irand )->clone();

		// allowing duplication
		//core::io::silent::SilentStructOP new_struct = library_central_.get_struct_random()->clone();
		// Add nuse
		core::Size nuse = new_struct->get_energy("nuse")+1;
		library_central_.get_struct( irand )->add_energy( "nuse", nuse );
		resultfeedback->decoys().add( new_struct );
	}

	send_MPI_workunit( resultfeedback, dest_rank );
}

void
MPI_Refinement::send_random_library_struct( core::Size dest_rank, core::Size newssid ) const {
	if ( library_central_.size() == 0 ) {
		TR.Error << "Have no structure to send" << std::endl;
		return;
	}

	// now fabricate a return WU
	WorkUnit_SilentStructStoreOP resultfeedback( new WorkUnit_SilentStructStore( ) );
	resultfeedback->set_wu_type( "resultfeedback" );
	core::io::silent::SilentStructOP new_struct = library_central_.get_struct_random()->clone();
	new_struct->add_energy("ssid", newssid);   // overwrite the ssid
	resultfeedback->decoys().add( new_struct );
	send_MPI_workunit( resultfeedback, dest_rank );
}

core::Real
MPI_Refinement::score( const core::io::silent::SilentStruct &ss ) const{
	// by deafult, the rosetta energy is the objective function
	return ss.get_energy("score");
}

void
MPI_Refinement::retag_library( protocols::wum::SilentStructStore &store,
	std::string const prefix ) const {

	core::Size i( 0 );
	for ( SilentStructStore::const_iterator it = store.begin();
			it != store.end(); ++it, ++i ) {
		std::stringstream newname("");
		newname << prefix << "_" << i;
		TR << "Retag: " << std::setw(15) << (*it)->decoy_tag();
		TR << " -> " << std::setw(15) << newname.str() << std::endl;
		(*it)->set_decoy_tag( newname.str() );
	}
}

std::string
MPI_Refinement::format_silent_struct( const core::io::silent::SilentStruct &ss ) const {
	std::stringstream sstream;
	//sstream << "["  <<I(4, ss.get_energy("ssid") )
	//    << " "  << std::setw(3) << ss.get_energy("usid");

	sstream << "[" << ObjexxFCL::format::A(15, ss.decoy_tag() );

	sstream << fobj_->formatted_objs_values( ss );

	sstream << " | " << ObjexxFCL::format::I(3, ss.get_energy("frontier") )
		//<< " | " << ObjexxFCL::format::I(3, ss.get_energy("iter") )
		<< " | " << ObjexxFCL::format::I(3, ss.get_energy("round") )
		//<< " |"    << ObjexxFCL::format::I(5,   time(NULL) - ss.get_energy("ltime") )
		//<< " |" << ObjexxFCL::format::I(3, ss.get_energy("master") )
		<< " | " << ObjexxFCL::format::I(2, ss.get_energy("nuse"));

	if ( native_given_ ) {
		sstream << " | " << ObjexxFCL::format::F(6,2, ss.get_energy("gdttm"))
			<< " | " << ObjexxFCL::format::F(6,2, ss.get_energy("gdtha"));
	}
	sstream  << " | " << ObjexxFCL::format::F(6,2, ss.get_energy("gdtha_i"));

	sstream << " ]";
	//sstream << "  " << ss.get_string_value("added");

	return sstream.str();
}

core::Real
MPI_Refinement::score( const core::io::silent::SilentStructOP &ss ) const {
	return score( *ss );
}

std::string
MPI_Refinement::format_silent_struct( const core::io::silent::SilentStructOP &ss ) const {
	return format_silent_struct( *ss );
}

void
MPI_Refinement::report_time( ) const {
	core::Real endtime = time(nullptr);
	core::Real dt = endtime - starttime_;
	core::Real dm = int(dt/60.0);
	core::Real ds = dt - dm*60.0;
	TR << dm << " min" << ds << " sec." << std::endl;
}

} // namespace mpi_refinement
} // namespace protocols
