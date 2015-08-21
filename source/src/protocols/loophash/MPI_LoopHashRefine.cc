// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/loophash/MPI_LoopHashRefine.cc
/// @brief
/// @author Mike Tyka
#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/loophash/MPI_LoopHashRefine.hh>
#include <protocols/loophash/WorkUnit_LoopHash.hh>
#include <protocols/relax/WorkUnit_BatchRelax.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <core/pose/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <numeric/random/random.hh>

// C/C++ headers
#ifndef _WIN32 // REQUIRED FOR WINDOWS
#endif

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <utility/vector1.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace loophash {

using namespace protocols::wum;

static thread_local basic::Tracer TR( "MPI.LHR" );


MPI_LoopHashRefine::MPI_LoopHashRefine( char machine_letter ):
	MPI_WorkUnitManager( machine_letter ),
	max_lib_size_(50),
	save_state_interval_(300),
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
MPI_LoopHashRefine::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	max_lib_size_   = option[ OptionKeys::lh::max_lib_size ]();
	save_state_interval_        = option[ OptionKeys::lh::mpi_save_state_interval ]();
	mpi_feedback_               = option[ OptionKeys::lh::mpi_feedback ]();
	mpi_metropolis_temp_        = option[ OptionKeys::lh::mpi_metropolis_temp ]();
	rms_limit_                  = option[ OptionKeys::lh::rms_limit ]();
	objective_function_         = option[ OptionKeys::lh::objective_function ]();
	mpi_resume_ = "";
	if ( option[ OptionKeys::lh::mpi_resume ].user() ) {
		mpi_resume_ = option[ OptionKeys::lh::mpi_resume ]();
	}
	jobname_ = option[ OptionKeys::lh::jobname ]();

	// Make ident string:
	// Make a medley of a string and the PID
	ident_string_ =  option[ lh::jobname ];
	TR << "IDENT: " << ident_string_ << std::endl;

	// make sure the state saves are randomly staggered - otherwise all the masters dump several hundred megs at once!
	last_save_state_ = time(NULL)  + core::Size( numeric::random::rg().uniform()  * (core::Real) save_state_interval_) ;
	TR << "Interlace dumps: " << last_save_state_ << "  " << time(NULL)  << "  " << last_save_state_ - time(NULL) << "  " << save_state_interval_ << "  " << std::endl;
}

void
MPI_LoopHashRefine::load_structures_from_cmdline_into_library( core::Size structure_read_offset )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	start_timer( TIMING_IO_READ  );
	TR << "Reading in structures..." << std::endl;

	core::chemical::ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

	core::Size count = 0;
	SilentStructStore temp_lib;


	// Rescore input structures for fair comparison with whatever comes back from relax calls
	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::get_score_function();
	core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn );


	while ( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );

		if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
			core::scoring::constraints::ConstraintSetOP cstset( core::scoring::constraints::ConstraintIO::get_instance()->read_constraints( core::scoring::constraints::get_cst_fa_file_option(), core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet ), pose  ) );
			TR << "Read constraints: " << ( cstset->has_constraints() ? "YES" : "NONE" ) << std::endl;
			pose.constraint_set( cstset );
		}

		(*scorefxn)(pose);

		core::pose::set_ss_from_phipsi( pose );
		core::io::silent::SilentStructOP ss = option[ OptionKeys::lh::bss]() ?
			core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary") :
			core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		ss->fill_struct( pose, "empty_tag" );
		ss->add_energy( "censcore", 0 );   // centroidscore of last time this structure was in centroid form
		ss->add_energy( "extra_score", 0 );   // extra_score is used if a different scoring scheme is used to select structures then to actually repack/relax them
		ss->add_energy( "combined_score", 0 );   // extra_score + score
		ss->add_energy( "parent_score", 0 );  // fa score of parent
		ss->add_energy( "lhcount", 0 );   // lh count is the number of times this structure has been used to generate new loophashed structures
		ss->add_energy( "round", 0 );     // round is the number of consecutive loophashes this structure has gone through
		ss->add_energy( "expire", 0 );    // how many times has this structure been retired to the temporor ? >1 means it's been sent out again
		if (  ss->get_string_value("usid") == "" ) {
			ss->add_string_value( "usid", "_" + string_of(count)); // A unique structure id (usid) that will identify the structure until it is modified.
		}
		ss->add_string_value( "husid", ss->get_string_value("usid") ); // history of usids
		ss->add_energy( "state", 0 );     // state: 0 init, 1 loophashed, 2 relaxed
		ss->add_energy( "ltime", time(NULL) ); // time when it became active. This is used to expire structures that havn't changed for a while.
		ss->add_energy( "master", mpi_rank() ); // what master is it curretnly on ? (volatile)
		ss->add_energy( "emperor_count", 0 ); // what master is it curretnly on ? (volatile)
		temp_lib.add( ss );
		count ++;
	}

	if ( temp_lib.size() == 0 ) {
		TR.Error << "Error reading starting structures: 0 valid structures read. Check your input files and parameters" << std::endl;
		utility_exit_with_message( "Error reading starting structures: 0 valid structures read. Check your input files and parameters" );
	}

	TR << "Loaded " << temp_lib.size() << " starting structures" << std::endl;

	core::Size position = structure_read_offset % temp_lib.size();
	count=1;
	while ( library_central_.size() < max_lib_size_ ) {
		core::io::silent::SilentStructOP ss;
		ss = temp_lib.get_struct( position )->clone();
		ss->add_energy( "ssid", count );
		library_central_.add( ss );
		position++;
		if ( position>=temp_lib.size() ) position = 0;
		count++;
	}

	TR << "Added " << library_central_.size() << " structures to library " << std::endl;

	runtime_assert(  library_central_.size() == max_lib_size_ );
	start_timer( TIMING_CPU );
}

void
MPI_LoopHashRefine::save_state(std::string prefix ){
	start_timer( TIMING_IO_WRITE );
	long starttime = time(NULL);
	write_queues_to_file( prefix + "." + string_of(mpi_rank()) );
	library_central_.serialize_to_file( prefix + "." + string_of(mpi_rank())+ ".lib.library_central" );
	long endtime = time(NULL);
	TR << "Saved state: " << endtime - starttime << "s " << inbound().size() << " + " << outbound().size() << " + " <<  library_central_.size() << " + " << ( inbound().size() + outbound().size() + library_central_.size() ) << std::endl;
	start_timer( TIMING_CPU );
}

void
MPI_LoopHashRefine::save_state_auto(){
	if ( (core::Size)(last_save_state_ + save_state_interval_ ) < (core::Size)time(NULL) ) {
		TR << "Saving state.. " << std::endl;
		last_save_state_ = time(NULL);
		save_state( ident_string_ );
	}
}

void
MPI_LoopHashRefine::load_state(std::string prefix ){
	start_timer( TIMING_IO_READ );
	inbound().clear();
	outbound().clear();
	read_queues_from_file( prefix + "." + string_of(mpi_rank()) );
	library_central_.clear();
	library_central_.read_from_file( prefix +  "." + string_of(mpi_rank()) + ".lib.library_central" );
	start_timer( TIMING_CPU );
}

void
MPI_LoopHashRefine::print_stats(){
	static int lasttime = 0;
	if ( (time(NULL) - lasttime) < 300 ) return;
	lasttime = time(NULL);

	TR << "STATL: "
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
MPI_LoopHashRefine::add_structure_to_library( core::io::silent::SilentStruct &pss, std::string add_algorithm ){
	// reset the lhcount to 0
	pss.add_energy( "lhcount", 0 );
	pss.add_energy( "ltime", time(NULL) );

	bool result = false;

	// default algorithm is set by class wide parameter
	if ( add_algorithm == "" ) add_algorithm = mpi_feedback_;
	if ( add_algorithm == "no" )             result = false;
	else if ( add_algorithm == "add_n_limit" )    result = add_structure_to_library_direct( pss );
	else if ( add_algorithm == "add_n_replace" )  result = add_structure_to_library_add_n_replace( pss );
	else if ( add_algorithm == "single_replace" ) result = add_structure_to_library_single_replace( pss );
	else {
		utility_exit_with_message( "FATAL ERROR:  Unknown adding algorithm: '" + add_algorithm + "'" );
	}
	return result;
}

bool
MPI_LoopHashRefine::add_structure_to_library_direct( core::io::silent::SilentStruct &pss )
{
	library_central_.add( pss );
	return true;
}

bool
MPI_LoopHashRefine::add_structure_to_library_add_n_replace( core::io::silent::SilentStruct &pss )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Real new_struct_score = objective_function(pss);
	pss.add_energy( "lhcount", 0 );
	start_timer( TIMING_CPU );
	TRDEBUG << "Checking for zero size .. " << std::endl;
	if ( library_central_.size() == 0 ) {

		add_structure_to_library_direct( pss );
		return true;
	}

	library_central_.sort_by();
	// if energy is wrse then worst member - ignore structure
	if ( new_struct_score > objective_function(library_central_.store().back()) ) {
		TR << "Ignoring struc: " << new_struct_score << " > " << format_silent_struct(library_central_.store().back()) << std::endl;
		return false;
	}

	// now find the closest RMS
	core::Real closest_rms = 1000000;
	SilentStructStore::iterator closest_struct;

	for ( SilentStructStore::iterator jt =  library_central_.begin(),
			end = library_central_.end(); jt != end; ++jt ) {
		core::Real the_rms;
		// downcast
		if ( option[ OptionKeys::lh::bss]() ) {
			core::io::silent::BinarySilentStruct *jt_pss = dynamic_cast < core::io::silent::BinarySilentStruct * > ( &(*(*jt)) );
			core::io::silent::BinarySilentStruct *ss = dynamic_cast < core::io::silent::BinarySilentStruct * > ( &(pss) );
			if ( jt_pss == NULL || ss == NULL ) utility_exit_with_message( "FATAL ERROR:  This code only runs with Binary Protein SilentStructs " );
			the_rms = ss->CA_rmsd( *jt_pss );
		} else {
			core::io::silent::ProteinSilentStruct *jt_pss = dynamic_cast < core::io::silent::ProteinSilentStruct * > ( &(*(*jt)) );
			core::io::silent::ProteinSilentStruct *ss = dynamic_cast < core::io::silent::ProteinSilentStruct * > ( &(pss) );
			if ( jt_pss == NULL || ss == NULL ) utility_exit_with_message( "FATAL ERROR:  This code only runs with Protein SilentStructs " );
			the_rms = ss->CA_rmsd( *jt_pss );
		}
		TRDEBUG << "The rms: " << the_rms << std::endl;
		if ( the_rms < closest_rms ) {
			//TR << "Found better: " << the_rms;
			closest_rms = the_rms;
			closest_struct = jt;
		}
	}

	core::Real rms_time = start_timer( TIMING_CPU );
	TR << "RMS_Time: " << rms_time << std::endl;

	if ( closest_rms < rms_limit_ ) {
		// replace if lower in energy
		core::Real energy_old = objective_function(*closest_struct);
		TR << "Enew vd Eold " << new_struct_score << " " << energy_old << std::endl;
		if ( new_struct_score < energy_old ) {
			pss.add_energy( "lhcount", 0 );
			*(*closest_struct) = pss;
			return true;
		}
	} else {
		// add
		pss.add_energy( "lhcount", 0 );
		library_central_.add( pss );
		library_central_.sort_by();
		return true;
	}
	TR << "No structure added" << std::endl;
	return false;
}

bool
MPI_LoopHashRefine::add_structure_to_library_single_replace( core::io::silent::SilentStruct &pss )
{
	core::Real ssid = pss.get_energy("ssid");

	// now find the library structure with the same ssid
	find_SilentStructOPs predic("ssid", ssid);
	SilentStructStore::iterator ssid_match = std::find_if( library_central_.begin(), library_central_.end(), predic );

	if ( ssid_match == library_central_.end() ) {
		TR << "ssid expired: " + ObjexxFCL::string_of( ssid ) + " - rejected structure" << std::endl;
		print_library();
		return false;
	}


	bool replace_it = false;
	// to get library replacement you must either: be of a more advanced round or have lower energy
	if ( pss.get_energy("round") > (*ssid_match)->get_energy("round") ) replace_it = true;
	else {
		if ( objective_function(pss) < objective_function(*ssid_match) ) replace_it = true;
	}

	if ( replace_it ) {
		core::Real new_energy = objective_function(pss);
		core::Real old_energy = objective_function(*ssid_match);

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
			TR << "ReplacingACC: " << format_silent_struct( *ssid_match) << " with " << format_silent_struct( pss ) << "  " << energy_diff_T << std::endl;
			*ssid_match = pss.get_self_ptr();
			return true;
		} else {
			TR << "ReplacingREJ: " << format_silent_struct( *ssid_match) << " with " << format_silent_struct( pss ) << "  " << energy_diff_T << std::endl;
		}
	}

	return false;
}

void
MPI_LoopHashRefine::print_library()
{
	TR << "N:" << library_central_.store().size() << "  <ssid>  <score> <censc>  <rms>  <round>  <time>  <master>  <lhcount> ----" << std::endl;
	for ( SilentStructStore::const_iterator it = library_central_.begin(); it != library_central_.end(); ++ it ) {
		TR << "LIB: " << format_silent_struct( *it ) << std::endl;
	}
}

// iterate through the structure store and add strucutres to the central library according to the algorithm specified or the default algorithm
bool
MPI_LoopHashRefine::add_structures_to_library( SilentStructStore &new_structs, std::string add_algorithm )
{
	bool result = false;

	for ( SilentStructStore::const_iterator it = new_structs.begin();
			it != new_structs.end(); ++it ) {
		runtime_assert( *it != 0 );
		TR << "Add structure... " << format_silent_struct( *it ) << std::endl;
		bool local_result = add_structure_to_library( *(*it), add_algorithm );
		result = result || local_result;
	}

	limit_library();
	print_library();
	return result;
}

void
MPI_LoopHashRefine::limit_library()
{
	// now shave off the worst structres
	library_central_.sort_by();
	while ( library_central_.size() > max_lib_size_ ) {
		library_central_.store().pop_back();
	}
}

void
MPI_LoopHashRefine::dump_structures( const SilentStructStore &new_structs, bool score_only ) const {
	start_timer( TIMING_IO_WRITE  );
	core::io::silent::SilentFileData sfd;
	std::string filename = jobname_ + "." + string_of( mpi_rank() ) + ".out";
	for ( SilentStructStore::const_iterator it = new_structs.begin();
			it != new_structs.end(); ++it ) {
		// don't add round 0 structures (they have conflicting silent struct headers and
		// mess up columns in the output
		if ( (*it)->get_energy("round") > 0 ) {
			sfd.write_silent_struct( *(*it), filename, score_only );
		} else {
			TR << "Trying to dump structure of round == 0. Refusing plainly. " << std::endl;
		}
	}
	core::Real write_time = start_timer( TIMING_CPU );
	TR << "Write time: " << write_time << std::endl;
}

void
MPI_LoopHashRefine::send_random_library_struct( core::Size dest_rank, core::Size newssid ) const {
	if ( library_central_.size() == 0 ) {
		TR.Error << "ERROR: Havce no structure to send" << std::endl;
		return;
	}

	// now fabricate a return WU
	WorkUnit_SilentStructStoreOP resultpack( new WorkUnit_SilentStructStore( ) );
	resultpack->set_wu_type( "resultpack" );
	core::io::silent::SilentStructOP new_struct = library_central_.get_struct_random()->clone();
	new_struct->add_energy("ssid", newssid);   // overwrite the ssid
	resultpack->decoys().add( new_struct );
	send_MPI_workunit( resultpack, dest_rank );
}

core::Real
MPI_LoopHashRefine::objective_function( const core::io::silent::SilentStruct &ss ) const{
	// by deafult, the rosetta energy is the objective function
	return ss.get_energy(objective_function_);
}

core::Real
MPI_LoopHashRefine::score( const core::io::silent::SilentStruct &ss ) const{
	// by deafult, the rosetta energy is the objective function
	return ss.get_energy("score");
}

std::string
MPI_LoopHashRefine::format_silent_struct( const core::io::silent::SilentStruct &ss ) const {
	std::stringstream sstream;
	sstream << "["<<I(4,    ss.get_energy("ssid") )
		<< " "     <<        ss.get_string_value("usid")
		<< " |"    << F(8,1, objective_function( ss ) )
		<< " |"    << F(8,1, ss.get_energy("score"))
		<< " |"    << F(8,1, ss.get_energy("censcore"))
		<< " |"    << F(5,1, ss.get_energy("rms") )
		<< " |"    << I(3,   ss.get_energy("round") )
		<< " |"    << I(5,   time(NULL) - ss.get_energy("ltime") )
		<< " |"    << I(3,   ss.get_energy("master") )
		<< " |"    << I(1,   ss.get_energy("lhcount")) << "]";
	return sstream.str();
}

core::Real
MPI_LoopHashRefine::objective_function( const core::io::silent::SilentStructOP &ss ) const {
	return objective_function( *ss );
}

core::Real
MPI_LoopHashRefine::score( const core::io::silent::SilentStructOP &ss ) const {
	return score( *ss );
}

std::string
MPI_LoopHashRefine::format_silent_struct( const core::io::silent::SilentStructOP &ss ) const {
	return format_silent_struct( *ss );
}

} // namespace loophash
} // namespace protocols
