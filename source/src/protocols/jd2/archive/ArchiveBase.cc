// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/jd2/MPIFileBufJobDistributor.cc
/// @brief  implementation of MPIFileBufJobDistributor
/// @author Oliver Lange olange@u.washington.edu

// Unit headers
#include <protocols/jd2/archive/ArchiveBase.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>

#include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>

// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/io/raw_data/DisulfideFile.hh>
// AUTO-REMOVED #include <core/scoring/ResidualDipolarCoupling.hh>



// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>

#include <core/chemical/ChemicalManager.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

//C++
// AUTO-REMOVED #include <iterator>

static basic::Tracer tr("protocols.jd2.Archive");
using basic::mem_tr;
// Utility headers
#include <basic/options/option_macros.hh>
#include <basic/prof.hh>

#include <utility/vector1.hh>
#include <boost/algorithm/string/erase.hpp>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


//
OPT_1GRP_KEY( Integer, iterative, pool_size )

bool protocols::jd2::archive::ArchiveBase::options_registered_( false );

using namespace basic::options;
using namespace basic::options::OptionKeys;
//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::jd2::archive::ArchiveBase::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		NEW_OPT( iterative::pool_size, "number of structures maintained in pool", 500 );
	}
}


namespace protocols {
namespace jd2 {
namespace archive {

/// @details Auto-generated virtual destructor
AbstractArchiveBase::~AbstractArchiveBase() {}
// using namespace basic::options;
// using namespace basic::options::OptionKeys;
using namespace core;

std::string const ArchiveBase::TAG_IN_FILE( "tag_in_file" );
std::string const ArchiveBase::SOURCE_FILE( "source_file" );

ArchiveBase::ArchiveBase( ArchiveManagerAP ptr ) :
	AbstractArchiveBase( ptr ),
	max_nstruct_( option[ iterative::pool_size ]() ),
	accepts_since_last_batch_( 0 ),
	total_accepts_( 0 ),
	proposed_since_last_batch_( 0 ),
	total_proposed_( 0 )
{

	runtime_assert( options_registered_ );

	//for debugging --- make archive smaller
	if ( ( option[ run::test_cycles ]() || option[ run::dry_run ]() ) && !option[ run::memory_test_cycles ]() ) {
		max_nstruct_ = 20;
	}

	nstruct_ = max_nstruct_;

	floating_acceptance_ratio_ = 1.0;
	min_structures_for_acceptance_statistics_ = max_nstruct_;
}

ArchiveBase::~ArchiveBase() {}

///@brief count the structure for the acceptance statistics
///   only count if not from expired batch
void ArchiveBase::count_structure( Batch const& batch, bool accepted ) {
	if ( still_interested( batch ) ) {
		accepts_since_last_batch_ += accepted ? 1 : 0;
		++proposed_since_last_batch_;

		//count with floating-average
		core::Real const inv_A( 1.0/min_structures_for_acceptance_statistics_ );
		if ( acceptance_history_.size() > min_structures_for_acceptance_statistics_ ) {
			if ( acceptance_history_.front() ) floating_acceptance_ratio_ -= inv_A;
			acceptance_history_.pop_front();
		} else { //assume before history we had full-acceptance ( in line with initializing floating_acc_rat = 1 ).
			floating_acceptance_ratio_ -= inv_A;
		}
		if ( accepted ) floating_acceptance_ratio_ += inv_A;
		acceptance_history_.push_back( accepted );
	}
}

///@brief count the structure for the acceptance statistics
///   only count if not from expired batch
void ArchiveBase::count_removed_structures( core::Size n_removed ) {
	if ( accepts_since_last_batch() > n_removed ) accepts_since_last_batch()-=n_removed;
	else accepts_since_last_batch() = 0;

	//count with floating-average
	core::Real const inv_A( 1.0/min_structures_for_acceptance_statistics_ );

	while ( n_removed-- > 0 && acceptance_history_.size() ) {
		if ( acceptance_history_.front() ) floating_acceptance_ratio_ -= inv_A;
		acceptance_history_.pop_front();
		acceptance_history_.push_back( false );
	}
	nstruct_ = decoys().size();
}

void ArchiveBase::init_from_decoy_set( core::io::silent::SilentFileData const& sfd ) {
	tr.Trace << "ArchiveBase::init_from_decoy_set" << std::endl;
	//make bogus batch that contains init-file
	Batch init_batch( 0 );
	using namespace core::io::silent;
	using namespace core::chemical;

	//transform to centroid
	SilentFileData centroid_sfd;
	for ( SilentFileData::const_iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
		std::string tag = it->decoy_tag();
		pose::Pose pose;
		(*it)->fill_pose( pose, *ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ) );
		SilentStructOP pss=(*it)->clone();
		pss->fill_struct( pose, tag );
		centroid_sfd.add_structure( pss );
	}
	read_structures( centroid_sfd, init_batch );
}


void ArchiveBase::read_structures( core::io::silent::SilentFileData& sfd, Batch const& batch ) {
	using namespace core;
	using namespace io::silent;

	tr.Debug << "read structures returned for " << batch.batch() << std::endl;
	Size ct( batch.decoys_returned() );
	Size accepted_ct( 0 );
	for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
		std::string tag = it->decoy_tag();
		it->set_decoy_tag( batch.batch()+"_"+ObjexxFCL::lead_zero_string_of( ++ct, 6 ) );

		(*it)->add_comment( TAG_IN_FILE, tag );
		(*it)->add_comment( SOURCE_FILE, batch.silent_out() );

		bool accept = add_structure( *it, batch );
		count_structure( batch, accept );
		accepted_ct += accept ? 1 : 0;
	}

	mem_tr << " ArchiveBase::read_structures now: " << decoys().size() << " decoys "<< std::endl;

	tr.Debug << "...done reading --- " << accepted_ct << " structures were accepted into the archive" << std::endl;
	tr.Trace << "now there are " << decoys().size() << " decoys in archive " << std::endl;
}

bool ArchiveBase::add_structure( core::io::silent::SilentStructOP new_decoy, Batch const& ) {
	add_structure_at_position( decoys_.begin(), new_decoy );
	return true;
}

void ArchiveBase::add_structure_at_position( SilentStructs::iterator iss, core::io::silent::SilentStructOP new_decoy ) {
	decoys_.insert( iss, new_decoy );
	if ( decoys_.size() > max_nstruct_ ) { //take all decoys until full
		decoys_.pop_back();
	}
}

void ArchiveBase::save_to_file( std::string suffix ) {
	std::string const dirname( name() + suffix );
	std::string const filename ( dirname + "/decoys.out" );
	std::string const backup_filename ( dirname + "/decoys.out.backup" );
	std::string const tmp_filename ( dirname + "/tmp_decoys.out" );
	using namespace core::io::silent;

	//don't write empty file
	if ( !decoys().size() ) return;

	utility::file::create_directory( dirname );

	SilentFileData sfd;

	//handle output myself... so it keeps the order of decoys.
	utility::io::ozstream output( tmp_filename );
	if ( decoys_.begin() != decoys_.end() ) (*decoys_.begin())->print_header( output );

	for ( SilentStructs::const_iterator it = decoys_.begin(); it != decoys_.end(); ++it ) {
		//		sfd.add_structure( *it ); //only add OP to sfd
		//	sfd.write_silent_struct( *it, tmp_filename );
		sfd.write_silent_struct( **it, output );
	}
	//	sfd.write_all( tmp_filename );

	//rename to final
	//delete old file
	if ( utility::file::file_exists( filename ) ) {
		rename( filename.c_str(), backup_filename.c_str() );
	}
	rename( tmp_filename.c_str(), filename.c_str() );

	utility::io::ozstream status( dirname+"/STATUS" );
	basic::show_time( tr,  "save "+name()+" status");
	save_status( status );
}

bool ArchiveBase::restore_from_file() {
	std::string const& dirname( name() );
	std::string const filename ( dirname + "/decoys.out" );
	bool b_have_restored( false );
	using namespace core::io::silent;
	mem_tr << "ArchiveBase::restore_from_file..." << std::endl;
	decoys_.clear();
	if ( utility::file::file_exists( filename ) ) {
		SilentFileData sfd;
		if ( !sfd.read_file( filename ) ) throw ( utility::excn::EXCN_BadInput( "problem reading silent file"+filename ) );
		for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
			decoys_.push_back( *it );
		}
	}
	mem_tr << "ArchiveBase::restored " << decoys_.size() << " decoys" << std::endl;
	utility::io::izstream status( dirname+"/STATUS" );
	if ( status.good() ) {
		restore_status( status );
		b_have_restored = true;
	} else {
		tr.Info << name() << ": no archive status found... " << std::endl;
	}

	mem_tr << "ArchiveBase::restored_from_file" << std::endl;
	return b_have_restored;
}

void ArchiveBase::restore_status( std::istream& is ) {
	std::string line;
	getline( is, line );
	is >> total_accepts_ >> accepts_since_last_batch_ >> total_proposed_ >> proposed_since_last_batch_;
	std::string tag;
	is >> tag;
	runtime_assert( tag == "acceptance_history:" );
	is >> floating_acceptance_ratio_;
	is >> tag;
	while ( tag != "END_AH" ) {
		runtime_assert( tag == "AH:" );
		getline( is, line ); //read rest of line
		std::istringstream line_stream( line );
		typedef std::istream_iterator< AcceptHistoryQueue::value_type > istream_iterator;
		std::copy(  istream_iterator( line_stream ), istream_iterator(), std::back_inserter( acceptance_history_ ) );
		is >> tag;
	}
	mem_tr << "ArchiveBase::restored_status" << std::endl;
}

void ArchiveBase::save_status( std::ostream& os ) const {
	using namespace ObjexxFCL::format;
	os << "total_accepts accepts_during_stage total_proposed proposed_during_stage\n"
		 << RJ( 14, total_accepts_+accepts_since_last_batch_ ) << RJ( 25, accepts_since_last_batch_ )
		 << RJ( 15, total_proposed_+proposed_since_last_batch_ ) << RJ( 30, proposed_since_last_batch_) << std::endl;
	os << "acceptance_history: " << floating_acceptance_ratio_ << "\nAH: ";
	Size const cols( 50 );
	Size ct( cols );
	for ( AcceptHistoryQueue::const_iterator it = acceptance_history_.begin(); it != acceptance_history_.end(); ++it, --ct ) {
		os << *it << " ";
		if ( ct <= 1 ) { os << "\nAH: "; ct = cols; }
	}
	os << "\nEND_AH" << std::endl;
}

void DebugArchive::generate_batch() {
	using namespace core::io::silent;
	++ct_batches_;
	tr.Debug << "generate batch number " << ct_batches_ << std::endl;
	if ( ct_batches_ <= 3 || decoys().size()==0 ) {

		//always call start_new_batch to generate a new batch
		Batch& batch( manager().start_new_batch() );

		utility::io::ozstream flags( batch.flag_file() );
		if ( ct_batches_ == 1 && make_mistake_ )	flags << "-abinitio::sskip_stages 1 2" << std::endl;
		flags << "-abinitio::skip_convergence_check" << std::endl;
		flags.close();

		//always finish your batch-generation with "finalize_batch"
		try {
			manager().finalize_batch( batch );
		} catch ( EXCN_Archive& excn ) {
			--ct_batches_;
			make_mistake_ = false; //don't make this mistake again
		}

	} else {
		SilentStructOPs start_decoys;
		std::copy( decoys().begin(), decoys().end(), std::back_inserter( start_decoys ) );

		Batch& batch( manager().start_new_batch( start_decoys ) );

		utility::io::ozstream broker( batch.broker_file() );
		broker << "USE_INPUT_POSE\n\nCLAIMER StartStructClaimer\nEND_CLAIMER\n" << std::endl;
		broker.close();

		manager().finalize_batch( batch );
	}
}

void DebugArchive::score( core::pose::Pose& pose ) const {
	runtime_assert( !pose.is_fullatom() );
	runtime_assert( cen_score_ );
	(*cen_score_)( pose );
}

DebugArchive::DebugArchive( ArchiveManagerAP ptr ) :
	ArchiveBase( ptr ),
	ct_batches_( 0 ),
	cen_score_( NULL ),
	make_mistake_( true )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	if ( !cen_score_ ) cen_score_ = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage4_patch ]() );
}

void DebugArchive::save_status( std::ostream& out ) const {
	ArchiveBase::save_status( out );
	out << "nr_batches_generated: "<< ct_batches_ << std::endl;
}

void DebugArchive::restore_status( std::istream& in ) {
	ArchiveBase::restore_status( in );
	if ( in.good() ) {
		std::string tag;
		in >> tag >> ct_batches_;
	}
}

bool DebugArchive::add_structure( core::io::silent::SilentStructOP decoy, Batch const& ) {
	if ( decoys().size() < 100 ) {
		decoys().push_back( decoy ); //of course this can't remain as simple as this.
	} else {
		decoys().front() = decoy;
	}
	return true;
}

}//archive
}//jd2
}//protocols
