// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

static thread_local basic::Tracer tr( "protocols.jd2.Archive" );
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
	SilentFileData alternative_sfd;
	read_structures( centroid_sfd, alternative_sfd, init_batch );
}

void ArchiveBase::erase_decoy( std::string const&  tag ) {
	SilentStructs::iterator it;
	for ( it = decoys_.begin(); it != decoys_.end(); ++it ) {
		if ( (*it)->decoy_tag() == tag ) break;
	}
	decoys_.erase( it );
}

void ArchiveBase::read_structures(
		core::io::silent::SilentFileData& sfd,
		core::io::silent::SilentFileData& alternative_decoys,
		Batch const& batch
) {
	using namespace core;
	using namespace io::silent;

	tr.Debug << "read structures returned for " << batch.batch() << std::endl;
	Size ct( batch.decoys_returned()-sfd.size() );
	tr.Debug << "first count for new decoys : " << ct << " in batch " << batch.id() << std::endl;
	Size accepted_ct( 0 );
	for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
		std::string tag = it->decoy_tag();
		it->set_decoy_tag( batch.batch()+"_"+ObjexxFCL::lead_zero_string_of( ++ct, 6 ) );

		(*it)->add_comment( TAG_IN_FILE, tag );
		(*it)->add_comment( SOURCE_FILE, batch.silent_out() );

		core::io::silent::SilentStructOP alternative_decoy = NULL;
		SilentFileData::iterator alternative_it = alternative_decoys.get_iterator_for_tag( tag );
		if ( alternative_it != alternative_decoys.end() ) {
			alternative_decoy = *alternative_it;
			alternative_decoy->set_decoy_tag( (*it)->decoy_tag() );
		}
		if ( batch.intermediate_structs() ) {
			if ( !alternative_decoy ) {
				std::cerr << "no alternative decoy found for " << tag << " in batch " << batch.id() << std::endl;
				//this might be for diversity pool bail-outs if they happen befor stage2 ...
				alternative_decoy = (*it)->clone();
			}
		}
		bool accept = add_structure( *it, alternative_decoy, batch );
		count_structure( batch, accept );
		accepted_ct += accept ? 1 : 0;
	}

	mem_tr << " ArchiveBase::read_structures now: " << decoys().size() << " decoys "<< std::endl;

	tr.Debug << "...done reading --- " << accepted_ct << " structures were accepted into the archive" << std::endl;
	tr.Trace << "now there are " << decoys().size() << " decoys in archive " << std::endl;
}

bool ArchiveBase::add_structure(
	core::io::silent::SilentStructOP new_decoy,
	core::io::silent::SilentStructOP alternative_decoy,
	Batch const&
) {
	add_structure_at_position( decoys_.begin(), new_decoy, alternative_decoy );
	return true;
}

void ArchiveBase::add_structure_at_position (
  SilentStructs::iterator iss,
	core::io::silent::SilentStructOP new_decoy,
	core::io::silent::SilentStructOP /*alternative_decoy*/
) {
	decoys_.insert( iss, new_decoy );
	if ( decoys_.size() > max_nstruct_ ) { //take all decoys until full
		decoys_.pop_back();
	}
}

void ArchiveBase::save_decoys( std::string const& dirname, std::string const& name, SilentStructs const& decoys ) {
	utility::file::create_directory( dirname );
	std::string const filename ( dirname + "/"+name+".out" );
	std::string const backup_filename ( dirname + "/"+name+".out.backup" );
	std::string const tmp_filename ( dirname + "/tmp_"+name+".out" );
	using namespace core::io::silent;
	SilentFileData sfd;

	//handle output myself... so it keeps the order of decoys.
	utility::io::ozstream output( tmp_filename );
	if ( decoys.begin() != decoys.end() ) (*decoys.begin())->print_header( output );

	for ( SilentStructs::const_iterator it = decoys.begin(); it != decoys.end(); ++it ) {
		sfd.write_silent_struct( **it, output );
	}

	//rename to final
	//delete old file
	if ( utility::file::file_exists( filename ) ) {
		rename( filename.c_str(), backup_filename.c_str() );
	}
	rename( tmp_filename.c_str(), filename.c_str() );
}


void ArchiveBase::save_to_file( std::string suffix ) {
	//don't write empty file
	if ( !decoys().size() ) return;

	std::string const dirname( name() + suffix );
	save_decoys( dirname, "decoys", decoys_ );

	//save status
	utility::io::ozstream status( dirname+"/STATUS" );
	basic::show_time( tr,  "save "+name()+" status");
	save_status( status );
}

bool ArchiveBase::restore_from_file() {
	bool b_have_restored( false );
	std::string const& dirname( name() );
	std::string const filename ( dirname + "/decoys.out" );
	load_decoys( filename, decoys_ );
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

void ArchiveBase::load_decoys( std::string const& filename, SilentStructs& decoys ) {
	using namespace core::io::silent;
	mem_tr << "ArchiveBase::restore_from_file..." << std::endl;
	decoys.clear();
	if ( utility::file::file_exists( filename ) ) {
		SilentFileData sfd;
		if ( !sfd.read_file( filename ) ) throw ( utility::excn::EXCN_BadInput( "problem reading silent file"+filename ) );
		for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
			decoys.push_back( *it );
		}
	}
	mem_tr << "ArchiveBase::restored " << decoys.size() << " decoys" << std::endl;
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


}//archive
}//jd2
}//protocols
