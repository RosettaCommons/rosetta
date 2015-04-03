// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum/WorkUnitBase.cc
/// @brief
/// @author Mike Tyka


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <string>

#include <utility/vector1.hh>


#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
	#ifndef WIN_PYROSETTA
	    #include <windows.h>
	#endif
#endif

namespace protocols {
namespace wum {


static thread_local basic::Tracer TR( "WorkUnitBase" );

WorkUnitBase::WorkUnitBase ( )
{
	TR.Debug << "Setting WorkUnitBaseType" << std::endl;
	header.id_ = 0;
	header.unixtime_creation_ = time(NULL);
	header.unixtime_start_ = 0;
	header.unixtime_stop_ = 0;
	header.extra_data_1_=0;
	header.extra_data_2_=0;
	header.extra_data_3_=0;
	header.extra_data_4_=0;
	last_received_from_ = 0;
	set_wu_type("uninitialized");
	set_options("uninitialized");
	// make a unique WU identifier.
	create_unique_id();
}

void WorkUnitBase::add_blacklist( int mpi_rank ) {
		blacklist_.push_back( mpi_rank );
}

void WorkUnitBase::clear_blacklist() {
		blacklist_.clear();
}

bool WorkUnitBase::in_blacklist( int mpi_rank ) {
		if( blacklist_.size() == 0 ) return false;
		std::vector< int >::iterator i;
		i = find(blacklist_.begin(), blacklist_.end(), mpi_rank );
		if( i != blacklist_.end() ) return true;
		return false;
}

void WorkUnitBase::raw_data_load( const unsigned char * raw_data_ptr, unsigned int size )
{
	TR.Debug << "Extracting header information:" << sizeof( WU_Header ) << "  " << size << std::endl;
	WU_Header *tgtheader = (WU_Header*) raw_data_ptr;
	header = *tgtheader;

	TR.Debug << "Extracting data information " << std::endl;
	// make sure the last byte of the data block is a 0 ( it should be a 0 terminated string)
	if( raw_data_ptr[size-1] != 0){
		TR.Error << "ERROR: cannot load data - terminal zero not found!" << std::endl;
		serial_data() = "";
		return;
	}

	std::string strdata( & (((char*)raw_data_ptr) [ sizeof( WorkUnitBase::WU_Header ) ]) );
	serial_data() = strdata;
}


void WorkUnitBase::run(){
	TR.Debug << "WorkUnitBase was called." << std::endl;
}

void
WorkUnitBase::print( std::ostream & out, bool verbose ) const {
	if( verbose ){
		out << "WU.wu_type_:   "  << header.wu_type_  			<< std::endl;
		out << "WU_id:         " << header.id_    << std::endl;
		out << "WU_time_create:" << header.unixtime_creation_<< std::endl;
		out << "WU_time_start: " << header.unixtime_start_<< std::endl;
		out << "WU_time_stop:  " << header.unixtime_stop_<< std::endl;
		out << "WU_int1:       " << header.extra_data_1_<< std::endl;
		out << "WU_int2:       " << header.extra_data_2_<< std::endl;
		out << "WU_int3:       " << header.extra_data_3_<< std::endl;
		out << "WU_int4:       " << header.extra_data_4_<< std::endl;
		out << "WU_options:    " << header.options_ << std::endl;
		out << "WU_serial:     " << serial_data().substr(0,40) << " [...] " << std::endl;
	} else {
		out << "WU: id: "<< header.id_ << " " << header.wu_type_
		#ifndef __CYGWIN__ // Workaround for CygWin and GCC 4.5
			<< " start: " << std::max( -1, (int)header.unixtime_start_ - (int)header.unixtime_creation_ )
			<< " stop: "  << std::max( -1, (int)header.unixtime_stop_ -  (int)header.unixtime_stop_ )
		#else
			<< " start: " << max( -1, (int)header.unixtime_start_ - (int)header.unixtime_creation_ )
			<< " stop: "  << max( -1, (int)header.unixtime_stop_ -  (int)header.unixtime_stop_ )
		#endif
		<< " dat: "
		 << header.extra_data_1_<<" "
		 << header.extra_data_2_<<" "
		 << header.extra_data_3_<<" "
		 << header.extra_data_4_<<" "
		<< " opt: " << header.options_
		<< " srl: <" << serial_data().substr(0,40)
		<< ((serial_data().size() > 40) ? "[...]>" : ">") << std::endl;


	}

}

unsigned int
WorkUnitBase::raw_data_size() const
{
	TR.Trace << "RawData Sizes: " <<  sizeof( WorkUnitBase::WU_Header ) << "  " << serial_data().size() << std::endl;
	return sizeof( WorkUnitBase::WU_Header ) + serial_data().size() + 1;
}

unsigned int
WorkUnitBase::raw_data_dump( unsigned char ** raw_data_ptr ) const
{
	unsigned int datasize = raw_data_size();
	(*raw_data_ptr) = new unsigned char [ datasize ];

	// first write the header:
	WU_Header *tgtheader = (WU_Header*) (*raw_data_ptr);
	*tgtheader = header;

	// now dump the contents of string
	unsigned char *tgtstring = &((*raw_data_ptr)[sizeof( WorkUnitBase::WU_Header ) ]);
	memcpy( (void*) tgtstring, serial_data().c_str(), serial_data().size() + 1 );

	TR.Debug << "Target string formatted: " << serial_data().substr(0, 40 ) << "[...]" << std::endl;

	return datasize;
}


void
WorkUnitBase::set_wu_type( const std::string &text ){
	#ifndef __CYGWIN__ // Workaround for CygWin and GCC 4.5
		unsigned int length = std::min( (int)text.length(), int(128) );
	#else
		unsigned int length = min( (int)text.length(), int(128) );
	#endif
	if( length == 0 ) return;
	strcpy( &header.wu_type_[0], text.c_str() );
}

std::string WorkUnitBase::get_wu_type() const {
	return std::string( &header.wu_type_[0] );
}


void
WorkUnitBase::set_options( const std::string &text ){
	#ifndef __CYGWIN__ // Workaround for CygWin and GCC 4.5
		unsigned int length = std::min( (int)text.length(), (int)128 );
	#else
		unsigned int length = min( (int)text.length(), (int)128 );
	#endif
	if( length == 0 ) return;
	strcpy( &header.options_[0], text.c_str() );
}

std::string WorkUnitBase::get_options() const {
	return std::string( &header.options_[0] );
}


void WorkUnitBase::set_run_start(){
	header.unixtime_start_ = time(NULL);
}

void WorkUnitBase::set_run_stop(){
	header.unixtime_stop_ = time(NULL);
}

core::Size WorkUnitBase::get_run_time(){
	return header.unixtime_stop_ - header.unixtime_start_;
}


void WorkUnit_Wait::run(){
	//TR << "Waiting for " << header.extra_data_1_ << std::endl;
#ifdef _WIN32
	#ifndef WIN_PYROSETTA
		Sleep( header.extra_data_1_ * 1000 );
	#endif
#else
	sleep( header.extra_data_1_ );
#endif
}


// @brief write decoys into serial data store overwritinge whatever was there before. It b// @brief Read decoys from serial data store. Overwrite what's in the SilentStruct store. asically syncs the silent struct store with the derial data
void
WorkUnit_SilentStructStore::serialize()
{
  decoys_.serialize( serial_data() );
}

// @brief Read decoys from serial data store. Overwrite what's in the SilentStruct store.
void
WorkUnit_SilentStructStore::deserialize()
{
  decoys_.clear();
  decoys_.read_from_string( serial_data() );
}


WorkUnit_MoverWrapper::WorkUnit_MoverWrapper( protocols::moves::MoverOP the_mover ):
	WorkUnit_SilentStructStore(),
	the_mover_(the_mover )
{
	// Figure out mover
	set_defaults();
}

void
WorkUnit_MoverWrapper::set_defaults()
{
}

void
WorkUnit_MoverWrapper::run(){
  using namespace core::pose;
	TR.Debug << "Executing mover wrapper..." << std::endl;

	// first figure what mover we should be using
	SilentStructStore result_store;
	if( decoys().size() == 0 ){
		TR.Error << "ERROR: WU did not contain any decoys!" << std::endl;
	}else{
		TR.Debug << "Applying the mover .. " << std::endl;
		for( SilentStructStore::const_iterator it = decoys().begin() ; it != decoys().end(); ++it ){
			Pose pose;
			runtime_assert(*it != 0);
			(*it)->fill_pose( pose );
			the_mover_->apply( pose );
			result_store.add( pose );
		}
	}
  decoys().clear();
  decoys().add( result_store );
}


}
}

