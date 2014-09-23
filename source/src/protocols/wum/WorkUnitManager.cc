// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopHashMap.cc
/// @brief
/// @author Mike Tyka

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/WorkUnitBase.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/util.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <utility/assert.hh>
// AUTO-REMOVED #include <ios>
#include <iostream>
#include <fstream>

//Auto Headers
#include <utility/vector1.hh>
#include <boost/function.hpp>


namespace protocols {
namespace wum {

static thread_local basic::Tracer TR( "WorkUnitManager" );

// 32 bit recognition integer to make sure we're infact about to read/write a WU to disk etc..
const unsigned int WUB_magic_header_integer = 0xAF34B14C;





WorkUnitBaseOP &
WorkUnitQueue::next()
{
	return *(wus_.begin());
}

WorkUnitBaseOP
WorkUnitQueue::pop_next()
{
	runtime_assert( size() != 0 );
	WorkUnitBaseOP tmp = next();
	wus_.pop_front();
	return tmp;
}

WorkUnitQueue::iterator WorkUnitQueue::erase( iterator i ) {
	return wus_.erase( i );
}


void WorkUnitQueue_Swapped::add( WorkUnitBaseOP new_wu )
{
	runtime_assert( n_swap_total_ >= n_swap_dead_ );

	// if either we have to many structures in memory OR we have a running file swap already,
	// add to swap pile, not to in-memory pile
	if( ( wus_.size() > memory_limit_) ||
      ( (n_swap_total_ - n_swap_dead_ ) > 0 ) ){
		add_to_swap( new_wu );
	} else {
		// or call normal parent version
		WorkUnitQueue::add( new_wu );
	}
}



void WorkUnitQueue_Swapped::add_to_swap( WorkUnitBaseOP new_wu ){
	swap_buffer_.add( new_wu );

	if( swap_buffer_.size() > max_swap_buffer_size_ ){
		// drain buffer to disk swap
		std::ofstream ofout( swap_file_.c_str() , std::ios::app | std::ios::binary );
		wum_->write_queue( swap_buffer_, ofout );
		ofout.close();

		// now empty the buffer, ready to take the next bunch of structures
		swap_buffer_.clear();
	}

}




void  WorkUnitManager::register_work_units( const protocols::wum::WorkUnitList &work_unit_list ){
	work_unit_list_.merge( work_unit_list );
}


void WorkUnitManager::write_queues_to_file( const std::string& prefix ) const {
	std::ofstream ofout( std::string(prefix + ".outbound.queue").c_str() , std::ios::binary );
	write_queue( outbound() , ofout );
	ofout.close();
	std::ofstream ofin( std::string(prefix + ".inbound.queue").c_str() , std::ios::binary );
	write_queue( inbound() , ofin );
	ofin.close();
}

void WorkUnitManager::read_queues_from_file( const std::string& prefix )  {
	std::ifstream ifout( std::string(prefix + ".outbound.queue").c_str() , std::ios::binary );
	TR << "Reading outbound queue..." << std::string(prefix + ".outbound.queue") << std::endl;
	read_queue( outbound() , ifout );
	ifout.close();
	std::ifstream ifin( std::string(prefix + ".inbound.queue").c_str() , std::ios::binary );
	TR << "Reading inbound queue..." << std::string(prefix + ".inbound.queue") << std::endl;
	read_queue( inbound() , ifin );
	ifin.close();
}


void WorkUnitManager::read_queue( WorkUnitQueue &the_queue, std::istream &fin ){
	core::Size count=0;
	while( !fin.eof() ){
		TR.Debug << "Read: " << count << std::endl;
		WorkUnitBaseOP new_wu;
		if( !read_work_unit( new_wu, fin ) ) break;
		the_queue.push_back( new_wu );
		count++;
	}
}


void WorkUnitManager::write_queue( const WorkUnitQueue &the_queue, std::ostream &out ) const {
	for( WorkUnitQueue::const_iterator it = the_queue.begin();
				it != the_queue.end(); ++it )
	{
		write_work_unit( *it, out );
	}
}



void WorkUnitManager::write_work_unit( const WorkUnitBaseOP& MPI_ONLY(wu), std::ostream& MPI_ONLY( out ) ) const {
	#ifdef USEMPI
	// serialize data
	double time1=MPI_Wtime();
	wu->serialize();
	double time2=MPI_Wtime();
	// now send data
	int size_of_raw_data;
	unsigned char * raw_data_ptr=NULL;
	size_of_raw_data = wu->raw_data_dump( &raw_data_ptr );
	TR.Debug << "Writing workunit .. " << std::endl;
	out.write( (char*) &WUB_magic_header_integer, 4 );
	out.write( (char*) &size_of_raw_data, 4 );
	out.write( (char*)raw_data_ptr, size_of_raw_data );
	TR.Debug << "  Wrote data.. " << std::endl;
	delete [] raw_data_ptr;
	TR.Debug << "  Deleted temp data.. " << std::endl;
	wu->clear_serial_data();
	double time3=MPI_Wtime();
	TR.Debug << "S: " << time3-time2 << "  " << time2-time1 << "  " << std::endl;
	#endif
}


bool WorkUnitManager::read_work_unit( WorkUnitBaseOP &qualified_wu,  std::istream &in ){
	unsigned int size_of_raw_data=0;
	unsigned char * raw_data_ptr=NULL;
	TR.Debug << "Reading a workunit..."  << std::endl;

	unsigned int my_WUB_magic_header_integer=0;
	// Read magic 32 bit int
	in.read( (char*) &my_WUB_magic_header_integer, 4 );
	if( in.eof() ){
		TR.Debug << "EOF" << std::endl;
		return false;
	}
	if(my_WUB_magic_header_integer != WUB_magic_header_integer){
		TR.Error << "Magic Integer in file: " << my_WUB_magic_header_integer << " != " << WUB_magic_header_integer << std::endl;
		TR.Error << "ERROR Reading in WorkUnit from stream - Magic integer does not match. " << std::endl;
		std::cerr << "Magic Integer in file: " << my_WUB_magic_header_integer << " != " << WUB_magic_header_integer << std::endl;
		utility_exit_with_message( "ERROR Reading in WorkUnit from stream - Magic integer does not match. " );
	}

	in.read( (char*)&size_of_raw_data, 4 );
	if( size_of_raw_data > (1024*1024*1024) ){
		TR.Error << "  Data corruption ? WorkUnitManager::read_work_unit found workunit with memory requirement > 1GB " << std::endl;
	}

	TR.Debug << "  READ WU: BLOCKSIZE: " << size_of_raw_data << std::endl;
	raw_data_ptr = new unsigned char [size_of_raw_data];

	in.read( (char*)raw_data_ptr, (std::streamsize) size_of_raw_data );

	if( raw_data_ptr[size_of_raw_data-1] != 0){
		utility_exit_with_message( "  ERROR: cannot load data - terminal zero not found!" );
		return false;
	}
	raw_data_ptr[size_of_raw_data-1] = 0;
	TR.Debug << "  READ WU: Data: " << std::endl;

	WorkUnitBaseOP wu = new WorkUnitBase;
  runtime_assert( wu != 0 );
	wu->raw_data_load( raw_data_ptr, size_of_raw_data );
	delete [] raw_data_ptr;

  // Here at this point we have a WorkUnitBaseOP to a workUnitBase.
  // Now we need to interpret the id field and upcast or somehow otherwise
  // create the right type of work unit such that the polymorphic code
  // for the interpretation of the serial data can take place.

	qualified_wu = work_unit_list().get_work_unit( *wu )->clone();
  runtime_assert( qualified_wu != 0 );
	// cope over data (the header and the serial data)
	(*qualified_wu) = (*wu);

	TR.Debug << "  Received: " << std::endl;
	if( TR.Debug.visible() ) qualified_wu->print( TR );

	qualified_wu->deserialize( );
	qualified_wu->clear_serial_data();

	TR.Debug << "DONE Receiving" << std::endl;
	return true;
}


core::Size
WorkUnitQueue::mem_foot_print() const {
	core::Size n_structs;
	core::Size structs_memory;
	core::Size WU_memory;
	mem_stats( n_structs, structs_memory, WU_memory);
	return WU_memory + structs_memory;
}

void
WorkUnitQueue::mem_stats(
	core::Size &n_structs,
	core::Size &structs_memory,
	core::Size &WU_memory
) const {
	n_structs=0;
	structs_memory=0;
	WU_memory=0;

	for( const_iterator it = begin(); it != end(); it++ ){
		WU_memory += (*it)->mem_footprint();
		WorkUnitBaseOP wu_op = *it;
		WorkUnit_SilentStructStoreOP structure_wu = dynamic_cast<  WorkUnit_SilentStructStore * > ( wu_op() );
		if ( structure_wu.get() == NULL ) continue;
		SilentStructStore &decoys = structure_wu->decoys();
		n_structs += structure_wu->decoys().size();
		for( SilentStructStore::iterator jt =  decoys.begin();
				jt != decoys.end(); jt ++ )
		{
			structs_memory += (*jt)->mem_footprint();
		}
	}


}



}
}

