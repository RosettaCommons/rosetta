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

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <utility/assert.hh>
#include <ios>
#include <iostream>
#include <fstream>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStructCreator.fwd.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverList.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/wum/SilentStructStore.fwd.hh>
#include <protocols/wum/WorkUnitBase.fwd.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace protocols {
namespace wum {

static basic::Tracer TR("WorkUnitManager");

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
  runtime_assert( wu );
	wu->raw_data_load( raw_data_ptr, size_of_raw_data );
	delete [] raw_data_ptr;

  // Here at this point we have a WorkUnitBaseOP to a workUnitBase.
  // Now we need to interpret the id field and upcast or somehow otherwise
  // create the right type of work unit such that the polymorphic code
  // for the interpretation of the serial data can take place.

	qualified_wu = work_unit_list().get_work_unit( *wu )->clone();
  runtime_assert( qualified_wu );
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
		WorkUnitBase*  wu_p = &(*wu_op);
		WorkUnit_SilentStructStore* structure_wu = dynamic_cast<  WorkUnit_SilentStructStore * > ( wu_p );
		if ( structure_wu == NULL ) continue;
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

