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
// AUTO-REMOVED #include <iostream>
#include <fstream>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/Adduct.fwd.hh>
// AUTO-REMOVED #include <core/chemical/Adduct.hh>
// AUTO-REMOVED #include <core/chemical/AtomICoor.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomICoor.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ElementSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueConnection.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueConnection.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/ICoorOrbitalData.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/OrbitalType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/sdf/MolData.fwd.hh>
// AUTO-REMOVED #include <core/chemical/sdf/MolData.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/conformation/RotamerSetBase.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/XYZEvent.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedStubID.fwd.hh>
// AUTO-REMOVED #include <core/id/TorsionID.fwd.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/PoseInputStream.hh>
// AUTO-REMOVED #include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructCreator.fwd.hh>
#include <core/io/silent/silent.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/DomainMap.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Jump.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MinimizerMapBase.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Stub.fwd.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.fwd.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pose/datacache/ObserverCache.fwd.hh>
// AUTO-REMOVED #include <core/pose/metrics/PoseMetricContainer.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/ConformationEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/DestructionEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/EnergyEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/GeneralEvent.fwd.hh>
// AUTO-REMOVED #include <core/scoring/Energies.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/LREnergyContainer.fwd.hh>
// AUTO-REMOVED #include <core/scoring/MinimizationGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionInfo.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethod.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
// AUTO-REMOVED #include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverList.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/wum/SilentStructStore.fwd.hh>
#include <protocols/wum/WorkUnitBase.fwd.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
// AUTO-REMOVED #include <utility/factory/WidgetRegistrator.hh>
// AUTO-REMOVED #include <utility/file/FileName.fwd.hh>
// AUTO-REMOVED #include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/PathName.fwd.hh>
// AUTO-REMOVED #include <utility/file/PathName.hh>
// AUTO-REMOVED #include <utility/keys/AutoKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/AutoKey.hh>
// AUTO-REMOVED #include <utility/keys/Key.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key.hh>
// AUTO-REMOVED #include <utility/keys/Key2Tuple.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key2Tuple.hh>
// AUTO-REMOVED #include <utility/keys/Key3Tuple.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key3Tuple.hh>
// AUTO-REMOVED #include <utility/keys/Key4Tuple.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key4Tuple.hh>
// AUTO-REMOVED #include <utility/keys/KeyLess.fwd.hh>
// AUTO-REMOVED #include <utility/keys/KeyLookup.fwd.hh>
// AUTO-REMOVED #include <utility/keys/KeyLookup.hh>
// AUTO-REMOVED #include <utility/keys/NoClient.fwd.hh>
// AUTO-REMOVED #include <utility/keys/NoClient.hh>
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.fwd.hh>
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.hh>
// AUTO-REMOVED #include <utility/keys/UserKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/VariantKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/VariantKey.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.hh>
// AUTO-REMOVED #include <utility/options/FileOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileOption.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.hh>
// AUTO-REMOVED #include <utility/options/IntegerVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerVectorOption.hh>
// AUTO-REMOVED #include <utility/options/Option.fwd.hh>
// AUTO-REMOVED #include <utility/options/Option.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.fwd.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.hh>
// AUTO-REMOVED #include <utility/options/PathOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathOption.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.hh>
// AUTO-REMOVED #include <utility/options/RealOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealOption.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.hh>
// AUTO-REMOVED #include <utility/options/StringOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringOption.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.hh>
// AUTO-REMOVED #include <utility/options/mpi_stderr.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/VectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/VectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/signals/BufferedSignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/BufferedSignalHub.hh>
// AUTO-REMOVED #include <utility/signals/Link.fwd.hh>
// AUTO-REMOVED #include <utility/signals/Link.hh>
// AUTO-REMOVED #include <utility/signals/LinkUnit.fwd.hh>
// AUTO-REMOVED #include <utility/signals/LinkUnit.hh>
// AUTO-REMOVED #include <utility/signals/SignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <numeric/numeric.functions.hh>
// AUTO-REMOVED #include <numeric/sphericalVector.fwd.hh>
// AUTO-REMOVED #include <numeric/trig.functions.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.fwd.hh>
// AUTO-REMOVED #include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
// AUTO-REMOVED #include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/TypeTraits.hh>
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
// AUTO-REMOVED #include <cmath>
#include <cstddef>
// AUTO-REMOVED #include <cstdlib>
// AUTO-REMOVED #include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <iostream>
#include <ostream>
// AUTO-REMOVED #include <set>
#include <sstream>
#include <string>
// AUTO-REMOVED #include <utility>
#include <vector>
// AUTO-REMOVED #include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.fwd.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <boost/bind.hpp>
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

