// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#if (defined _WIN32) && (!defined WIN_PYROSETTA)
#define ZLIB_WINAPI  // REQUIRED FOR WINDOWS
#endif

#include <protocols/frags/TorsionFragment.hh>

// Rosetta Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>

// Utility Headers
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.frags.TorsionFragment");

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <map>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <protocols/frags/TorsionFragment.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
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
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>


// C++ Headers
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// #include <fstream>
// #include <sstream>

namespace protocols {
namespace frags {

using namespace core;

///\brief insert this piece of fragment to a pose at position "begin"
///
///call pose.set_torsion which maps TorsionID to DOF_ID, so it is safe
///to use even if the folding direction is not standard as N to C.
void
TorsionFragment::insert( pose::Pose & pose, Size const begin ) const
{
	std::cerr << "ERROR! USING **DEPRECATED** FRAGMENTS. PLEASE SWITCH TO src/core/fragments/* FRAGMENTS!" << std::endl;
	for ( Size i=1; i<= size(); ++i ) {
		utility::vector1< Real > const & bb( torsions_[i] );
		Size const seqpos( begin + i - 1 );
		for ( Size j=1; j<= bb.size(); ++j ) {
			pose.set_torsion( id::TorsionID( seqpos, id::BB, j ), bb[j] );
		}
		pose.set_secstruct( seqpos, secstruct_[i] );
	}
}
//////////////////////////////////////////////////////////////////////////////////
///\brief initialize fragment data from a classic Rosetta fragment library
///
///frag_size is the size of fragment in the library and nbb is the number of backbone
///torision angles stored in each fragment, for protein, this will be 3. Return false if
///there is a reading error and all the data which have been stored will be erased.
///
bool
TorsionFragmentLibrary::read_file(
	std::string const filename,
	Size const frag_size,
	Size const nbb
)
{
	utility::io::izstream data ( filename );
	std::cerr << "ERROR! USING **DEPRECATED** FRAGMENTS. PLEASE SWITCH TO src/core/fragments/* FRAGMENTS!" << std::endl;
	if ( !data ) {
		std::cerr << "Cannot open " + ObjexxFCL::string_of(frag_size) + "-mer fragment library file: " + filename + "\n" ;
		return false;
	}
	std::string line;
	while( getline( data, line ) ) {
		std::istringstream line_stream( line );
		std::string tag1, tag2;
		Size position, neighbors;
		line_stream >> tag1 >> position >> tag2 >> neighbors;
		if (line_stream.fail() || tag1 != "position:" || tag2 != "neighbors:" ) {
			std::cerr << " format errors in fragment library file: " << line << std::endl;
			resize(0);
			return false;
		}
		resize( position );
		getline(data, line); // skip blank line
		SingleResidueTorsionFragmentLibrary & current_position_library( (*this)[position] );
		for ( Size i=1; i<=neighbors; ++i ) {
			TorsionFragmentOP fragment = new TorsionFragment( frag_size, nbb );
			// read lines within each fragment
			std::string last_pdb; char last_chain('0'); Size last_seqpos(0);
			for ( Size j=1; j<=frag_size; ++j ) {
				getline(data, line);
				std::istringstream line_stream( line );
				std::string pdb; char chain, aa, secstruct; Size seqpos;
				line_stream >> pdb >> chain >> seqpos >> aa >> secstruct;
				if ( j == 1 ) {
					last_pdb = pdb;
					last_chain = chain;
					last_seqpos = seqpos;
				} else {
					if (last_pdb != pdb  || last_chain != chain || last_seqpos != (seqpos-j+1) ) {
						std::cerr << "fragment reading error -- pdb, chain and seqpos mismatch\n"
										 << "position: " << position  << "; neighbor: " << i  << "; line: " << j << std::endl;
						resize(0);
						return false;
					}
				}
				fragment->set_secstruct(j,secstruct);
				for ( Size k=1; k<=nbb; ++k ) {
					Real torsion;
					line_stream >> torsion;
					fragment->set_torsion(j,k,torsion);
				} // bb torsion
			} // fragment
			// append this fragment to SingleResidueFragmentLibrary
			current_position_library.append_fragment( fragment );
			getline(data,line); // skip blank line
		} // SingleResidueFragmentLibrary
	} // TorsionFragmentLibrary
	TR.Info << "read succesfully " << frag_size << "-mer fragment library file: " << filename << std::endl;
	return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////
///\brief extract a fragment library with smaller fragment size from the one with larger lize
///
/// for example, 1-mer library can ben extracted from 3-mer library. This function is set up
/// in a general way that both sizes of smaller fragment and larger fragment can be flexibly
/// specified by my_size and src_size. Also, the source library does not have to contain
/// fragment data for every residue position. For example, for loop modeling, we only need
/// fragment for loop segment and after this function is called, smaller library is created
/// only for regions in which there is data in the larger library. Lastly, smaller fragment
/// is aligned in the center of larger fragment to have data extracted.
bool
TorsionFragmentLibrary::derive_from_src_lib(
	Size my_size,
	Size src_size,
	TorsionFragmentLibraryCOP src_lib_op
)
{
	// can only go from larger size to smaller size
	runtime_assert( my_size < src_size );
	// get continuous segment in src_lib_op
	std::map< Size, Size > seg_map;
	Size prev_nbrs(0);
	Size seg_start(0);
	for ( Size i = 1; i <= src_lib_op->size(); ++i ) {
		Size current_nbrs( (*src_lib_op)[i].size() );
		if ( !prev_nbrs && current_nbrs ) {
			runtime_assert( !seg_start );
			seg_start = i;
		} else if ( prev_nbrs && !current_nbrs ) {
			runtime_assert( seg_start && seg_start <i );
			seg_map.insert( std::make_pair( seg_start, i-1 ) );
			seg_start = 0;
		} else if ( i == src_lib_op->size() && prev_nbrs && current_nbrs ) {
			runtime_assert( seg_start && seg_start < i );
			seg_map.insert( std::make_pair( seg_start, i ) );
			seg_start = 0;
		}
		prev_nbrs = current_nbrs;
	}
	runtime_assert( !seg_start );

	if ( !seg_map.size() ) {
		std::cerr << "Warning: source fragment library does not have any data!" << std::endl;
		return false;
	}

	// resize myself lib to correct size
	Size const size_diff(src_size - my_size);
	resize( src_lib_op->size() + size_diff );
	// find offset by which shorter frag is aligned to longer frag
	Size const offset( size_diff%2 == 0 ? size_diff/2 : (size_diff-1)/2 );

	// loop through each segment and derive fragments
	for ( std::map<Size, Size>::const_iterator it = seg_map.begin(),
					it_end = seg_map.end(); it != it_end; it++ ) {
		Size const seg_begin(it->first);
		Size const seg_end(it->second);
		Size lib_index, copy_start;
		for ( Size my_index = seg_begin, my_end = seg_end+size_diff; my_index <= my_end; my_index++ ) {
			// get hold of current residue frag lib
			SingleResidueTorsionFragmentLibraryOP my_residue_lib( fragments_[my_index] );
			// figure out which subset of fragment data we should copy
			if ( my_index < seg_begin + offset ) {
			 	lib_index = seg_begin;
				copy_start = 1 + my_index - seg_begin;
			} else if ( my_index > seg_end + offset ) {
				lib_index = seg_end;
				copy_start = 1 + my_index - seg_end;
			}	else {
				lib_index = my_index - offset;
				copy_start = 1 + offset;
			}
			// do actual copying from src residue_fragment_library
			SingleResidueTorsionFragmentLibrary const & src_residue_lib( (*src_lib_op)[lib_index] );
			for ( Size i = 1; i <= src_residue_lib.size(); ++i ) {
				TorsionFragment const & src_fragment( src_residue_lib[i] );
				Size const nbb( src_fragment.nbb() );
				TorsionFragmentOP my_fragment( new TorsionFragment( my_size, nbb ) );
				for ( Size j = 1, jj = copy_start; j <= my_size; ++j, ++jj ) {
					for ( Size k = 1; k <= nbb; ++k ) {
						my_fragment->set_torsion( j, k, src_fragment.get_torsion(jj, k) );
					} // each torsion
					my_fragment->set_secstruct( j, src_fragment.get_secstruct(jj) );
				} // each frag_pos
				my_residue_lib->append_fragment( my_fragment );
			} // each neighbor
		} // each residue position
	} //seg_map

	// successful
	return true;
}

/// @details  Show size info about library
void
TorsionFragmentLibrary::print( std::ostream & os ) const
{
	os << "TorsionFragmentLibrary:: size= " << size() << std::endl;
	for ( Size i=1; i<= size(); ++i ) {
		os << "window: " << i << " depth: " << fragments_[i]->size() << std::endl;
	}
}


} // ns frags
} // ns protocols

