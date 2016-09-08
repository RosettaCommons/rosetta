// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief general functions for generating typical kinds of Loops sets.
/// @author ashworth

#include <protocols/loops/make_loops.hh>
#include <protocols/loops/Loops.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>

#include <algorithm>

#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>

//Auto Headers


namespace protocols {
namespace loops {

using namespace core;
using namespace chemical;
using utility::vector1;

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static THREAD_LOCAL basic::Tracer TR( "protocols.loops.make_loops", t_info );

void add_loop(
	Size seg_begin,
	Size seg_end,
	pose::Pose const & pose,
	loops::Loops & loops
)
{
	if ( seg_begin >= seg_end || pose.chain( seg_begin ) != pose.chain( seg_end ) ) {
		TR(t_info) << "WARNING: skipping segment with illegal beginning/ending: ";
		if ( pose.pdb_info() ) {
			pose::PDBInfo const & pdb( *pose.pdb_info() );
			TR(t_info) << pdb.chain(seg_begin) << "." << pdb.number(seg_begin) << "."
				<< pose.residue_type(seg_begin).name3() << "/"
				<< pdb.chain(seg_end) << "." << pdb.number(seg_end) << "."
				<< pose.residue_type(seg_end).name3() << std::endl;
		} else {
			TR(t_info) << pose.chain(seg_begin) << "." << seg_begin << "."
				<< pose.residue_type(seg_begin).name3() << "/"
				<< pose.chain(seg_end) << "." << seg_end << "."
				<< pose.residue_type(seg_end).name3() << std::endl;
		}
		return;
	}
	Size const length( seg_end - seg_begin );
	// random cutpoint - needs to be "seg_start < cut < seg_end"
	Size cut(0), safety(0);
	// ensure (cutpoint+1) is not a proline, as this will cause problems (?)
	while ( safety < 100 && ( cut == 0 || pose.residue_type( cut+1 ).aa() == aa_pro ) ) {
		cut = seg_begin + static_cast< Size >( numeric::random::uniform() * ( length ) );
		++safety;
	}
	runtime_assert( pose.residue_type( cut+1 ).aa() != aa_pro );

	if ( pose.pdb_info() ) {
		pose::PDBInfo const & pdb( *pose.pdb_info() );
		TR(t_info) << "adding segment: " << pdb.number( seg_begin ) << "-(" << pdb.number( cut ) << ")-"
			<< pdb.number( seg_end ) << " Chain " << pdb.chain( seg_end ) << std::endl;
	} else {
		TR(t_info) << "adding segment: " << seg_begin << "-(" << cut << ")-" << seg_end <<
			" Chain " << pose.chain( seg_end ) << std::endl;
	}

	loops.add_loop( seg_begin, seg_end, cut, 0, false );
}

/// @brief add a set of loops 'built' around the provided residue indices
/// @details A maximum of gapsize residues are allowed between specified residues for any given loop, loop ends are extended by extend residues, and chain discontinuity starts a new loop
/// @author ashworth
void loops_around_residues(
	loops::Loops & loops,
	pose::Pose const & pose,
	vector1< Size > const & residue_indices,
	Size gapsize /* = 6 */,
	Size extend /* = 2 */
)
{
	if ( residue_indices.empty() ) {
		TR(t_info) << "WARNING: no residues provided--can not define any loops" << std::endl;
		return;
	}
	// loops should be empty, in order to guarantee no overlap
	runtime_assert( loops.num_loop() == 0 );
	// loop ends should be extended by at least 1 residue on each end
	runtime_assert( extend >= 1 );
	core::conformation::Conformation const & conf( pose.conformation() );
	Size seg_begin( residue_indices.front() ), seg_end( residue_indices.front() );
	int last_chain( pose.chain( residue_indices.front() ) );
	Size const nres( pose.size() );
	for ( auto index( residue_indices.begin() ); index != residue_indices.end(); ++index ) {
		runtime_assert( *index > 0 );
		runtime_assert( *index <= nres );
		int const chain( pose.chain(*index) );
		Size const chain_begin( chain > 0 ? conf.chain_begin( chain ) : 1 ),
			chain_end( chain > 0 ? conf.chain_end( chain ) : nres );
		// subtracting from unsigned ints is dangerous!
		Size begin(*index);
		for ( Size i(1); i <= extend; ++i ) {
			if ( begin == chain_begin ) break;
			begin -= 1;
		}
		Size const end( std::min( chain_end, *index + extend ) );
		TR(t_trace) << *index << " ch " << chain << " chbgn " << chain_begin << " chend " << chain_end << " bgn " << begin << " end " << end << std::endl;
		// first residue always opens a segment
		if ( index == residue_indices.begin() ) {
			seg_begin = begin;
			// residue closes open segment and opens a new one only if gapsize is exceeded, or new chain
		} else if ( begin > seg_end + gapsize || chain != last_chain ) {
			// note: adds loop ending on previous residue, not current one
			add_loop( seg_begin, seg_end, pose, loops );
			seg_begin = begin;
		}
		seg_end = end;
		last_chain = chain;
		// last residue always closes last segment
		if ( index + 1 == residue_indices.end() ) add_loop( seg_begin, seg_end, pose, loops );
	}
}

} //namespace loops
} //namespace protocols

