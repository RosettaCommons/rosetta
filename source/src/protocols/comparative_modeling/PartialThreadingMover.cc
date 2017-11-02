// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/PartialThreadingMover.hh
/// @brief
/// @author James Thompson
/// @author fpd some updates

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/PartialThreadingMover.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/util.hh>

#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/pack/optimizeH.hh>
#include <core/pack/pack_missing_sidechains.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// C++ headers
#include <string>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace comparative_modeling {

static THREAD_LOCAL basic::Tracer tr( "protocols.comparative_modeling.threading" );

PartialThreadingMover::PartialThreadingMover(
	core::sequence::SequenceAlignment const & align,
	core::pose::Pose const & template_pose
) : template_pose_( template_pose ), align_( align ) {}

void PartialThreadingMover::apply( core::pose::Pose & query_pose ) {
	using namespace core::scoring;
	using core::Size;
	using basic::Tracer;
	using core::id::SequenceMapping;

	// alignment
	SequenceMapping query_to_pdbseq = get_qt_mapping_general(query_pose, align_, template_pose_, 1, 2 );

	//fpd update PDBinfo to have correct residue numbering
	utility::vector1< int > pdb_numbering;
	utility::vector1< char > pdb_chains;

	//////////
	// 1) delete all unaligned residues
	tr.Debug << "current sequence is " << query_pose.sequence() << std::endl;
	for ( Size resi = query_pose.size(); resi >= 1; --resi ) {
		Size const t_resi = query_to_pdbseq[ resi ];

		if ( t_resi == 0 && query_pose.size() > 1 ) {
			query_pose.conformation().delete_residue_slow(resi);
		} else {
			pdb_numbering.push_back( resi );
			pdb_chains.push_back( 'A' );
		}
	} // for resi

	std::reverse(pdb_numbering.begin(), pdb_numbering.end());
	core::pose::PDBInfoOP new_pdb_info( new core::pose::PDBInfo(query_pose,true) );
	new_pdb_info->set_numbering( pdb_numbering );
	new_pdb_info->set_chains( pdb_chains );
	query_pose.pdb_info( new_pdb_info );
	query_pose.pdb_info()->obsolete( false );

	//////////
	// 2) update coords
	std::string const template_id( utility::file_basename( template_pose_.pdb_info()->name() ) );
	core::pose::add_score_line_string( query_pose, "template", template_id );

	core::id::AtomID_Mask missing( true );
	core::pose::initialize_atomid_map( missing, query_pose ); // used for repacking atoms

	utility::vector1< core::id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;

	for ( Size resi = 1; resi <= query_pose.size(); resi++ ) {
		Size const t_resi = query_to_pdbseq[ pdb_numbering[resi] ];

		// skip this residue if we're not aligned
		if ( t_resi == 0 ) continue;
		if ( t_resi > template_pose_.size() ) continue;

		for ( Size atomj = 1; atomj <= query_pose.residue(resi).natoms(); ++atomj ) {
			std::string atom_name( query_pose.residue(resi).atom_name( atomj ));

			// unless match, don't copy BB
			if ( !query_pose.residue(resi).atom_is_backbone(atomj) &&
					query_pose.residue(resi).aa() != template_pose_.residue(t_resi).aa() ) continue;
			if ( !template_pose_.residue_type(t_resi).has( atom_name ) ) continue;

			// match!
			core::id::AtomID atm_id( atomj, resi );
			missing[ atm_id ] = false;
			atm_ids.push_back( atm_id );
			atm_xyzs.push_back( template_pose_.residue(t_resi).xyz( atom_name ) );
		} // for atom_i
	}

	query_pose.batch_set_xyz( atm_ids, atm_xyzs );
	core::conformation::idealize_position( query_pose.size(), query_pose.conformation() );

	//////////
	// 3) repack all missing atoms & idealize H
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	core::pack::pack_missing_sidechains( query_pose, missing );

	for ( core::Size ii = 1; ii <= query_pose.size(); ++ii ) {
		core::conformation::ResidueOP iires = query_pose.residue( ii ).clone();
		core::conformation::idealize_hydrogens( *iires, query_pose.conformation() );
		query_pose.replace_residue( ii, *iires, false );
	}
	core::pack::optimize_H_and_notify( query_pose, missing );
} // apply


std::string
PartialThreadingMover::get_name() const {
	return "PartialThreadingMover";
}

} // comparative_modeling
} // protocols
