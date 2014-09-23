// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/ExtraThreadingMover.cc
/// @brief method declarations for ExtraThreadingMover.
/// @author James Thompson

// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/ExtraThreadingMover.hh>

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// AUTO-REMOVED #include <core/conformation/Residue.functions.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>

#include <core/sequence/SWAligner.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>
#include <core/scoring/rms_util.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>

// option key includes
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace protocols {
namespace comparative_modeling {

static thread_local basic::Tracer tr( "protocols.comparative_modeling.extra_threading" );

ExtraThreadingMover::ExtraThreadingMover(
	core::sequence::SequenceAlignment const & align,
	core::pose::Pose const & template_pose,
	utility::vector1< core::Size > const & residue_selection
) :
	protocols::moves::Mover( "ExtraThreadingMover" ),
	query_index_( 1 ),
	template_index_( 2 ),
	template_pose_( template_pose ),
	align_( align ),
	residue_selection_(residue_selection)
{}

/// @brief Returns the index of the query sequence in SequenceAlignment
/// object.
core::Size ExtraThreadingMover::query_index() const {
	return query_index_;
}

/// @brief Returns the index of the template sequence in SequenceAlignment
/// object.
core::Size ExtraThreadingMover::template_index() const {
	return template_index_;
}

/// @brief Returns the SequenceAlignment object used in threading.
core::sequence::SequenceAlignment ExtraThreadingMover::alignment() {
	return align_;
}

core::id::SequenceMapping ExtraThreadingMover::get_qt_mapping(
	core::pose::Pose const & query_pose
) const {
	// put checks here to make sure that the template pose and the template
	// pose in the alignment file match up!
	using namespace core::id;
	using namespace core::sequence;

	SequenceOP query_sequence( new Sequence(
			align_.sequence( 1 )->ungapped_sequence(),
			align_.sequence( 1 )->id(),
			align_.sequence( 1 )->start()
		) );

	SequenceOP aligned_template(
		align_.sequence(template_index_)->clone()
	);

	SequenceOP t_align_seq( new Sequence(
			aligned_template->ungapped_sequence(),
			aligned_template->id() + "_align_seq",
			aligned_template->start()
		) );

	SequenceOP t_pdb_seq( new Sequence (
			template_pose_.sequence(),
			aligned_template->id() + "_pdb_seq",
			1
		) );

	// construct an intermediate alignment of the sequence from the alignment
	// to the sequence in the PDB file.
	SWAligner sw_align;
	ScoringSchemeOP ss( new SimpleScoringScheme( 120, 0, -100, 0 ) );

	tr.Debug << "query sequence         : " << query_pose.sequence() << std::endl;
	tr.Debug << "query sequence         : " << (*query_sequence) << std::endl;
	tr.Debug << "aligned_template       : " << (*aligned_template) << std::endl;
	tr.Debug << "template_sequence (aln): " << (*t_align_seq) << std::endl;
	tr.Debug << "template_sequence (pdb): " << (*t_pdb_seq) << std::endl;

	SequenceAlignment intermediate = sw_align.align( t_align_seq, t_pdb_seq, ss );

	if ( intermediate.identities() != intermediate.length() ) {
		tr.Warning << "Error: potential mismatch between sequence from alignment ";
		tr.Warning << " and sequence from PDB!" << std::endl;
		tr.Warning << "alignment: " << std::endl << intermediate
			<< std::endl;
	}

	SequenceMapping query_to_fullseq = align_.sequence_mapping(
		query_index_, template_index_
	);
	tr.Debug << "Query:    " << *align_.sequence( query_index_ ) << std::endl;
	tr.Debug << "Template: " << *align_.sequence( template_index_ ) << std::endl;
	tr.Debug << "Original Mapping:" <<  query_index_ << "-->" << template_index_
		<<  std::endl;
	query_to_fullseq.show( tr.Debug );

	SequenceMapping intermed_map = intermediate.sequence_mapping( 1, 2 );

	// final mapping is the mapping from query to the template PDB sequence,
	// rather then the direct template sequence.
	SequenceMapping query_to_pdbseq = core::sequence::transitive_map(
		query_to_fullseq, intermed_map
	);
	tr.Debug << "Transitive Map" << std::endl;
	query_to_pdbseq.show( tr.Debug );

	return query_to_pdbseq;
}

void ExtraThreadingMover::apply(
	core::pose::Pose & query_pose
) {
	using core::Real;
	using core::Size;
	using std::string;
	using utility::vector1;
	using core::id::AtomID;
	using core::id::AtomID_Map;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::id::SequenceMapping mapping = get_qt_mapping(query_pose);
	tr.Debug  << "mapping from " << query_index_ << " to " << template_index_
		<< std::endl;
	mapping.show( tr.Debug );

	std::string const template_id(utility::file_basename(template_pose_.pdb_info()->name()));
	core::pose::add_score_line_string( query_pose, "extra_template", template_id );

	// superimpose query onto template
	AtomID_Map< AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, query_pose, core::id::BOGUS_ATOM_ID );

	for ( Size ii = 1; ii <= query_pose.total_residue(); ++ii ) {
		Size const templ_ii( mapping[ii] );
		if ( templ_ii == 0 ) continue;
		if ( !query_pose.residue(ii).has("CA") || !template_pose_.residue(templ_ii).has("CA") ) continue;

		AtomID const id1( query_pose.residue(ii).atom_index("CA"), ii );
		AtomID const id2( template_pose_.residue(templ_ii).atom_index("CA"), templ_ii );
		atom_map.set( id1, id2 );
	}

	using core::scoring::superimpose_pose;
	superimpose_pose( query_pose, template_pose_, atom_map );

	// Iterate over residue selection, copy relevant residues.
	// Make sure that fold-tree is coherent by checking for bonded_partners
	using core::conformation::Residue;
 	typedef vector1< Size >::const_iterator iter;
	std::map< Size, Size > template_to_query_res; // keeps track of residues that we've added to query_pose.
	Size const query_jump_anchor( query_pose.total_residue() ); // maybe do this with a VRT residue later.
	//Size current_chain( query_pose.residue(query_jump_anchor).chain() );

	for ( iter it = residue_selection_.begin(), end = residue_selection_.end(); it != end; ++it ) {
		Residue const & template_res(template_pose_.residue(*it));
		Residue new_rsd(template_res);

		bool append_by_jump(true);
		// Check for chemical bond to a residue that we've already added to query_pose.
		// If that bond exists, add this residue by bond.
		typedef std::map< Size, Size >::const_iterator map_iter;
		for ( map_iter m_iter = template_to_query_res.begin(), m_end = template_to_query_res.end(); m_iter != m_end; ++m_iter ) {
			Residue const & partner_res( template_pose_.residue( m_iter->first ) );
			if ( template_res.is_bonded( partner_res ) ) {
				// append by bond to analogous residue in query
				query_pose.conformation().append_residue_by_bond(new_rsd);
				append_by_jump = false;
			}
		}

		if ( append_by_jump ) {
			query_pose.conformation().append_residue_by_jump(new_rsd, query_jump_anchor, "", "", true);
		}

		// keep track of what extra residues have been added to query_pose.
		template_to_query_res[*it] = query_pose.total_residue();
	}

	tr.flush();
} // apply

std::string
ExtraThreadingMover::get_name() const {
	return "ExtraThreadingMover";
}

} // comparative_modeling
} // protocols
