// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/ExtraThreadingMover.cc
/// @brief method declarations for ExtraThreadingMover.
/// @author James Thompson

#include <protocols/comparative_modeling/ExtraThreadingMover.hh>
#include <protocols/comparative_modeling/util.hh>

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/Residue.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>

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

#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace protocols {
namespace comparative_modeling {

static THREAD_LOCAL basic::Tracer tr( "protocols.comparative_modeling.extra_threading" );

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
	return get_qt_mapping_general( query_pose, align_, template_pose_, query_index_, template_index_ );
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

	for ( Size ii = 1; ii <= query_pose.size(); ++ii ) {
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
	Size const query_jump_anchor( query_pose.size() ); // maybe do this with a VRT residue later.
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
		template_to_query_res[*it] = query_pose.size();
	}

	tr.flush();
} // apply

std::string
ExtraThreadingMover::get_name() const {
	return "ExtraThreadingMover";
}

} // comparative_modeling
} // protocols
