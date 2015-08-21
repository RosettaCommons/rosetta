// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/PartialThreadingMover.hh
/// @brief
/// @author James Thompson

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/PartialThreadingMover.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>

#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/PDBInfo.hh>

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

PartialThreadingMover::PartialThreadingMover(
	core::sequence::SequenceAlignment const & align,
	core::pose::Pose const & template_pose
) :
	ThreadingMover( align, template_pose )
{}

void PartialThreadingMover::apply(
	core::pose::Pose & query_pose
) {
	using core::Size;
	using basic::Tracer;
	using core::id::SequenceMapping;

	static Tracer tr("protocols.comparative_modeling.partial_threading");
	build_loops(false);
	repack_query(true);
	ThreadingMover::apply(query_pose);

	SequenceMapping query_to_pdbseq = get_qt_mapping(query_pose);

	//fpd update PDBinfo to have correct residue numbering
	utility::vector1< int > pdb_numbering;
	utility::vector1< char > pdb_chains;

	// iterate backwards as we change downstream sequence numbering with deletions
	tr.Debug << "current sequence is " << query_pose.sequence() << std::endl;
	for ( Size resi = query_pose.total_residue(); resi >= 1; --resi ) {
		Size const t_resi = query_to_pdbseq[ resi ];

		if ( t_resi == 0 && query_pose.total_residue() > 1 ) {
			query_pose.conformation().delete_residue_slow(resi);
		} else {
			pdb_numbering.push_back( resi );
			pdb_chains.push_back( 'A' );
		}
	} // for resi
	tr.Debug << "final sequence is " << query_pose.sequence() << std::endl;

	std::reverse(pdb_numbering.begin(), pdb_numbering.end());
	core::pose::PDBInfoOP new_pdb_info( new core::pose::PDBInfo(query_pose,true) );
	new_pdb_info->set_numbering( pdb_numbering );
	new_pdb_info->set_chains( pdb_chains );
	query_pose.pdb_info( new_pdb_info );
	query_pose.pdb_info()->obsolete( false );

	tr.flush_all_channels();
} // apply

void PartialThreadingMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /* data */,
	protocols::filters::Filters_map const & /* filters */,
	protocols::moves::Movers_map const & /* movers */,
	core::pose::Pose const & /* pose */
) {
	// need to provide aln_fn, template Pose somehow
	runtime_assert( tag->hasOption("aln_fn") );
	runtime_assert( tag->hasOption("aln_id") );
	runtime_assert( tag->hasOption("template_pdb_fn") );

	using std::string;
	using utility::vector1;
	using core::pose::Pose;
	using core::import_pose::pose_from_pdb;
	using core::sequence::read_aln;
	using core::sequence::SequenceAlignment;

	string const template_pdb_fn(
		tag->getOption< string >("template_pdb_fn")
	);
	Pose template_pose;
	core::import_pose::pose_from_pdb(template_pose,template_pdb_fn);
	ThreadingMover::template_pose(template_pose);

	string const aln_fn( tag->getOption< string >("aln_fn") );
	string const aln_id( tag->getOption< string >("aln_id") );
	string aln_format("grishin");
	if ( tag->hasOption("aln_format") ) {
		aln_format = tag->getOption< string >("aln_format");
	}

	vector1< SequenceAlignment > alns = read_aln( aln_fn, aln_format );
	typedef vector1< SequenceAlignment >::const_iterator iter;
	bool found_aln(false);
	for ( iter it = alns.begin(), end = alns.end(); it != end; ++it ) {
		if ( it->alignment_id() == aln_id ) {
			found_aln = true;
			ThreadingMover::alignment(*it);
		}
	}

	if ( !found_aln ) {
		string const msg(
			"Error: couldn't find aln with id " + aln_id +
			" in aln_file " + aln_fn + "!"
		);
		utility_exit_with_message("Error!");
	}
}

std::string
PartialThreadingMover::get_name() const {
	return "PartialThreadingMover";
}

} // comparative_modeling
} // protocols
