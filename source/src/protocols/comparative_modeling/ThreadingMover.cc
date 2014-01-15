// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/ThreadingMover.cc
/// @brief method declarations for ThreadingMover.
/// @author James Thompson

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/LoopMoverFactory.hh>

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/Residue.functions.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/fragment/FragSet.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/NamedAtomID.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>

#include <core/pack/optimizeH.hh>
#include <core/pack/pack_missing_sidechains.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <core/util/SwitchResidueTypeSet.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace protocols {
namespace comparative_modeling {

static basic::Tracer tr("protocols.comparative_modeling.threading");

// Empty constructor
ThreadingMover::ThreadingMover() : protocols::moves::Mover("ThreadingMover"),
		query_index_(1),
		template_index_(2),
		build_query_loops_(true),
		repack_query_(true),
		randomize_loop_coords_(false),
		min_loop_size_(3)
{}

// Copy constructor
ThreadingMover::ThreadingMover(ThreadingMover const & object_to_copy) : protocols::moves::Mover(object_to_copy),
		query_index_(object_to_copy.query_index_),
		template_index_(object_to_copy.template_index_),
		template_pose_(object_to_copy.template_pose_),
		align_(object_to_copy.align_),
		build_query_loops_(object_to_copy.build_query_loops_),
		repack_query_(object_to_copy.repack_query_),
		randomize_loop_coords_(object_to_copy.randomize_loop_coords_),
		min_loop_size_(object_to_copy.min_loop_size_)
{}

ThreadingMover::ThreadingMover(
	core::sequence::SequenceAlignment const & align,
	core::pose::Pose const & template_pose
) :
	protocols::moves::Mover( "ThreadingMover" ),
	query_index_( 1 ),
	template_index_( 2 ),
	template_pose_( template_pose ),
	align_( align ),
	build_query_loops_( true ),
	repack_query_( true ),
	randomize_loop_coords_( false ),
	min_loop_size_( 3 )
{}

/// @brief Returns the index of the query sequence in SequenceAlignment
/// object.
core::Size ThreadingMover::query_index() const {
	return query_index_;
}

/// @brief Returns the index of the template sequence in SequenceAlignment
/// object.
core::Size ThreadingMover::template_index() const {
	return template_index_;
}

/// @brief Returns the SequenceAlignment object used in threading.
core::sequence::SequenceAlignment ThreadingMover::alignment() {
	return align_;
}

/// @brief Sets the index of the query sequence in SequenceAlignment object.
void ThreadingMover::query_index( core::Size new_index ) {
	query_index_ = new_index;
}

/// @brief Sets the index of the template sequence in SequenceAlignment
/// object.
void ThreadingMover::template_index( core::Size new_index ) {
	template_index_ = new_index;
}

/// @brief Sets the SequenceAlignment associated with this object.
void ThreadingMover::alignment( core::sequence::SequenceAlignment new_align ) {
	align_ = new_align;
}

void ThreadingMover::template_pose( core::pose::Pose template_pose ) {
	template_pose_ = template_pose;
}

//boolean setters
void ThreadingMover::build_loops( bool setting ) {
	build_query_loops_ = setting;
}

void ThreadingMover::randomize_loop_coords( bool setting ) {
	randomize_loop_coords_ = setting;
}

void ThreadingMover::repack_query( bool setting ) {
	repack_query_ = setting;
}

//boolean getters
bool ThreadingMover::build_loops() const {
	return build_query_loops_;
}

bool ThreadingMover::repack_query() const {
	return repack_query_;
}

bool ThreadingMover::randomize_loop_coords() {
	return randomize_loop_coords_;
}

void ThreadingMover::min_loop_size( core::Size const new_size ) {
	min_loop_size_ = new_size;
}

core::Size ThreadingMover::min_loop_size() const {
	return min_loop_size_;
}

utility::vector1< core::fragment::FragSetOP > ThreadingMover::frag_libs() const {
    return frag_libs_;
}

void ThreadingMover::frag_libs(
   utility::vector1< core::fragment::FragSetOP > new_libs
) {
   frag_libs_ = new_libs;
}

core::id::SequenceMapping ThreadingMover::get_qt_mapping(
	core::pose::Pose const & query_pose
) const {
	// put checks here to make sure that the template pose and the template
	// pose in the alignment file match up!
	using namespace core::id;
	using namespace core::sequence;

	SequenceOP query_sequence(
		new Sequence(
			align_.sequence( 1 )->ungapped_sequence(),
			align_.sequence( 1 )->id(),
			align_.sequence( 1 )->start()
		)
	);

	SequenceOP aligned_template(
		align_.sequence(template_index_)->clone()
	);

	SequenceOP t_align_seq(
		new Sequence(
			aligned_template->ungapped_sequence(),
			aligned_template->id() + "_align_seq",
			aligned_template->start()
		)
	);

	SequenceOP t_pdb_seq(
		new Sequence (
			template_pose_.sequence(),
			aligned_template->id() + "_pdb_seq",
			1
		)
	);

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

void ThreadingMover::apply(
	core::pose::Pose & query_pose
) {
	using core::Real;
	using core::Size;
	using std::string;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::id::SequenceMapping query_to_pdbseq = get_qt_mapping(query_pose);
	tr.Debug  << "mapping from " << query_index_ << " to " << template_index_
		<< std::endl;
	query_to_pdbseq.show( tr.Debug );

	std::string const template_id(
		utility::file_basename( template_pose_.pdb_info()->name() )
	);
	core::pose::add_score_line_string( query_pose, "template", template_id );

	core::Size n_copied( 0 );

	core::id::AtomID_Mask missing( true );
	core::pose::initialize_atomid_map( missing, query_pose ); // used for repacking atoms
	for ( Size resi = 1; resi <= query_pose.total_residue(); resi++ ) {
		Size const t_resi = query_to_pdbseq[ resi ];

		// skip this residue if we're not aligned
		if ( t_resi == 0 ) {
			continue;
		}

		if ( t_resi > template_pose_.total_residue() ) {
			tr.Error << "Error: don't have residue " << t_resi
				<< " in template_pose!" << std::endl;

			tr.Error << "template_pose.total_residue() == "
				<< template_pose_.total_residue() << std::endl;
			continue;
		}

		tr.Debug << "copying residue in query " << resi << " from template residue "
			<< t_resi << std::endl;
		for ( Size atomj = 1; atomj <= query_pose.residue(resi).natoms(); ++atomj ) {
			std::string const atom_name( query_pose.residue(resi).atom_name( atomj ));
			std::string t_atom_name = atom_name;

			//determine if we copy this atom: backbone and CB always, sidechain if identical aa
			// core::id::AtomID cb_id(
			// 	core::id::NamedAtomID( "CB", resi ),
			// 	query_pose,
			// 	false /*don't raise exception we check explicitly*/
			// );
			if ( !query_pose.residue(resi).atom_is_backbone(atomj)
				&& query_pose.residue(resi).name1() != template_pose_.residue(t_resi).name1()
				&& query_pose.residue(resi).atom_index("CB") != atomj
			) continue;

			// fix OXT/O ambiguity in template
			if ( !template_pose_.residue_type(t_resi).has( atom_name ) ) {
				if ( template_pose_.residue_type(t_resi).has("OXT") && atom_name == " O  " )
					t_atom_name = " OXT";
				else if ( template_pose_.residue_type(t_resi).has("O") && atom_name == " OXT" )
					t_atom_name = " O  ";
			}

			core::id::NamedAtomID query_id   ( atom_name, resi   );
			core::id::NamedAtomID template_id( t_atom_name, t_resi );

			// check to make sure that both poses have this atom_name
			if ( !query_pose.residue_type(resi).has( atom_name ) ) {
				tr.Warning 	<< "skipping atom,position " << atom_name << "," << resi
					<< " because query doesn't have atom " << atom_name << "."
					<< std::endl;
				continue;
			}
			if ( !template_pose_.residue_type(t_resi).has( t_atom_name ) ) {
				tr.Warning 	<< "skipping atom,position " << atom_name << "," << resi
					<< " because template doesn't have atom " << t_atom_name << "."
					<< std::endl;
				continue;
			}

			missing[ core::id::AtomID( atomj, resi ) ] = false;
			query_pose.set_xyz( query_id, template_pose_.xyz( template_id ) );
		} // for atom_i

		//fpd idealize all H
		core::conformation::Residue res_to_fix = query_pose.residue(resi);
		core::conformation::idealize_hydrogens( res_to_fix, query_pose.conformation() );
		query_pose.replace_residue( resi, res_to_fix, false );

		++n_copied;
	} // for resi

	tr.Debug << "Built threading model for sequence "
		<< query_pose.sequence() << std::endl;
	tr.Debug	<< "Copied " << n_copied << " / "
		<< query_pose.total_residue() << " from "
		<< align_.sequence(template_index_)->id() << std::endl;

	if ( randomize_loop_coords() ) {
		for ( core::Size ii = 1; ii <= query_pose.total_residue(); ++ii ) {
			core::conformation::ResidueOP iires = query_pose.residue( ii ).clone();
			core::conformation::idealize_hydrogens( *iires, query_pose.conformation() );
			query_pose.replace_residue( ii, *iires, false );
		}

		// randomize missing atoms
		randomize_selected_atoms(query_pose,missing);
	} // if randomize_loop_coords()

	if ( build_loops() ) {
		using protocols::loops::Loops;
		tr.Debug << "building query loops." << std::endl;
		loops::LoopsOP query_loops = loops_from_alignment(
			 query_pose.total_residue(), align_, min_loop_size()
		);
		query_loops->choose_cutpoints( query_pose );
		tr.Debug << query_loops << std::endl;

		if ( query_loops->size() > 0 ) {
			// switch to centroid ResidueTypeSet for loop remodeling
			std::string const orig_rsd_set_name(
				query_pose.residue_type(1).residue_type_set().name()
			);

			using core::util::switch_to_residue_type_set;
			core::util::switch_to_residue_type_set( query_pose, core::chemical::CENTROID );

			loops::loop_mover::LoopMoverOP loop_mover = protocols::loops::LoopMoverFactory::get_instance()->create_loop_mover(
				option[ cm::loop_mover ](), query_loops
			);
			for ( Size ii = 1; ii <= frag_libs().size(); ++ii ) {
				loop_mover->add_fragments( frag_libs()[ii] );
			}
			loop_mover->apply( query_pose );

			// switch back to original ResidueTypeSet after loop modeling
			if ( orig_rsd_set_name != core::chemical::CENTROID ) {
				core::util::switch_to_residue_type_set(
					query_pose,
					orig_rsd_set_name
				);
			}
		} else {
			 tr.Warning << "No loops found!" << std::endl;
		}
	} // build_query_loops

	// steal side chains
	//std::cout << "copying sidechains" << std::endl;
	//StealSideChainsMover sc_mover( template_pose_, query_to_pdbseq );
	//sc_mover.apply( query_pose );
	//std::cout << "finished copying sidechains" << std::endl;

	// repack the structure if specified by the user
	if ( repack_query() ) {
		using namespace core::scoring;
		ScoreFunctionOP scorefxn( getScoreFunction() );

		// repack missing sidechains
		core::id::AtomID_Mask missing( true );
		core::pose::initialize_atomid_map( missing, query_pose );
		tr.Debug << "repacking residues on pose with ScoreFunction: " << std::endl;
		scorefxn->show( tr.Debug );
		tr.Debug << std::endl << std::endl;
		core::pack::pack_missing_sidechains( query_pose, missing );

		tr.Debug << "setting up ideal hydrogen geometry on all residues."
			<< std::endl;
		for ( core::Size ii = 1; ii <= query_pose.total_residue(); ++ii ) {
			core::conformation::ResidueOP iires = query_pose.residue( ii ).clone();
			core::conformation::idealize_hydrogens( *iires, query_pose.conformation() );
			query_pose.replace_residue( ii, *iires, false );
		}

		tr.Debug << "optimizing hydrogen placement with the packer."
			<< std::endl;
		core::pack::optimize_H_and_notify( query_pose, missing );

		scorefxn->set_weight( core::scoring::peptide_bond, 1.0 );
		(*scorefxn)(query_pose);
		scorefxn->show( tr.Debug, query_pose );
	} // repack_query

	tr.flush();
} // apply

std::string
ThreadingMover::get_name() const {
	return "ThreadingMover";
}

protocols::moves::MoverOP
ThreadingMover::clone() const
{
	return new ThreadingMover(*this);
}

protocols::moves::MoverOP
ThreadingMover::fresh_instance() const
{
	return new ThreadingMover();
}

} // comparative_modeling
} // protocols
