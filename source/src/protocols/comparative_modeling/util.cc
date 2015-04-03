// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/comparative_modeling/util.cc
/// @brief set of utilities used in comparative modeling of protein structures
/// @author James Thompson

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

// Symmetry
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/AlignerFactory.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/AlignmentSet.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/LoopMoverFactory.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/sequence/Aligner.hh>
#include <core/sequence/ScoringScheme.hh>

#include <numeric/random/random.hh>

#include <core/import_pose/import_pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>

namespace protocols {
namespace comparative_modeling {

#define NO_LOOP_SIZE_CST 0

using utility::vector1;
using std::string;

static thread_local basic::Tracer tr( "protocols.comparative_modeling.util" );

/// @detail The premise underlying this tortuous method is simple--
/// identify aligned/unaligned regions in a sequence alignment with
/// the constraint that each region has a certain minimum length.
///
/// The current implementation achieves this goal in a roundabout
/// manner by making use of existing, less specialized utility
/// functions.
void bounded_loops_from_alignment(
    const core::Size num_residues,
    const core::Size min_size,
    const core::sequence::SequenceAlignment& alignment,
    protocols::loops::LoopsOP & unaligned_regions) {

  using core::Size;
  using core::id::SequenceMapping;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  assert(unaligned_regions);

  const Size query_idx = 1;
  const Size templ_idx = 2;
  SequenceMapping mapping(alignment.sequence_mapping(query_idx, templ_idx));

  vector1<Size> unaligned_residues;
  for (Size resi = 1; resi <= num_residues; resi++) {
    Size t_resi = mapping[resi];

    bool const gap_exists(
      t_resi == 0 || // query residue maps to a gap
      (resi > 1 && mapping[ resi - 1 ] != t_resi - 1) ||            // last residue was gapped
      (resi < num_residues && mapping[ resi + 1 ] != t_resi + 1));  // next residue is gapped

    if (gap_exists) unaligned_residues.push_back( resi );
  }

  // Ensure the unaligned regions meet size constraints.
  // Aligned regions are incorrect at this point.
  protocols::loops::LoopsOP unaligned_ok = pick_loops_unaligned(num_residues, unaligned_residues, min_size);

  // Ensure the aligned regions meet size constraints.
  unaligned_residues.clear();
  for (Loops::const_iterator i = unaligned_ok->begin(); i != unaligned_ok->end(); ++i) {
    const Loop& loop = *i;
    for (Size j = loop.start(); j <= loop.stop(); ++j) {
      unaligned_residues.push_back(j);
    }
  }

  vector1<Size> bounded_unaligned_residues(unaligned_residues);
  for (Size i = 2; i <= unaligned_residues.size(); ++i) {
    Size prev_residue = unaligned_residues[i - 1];
    Size curr_residue = unaligned_residues[i];

    // Length of the unaligned region is (curr - 1) - (prev + 1) + 1,
    Size delta = curr_residue - prev_residue - 1;
    if (delta == 0 || delta >= min_size)
      continue;

    for (Size j = (prev_residue + 1); j <= (curr_residue - 1); ++j)
      bounded_unaligned_residues.push_back(j);
  }
  std::sort(bounded_unaligned_residues.begin(), bounded_unaligned_residues.end());

  // Retrieve loops without affecting unaligned region length
  unaligned_regions = pick_loops_unaligned(num_residues, bounded_unaligned_residues, NO_LOOP_SIZE_CST);
}

protocols::loops::LoopsOP loops_from_alignment(
	core::Size nres,
	core::sequence::SequenceAlignment const & aln,
	core::Size const min_loop_size
) {
	using core::Size;

	Size const query_idx( 1 );
	Size const templ_idx( 2 );
	core::id::SequenceMapping mapping_(
		aln.sequence_mapping( query_idx, templ_idx )
	);
	tr.Debug << "called loops_from_alignment with arguments:" << std::endl;
	tr.Debug << nres << std::endl;
	tr.Debug << aln << std::endl;
	tr.Debug << min_loop_size << std::endl;
	mapping_.show( tr.Debug );
	vector1< core::Size > unaligned_residues;
	for ( Size resi = 1; resi <= nres; resi++ ) {
		Size t_resi = mapping_[ resi ];

		// gap checks
		bool const gap_exists(
			t_resi == 0 || // query residue maps to a gap
			( resi > 1    && mapping_[ resi - 1 ] != t_resi - 1 ) || // last residue was gapped
			( resi < nres && mapping_[ resi + 1 ] != t_resi + 1 ) // next residue is gapped
		);
		if ( gap_exists ) unaligned_residues.push_back( resi );
	}

	tr.flush_all_channels();

	return pick_loops_unaligned(
		nres,
		unaligned_residues,
		min_loop_size
	);
}


//fpd  build a loopfile from the intersection of loops from multiple aln files
protocols::loops::LoopsOP loops_from_transitive_alignments(
	core::Size nres1,
	core::sequence::SequenceAlignment const & aln1,
	core::Size nres2,
	core::sequence::SequenceAlignment const & aln2,
	core::Size const min_loop_size
) {
	using core::Size;

	Size const query_idx( 1 );
	Size const templ_idx( 2 );
	core::id::SequenceMapping mapping1_( aln1.sequence_mapping( query_idx, templ_idx ) );
	core::id::SequenceMapping mapping2_( aln2.sequence_mapping( query_idx, templ_idx ) );
	tr.Debug << "called loops_from_multiple_alignments with arguments:" << std::endl;
	tr.Debug << nres1 << std::endl;
	tr.Debug << aln1 << std::endl;
	tr.Debug << nres1 << std::endl;
	tr.Debug << aln2 << std::endl;
	tr.Debug << min_loop_size << std::endl;
	mapping1_.show( tr.Debug );
	mapping2_.show( tr.Debug );
	vector1< core::Size > unaligned_residues;
	for ( Size resi = 1; resi <= nres1; resi++ ) {
		Size t_resi1 = mapping1_[ resi ];

		// gap checks
		//fpd  First check to see if there is a gap in the alignment
		bool gap_exists =
			t_resi1 == 0 || // query residue maps to a gap (aln1)
			( resi > 1    && mapping1_[ resi - 1 ] != t_resi1 - 1 ) || // last residue was gapped
			( resi < nres1 && mapping1_[ resi + 1 ] != t_resi1 + 1 ); // next residue is gapped

		//fpd Now check if there is this maps to a part of the template sequence
		//fpd  that is missing in the input template PDB
		if (!gap_exists) {
			Size t_resi2 = mapping2_[ t_resi1 ];
			gap_exists = t_resi2 == 0 || // query residue maps to a gap (aln2)
				( t_resi1 > 1    && mapping2_[ t_resi1 - 1 ] != t_resi2 - 1 ) || // last residue was gapped
				( t_resi1 < nres2 && mapping2_[ t_resi1 + 1 ] != t_resi2 + 1 ); // next residue is gapped
		}

		if ( gap_exists )
			unaligned_residues.push_back( resi );
	}

	tr.flush_all_channels();

	return pick_loops_unaligned(
		nres1,
		unaligned_residues,
		min_loop_size
	);
}

protocols::loops::LoopsOP pick_loops_unaligned(
	core::Size nres,
	utility::vector1< core::Size > const & unaligned_residues,
	core::Size min_loop_size
) {
	typedef core::Size Size;

	protocols::loops::LoopsOP query_loops( new protocols::loops::Loops() );
	if ( unaligned_residues.size() == 0 ) {
		tr.Warning << "No unaligned residues, no loops found." << std::endl;
		return query_loops;
	}

	Size loop_stop ( *unaligned_residues.begin() );
	Size loop_start( *unaligned_residues.begin() );

	for ( vector1< Size >::const_iterator it = unaligned_residues.begin(),
			next = it + 1,
			end  = unaligned_residues.end();
			next != end; ++it, ++next
	) {
		tr.Debug << "residue " << *it << " is unaligned." << std::endl;
		if ( *next - *it > 1 ) {
			// add loop
			loop_stop = *it;
			while ( (loop_stop - loop_start + 1) < min_loop_size ) {
				if ( loop_stop < nres )
					++loop_stop;
				if ( loop_start > 1 && (loop_stop - loop_start + 1) < min_loop_size )
					--loop_start;
			}
			tr.Debug << "adding loop from " << loop_start << " to " << loop_stop
				<< std::endl;
			protocols::loops::Loop loop( loop_start, loop_stop, 0, 0, false );
			query_loops->add_loop( loop, 1 );

			loop_start = *next;
		}
	}

	loop_stop = ( *(unaligned_residues.end() - 1) );

	while ( (loop_stop - loop_start + 1) < min_loop_size ) {
		if ( loop_stop < nres ) ++loop_stop;
		if ( loop_start > 1 ) --loop_start;
	}
	tr.Debug << "adding loop from " << loop_start << " to " << loop_stop
		<< std::endl;
	protocols::loops::Loop loop( loop_start, loop_stop, 0, 0, false );
	query_loops->add_loop( loop , 1 );

	tr.flush_all_channels();

	return query_loops;
} // pick_loops

protocols::loops::LoopsOP pick_loops_chainbreak(
	core::pose::Pose & query_pose,
	core::Size min_loop_size
) {
	typedef core::Size Size;

	core::Real const chainbreak_cutoff( 4.0 );
	core::Size nres = query_pose.total_residue();

	//fpd symm
	if ( core::pose::symmetry::is_symmetric(query_pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( query_pose.conformation()) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}

	vector1< Size > residues_near_chainbreak;
	for ( Size i = 1; i <= nres - 1; ++i ) {
		if ( query_pose.residue_type(i).is_protein() &&  query_pose.residue_type(i+1).is_protein()) {
			core::Real dist = query_pose.residue(i).xyz("CA").distance(
				query_pose.residue(i+1).xyz("CA")
			);
			//std::cout << "dist(" << i << "," << i+1 << ") = " << dist << std::endl;
			if ( dist > chainbreak_cutoff ) {
				residues_near_chainbreak.push_back( i );
			}
		}
	} // for ( Size i )

	if ( residues_near_chainbreak.size() == 0 ) {
		tr.Warning << "No chainbreaks found, so not picking any loops!"
			<< std::endl;
	}

	tr.flush();

	return pick_loops_unaligned(
		query_pose.total_residue(),
		residues_near_chainbreak,
		min_loop_size
	);
} // pick_loops

void rebuild_loops_until_closed(
	core::pose::Pose & query_pose,
	core::Size const min_loop_size,
	core::Size const max_rebuild,
	std::string const & loop_mover_name
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// switch to centroid ResidueTypeSet for loop remodeling
	protocols::loops::LoopsOP my_loops = pick_loops_chainbreak(
		query_pose,
		min_loop_size
	);

	if ( my_loops->size() == 0 ) {
		tr.Debug << "no loops found." << std::endl;
		return;
	}

	std::string const orig_rsd_set_name(
		query_pose.residue_type(1).residue_type_set().name()
	);
	core::util::switch_to_residue_type_set( query_pose, core::chemical::CENTROID );

	bool closed( false );
	for ( core::Size iter = 1; !closed && iter <= max_rebuild; iter++ ) {
		loops::loop_mover::LoopMoverOP loop_mover = protocols::loops::LoopMoverFactory::get_instance()->create_loop_mover(
			loop_mover_name, my_loops
		);
		loop_mover->apply( query_pose );

		my_loops = pick_loops_chainbreak(
			query_pose,
			min_loop_size
		);

		if ( my_loops->size() == 0 ) {
			tr.Debug << "closed loops on iteration " << iter << " ." << std::endl;
			closed = true;
		}
	}

	tr.flush();

	core::util::switch_to_residue_type_set( query_pose, orig_rsd_set_name );
} // rebuild_loops_until_closed

void steal_ligands(
	core::pose::Pose & dest_pose,
	core::pose::Pose const & source_pose_in,
	core::id::NamedAtomID const anchor_atom_dest,
	core::id::NamedAtomID const anchor_atom_source,
	utility::vector1< core::id::NamedAtomID > const ligand_indices
) {
	using core::Size;
	using utility::vector1;

	// add some runtime asserts here!
	if ( !anchor_atom_dest.valid() ) {
		tr.Error << "Error: can't place ligands. "
			<< "Destination anchor atom is not valid!" << anchor_atom_dest
			<< std::endl;
		return;
	}
	if ( !anchor_atom_source.valid() ) {
		tr.Error << "Error: can't place ligands. "
			<< "Source anchor atom is not valid! (" << anchor_atom_source << ")"
			<< std::endl;

		return;
	}

	// create a copy to avoid modifying original
	core::pose::Pose source_pose = source_pose_in;

	// set up FoldTree for source_pose that has the jump orientation that we want
	core::kinematics::FoldTree new_fold_tree;
	core::Size old_fold_tree_end(
		source_pose.total_residue() - ligand_indices.size()
	); // stupid assumption!
	new_fold_tree.add_edge(
		1,
		anchor_atom_source.rsd(),
		core::kinematics::Edge::PEPTIDE
	);
	new_fold_tree.add_edge(
		anchor_atom_source.rsd(),
		old_fold_tree_end,
		core::kinematics::Edge::PEPTIDE
	);
	tr.Debug << "adding ligand residues to fold-tree" << std::endl;

	// add edges from anchor to ligand residues
	for ( Size jj = 1; jj <= ligand_indices.size(); ++jj ) {
		tr.Error << "adding " << jj << std::endl;
		core::kinematics::Edge out_edge(
			anchor_atom_source.rsd(), // start
			ligand_indices[jj].rsd(), // stop
			static_cast< int > (jj), // label
			anchor_atom_source.atom(), // start_atom
			ligand_indices[jj].atom(), // stop_atom
			false // bKeepStubInResidue
		);
		tr.Error << out_edge << std::endl;
		new_fold_tree.add_edge( out_edge );
	}

	tr.Error << source_pose.fold_tree();
	source_pose.fold_tree( new_fold_tree );

	// copy the residues from the source_pose into the dest_pose
	// using the jump geometry defined above.
	for ( Size jj = 1; jj <= ligand_indices.size(); ++jj ) {
		dest_pose.append_residue_by_jump(
			source_pose.residue( ligand_indices[jj].rsd() ),
			anchor_atom_dest.rsd(),
			anchor_atom_dest.atom(),
			ligand_indices[jj].atom()
		);
		dest_pose.set_jump(
			static_cast< int > (jj),
			source_pose.jump( static_cast< int > (jj) )
		);
	}

	tr.flush();
} // steal_ligands

void initialize_ss( core::pose::Pose & pose ) {
	using namespace core::pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	bool psipred_ss2_ok = loops::set_secstruct_from_psipred_ss2( pose );
	if ( !psipred_ss2_ok ) {
		std::string dssp_name( option[ in::file::dssp ]().name() );
		bool dssp_ok = loops::set_secstruct_from_dssp(pose, dssp_name);
		if ( !dssp_ok ) {
			set_ss_from_phipsi( pose );
		}
	}
}

utility::vector1< core::pose::Pose >
templates_from_cmd_line() {
	using std::string;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	vector1< string > template_pdb_fns(
		option[ in::file::template_pdb ]()
	);
	vector1< core::pose::Pose > template_poses
		= core::import_pose::poses_from_pdbs( template_pdb_fns );

	return template_poses;
}

bool loops_are_closed( core::pose::Pose & pose ) {
	return ( pick_loops_chainbreak(pose, 3)->size() == 0 ) ;
}

std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::file::file_exists;
	using core::import_pose::pose_from_pdb;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
	map< string, Pose > poses;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		if ( file_exists(*it) ) {
			Pose pose;
			core::import_pose::pose_from_pdb( pose, *rsd_set, *it );
			string name = utility::file_basename( *it );
			name = name.substr( 0, 5 );
			poses[name] = pose;
		}
	}

	return poses;
}

AlignmentSet
alignments_from_cmd_line() {
	using core::Real;
	using std::string;
	using utility::vector1;
	using utility::file::FileName;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::sequence;

	// options set up
	FileName fn1( option[ in::file::pssm ]()[1] );
	FileName fn2( option[ in::file::pssm ]()[2] );
	string const aligner_type( option[ cm::aligner ]() );
	string const seq_score( option[ cm::seq_score ]()[1] );
	Real const min_gap_open( option[ cm::min_gap_open ]() );
	Real const max_gap_open( option[ cm::max_gap_open ]() );
	Real const min_gap_extend( option[ cm::min_gap_extend ]() );
	Real const max_gap_extend( option[ cm::max_gap_extend ]() );
	Real const step_size( 0.5 ); // maybe make this an option?

	runtime_assert( min_gap_open <= max_gap_open );
	runtime_assert( min_gap_extend <= max_gap_extend );

	// setup objects
	ScoringSchemeFactory ssf;
	AlignerOP aligner( AlignerFactory::get_aligner( aligner_type ) );
	ScoringSchemeOP ss( ssf.get_scoring_scheme( seq_score ) );

	SequenceProfileOP prof1( new SequenceProfile );
	prof1->read_from_file( fn1 );
	prof1->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()

	SequenceProfileOP prof2( new SequenceProfile );
	prof2->read_from_file( fn2 );
	prof2->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()

	// eliminate leading paths from prof1 and prof2
	prof1->id( FileName( prof1->id() ).base() );
	prof2->id( FileName( prof2->id() ).base() );

	AlignmentSet set;
	for ( Real o = min_gap_open; o <= max_gap_open; o += step_size ) {
		for ( Real e = min_gap_extend; e <= max_gap_extend;
					e += step_size
		) {
			ss->gap_open  ( o );
			ss->gap_extend( e );

			SequenceAlignment align = aligner->align( prof1, prof2, ss );
			set.insert( align );
		} // g_extend
	} // g_open

	// add i/o of alignments from files here

	return set;
} // alignments_from_cmd_line


void randomize_selected_atoms(
	core::pose::Pose & query_pose,
	core::id::AtomID_Mask const & selected
) {
	using core::Size;
	for ( Size pos = 1; pos <= query_pose.total_residue(); ++pos ) {
		Size atomj( 1 );
		for ( core::id::AtomID_Mask::AtomMap::const_iterator
				it = selected[ pos ].begin(), eit = selected[ pos ].end(); it != eit;
				++it, ++atomj
		) {

			if ( query_pose.residue( pos ).atom_is_hydrogen( atomj ) ) continue;
			if ( *it ) { //entry is missng == true
				core::Vector ai(
					900.000 + numeric::random::rg().uniform()*100.000,
					900.000 + numeric::random::rg().uniform()*100.000,
					900.000 + numeric::random::rg().uniform()*100.000
				);
				query_pose.set_xyz( core::id::AtomID( atomj, pos ), ai );
				//now randomize also attached hydrogens
				for ( Size atom_nr = query_pose.residue( pos ).attached_H_begin( atomj );
							atom_nr <= query_pose.residue( pos ).attached_H_end( atomj ); ++atom_nr ) {
					core::Vector ai(
							900.000 + numeric::random::rg().uniform()*100.000,
							900.000 + numeric::random::rg().uniform()*100.000,
							900.000 + numeric::random::rg().uniform()*100.000
					);
					query_pose.set_xyz( core::id::AtomID( atom_nr, pos ), ai );
				}
			}
		}
	} // for selected atoms
} // randomize_selected_atoms

} // comparative_modeling
} // protocols
