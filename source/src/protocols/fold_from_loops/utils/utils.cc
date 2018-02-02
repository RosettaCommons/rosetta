// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utils.cc
/// @brief  Helper functions for nub initio
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

// Unit headers
#include <protocols/fold_from_loops/utils/utils.hh>
#include <protocols/fold_from_loops/selectors/CutpointResidueSelector.hh>
#include <protocols/fold_from_loops/selectors/ConstraintResidueSelector.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <iterator>
#include <algorithm>

namespace protocols {
namespace fold_from_loops {
namespace utils {

static basic::Tracer TR( "protocols.fold_from_loops.utils", basic::t_trace );

/// @brief It will take each range that is not the first or the last one
/// and split it into two, following the logic of ```find_cutpoint_from_secondary_structure```
/// Input Ranges are not modified, a new ResidueRanges is returned.
/// WARNING1: When splitting a range of 1 residue, it splits into 1 range of 1 and another from 0 to 0, so
/// that it is easy to filter afterwards ( for the purpose of NubInitio, this means adding and empty pose
/// to the vector, for example ).
/// WARNING2: If the first ResidueRange does not start with 1, an empty_range (0,0) is set at the begining too.
/// WARNING3: If the final number of ranges is odds, we will assume that the last insertion is C-terminal and
/// add a final empty range at the end.
core::select::residue_selector::ResidueRanges
split_mid_ranges(
	std::string structure,
	core::select::residue_selector::ResidueRanges const & ranges )
{
	core::select::residue_selector::ResidueRanges new_ranges;
	core::select::residue_selector::ResidueRange empty_range( 0, 0 );
	core::Size start = 2;
	if ( ranges[1].start() != 1 ) { // Means that the first insertion is at N-term
		new_ranges.push_back( empty_range );
		start = 1;
	} else {
		new_ranges.push_back( ranges[1] );
	}
	for ( core::Size i = start; i < ranges.size(); ++i ) {
		core::select::residue_selector::ResidueRange range = ranges[i];
		std::string wstruct = structure.substr( range.start() - 1, range.stop() - range.start() + 1 );
		TR << "STR: " << wstruct << std::endl;
		TR << "FULL:  " << range.to_string() << std::endl;
		if ( wstruct.size() == 1 ) {
			new_ranges.push_back( range );
			new_ranges.push_back( empty_range );
		} else {
			core::Size cutpoint = find_cutpoint_from_secondary_structure( wstruct ) + range.start();
			core::select::residue_selector::ResidueRange range1( range.start(), cutpoint );
			TR << "SPLIT: " << range1.to_string() << std::endl;
			new_ranges.push_back( range1 );
			core::select::residue_selector::ResidueRange range2( cutpoint + 1, range.stop() );
			TR << "SPLIT: " << range2.to_string() << std::endl;
			new_ranges.push_back( range2 );
		}
	}
	new_ranges.push_back( ranges[ ranges.size() ] );
	return new_ranges;
}

/// @brief Given a Secondary Structure definition as a string, it finds what
/// it considers the optimal cutpoint.
/// This means that, if there is no loop 'L' segment it will return the middle
/// point as a residue; otherwise it will try to split the larger loop 'L' region.
core::Size
find_cutpoint_from_secondary_structure( std::string structure )
{
	// If there is only one or two residues, we return the firts.
	if ( structure.size() < 2 ) {
		return 0;
	}

	core::select::residue_selector::ResidueSubset boolstruct;
	for ( auto const & c : structure ) { boolstruct.push_back( c == 'L'); }
	core::select::residue_selector::ResidueRanges ranges( boolstruct );
	core::Size midpoint = structure.size() / 2;

	TR << "STR:  " << structure << std::endl;
	TR << "BOOL: " << boolstruct << std::endl;

	// If there is no 'L' secondary structure, just return the middle of the string.
	if ( ranges.size() == 0 ) {
		return midpoint;
	} else if ( ranges.size() <= 2 ) {
		// Check out when we have only one or two 'L' patches
		core::Size i = 1;
		if ( ranges.size() > 1 ) {
			if ( ( ranges[2].stop() - ranges[2].start() + 1 ) > ( ranges[1].stop() - ranges[1].start() + 1 ) ) {
				i = 2;
			}
		}
		return ( ranges[i].start() + ( ranges[i].stop() - ranges[i].start() + 1 ) / 2 ) -  2;
	} else {
		// When we have more than two 'L' patches
		std::vector<int> sizes;
		for ( core::Size i = 2; i < ranges.size() ; ++i ) {
			sizes.push_back( ranges[i].stop() - ranges[i].start() + 1 );
		}
		// +2: +1 for vector1 count, +1 because position 1 is not there
		core::Size max = std::distance( sizes.begin(), std::max_element(sizes.begin(), sizes.end() ) ) + 2 ;
		return ( ranges[max].start() + ( ranges[max].stop() - ranges[max].start() + 1 ) / 2 ) -  2;
	}

	return 1;
}

/// @brief To a given scaffold it attaches a unfolded version of a pose to the N-terminal and another to the C-terminal
/// Generates a FoldTree that goes from the middle of 'scaffold' towards both directions
void
attach_n_and_c_unfolded_poses_to_pose( core::pose::Pose const & n_insert, core::pose::Pose & scaffold, core::pose::Pose const & c_insert )
{
	using namespace core::chemical;
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ));
	core::Size kmidpoint = scaffold.size() / 2;
	if ( kmidpoint == 0 ) {
		kmidpoint = 1;
	}
	core::Size key = n_insert.size() + kmidpoint;
	TR << "inserting n_pose: " << n_insert.secstruct() << std::endl;
	attach_unfolded_pose_to_pose_n_term( n_insert, scaffold, rsd_set );
	// scaffold.pdb_info()->copy(*(n_insert.pdb_info()), 1, n_insert.size(), 1);
	// core::Size current_length = scaffold.size();
	TR << "inserting n_pose: " << c_insert.secstruct() << std::endl;
	attach_unfolded_pose_to_pose_c_term( c_insert, scaffold, rsd_set );
	// scaffold.pdb_info()->copy(*(c_insert.pdb_info()), 1, c_insert.size(), current_length + 1);
	scaffold.pdb_info()->obsolete(false);
	core::kinematics::FoldTree fold_tree;
	fold_tree.add_edge( key, 1, core::kinematics::Edge::PEPTIDE );
	fold_tree.add_edge( key, scaffold.size(), core::kinematics::Edge::PEPTIDE );
	scaffold.fold_tree( fold_tree );
	TR << scaffold.sequence() << std::endl;
	TR << scaffold.secstruct() << std::endl;
	TR << fold_tree << std::endl;
}

/// @brief adds the 'insert' unfolded pose to the N-terminal of the provided pose
void
attach_unfolded_pose_to_pose_n_term( core::pose::Pose const & insert, core::pose::Pose & scaffold, core::chemical::ResidueTypeSetCOP rsd_set )
{
	if ( insert.size() > 0 ) {
		using namespace core::conformation;
		core::pose::Pose insert_local( insert );
		core::pose::remove_variant_type_from_pose_residue( insert_local, core::chemical::UPPER_TERMINUS_VARIANT, insert_local.size() );
		for ( core::Size i = insert_local.size(); i >= 1; --i ) {
			TR << "Prepending from insert pose residue: " << std::to_string( i ) << " " << insert.residue(i).name3() << std::endl;
			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( insert.residue(i).name1() ) );
			ResidueOP new_res( ResidueFactory::create_residue( *new_rsd_type ) );
			scaffold.conformation().safely_prepend_polymer_residue_before_seqpos( *new_res, 1, true );
			scaffold.set_secstruct( 1, insert.secstruct(i) );
			scaffold.set_omega( 1, 180.0 );
		}
		scaffold.pdb_info()->copy(*(insert.pdb_info()), 1, insert.size(), 1);
	} else {
		TR << "Insertion pose size is 0! Skip!" << std::endl;
	}
}

/// @brief adds the 'insert' unfolded pose to the C-terminal of the provided pose
void
attach_unfolded_pose_to_pose_c_term( core::pose::Pose const & insert, core::pose::Pose & scaffold, core::chemical::ResidueTypeSetCOP rsd_set )
{
	if ( insert.size() > 0 ) {
		using namespace core::conformation;
		core::pose::Pose insert_local( insert );
		core::pose::remove_variant_type_from_pose_residue( insert_local, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
		core::Size current_length = scaffold.size();
		for ( core::Size i = 1; i <= insert_local.size(); ++i ) {
			TR << "Appending from insert pose residue: " << std::to_string( i ) << " " << insert.residue(i).name3() << std::endl;
			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( insert.residue(i).name1() ) );
			ResidueOP new_res( ResidueFactory::create_residue( *new_rsd_type ) );
			scaffold.conformation().safely_append_polymer_residue_after_seqpos( *new_res, scaffold.size(), true );
			scaffold.set_secstruct( scaffold.size() - 1, insert.secstruct(i) );
			scaffold.set_omega( scaffold.size() - 1, 180.0 );
		}
		scaffold.pdb_info()->copy(*(insert.pdb_info()), 1, insert.size(), current_length + 1);
	} else {
		TR << "Insertion pose size is 0! Skip!" << std::endl;
	}
}

/// @brief Basically works as 'core::pose::append_pose_to_pose' but the FoldTree is build by joining by jump the
/// root of 'insert' FoldTree to that of 'scaffold' FoldTree.
void
append_pose_to_pose_keep_fold_tree( core::pose::Pose & scaffold, core::pose::Pose const & insert, bool new_chain )
{
	if ( scaffold.empty() ) {
		scaffold = insert;
	} else {
		core::Size current_length = scaffold.size();
		core::kinematics::FoldTree keeptree  = scaffold.fold_tree();
		keeptree.insert_fold_tree_by_jump( insert.fold_tree(), scaffold.size() + 1, keeptree.root(), keeptree.root() );
		core::pose::append_pose_to_pose( scaffold, insert, new_chain );
		scaffold.fold_tree( keeptree );
		scaffold.pdb_info()->copy(*(insert.pdb_info()), 1, insert.size(), current_length + 1);
		scaffold.pdb_info()->obsolete(false);
		std::string sse = insert.secstruct();
		for ( core::Size sr = current_length + 1, ss=0; sr <= scaffold.size(); ++sr, ++ss ) {
			scaffold.set_secstruct( sr, sse[ss] );
		}
		TR << "FINAL TREE "<< scaffold.fold_tree() << std::endl;
	}
}

/// @brief Creates a sequence mapping between two proteins assuming that the False selections on a ResidueSubset
/// marks the residues that are the same for both of them. Assumes the same number of non-selected patches for each pose.
core::id::SequenceMapping
map_by_residue_subsets(
	core::pose::Pose const & p1,
	core::select::residue_selector::ResidueSubset const & r1,
	core::pose::Pose const & p2,
	core::select::residue_selector::ResidueSubset const & r2 )
{
	core::id::SequenceMapping seqmap;
	if ( not core::select::residue_selector::are_selections_equal( r1, r2 ) ) {
		core::sequence::SequenceOP refseq( new core::sequence::Sequence( p1.sequence(), "refseq" ) );
		core::sequence::SequenceOP newseq( new core::sequence::Sequence( p2.sequence(), "newseq" ) );

		core::select::residue_selector::ResidueRanges R1( r1 );
		core::select::residue_selector::ResidueRanges R2( r2 );

		runtime_assert_msg( R1.size() == R2.size(), "The number of putative different patches should be the same.");

		TR << "r1" << r1 << std::endl;
		TR << "r2" << r2 << std::endl;
		core::Size refseq_ins = 0, newseq_ins = 0;
		for ( core::Size i = 1; i <= R1.size(); ++i ) {
			TR << "R1 " << std::to_string(i) << " " << R1[i].to_string() << std::endl;
			TR << "R2 " << std::to_string(i) << " " << R2[i].to_string() << std::endl;
			core::Size ref_size = R1[i].stop() - R1[i].start() + 1;
			core::Size new_size = R2[i].stop() - R2[i].start() + 1;
			if ( ref_size < new_size ) {
				core::Size mid = R1[i].start() + ( ref_size / 2 );
				for ( core::Size j = 0; j < (new_size - ref_size); ++j ) {
					refseq->insert_gap( mid + refseq_ins );
				}
				refseq_ins = new_size - ref_size;
			} else if ( ref_size > new_size ) {
				core::Size mid = R2[i].start() + ( new_size / 2 );
				for ( core::Size j = 0; j < (ref_size - new_size); ++j ) {
					newseq->insert_gap( mid + newseq_ins );
				}
				newseq_ins = ref_size - new_size;
			}
		}

		core::sequence::SequenceAlignment seqal = core::sequence::SequenceAlignment();
		seqal.add_sequence( refseq );
		seqal.add_sequence( newseq );
		seqmap = seqal.sequence_mapping( 1, 2 );
	} else {
		TR << "Global Size and disposition of template residues has not changed. Mapping is identity" << std::endl;
		seqmap = core::id::SequenceMapping().identity( p2.size() );
	}
	TR << seqmap.to_string() << std::endl;
	return seqmap;
}

void
report_unfolded( core::pose::Pose const & pose, core::kinematics::MoveMapOP movemap )
{
	if ( TR.visible() ) {
		using namespace core::select::residue_selector;

		std::string sequence  = pose.sequence();
		std::string structure = pose.secstruct();
		std::set< std::string > labels;
		for ( core::Size i = 1; i <= sequence.size(); i++ ) {
			utility::vector1< std::string > rlabs = pose.pdb_info()->get_reslabels( i );
			for ( auto rlab : rlabs ) {
				labels.emplace( rlab );
			}
		}
		std::map< std::string, std::string > strlabels;
		for ( auto label : labels ) {
			strlabels[label] = "";
		}
		for ( core::Size i = 1; i <= sequence.size(); i++ ) {
			utility::vector1< std::string > rlabs = pose.pdb_info()->get_reslabels( i );
			for ( auto label : labels ) {
				if ( std::find(rlabs.begin(), rlabs.end(), label) != rlabs.end() ) {
					strlabels[label].append( "*" );
				} else {
					strlabels[label].append( " " );
				}
			}
		}

		ResidueSelectorOP chainbreaks( new selectors::CutpointResidueSelector );
		std::string chainbreaks_repr = represent_residue_selector( chainbreaks->apply( pose ) );
		bool chainbreaks_show = chainbreaks_repr.find("*")!=std::string::npos;
		ResidueSelectorOP constraints( new selectors::ConstraintResidueSelector );
		std::string constraint_repr = represent_residue_selector( constraints->apply( pose ) );
		bool constraint_show = constraint_repr.find("*")!=std::string::npos;

		std::string bb  = "";
		std::string chi = "";
		for ( core::Size i = 1; i <= sequence.size(); i++ ) {
			if ( movemap->get_bb ( i ) ) {
				bb.append( "*" );
			} else {
				bb.append( " " );
			}
			if ( movemap->get_chi ( i ) ) {
				chi.append( "*" );
			} else {
				chi.append( " " );
			}
		}

		core::Size margin = 15;
		TR << std::left << std::setw( margin ) << "SEQUENCE";
		TR << sequence << std::endl;
		TR << std::left << std::setw( margin ) << "STRUCTURE";
		TR << structure << std::endl;
		for ( auto label : labels ) {
			TR << std::left << std::setw( margin ) << label;
			TR << strlabels[label] << std::endl;
		}
		if ( constraint_show ) {
			TR << std::left << std::setw( margin ) << "CONSTRAINTS";
			TR << constraint_repr << std::endl;
		}
		if ( chainbreaks_show ) {
			TR << std::left << std::setw( margin ) << "CUTPOINT";
			TR << chainbreaks_repr << std::endl;
		}
		TR << std::left << std::setw( margin ) << "MVMP_BB";
		TR << bb << std::endl;
		TR << std::left << std::setw( margin ) << "MVMP_CHI";
		TR << chi << std::endl;
		TR << pose.fold_tree() << std::endl;
		TR << pose.annotated_sequence() << std::endl;
		pose.pdb_info()->show(TR);
	}
}

}
}
}
