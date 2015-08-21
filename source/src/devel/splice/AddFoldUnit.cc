// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AddFoldUnitMover.cc
/// @brief

// Unit headers
#include <devel/splice/AddFoldUnit.hh>
#include <devel/splice/Splice.hh>
#include <devel/splice/AddFoldUnitCreator.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR("devel.splice.AddFoldUnitMover");

#include <utility/tag/Tag.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <string>
#include <utility/string_util.hh>
//#include <sstream>
#include <core/pose/util.hh>
#include <utility/io/izstream.hh>
#include <boost/foreach.hpp>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <numeric/random/random.hh>
#include <boost/unordered/unordered_map.hpp>

namespace devel {
namespace splice {
using namespace protocols;
using namespace core;
using namespace std;
using utility::vector1;

/// @brief return a subset of entries that match all of the selection criteria
vector1< Size >
FoldUnitUtils::entry_subset() const{
	typedef boost::unordered_multimap< Size, Size > UM;
	utility::vector1< Size > bb_dof_entries;
	bb_dof_entries.clear();
	pair< UM::const_iterator, UM::const_iterator > n_entries( entry_pairs_quick_access_N_C_.cbegin(), entry_pairs_quick_access_N_C_.cbegin() );
	pair< UM::const_iterator, UM::const_iterator > c_entries( entry_pairs_quick_access_C_N_.cbegin(), entry_pairs_quick_access_C_N_.cbegin()  );
	if ( n_term_entry() ) {
		n_entries = entry_pairs_quick_access_N_C_.equal_range( n_term_entry() );
	}
	if ( c_term_entry() ) {
		c_entries = entry_pairs_quick_access_C_N_.equal_range( c_term_entry() );
	}

	if ( n_term_entry() && !c_term_entry() ) {
		for ( UM::const_iterator itn = n_entries.first; itn != n_entries.second; ++itn ) {
			bb_dof_entries.push_back( itn->second );
		}
		TR<<"Found "<<bb_dof_entries.size()<<" entries that meet restrictions."<<std::endl;
		return bb_dof_entries;
	}

	// for( UM::const_iterator it = entry_pairs_quick_access_C_N_.begin(); it != entry_pairs_quick_access_C_N_.end(); ++it )
	//  TR<<it->first<<" "<<it->second<<std::endl;
	TR<<"c_term_entry: "<<c_term_entry()<<std::endl;
	TR<<"n_term_entry: "<<n_term_entry()<<std::endl;
	if ( c_term_entry() && !n_term_entry() ) {
		for ( UM::const_iterator itc = c_entries.first; itc != c_entries.second; ++itc ) {
			bb_dof_entries.push_back( itc->second );
		}
		TR<<"Found "<<bb_dof_entries.size()<<" entries that meet restrictions."<<std::endl;
		return bb_dof_entries;
	}

	if ( !n_term_entry() && !c_term_entry() ) {
		for ( UM::const_iterator it = entry_pairs_quick_access_N_C_.cbegin(); it != entry_pairs_quick_access_N_C_.cend(); ++it ) {
			bb_dof_entries.push_back( it->first );
		}
		sort( bb_dof_entries.begin(), bb_dof_entries.end() );
		vector1< Size >::iterator last = unique( bb_dof_entries.begin(), bb_dof_entries.end() );
		bb_dof_entries.erase( last, bb_dof_entries.end());
		TR<<"Found "<<bb_dof_entries.size()<<" entries that meet restrictions."<<std::endl;
		return bb_dof_entries;
	}

	for ( UM::const_iterator itn = n_entries.first; itn != n_entries.second; ++itn ) {
		for ( UM::const_iterator itc = c_entries.first; itc != c_entries.second; ++itc ) {
			if ( itn->second == itc->second ) {
				bb_dof_entries.push_back( itn->second );
				break;
			}
		}
	}
	TR<<"Found "<<bb_dof_entries.size()<<" entries that meet restrictions."<<std::endl;
	return( bb_dof_entries );
}


/// @brief return a subset of entries that match all of the selection criteria
vector1< Size >
FoldUnitUtils::entry_subset_slow() const{
	utility::vector1< Size > bb_dof_entries;
	bb_dof_entries.clear();
	for ( Size bb_dof_entry = 1; bb_dof_entry <= bbdofs().size(); ++bb_dof_entry ) {
		if ( legal_bbdofs_[ bb_dof_entry ] &&
				bbdofs()[ bb_dof_entry ].size() <= max_length() ) {
			vector1< pair< Size, Size > >::const_iterator it;
			Size entry_pairs_idx( 0 );
			if ( n_term_entry() ) {
				if ( !fragment_compatibility_check( n_term_entry(), bb_dof_entry, it ) ) {
					continue;
				}
				entry_pairs_idx = it - entry_pairs_.begin() + 1;
				if ( overlap_length_[ entry_pairs_idx ] < min_overlap() || overlap_rmsd_[ entry_pairs_idx ] > max_rmsd() ) {
					continue;
				}
			}
			if ( c_term_entry() ) {
				if ( !fragment_compatibility_check( bb_dof_entry, c_term_entry(), it ) ) {
					continue;
				}
				entry_pairs_idx = it - entry_pairs_.begin() + 1;
				if ( overlap_length_[ entry_pairs_idx ] < min_overlap() || overlap_rmsd_[ entry_pairs_idx ] > max_rmsd() ) {
					continue;
				}
			}
			bb_dof_entries.push_back( bb_dof_entry );
		}
	}
	TR<<"Found "<<bb_dof_entries.size()<<" entries that meet restrictions."<<std::endl;
	return( bb_dof_entries );
}

bool
FoldUnitUtils::read_dbase(){
	utility::io::izstream data( fragment_dbase_ );
	runtime_assert( data );
	string line;
	while ( getline(data, line) ) {
		if ( line.length() == 0 ) {
			TR << "Fold bb database file empty or corrupted. Not loading." << std::endl;
			return false;
		}
		istringstream line_stream(line);
		string pdb, dssp, sequence;
		Size from_res, to_res, entry/*dummy variable needed for external scripts but not in Rosetta since all of the entries are in vectors*/;
		line_stream >> entry >> pdb >> from_res >> to_res >> dssp >> sequence;
		//TR<<"entry: "<<entry<<" pdb: "<<pdb<<" from_res: "<<from_res<<" to_res: "<<to_res<<" dssp: "<<dssp<< " sequence: "<<sequence<<std::endl;
		ResidueBBDofs bbdof;
		if ( dssp.length() != sequence.length() ) {
			//    TR<<"Illegal entry"<<std::endl;
			legal_bbdofs_.push_back( false );
		} else {
			legal_bbdofs_.push_back( true );
			bbdof.start_loop( from_res );
			bbdof.stop_loop( to_res );
			bbdof.source_pdb( pdb );
			bbdof.dssp( dssp );
			bbdof.aa_sequence( sequence );

			Size count = 0;
			while ( !line_stream.eof() ) {
				Real phi, psi, omega;
				line_stream >> phi >> psi >> omega;
				bbdof.push_back( BBDofs( count, phi, psi, omega, sequence.substr( count, 1 ) ) );
				++count;
			}
		}
		bbdofs_.push_back( bbdof );
		runtime_assert( bbdofs_.size() == entry );
	}
	data.close();
	TR<<"Read "<<bbdofs_.size()<<" torsion dbase entries"<<std::endl;
	Size count_illegal( 0 );
	TR<<"Illegal entries: ";
	for ( vector1< bool >::const_iterator b = legal_bbdofs_.begin(); b != legal_bbdofs_.end(); ++b ) {
		if ( !*b ) {
			TR<<b - legal_bbdofs_.begin() + 1<<", ";
			++count_illegal;
		}
	}
	TR<<"\ntotal of "<<count_illegal<<" illegal entries."<<std::endl;

	/// read the junctions database specifying which pairs can be combined and how to combine them
	utility::io::izstream pairs( pair_dbase_ );
	runtime_assert( pairs );
	Size count( 0 );
	while ( getline( pairs, line ) ) {
		istringstream line_stream(line);
		Size fragi, fragj, overlap_length;
		Real overlap_rmsd;
		line_stream >> fragi >> fragj >> overlap_length >> overlap_rmsd;
		// SJF 19Oct14: entry_pairs_ is no longer used since it's slow. It is more efficient with memory though, so if needed, the next line should be reactivated.
		// entry_pairs_.push_back( pair< Size, Size >( fragi, fragj ) );
		if ( overlap_length >= min_overlap() && overlap_rmsd <= max_rmsd() ) { // the quick_access multimaps precompute conditions
			entry_pairs_quick_access_N_C_.insert( pair< Size, Size >( fragi, fragj ));
			entry_pairs_quick_access_C_N_.insert( pair< Size, Size >( fragj, fragi ));
			count++;
		}
		overlap_length_.push_back( overlap_length );
		overlap_rmsd_.push_back( overlap_rmsd );
	}
	pairs.close();
	TR<<"Read "<<entry_pairs_quick_access_N_C_.size()<<" pairs."<<std::endl;
	TR<<"count: "<<count<<std::endl;
	/// remove from considerations entries that are not within the pairs list.
	// the following is unnecessary with the new quick_access multimaps 16Oct14
	//  for( Size i = 1; i <= bbdofs_.size(); ++i ){
	//   vector1< pair< Size, Size > >::const_iterator it;
	//   bool found = false;
	//   for( it = entry_pairs_.begin(); it != entry_pairs_.end() && found == false; ++it ){
	//    if( it->first == i || it->second == i )
	//     found = true;
	//   }
	//   if( !found )
	//    legal_bbdofs_[ i ] = false;
	//  }
	return true;
}

bool
FoldUnitUtils::fragment_compatibility_check( core::Size const i, core::Size const j, vector1< pair< Size, Size > >::const_iterator & it, Real const max_rmsd/*=10.0*/ ) const{
	if ( i == 0 || j == 0 ) {
		return true;
	}
	// the following shouldn't happen. Let's not test it
	// if( !legal_bbdofs_[ i ] || !legal_bbdofs_[ j ] )
	//  return false;
	it = find( entry_pairs_.begin(), entry_pairs_.end(), pair< Size, Size >( i, j ));
	if ( it == entry_pairs_.end() ) {
		return false;
	}
	if ( overlap_rmsd_[ it - entry_pairs_.begin() + 1 ] > max_rmsd ) {
		return false;
	}
	return true;
}

/// add_fragment_to_pose and replace_fragment_in_pose simply set the entry in the correct segment and call pose_from_fragment_info
void
FoldUnitUtils::add_fragment_to_pose( core::pose::Pose & pose, PoseFragmentInfo & fragment_info, core::Size const entry, bool const c_term ) const {
	Size const pose_fragments = fragment_info.size();

	if ( c_term ) {
		fragment_info.add_pair( pose_fragments + 1, entry );
	} else { //n_term: insert the new entry at the top and shift the others down by one
		PoseFragmentInfo new_frag_info;
		new_frag_info[ 1 ] = entry;
		for ( Size seg = 1; seg <= fragment_info.size(); ++seg ) {
			new_frag_info.add_pair( seg + 1, fragment_info[ seg ] );
		}
		fragment_info = new_frag_info;
	}
	pose_from_fragment_info( pose, fragment_info );
	fragment_info.set_fragment_info_in_pose( pose );
}

void
FoldUnitUtils::replace_fragment_in_pose( core::pose::Pose & pose, PoseFragmentInfo & fragment_info, core::Size const entry, core::Size const fragment_num ) {
	fragment_info[ fragment_num ] = entry;
	pose_from_fragment_info( pose, fragment_info );
	fragment_info.set_fragment_info_in_pose( pose );
}

/// @brief this is the place where the real 'logic' happens. Using the PoseFragmentInfo object the method constructs a pose from scratch, eliminating
// the contents of the incoming pose in the process.
void
FoldUnitUtils::pose_from_fragment_info( core::pose::Pose & pose, PoseFragmentInfo const & pose_fragment_info ) const{
	/// generate the sequence accounting for overlaps between adjoining segments
	string seq( "" ); // seq comprises the concatenated sequence from the adjoining segments
	for ( core::Size i = 1; i <= pose_fragment_info.size(); ++i ) {
		core::Size const dbase_entry = pose_fragment_info[ i ];
		ResidueBBDofs dofs = bbdofs_[ dbase_entry ];
		Size overlap( 0 );
		if ( i < pose_fragment_info.size() ) {
			Size const p_idx = pair_index( pose_fragment_info[ i ], pose_fragment_info[ i + 1 ] );
			overlap = overlap_length()[ p_idx ];
		}

		seq += dofs.aa_sequence().substr( 0, dofs.aa_sequence().length() - overlap );
	}
	// impose the dihedral angles
	using namespace core::chemical;
	ResidueTypeSet const & residue_set( pose.total_residue() ? pose.residue( 1 ).residue_type_set() : *ChemicalManager::get_instance()->residue_type_set( CENTROID ) ); // residuetypeset is noncopyable
	core::pose::make_pose_from_sequence( pose, seq, residue_set );
	/// set bb dofs
	Size start( 1 );
	for ( core::Size i = 1; i <= pose_fragment_info.size(); ++i ) {
		core::Size const dbase_entry = pose_fragment_info[ i ];
		if ( i > 1 ) {
			Size const p_idx = pair_index( pose_fragment_info[ i - 1 ], pose_fragment_info[ i ] );
			start -= overlap_length()[ p_idx ];
		}
		ResidueBBDofs dofs = bbdofs_[ dbase_entry ];
		for ( Size resi = start; resi <= start + dofs.size() - 1; ++resi ) {
			pose.set_phi(   resi, dofs[ resi - start + 1 ].phi() );
			pose.set_psi(   resi, dofs[ resi - start + 1 ].psi() );
			pose.set_omega( resi, dofs[ resi - start + 1 ].omega() );
		}
		start += dofs.size();
	}
}

std::string
AddFoldUnitMoverCreator::keyname() const
{
	return AddFoldUnitMoverCreator::mover_name();
}

protocols::moves::MoverOP
AddFoldUnitMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddFoldUnitMover );
}

std::string
AddFoldUnitMoverCreator::mover_name()
{
	return "AddFoldUnit";
}


AddFoldUnitMover::AddFoldUnitMover(): moves::Mover("AddFoldUnit"),
	fragment_dbase_( "" ),
	max_length_( 40 ),
	min_overlap_( 5 ),
	max_rmsd_( 10 ),
	terminus_( 'b' ),
	replace_fragment_( false ),
	replace_fragment_segment_number_( 0 ),
	fold_unit_utils_( NULL ),
	max_segments_( 5 )
{
}

AddFoldUnitMover::~AddFoldUnitMover()
{
}

void
AddFoldUnitMover::apply( Pose & pose ){
	PoseFragmentInfo pfi;
	pfi.load_fragment_info_from_pose( pose );
	TR<<"Pose size: "<<pose.total_residue()<<std::endl;
	Size last_segment_entry( 0 ), first_segment_entry( 0 );
	if ( !replace_fragment() && pfi.size() >= max_segments() ) {
		TR<<"Added maximum number of segments. Not adding any more."<<std::endl;
		return;
	}
	if ( pfi.size() ) {
		last_segment_entry = pfi[ pfi.size() ];
		first_segment_entry = pfi.size() == 0 ? 0 : pfi[ 1 ];
	}
	fold_unit_utils_->n_term_entry( 0 );
	fold_unit_utils_->c_term_entry( 0 );

	fold_unit_utils_->max_length( max_length() );
	fold_unit_utils_->min_overlap( min_overlap() );
	fold_unit_utils_->max_rmsd( max_rmsd() );
	core::Size segment_to_replace( 0 );
	bool const c_terminus = ( (terminus_ == 'b') ? ( (int) (numeric::random::rg().uniform() * 2.0) == 1 ) : (terminus_ == 'c') );
	TR<<"c_terminus? "<<c_terminus<<std::endl;
	if ( replace_fragment() ) {
		segment_to_replace = replace_fragment_segment_number();
		if ( segment_to_replace == 0 ) {
			segment_to_replace = (int) (numeric::random::rg().uniform() * (double) pfi.size() ) + 1;
			TR<<"Out of "<<pfi.size()<<" segments, replacing segment "<<segment_to_replace<<std::endl;
		}
		fold_unit_utils_->n_term_entry( 0 );
		fold_unit_utils_->c_term_entry( 0 );

		if ( segment_to_replace > 1 ) {
			fold_unit_utils_->n_term_entry( pfi[ segment_to_replace - 1 ] );
		}
		if ( segment_to_replace < pfi.size() ) {
			fold_unit_utils_->c_term_entry( pfi[ segment_to_replace + 1 ] );
		}
	} else {
		fold_unit_utils_->n_term_entry( 0 );
		fold_unit_utils_->c_term_entry( 0 );

		if ( c_terminus ) {
			fold_unit_utils_->n_term_entry( last_segment_entry );
		} else {
			fold_unit_utils_->c_term_entry( first_segment_entry );
		}
	}
	utility::vector1< Size > const entries( fold_unit_utils_->entry_subset() );
	if ( entries.size() == 0 ) {
		TR<<"Found no entry that matches all selection criteria; failing."<<std::endl;
		return;
	}
	Size const rand_entry = entries[ (Size) ( numeric::random::rg().uniform() * (double) entries.size()) + 1 ];

	if ( replace_fragment() ) {
		fold_unit_utils_->replace_fragment_in_pose( pose, pfi, rand_entry, segment_to_replace );
	} else {
		fold_unit_utils_->add_fragment_to_pose( pose, pfi, rand_entry, c_terminus );
	}
}

std::string
AddFoldUnitMover::get_name() const {
	return AddFoldUnitMoverCreator::mover_name();
}

moves::MoverOP
AddFoldUnitMover::clone() const
{
	return moves::MoverOP( new AddFoldUnitMover( *this ) );
}

moves::MoverOP
AddFoldUnitMover::fresh_instance() const
{
	return moves::MoverOP( new AddFoldUnitMover );
}

void
AddFoldUnitMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	fragment_dbase( tag->getOption< string >( "fragment_dbase", "in/fold_from_frags.db" ));
	pair_dbase( tag->getOption< string >( "pair_dbase", "in/junctions.db" ));
	max_length( tag->getOption< Size >( "max_length", 40 ));
	max_segments( tag->getOption< Size > ("max_segments", 5 ) );
	min_overlap( tag->getOption< Size >( "min_overlap", 5 ) );
	max_rmsd( tag->getOption< Real > ( "max_rmsd", 10.0 ) );
	if ( !fold_unit_utils_ ) {
		fold_unit_utils_ = FoldUnitUtilsOP( new FoldUnitUtils );
	}

	fold_unit_utils_->max_length( max_length() );
	fold_unit_utils_->min_overlap( min_overlap() );
	fold_unit_utils_->max_rmsd( max_rmsd() );
	terminus( tag->getOption< char >( "terminus", 'b' ) );
	fold_unit_utils_->fragment_dbase( fragment_dbase() );
	fold_unit_utils_->pair_dbase( pair_dbase() );
	bool const success = fold_unit_utils_->read_dbase();
	TR<<"fold_unit_utils dbase read status: "<<success<<std::endl;
	runtime_assert( success );
	replace_fragment( tag->getOption< bool >( "replace_fragment", false ) );
	replace_fragment_segment_number( tag->getOption< Size >( "replace_fragment_segment_number", 0 ) );
	TR<<"fragment_dbase: "<<fragment_dbase()<<" max_length: "<<max_length()<<" max_segments: "<<max_segments()<<" terminus: "<<terminus_<<" min_overlap: "<<min_overlap()<<" max_rmsd: "<<max_rmsd()<<" replace_fragment: "<<replace_fragment()<<" replace_fragment_segment_number: "<<replace_fragment_segment_number()<<std::endl;
}

void
PoseFragmentInfo::load_fragment_info_from_pose( core::pose::Pose const & pose ){
	map<std::string, std::string> comments = core::pose::get_all_comments(pose);
	Size count = 0;
	fragment_map_.clear();
	for ( map< string, string >::const_iterator com = comments.begin(); com != comments.end(); ++com ) {
		if ( com->first.substr(0, 19 ) != "fragment_definition" ) {
			continue;
		}
		++count;
		fragment_map_[ count ] = utility::from_string< Size >(com->second, Size() );
	}
	TR<<"loaded "<<count<<" fragment definitions from pose"<<std::endl;
}

/// @brief first delete then write from scratch the fragment_definition comments to the pose
void
PoseFragmentInfo::set_fragment_info_in_pose( core::pose::Pose & pose )const{

	map<std::string, std::string> comments = core::pose::get_all_comments(pose);
	for ( map< string, string >::const_iterator com = comments.begin(); com != comments.end(); ++com ) {
		if ( com->first.substr(0, 19 ) == "fragment_definition" ) {
			core::pose::delete_comment( pose, com->first );
		}
	}

	Size count( 1 );
	for ( map< Size, Size >::const_iterator frag = fragment_map_.begin(); frag != fragment_map_.end(); ++frag ) {
		core::pose::add_comment( pose, "fragment_definition" + utility::to_string(count), utility::to_string( frag->second ) ); // fragment_definition# 1_14
		++count;
	}
	TR<<"Added "<<fragment_map_.size()<<" fragment comments to pose"<<std::endl;
}

/// @brief find the index within the entry_pairs_ array of pair i,j
Size
FoldUnitUtils::pair_index( Size const i, Size const j ) const{
	if ( i == 0 || j == 0 ) {
		return 0;
	}
	vector1< pair< Size, Size > >::const_iterator it = find( entry_pairs_.begin(), entry_pairs_.end(), pair< Size, Size >( i, j ) );
	return it - entry_pairs_.begin() + 1;
}


void
PoseFragmentInfo::fragment_start_end( FoldUnitUtils const fuu, core::Size const fragment, core::Size & start, core::Size & end ){
	runtime_assert( fragment_map_.size() >= fragment );
	start = 0; end = 0;
	Size prev_frag( 0 );
	Size overlap = 0;
	for ( Size frag = 1; frag <= fragment; ++frag ) {
		start = end + 1;
		Size const curr_frag = fragment_map_[ frag ];
		if ( prev_frag ) {
			overlap += fuu.overlap_length()[ fuu.pair_index( prev_frag, curr_frag ) ];
		}

		if ( fragment > 0 ) {
			end += fuu.bbdofs()[ curr_frag ].size();
		}
		prev_frag = curr_frag;
		TR<<"frag: "<<frag<<" start: "<<start<<" end: "<<end<<" cumulative overlap: "<<overlap<<std::endl;
	}
	start -= overlap; end -= overlap;
	TR<<"Fragment "<<fragment<<" start: "<<start<<" end: "<<end<<std::endl;
}

core::Size &
PoseFragmentInfo::operator[]( core::Size const s ) {
	runtime_assert( s );
	return fragment_map_[ s ];
}

core::Size
PoseFragmentInfo::operator[] ( core::Size const s ) const {
	if ( fragment_map_.size() < s ) {
		return (Size) 0;
	}
	return fragment_map_.at( s );
}

FoldUnitUtils::~FoldUnitUtils() {}

PoseFragmentInfo::~PoseFragmentInfo() {}

StartFreshMover::StartFreshMover() : Mover( "StartFresh" ), residue_type_set_( "centroid" ) {}
StartFreshMover::~StartFreshMover() {}

void
StartFreshMover::apply( core::pose::Pose & pose ){
	core::pose::Pose new_pose;

	core::util::switch_to_residue_type_set( new_pose, residue_type_set_ );
	pose = new_pose;
}

void
StartFreshMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	residue_type_set_ = tag->getOption< string >( "residue_type_set", "centroid" );
	TR<<"residue_type_set: "<<residue_type_set()<<std::endl;
}

moves::MoverOP
StartFreshMover::clone() const
{
	return moves::MoverOP( new StartFreshMover( *this ) );
}

moves::MoverOP
StartFreshMover::fresh_instance() const
{
	return moves::MoverOP( new StartFreshMover );
}

std::string
StartFreshMover::get_name() const {
	return StartFreshMoverCreator::mover_name();
}

std::string
StartFreshMoverCreator::keyname() const
{
	return StartFreshMoverCreator::mover_name();
}

protocols::moves::MoverOP
StartFreshMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StartFreshMover );
}

std::string
StartFreshMoverCreator::mover_name()
{
	return "StartFresh";
}
} // simple_moves
} // protocols

