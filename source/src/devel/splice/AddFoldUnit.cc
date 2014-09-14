// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
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
static basic::Tracer TR("devel.splice.AddFoldUnitMover");
#include <utility/tag/Tag.hh>
#include <core/util/SwitchResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <string>
#include <utility/string_util.hh>
//#include <sstream>
#include <core/pose/util.hh>
#include <utility/io/izstream.hh>
#include <boost/foreach.hpp>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>

namespace devel {
namespace splice {

using namespace protocols;
using namespace core;
using namespace std;
using utility::vector1;

bool
FoldUnitUtils::read_dbase(){
		utility::io::izstream data( fragment_dbase_ );
		runtime_assert( data );
		string line;
		while( getline(data, line) ){
		  if (line.length() == 0) {
		    TR << "Fold bb database file empty or corrupted. Not loading." << std::endl;
		    return false;
		  }
		  istringstream line_stream(line);
			string pdb, dssp, sequence;
			Size from_res, to_res;
			line_stream >> pdb >> from_res >> to_res >> dssp >> sequence;
			TR<<"pdb: "<<pdb<<" from_res: "<<from_res<<" to_res: "<<to_res<<" dssp: "<<dssp<< " sequence: "<<sequence<<std::endl;
			ResidueBBDofs bbdof;
			bbdof.start_loop( from_res );
			bbdof.stop_loop( to_res );
			bbdof.source_pdb( pdb );
			bbdof.dssp( dssp );
			bbdof.aa_sequence( sequence );

			Size count = 0;
			while( !line_stream.eof() ){
				Real phi, psi, omega;
				line_stream >> phi >> psi >> omega;
				bbdof.push_back( BBDofs( count, phi, psi, omega, sequence.substr( count, 1 ) ) );
				++count;
			}
			bbdofs_.push_back( bbdof );
			TR<<"read "<<count<<" dofs"<<std::endl;
		}
		data.close();

		utility::io::izstream pairs( fragment_dbase_ );
		runtime_assert( pairs );
		while( getline( pairs, line )){
		  istringstream line_stream(line);
			Size fragi, fragj, overlap_length, overlap_rmsd;
			line_stream >> fragi >> fragj >> overlap_length >> overlap_rmsd;
			entry_pairs_.push_back( pair< Size, Size >( fragi, fragj ) );
			overlap_length_.push_back( overlap_length );
			overlap_rmsd_.push_back( overlap_rmsd );
		}
		pairs.close();
		return true;
}

bool
FoldUnitUtils::fragment_compatibility_check( core::Size const i, core::Size const j, vector1< pair< Size, Size > >::const_iterator it, Real const max_rmsd/*=10.0*/ ) const{
	it = find( entry_pairs_.begin(), entry_pairs_.end(), pair< Size, Size >( i, j ));
	if( it == entry_pairs_.end() )
		return false;
	if( overlap_rmsd_[ it - entry_pairs_.begin() + 1 ] > max_rmsd )
		return false;
	return true;
}

/// @add fragment to pose, taking care of overlaps between subsequence fragments
void
FoldUnitUtils::add_fragment_to_pose( core::pose::Pose & pose, PoseFragmentInfo fragment_info, core::Size const entry, bool const c_term ) const {
	ResidueBBDofs dofs = bbdofs_[ entry ];
	Size const pose_fragments = fragment_info.size();
	Size overlap( 0 );
	if( pose_fragments ){ //check compatibility between fragments
		Size const last_fragment = fragment_info[ pose_fragments ];
		vector1< pair< Size, Size > >::const_iterator it;
		fragment_compatibility_check( last_fragment, entry, it);
		if( it == entry_pairs_.end() )
			TR<<"requested to add fragment "<<entry<<" after fragment "<<last_fragment<<" but they are incompatible!"<<std::endl;
		else{
			overlap = overlap_length_[ it - entry_pairs_.begin() + 1 ];
			TR<<"requested to add fragment "<<entry<<" after fragment "<<last_fragment<<" which are compatible. overlap is "<<overlap<<std::endl;
		}
	}

	using namespace core::chemical;
	using namespace core::conformation;

  ResidueTypeSet const & residue_set( pose.total_residue() ? pose.residue( 1 ).residue_type_set() : *ChemicalManager::get_instance()->residue_type_set( CENTROID ) ); // residuetypeset is noncopyable
	Size const start_pose_size = ( c_term ? pose.total_residue() : 0 );
	for( Size resi = 1 + overlap; resi <= dofs.size(); ++resi ){
		ResidueCOP new_res = ResidueFactory::create_residue( residue_set.name_map( name_from_aa( aa_from_oneletter_code( dofs.aa_sequence()[ resi ] ) ) ) );
		if( c_term )
  		pose.conformation().append_polymer_residue_after_seqpos( *new_res, pose.total_residue(), true/*build_ideal*/);
		else
  		pose.conformation().prepend_polymer_residue_before_seqpos( *new_res, 1, true/*build_ideal*/);
	}
/// separate sequence length change from dihedral changes, b/c omega is not defined for the last residue
	for( Size resi = 1; resi <= dofs.size() - overlap; ++resi ){
		pose.set_phi(   start_pose_size + resi, dofs[ start_pose_size + resi + overlap ].phi() );
		pose.set_psi(   start_pose_size + resi, dofs[ start_pose_size + resi + overlap ].psi() );
		pose.set_omega( start_pose_size + resi, dofs[ start_pose_size + resi + overlap ].omega() );
	}
}

void
FoldUnitUtils::replace_fragment_in_pose( core::pose::Pose & pose, PoseFragmentInfo fragment_info, core::Size const entry, core::Size const fragment_num ) {
	using namespace core::chemical;
	using namespace core::conformation;

  ResidueTypeSet const & residue_set( pose.residue( 1 ).residue_type_set() ); // residuetypeset is noncopyable
	ResidueBBDofs dofs = bbdofs_[ entry ];
	Size from_res, to_res;
	fragment_info.fragment_start_end( *this, fragment_num, from_res, to_res );
	int const residue_diff = dofs.size() - ( to_res - from_res + 1 );
	if( residue_diff < 0 ){
		for( int i = 1; i >= residue_diff; --i )
			pose.delete_polymer_residue( from_res );
	}
	else if( residue_diff > 0 ){
		ResidueCOP new_res = ResidueFactory::create_residue( residue_set.name_map( name_from_aa( aa_from_oneletter_code( 'A' ) ) ) );
		for( Size i = 1; i <= ( Size ) residue_diff; ++i )
			pose.append_polymer_residue_after_seqpos( *new_res, to_res, true );
	}

	Size overlap_n( 0 ), overlap_c( 0 );
	Size const prev_fragment = fragment_info[ fragment_num - 1 ];
	Size const next_fragment = fragment_info[ fragment_num + 1 ];
	vector1< pair< Size, Size > >::const_iterator it_n, it_c;
	fragment_compatibility_check( prev_fragment, entry, it_n );
	if( it_n == entry_pairs_.end() )
		TR<<"requested to combine fragments "<<prev_fragment<<" and "<<entry<<" but they're incompatible!"<<std::endl;
	else{
		overlap_n = overlap_length_[ it_n - entry_pairs_.begin() + 1 ];
		TR<<"requested to combine fragments "<< prev_fragment<<" and "<<entry<<" with overlap "<<overlap_n<<std::endl;
	}
	fragment_compatibility_check( entry, next_fragment, it_c );
	if( it_c == entry_pairs_.end() )
		TR<<"requested to combine fragments "<<entry<<" and "<<next_fragment<<" but they're incompatible!"<<std::endl;
	else{
		overlap_c = overlap_length_[ it_c - entry_pairs_.begin() + 1 ];
		TR<<"requested to combine fragments "<< entry<<" and "<<next_fragment<<" with overlap "<<overlap_c<<std::endl;
	}

	for( Size resi = from_res + overlap_n; resi <= to_res + residue_diff - overlap_c; ++resi ){
		Size const resi_in_dofs_array( resi - from_res + 1 );
		ResidueCOP new_res = ResidueFactory::create_residue( residue_set.name_map( name_from_aa( aa_from_oneletter_code( dofs.aa_sequence()[ resi_in_dofs_array ] ) ) ) );
		pose.replace_residue( resi, *new_res, false/*orient bb*/ );
		pose.set_phi(   resi, dofs[ resi_in_dofs_array ].phi() );
		pose.set_psi(   resi, dofs[ resi_in_dofs_array ].psi() );
		pose.set_omega( resi, dofs[ resi_in_dofs_array ].omega() );
	}
}

std::string
AddFoldUnitMoverCreator::keyname() const
{
	return AddFoldUnitMoverCreator::mover_name();
}

protocols::moves::MoverOP
AddFoldUnitMoverCreator::create_mover() const {
	return new AddFoldUnitMover;
}

std::string
AddFoldUnitMoverCreator::mover_name()
{
	return "AddFoldUnit";
}


AddFoldUnitMover::AddFoldUnitMover(): moves::Mover("AddFoldUnit"),
	fragment_dbase_( "" ),
	max_length_( 40 ),
	add_c_( true )
{
}

AddFoldUnitMover::~AddFoldUnitMover()
{
}

void
AddFoldUnitMover::apply( Pose & ){
	utility::vector1< Size > entries;
	for( Size entry = 1; entry <= fold_unit_utils_->bbdofs().size(); ++entry ){
		if( fold_unit_utils_->bbdofs()[ entry ].size() <= max_length() )
			entries.push_back( entry );
	}
	TR<<"Found "<<entries.size()<<" entries that max restrictions."<<std::endl;
}

std::string
AddFoldUnitMover::get_name() const {
	return AddFoldUnitMoverCreator::mover_name();
}

moves::MoverOP
AddFoldUnitMover::clone() const
{
	return new AddFoldUnitMover( *this );
}

moves::MoverOP
AddFoldUnitMover::fresh_instance() const
{
	return new AddFoldUnitMover;
}

void
AddFoldUnitMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	fragment_dbase( tag->getOption< string >( "fragment_dbase" ));
	max_length( tag->getOption< Size >( "max_length", 40 ));
	add_c( tag->getOption< bool >( "add_c", true ) );
	if( !fold_unit_utils_ ){
		fold_unit_utils_ = new FoldUnitUtils;
		fold_unit_utils_->fragment_dbase( fragment_dbase() );
		bool const success = fold_unit_utils_->read_dbase();
		TR<<"fold_unit_utils dbase read status: "<<success<<std::endl;
		runtime_assert( success );
	}
	TR<<"fragment_dbase: "<<fragment_dbase()<<" max_length: "<<max_length()<<" add_c: "<<add_c()<<std::endl;
}

void
PoseFragmentInfo::load_fragment_info_from_pose( core::pose::Pose const & pose ){
	 map<std::string, std::string> comments = core::pose::get_all_comments(pose);
	 Size count = 0;
	 for( map< string, string >::const_iterator com = comments.begin(); com != comments.end(); ++com ){
		 if( com->first.substr(0, 19 ) != "fragment_definition" )
			 continue;
		 utility::vector1< Size > first_second = utility::string_split( com->second, '_', Size() );
		 fragment_map_[ first_second[ 1 ] ] = first_second[ 2 ];
		 ++count;
	 }
	 TR<<"loaded "<<count<<" fragment definitions from pose"<<std::endl;
}

void
PoseFragmentInfo::set_fragment_info_in_pose( core::pose::Pose & pose )const{
	Size count( 1 );
	for( map< Size, Size >::const_iterator frag = fragment_map_.begin(); frag != fragment_map_.end(); ++frag ){
		core::pose::add_comment( pose, "fragment_definition" + utility::to_string(count), utility::to_string( frag->first ) + "_" + utility::to_string( frag->second ) ); // fragment_definition# 1_14
		++count;
	}
	TR<<"Added "<<fragment_map_.size()<<" fragment comments to pose"<<std::endl;
}

void
PoseFragmentInfo::fragment_start_end( FoldUnitUtils & fuu, core::Size const fragment, core::Size & start, core::Size & end ){
	runtime_assert( fragment_map_.size() >= fragment );
	start = 0; end = 0;
	for( Size frag = 1; frag <= fragment; ++frag ){
		start = end + 1;
		end += fuu.bbdofs()[ fragment_map_[ frag ] ].size();
	}
	TR<<"Fragment "<<fragment<<" start: "<<start<<" end: "<<end<<std::endl;
}

core::Size
PoseFragmentInfo::operator[]( core::Size const s ) { return fragment_map_[ s ]; }

FoldUnitUtils::~FoldUnitUtils() {}

PoseFragmentInfo::~PoseFragmentInfo() {}

StartFreshMover::StartFreshMover() : Mover( "StartFresh" ), residue_type_set_( "CENTROID" ) {}
StartFreshMover::~StartFreshMover() {}

void
StartFreshMover::apply( core::pose::Pose & pose ){
	core::pose::Pose new_pose;

	core::util::switch_to_residue_type_set( pose, residue_type_set_ );
	pose = new_pose;
}

void
StartFreshMover::parse_my_tag(
    utility::tag::TagCOP tag,
    basic::datacache::DataMap &,
    protocols::filters::Filters_map const &,
    protocols::moves::Movers_map const &,
    core::pose::Pose const & ){
	residue_type_set_ = tag->getOption< string >( "residue_type_set", "CENTROID" );
	TR<<"residue_type_set: "<<residue_type_set()<<std::endl;
}

moves::MoverOP
StartFreshMover::clone() const
{
	return new StartFreshMover( *this );
}

moves::MoverOP
StartFreshMover::fresh_instance() const
{
	return new StartFreshMover;
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
	return new StartFreshMover;
}

std::string
StartFreshMoverCreator::mover_name()
{
	return "StartFresh";
}
} // simple_moves
} // protocols

