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
			Size fragi, fragj;
			line_stream >> fragi >> fragj;
			entry_pairs_.push_back( pair< Size, Size >( fragi, fragj ) );
		}
		pairs.close();
		return true;
}

void
FoldUnitUtils::add_fragment_to_pose( core::pose::Pose & pose, core::Size const entry, bool const c_term ) const {
	using namespace core::chemical;
	using namespace core::conformation;

	ResidueBBDofs dofs = bbdofs_[ entry ];

  ResidueTypeSet const & residue_set( pose.residue( 1 ).residue_type_set() ); // residuetypeset is noncopyable
	Size const start_pose_size = ( c_term ? pose.total_residue() : 1 );
	for( Size resi = 1; resi <= dofs.size(); ++resi ){
		ResidueCOP new_res = ResidueFactory::create_residue( residue_set.name_map( name_from_aa( aa_from_oneletter_code( dofs.aa_sequence()[ resi ] ) ) ) );
		if( c_term )
  		pose.conformation().append_polymer_residue_after_seqpos( *new_res, pose.total_residue(), true/*build_ideal*/);
		else
  		pose.conformation().prepend_polymer_residue_before_seqpos( *new_res, 1, true/*build_ideal*/);
	}
/// separate sequence length change from dihedral changes, b/c omega is not defined for the last residue
	for( int resi = start_pose_size; resi <= (int) ( start_pose_size + dofs.size() ); ++resi ){
		pose.set_phi( resi, dofs[ resi ].phi() );
		pose.set_psi( resi, dofs[ resi ].psi() );
		pose.set_omega( resi, dofs[ resi ].omega() );
	}
}

void
FoldUnitUtils::replace_fragment_in_pose( core::pose::Pose & pose, core::Size const entry, core::Size const from_res, core::Size const to_res ) const {
	using namespace core::chemical;
	using namespace core::conformation;

  ResidueTypeSet const & residue_set( pose.residue( 1 ).residue_type_set() ); // residuetypeset is noncopyable
	ResidueBBDofs dofs = bbdofs_[ entry ];
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
	for( Size resi = from_res; resi <= to_res + residue_diff; ++resi ){
		ResidueCOP new_res = ResidueFactory::create_residue( residue_set.name_map( name_from_aa( aa_from_oneletter_code( dofs.aa_sequence()[ resi - from_res + 1 ] ) ) ) );
		pose.replace_residue( resi, *new_res, false/*orient bb*/ );
		pose.set_phi( resi, dofs[ resi ].phi() );
		pose.set_psi( resi, dofs[ resi ].psi() );
		pose.set_omega( resi, dofs[ resi ].omega() );
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

} // simple_moves
} // protocols

