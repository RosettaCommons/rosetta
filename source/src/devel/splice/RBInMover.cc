// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RBInMover.cc
/// @brief

// Unit headers
#include <devel/splice/RBInMover.hh>
#include <devel/splice/RBOutMover.hh>
#include <devel/splice/RBInMoverCreator.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("devel.splice.RBInMover");
#include <utility/tag/Tag.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <string>
//#include <sstream>
#include <core/pose/util.hh>
#include <algorithm>
#include <utility/io/izstream.hh>
#include <protocols/protein_interface_design/movers/SetAtomTree.hh>
#include <boost/foreach.hpp>

namespace devel {
namespace splice {

using namespace::protocols;

std::string
RBInMoverCreator::keyname() const
{
	return RBInMoverCreator::mover_name();
}

protocols::moves::MoverOP
RBInMoverCreator::create_mover() const {
	return new RBInMover;
}

std::string
RBInMoverCreator::mover_name()
{
	return "RBIn";
}

RBInMover::RBInMover(): moves::Mover("RBIn"),
    RB_dbase_(""),
    from_entry_( 1 ),
    to_entry_( 1 ),
    current_entry_( 1 ),
    randomize_( true ),
		checkpointing_file_( "" ),
		modify_foldtree_( true )
{
	jump_library_.clear();
}

RBInMover::~RBInMover()
{
}

void
RBInMover::jump_library( utility::vector1< core::kinematics::Jump > j ){ jump_library_ = j; }

utility::vector1< core::kinematics::Jump >
RBInMover::jump_library() const{ return jump_library_; }

bool
RBInMover::checkpoint() const{
	if( checkpointing_file_ == "" )
		return false;
//write to file
	std::ofstream data;
	data.open( checkpointing_file_.c_str(), std::ios::out );
	if( !data.good() )
		utility_exit_with_message( "Unable to open checkpointing file for writing: " + checkpointing_file() + "\n" );
	data<<current_entry_<<'\n';
	BOOST_FOREACH( core::kinematics::Jump const jump, jump_library() )
		data<<jump<<'\n';
	data.flush();
	TR<<"Finished writing checkpointing information to "<<checkpointing_file()<<std::endl;
	return true;
}

bool
RBInMover::checkpoint_recovery(){
    if( checkpointing_file_ == "" )
        return false;

    utility::io::izstream data( checkpointing_file_ );
    if( !data ){
        TR<<"Checkpointing file not yet written. Not recovering."<<std::endl;
        return false;
    }
    TR<<"Checkpointing file found. Recovering."<<std::endl;
    runtime_assert( jump_library_.size() == 0 ); /// if it's not 0, someone has already initialized the library and that doesn't make too much sense...
    std::string line;
    getline( data, line );
    std::istringstream data_stream( line );
    data_stream >> current_entry_;
    TR<<"Recovered entry: "<<current_entry_<<std::endl;
    while( getline( data, line ) ) {
        TR<<"pushing jump into library: "<<line<<std::endl;
        core::kinematics::Jump jump;
        std::istringstream jump_stream( line );
        jump_stream >> jump;
        jump_library_.push_back( jump );
    }
    TR<<"Read "<<jump_library_.size()<<" jumps into library"<<std::endl;
    return true;
    }

void
RBInMover::init(){
    bool const checkpoint_recover( checkpoint_recovery() );
    if( checkpoint_recover )
        return;
    bool const first_pass( jump_library_.size() == 0 );
    if( !first_pass )// already initialized; don't repeat
        return;

    TR<<"Initializing rigid-body jump library"<<std::endl;
    utility::io::izstream data( RB_dbase_ );
    if ( !data ) {
        TR << "cannot open rigid-body database " << RB_dbase_ << std::endl;
        utility_exit();
    }
    std::string line;
    core::Size count = 0;
    while( getline( data, line ) ) {
        TR<<"pushing jump into library: "<<line<<std::endl;
        core::kinematics::Jump jump;
        std::istringstream data_stream( line );
        data_stream >> jump;
        count++;
        if( count <= to_entry() && count >= from_entry() )
            jump_library_.push_back( jump );
    }
    if( randomize() ){
        TR<<"Randomizing the order in the jump_library"<<std::endl;
        std::random_shuffle( jump_library_.begin(), jump_library_.end() );
    }
    current_entry_ = 1;/// this should already be initialized, but just in case
    checkpoint(); /// first pass should also checkpoint to keep the order of entries
    TR<<"Done making jump library with "<<jump_library_.size()<<std::endl;
}

void
RBInMover::set_fold_tree( Pose & pose ) const{
    utility::vector1< std::pair< core::Size, core::Size > > const disulfs( find_disulfs_in_range( pose, 1, pose.conformation().chain_end( 1 ) ) );
    protocols::protein_interface_design::movers::SetAtomTree sat;
    sat.two_parts_chain1( true );
    sat.apply( pose );
    core::kinematics::FoldTree ft( pose.fold_tree() );
    TR<<"Fold tree before slide jump: "<<ft<<std::endl;
    TR<<"Setting slide jump from: "<<disulfs[ 1 ].second<<" to "<<disulfs[ 2 ].second<<std::endl;
    ft.slide_jump( 1, disulfs[ 1 ].second, disulfs[ 2 ].second );
    ft.reorder(1);
    TR<<"Fold tree after slide jump: "<<ft<<std::endl;
    pose.fold_tree( ft );
}


void
RBInMover::apply( Pose & pose ){
    init();
		if( modify_foldtree() )
			set_fold_tree( pose );
    TR<<"Setting jump now"<<std::endl;
    pose.set_jump( 1, jump_library_[ current_entry_ ] );

    // Update RBO comment
    //std::ostream RBMoveStream;
    //RBMoveStream<<jump_library_[ current_entry_ ];
    //std::string RBMoveString = RBMoveStream.str();
    //std::replace( RBMoveString.begin(), RBMoveString.end(), ' ', ',');
    //TR<<"Current RBO jump: "<<RBMoveString<<std::endl;
    //core::pose::add_comment(pose,"RBO ",RBMoveString);

    if( current_entry_ == jump_library_.size() ) // rewind to the beginning of the library
			current_entry_ = 1;
		else
			++current_entry_; // go up one entry
	checkpoint();
}

std::string
RBInMover::get_name() const {
	return RBInMoverCreator::mover_name();
}

moves::MoverOP
RBInMover::clone() const
{
	return new RBInMover( *this );
}

moves::MoverOP
RBInMover::fresh_instance() const
{
	return new RBInMover;
}

void
RBInMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
  RB_dbase( tag->getOption< std::string >( "rigid_body_dbase", "" ) );
  from_entry( tag->getOption< core::Size >( "from_entry", 1 ) );
  to_entry( tag->getOption< core::Size >( "to_entry", 1 ) );
  randomize( tag->getOption< bool >( "randomize", true  ) );
	checkpointing_file( tag->getOption< std::string >( "checkpointing_file", "" ) );
	modify_foldtree( tag->getOption< bool >( "modify_foldtree", true ) );

  runtime_assert( from_entry() <= to_entry() && from_entry() >= 1 );
  TR<<"RB_dbase: "<<RB_dbase()<<" from_entry: "<<from_entry()<<" to_entry: "<<to_entry()<<" randomize: "<<randomize()<<" checkpointing file: "<<checkpointing_file()<<" modify_foldtree: "<<modify_foldtree()<<std::endl;
}
} // simple_moves
} // protocols

