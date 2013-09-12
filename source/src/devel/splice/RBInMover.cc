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
#include <protocols/moves/DataMapObj.hh>
#include <protocols/moves/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.simple_moves.RBInMover");
#include <utility/tag/Tag.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <fstream>
#include <algorithm>
#include <utility/io/izstream.hh>
#include <protocols/protein_interface_design/movers/SetAtomTree.hh>

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
    randomize_( true )
{
}

RBInMover::~RBInMover()
{
}

void
RBInMover::jump_library( utility::vector1< core::kinematics::Jump > j ){ jump_library_ = j; }

utility::vector1< core::kinematics::Jump >
RBInMover::jump_library() const{ return jump_library_; }
 
void
RBInMover::init(){
   if( jump_library_.size() > 0 )
	return;
  utility::io::izstream data( RB_dbase_ );
  if ( !data ) {
    TR << "cannot open rigid-body database " << RB_dbase_ << std::endl;
    utility_exit();
  }
  std::string line;
  while( getline( data, line ) ) {
    TR<<"pushing jump into library: "<<line<<std::endl;
    core::kinematics::Jump jump;
    std::istringstream data_stream( line );
    data_stream >> jump;
    jump_library_.push_back( jump );
  }
  if( randomize() ){
	TR<<"Randomizing the order in the jump_library"<<std::endl;
	std::random_shuffle( jump_library_.begin(), jump_library_.end() );
  }
  TR<<"Done making jump library with "<<jump_library_.size()<<std::endl;
}

void
RBInMover::apply( Pose & pose )
{
    init();
    utility::vector1< std::pair< core::Size, core::Size > > const disulfs( find_disulfs_in_range( pose, 1, pose.conformation().chain_end( 1 ) ) );
    protocols::protein_interface_design::movers::SetAtomTree sat;
    sat.two_parts_chain1( true );
    sat.apply( pose );  
    core::kinematics::FoldTree ft( pose.fold_tree() );
    TR<<"Fold tree before slide jump: "<<ft<<std::endl;
    ft.slide_jump( 1, disulfs[ 1 ].second, disulfs[ 2 ].second );
    TR<<"Fold tree after slide jump: "<<ft<<std::endl;
    pose.fold_tree( ft );
    TR<<"Setting jump now"<<std::endl;
    pose.set_jump( 1, jump_library_[ current_entry_ ] );
    ++current_entry_;
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
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
  RB_dbase( tag->getOption< std::string >( "rigid_body_dbase", "" ) );
  from_entry( tag->getOption< core::Size >( "from_entry", 1 ) );
  to_entry( tag->getOption< core::Size >( "to_entry", 1 ) );
  randomize( tag->getOption< bool >( "randomize", true  ) );

  runtime_assert( from_entry() <= to_entry() && from_entry() >= 1 );
  TR<<"RB_dbase: "<<RB_dbase()<<" from_entry: "<<from_entry()<<" to_entry: "<<to_entry()<<" randomize: "<<randomize()<<std::endl;
}
} // simple_moves
} // protocols

