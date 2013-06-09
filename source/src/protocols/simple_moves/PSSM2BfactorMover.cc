// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PSSM2BfactorMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/PSSM2BfactorMover.hh>
#include <protocols/simple_moves/PSSM2BfactorMoverCreator.hh>
#include <protocols/moves/DataMapObj.hh>
#include <protocols/moves/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.simple_moves.PSSM2BfactorMover");
#include <utility/tag/Tag.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

namespace protocols {
namespace simple_moves {

std::string
PSSM2BfactorMoverCreator::keyname() const
{
	return PSSM2BfactorMoverCreator::mover_name();
}

protocols::moves::MoverOP
PSSM2BfactorMoverCreator::create_mover() const {
	return new PSSM2BfactorMover;
}

std::string
PSSM2BfactorMoverCreator::mover_name()
{
	return "PSSM2Bfactor";
}

PSSM2BfactorMover::PSSM2BfactorMover()
	: moves::Mover("PSSM2Bfactor"),
	offset_( 8.0 ),
	scaling_factor_( 3.0 )
{
}

void
PSSM2BfactorMover::apply( Pose & pose )
{
  using namespace core::scoring::constraints;
  using namespace core::sequence;

  std::map< char, core::Size > order;
  order['A'] = 1; order['C'] = 2; order['D'] = 3; order['E'] = 4; order['F'] = 5; order['G'] = 6; order['H'] = 7;
  order['I'] = 8; order['K'] = 9; order['L'] = 10; order['M'] = 11; order['N'] = 12; order['P'] = 13; order['Q'] = 14;
  order['R'] = 15; order['S'] = 16; order['T'] = 17; order['V'] = 18; order['W'] = 19; order['Y'] = 20;

  ConstraintCOPs constraints( pose.constraint_set()->get_all_constraints() );
  TR<<"Total number of constraints in pose: "<<constraints.size()<<std::endl;
  core::Size cst_num( 0 );
  foreach( ConstraintCOP const c, constraints ){
    if( c->type() == "SequenceProfile" ){
      SequenceProfileConstraintCOP seqprof_cst( dynamic_cast< SequenceProfileConstraint const * >( c() ) );
      runtime_assert( seqprof_cst );
      core::Size const seqpos( seqprof_cst->seqpos() );
      SequenceProfileCOP seqprof_pos( seqprof_cst->sequence_profile() );
      core::Real const PSSM_score( seqprof_pos->prof_row( seqpos )[ order[ pose.residue( seqpos ).name1() ] ] );
		  core::Real const transformed_score( ( PSSM_score + offset() ) * scaling_factor() );
		  TR<<"Position: "<<seqpos<<" pssm_val: "<<PSSM_score<<" transformed score: "<<transformed_score<<std::endl;
		  for( core::Size idx = 1; idx <= pose.residue( seqpos ).natoms(); ++idx ){
		 	  pose.pdb_info()->temperature( seqpos, idx, transformed_score );
		  }
      cst_num++;
    }//fi c->type()=="sequenceprofile"
  }//foreach
	TR<<"Read "<<cst_num<<" sequence constraints"<<std::endl;
}

std::string
PSSM2BfactorMover::get_name() const {
	return PSSM2BfactorMoverCreator::mover_name();
}

moves::MoverOP
PSSM2BfactorMover::clone() const
{
	return new PSSM2BfactorMover( *this );
}

moves::MoverOP
PSSM2BfactorMover::fresh_instance() const
{
	return new PSSM2BfactorMover;
}

void
PSSM2BfactorMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	offset( tag->getOption< core::Real >( "offset", 8.0 ) );
	scaling_factor( tag->getOption< core::Real >( "scaling_factor", 3.0 ));
	TR<<"PSSM2Bfactor sets offset: "<<offset()<<" scaling_factor: "<<scaling_factor()<<std::endl;
}
} // simple_moves
} // protocols
