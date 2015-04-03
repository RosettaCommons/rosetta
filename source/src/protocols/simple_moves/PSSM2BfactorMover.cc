// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.PSSM2BfactorMover" );
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <boost/foreach.hpp>
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
	return protocols::moves::MoverOP( new PSSM2BfactorMover );
}

std::string
PSSM2BfactorMoverCreator::mover_name()
{
	return "PSSM2Bfactor";
}

PSSM2BfactorMover::PSSM2BfactorMover()
: moves::Mover("PSSM2Bfactor"),
  min_value_( -1.0 ),
  max_value_( 5.0 )
{
}

PSSM2BfactorMover::PSSM2BfactorMover(core::Real const min_in , core::Real const max_in)
: moves::Mover("PSSM2Bfactor"),
  min_value_( min_in ),
  max_value_( max_in )
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

	core::pose::add_comment(pose,"PSSMSTYLE residueNo", " A C D E F G H I K L M N P Q R S T V W Y");

	BOOST_FOREACH( ConstraintCOP const c, constraints ){
		if( c->type() == "SequenceProfile" ){
			SequenceProfileConstraintCOP seqprof_cst( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::SequenceProfileConstraint const > ( c ) );
			runtime_assert( seqprof_cst != 0 );
			runtime_assert( seqprof_cst->profile_mapping() != 0 );
			core::Size const seqpos( seqprof_cst->seqpos() );
			TR<<"sepos="<<seqpos<<std::endl;
			core::id::SequenceMappingCOP SM = seqprof_cst->profile_mapping();
			core::id::SequenceMapping tempSM = *SM;
			if (SM==0) {
				TR << "asdfasd" << std::endl;
			}
			TR<<"seqpos_mapping:"<<tempSM[seqpos]<<std::endl;
			SequenceProfileCOP seqprof_pos( seqprof_cst->sequence_profile() );
			TR<<"the size of seqprof: "<<seqprof_pos->size()<<std::endl;
			core::Real const PSSM_score( seqprof_pos->prof_row( tempSM[seqpos] )[ order[ pose.residue( seqpos ).name1() ] ] );
			//Jun13 Gideon Lapidoth and Chris Norn added this so user can choose cut off PSSM values

			core::Real min = 0.0 ;// Minimum Bfactor level used for pymol coloring
			core::Real max = 50.0 ;// Maximum Bfactor level used for pymol coloring
			core::Real alpha = ( min - max ) / ( max_value_ -  min_value_ );
			core::Real beta = ( min * min_value_ - max * max_value_ ) / ( min_value_ - max_value_ ) ;

			core::Real tmpscore = PSSM_score;
			if ( tmpscore >= max_value_ ) {
				TR<<tmpscore<<" greater than "<<max_value_<<std::endl;
				tmpscore = max_value_;
			}
			else if ( tmpscore <= min_value_ ) {
				tmpscore = min_value_;
			}

			core::Real const transformed_score( tmpscore * alpha + beta );
			TR<<"Position: "<<seqpos<<" pssm_val: "<<PSSM_score<<" min val: "<<min_value_<<" max val: "<<max_value_<<" tmpscore: "<<tmpscore<<" transformed score: "<<transformed_score<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;

			std::ostringstream residuesSpecificPSSMScores;
			for (int i = 1; i<=20; ++i){
				residuesSpecificPSSMScores << seqprof_pos->prof_row( tempSM[seqpos] )[i] << " ";
			}

			std::ostringstream PSSM_resNo;
			PSSM_resNo << "PSSM " << seqpos << " ";
			core::pose::add_comment(pose,PSSM_resNo.str(),residuesSpecificPSSMScores.str());


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
	return moves::MoverOP( new PSSM2BfactorMover( *this ) );
}

moves::MoverOP
PSSM2BfactorMover::fresh_instance() const
{
	return moves::MoverOP( new PSSM2BfactorMover );
}

void
PSSM2BfactorMover::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	max_value( tag->getOption< core::Real >( "Value_for_blue", 5.0 ) );
	chain_num( tag->getOption< core::Real >( "chain_num", 1 ) );
	min_value( tag->getOption< core::Real >( "Value_for_red", -1.0 ));
	TR<<"PSSM2Bfactor sets Value_for_blue: "<<max_value_<<" Value_for_red: "<<min_value_<<std::endl;
}
} // simple_moves
} // protocols
