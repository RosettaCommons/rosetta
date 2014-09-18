// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MaxSeqSepConstraintSet.cc
/// @brief
/// @detailed
/// @author Oliver Lange
///


// Unit Headers
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>

// Package Headers
//#include <protocols/jumping/JumpSetup.fwd.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
/*
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyMap.fwd.hh>
*/
/*
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh> // SequenceMover
#include <protocols/moves/TrialMover.hh>
*/

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/numeric.functions.hh>

// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>


//// C++ headers
#include <cstdlib>
#include <string>
// AUTO-REMOVED #include <fstream>

#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.constraints_additional.MaxSeqSepConstraintSet", basic::t_info );

using core::scoring::constraints::ConstraintSet;
using core::scoring::constraints::ConstraintSetOP;
using core::kinematics::ShortestPathInFoldTree;
using core::kinematics::ShortestPathInFoldTreeOP;

namespace protocols {
namespace constraints_additional {

using namespace core;

/// @brief a ConstraintsSet whose constraints can be switched off, according to sequence separation in residues
/// between residue pair constraints.
MaxSeqSepConstraintSet::~MaxSeqSepConstraintSet() {}
MaxSeqSepConstraintSet::MaxSeqSepConstraintSet( ConstraintSet const & other, core::kinematics::FoldTree const&f ) :
  ConstraintSet( other )
{
  tr.Trace << f ;
	// if ( f.size() > 2 ) { //assuming simple fold-tree has eges 1 1 -1 and 1 nres -1
	shortest_path_ = new ShortestPathInFoldTree( f );
		//} else {
    //shortest_path_ = NULL;
		//  }
}

/// @copy constructor. Does nothing.
MaxSeqSepConstraintSet::MaxSeqSepConstraintSet( MaxSeqSepConstraintSet const &other )
  : ConstraintSet( other ),
    max_seq_sep_ ( other.max_seq_sep_ ),
    shortest_path_( other.shortest_path_ )
{	}

MaxSeqSepConstraintSet::MaxSeqSepConstraintSet( ConstraintSet const &other, ShortestPathInFoldTreeOP sp ) :
	ConstraintSet( other ),
	shortest_path_( sp )
{}

ConstraintSetOP MaxSeqSepConstraintSet::remapped_clone(
    core::pose::Pose const& src,
		core::pose::Pose const& dest,
		core::id::SequenceMappingCOP smap
) const {
	MaxSeqSepConstraintSetOP clone_ptr = new MaxSeqSepConstraintSet( *Parent::remapped_clone( src, dest, smap ), shortest_path_ );
	clone_ptr->set_max_seq_sep( max_seq_sep() );
	return clone_ptr;
}

/*void
MaxSeqSepConstraintSet::eval_atom_derivative_for_residue_pairs (
  id::AtomID const & atom_id,
  pose::Pose const & pose,
  scoring::ScoreFunction const &,
  scoring::EnergyMap const & weights,
  Vector & F1,
  Vector & F2
) const
{
	using scoring::constraints::ResidueConstraints;
	// residue pair constraints:
	Size const seqpos( atom_id.rsd() );
	for ( ResidueConstraints::const_iterator
			it= residue_pair_constraints_begin( seqpos ), ite = residue_pair_constraints_end( seqpos );
			it != ite; ++it ) {
		Size const seqpos2( it->first );
		if ( !too_far( seqpos2, seqpos ) ) {
			//				trDebug << "evaluate " << seqpos << " " << seqpos2;  it->second->show(trDebug); trDebug << "\n";
			it->second->eval_atom_derivative( atom_id, pose.conformation(), weights, F1, F2 );
		}
	}
}*/


///
void
MaxSeqSepConstraintSet::residue_pair_energy(
   Residue const & rsd1,
   Residue const & rsd2,
   Pose const & pose,
   core::scoring::ScoreFunction const & scorefxn,
   core::scoring::EnergyMap & emap
) const
{
  int const pos1( rsd1.seqpos() ), pos2 ( rsd2.seqpos() );
  //	if ( tr.Trace.visible() ) tr.Trace << "max_seq_sep(): " << max_seq_sep() << " check seqsep for residues " << rsd1.seqpos() << "  and  " << rsd2.seqpos() << "....";
  if ( too_far( pos1, pos2 ) ) {
    //if ( tr.Trace.visible() ) tr.Trace << "\n";
    return; //cast avoids warning
  }
  //	if ( tr.Trace.visible() ) tr.Trace << "evaluated\n";
  ConstraintSet::residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap );
}

void
MaxSeqSepConstraintSet::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & minmap,
	core::scoring::ResSingleMinimizationData const & res1_data_cache,
	core::scoring::ResSingleMinimizationData const & res2_data_cache,
	core::scoring::ResPairMinimizationData & respair_data_cache
) const
{
	if ( too_far( rsd1.seqpos(), rsd2.seqpos() ) ) return;
	Parent::setup_for_minimizing_for_residue_pair( rsd1, rsd2, pose, sfxn, minmap, res1_data_cache, res2_data_cache, respair_data_cache );
}


Size
MaxSeqSepConstraintSet::show_violations( std::ostream& out, pose::Pose& pose, Size verbose_level, Real threshold ) {
  out << " total constr: " << get_all_constraints().size() << "  ";
  out << " max_seq_sep: " << max_seq_sep() << " ";
  return Parent::show_violations( out, pose, verbose_level, threshold );
}

bool
MaxSeqSepConstraintSet::too_far( int const pos1, int const pos2 ) const {
  if ( shortest_path_ ) return max_seq_sep() < shortest_path_->dist( pos1, pos2 );
  return max_seq_sep() < (core::Size) std::abs( pos1 - pos2 );
}

Size
MaxSeqSepConstraintSet::largest_possible_sequence_sep( core::pose::Pose const& pose ) const {
  if ( shortest_path_ ) return shortest_path_->max_dist();
  return pose.total_residue();
}

/// Does *NOT* zero the emap values, just adds the additional contribution to the
/// existing emap energies (so can be called inside finalize_total_energies)
void
MaxSeqSepConstraintSet::eval_non_residue_pair_energy(
  core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::scoring::EnergyMap & emap
) const
{
	//non_residue_pair_constraints_.conformation_energy( pose.conformation(), sfxn.weights(), emap );
	core::scoring::func::ConformationXYZ const xyz_func( pose.conformation() );
	core::scoring::constraints::ConstraintCOPs const& csts( non_residue_pair_constraints().constraints() );
	for ( core::scoring::constraints::ConstraintCOPs::const_iterator it=csts.begin(), ite = csts.end(); it != ite; ++it ) {
		core::scoring::constraints::Constraint const & cst( **it );
		if ( cst.effective_sequence_separation( shortest_path() ) < max_seq_sep() ) cst.score( xyz_func, sfxn.weights(), emap );
	}
}




} //abinitio
} //protocols
