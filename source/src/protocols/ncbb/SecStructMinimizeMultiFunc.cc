// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// Unit headers
#include <protocols/ncbb/SecStructMinimizeMultiFunc.hh>
#include <protocols/ncbb/util.hh>

/// Package headers
#include <numeric/crick_equations/HelixParams.hh>

/// Project headers
#include <core/pose/Pose.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <numeric/model_quality/rms.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <core/optimization/Multifunc.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/util.tmpl.hh>

/// Utility headers
#include <cmath>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

static basic::Tracer TR("protocols.ncbb.SecStructMinimizeMultiFunc");

namespace protocols {
namespace ncbb {

using namespace core;
using namespace core::optimization;

SecStructMinimizeMultiFunc::~SecStructMinimizeMultiFunc() {}

SecStructMinimizeMultiFunc::SecStructMinimizeMultiFunc(
	core::pose::Pose &pose,
	scoring::ScoreFunction & scorefxn_in,
	core::optimization::MinimizerMap & min_map,
	std::string const alpha_beta_pattern,
	std::string const dihedral_pattern
) :
	pose_( pose ),
	scorefxn_( scorefxn_in ),
	min_map_( min_map ),
	pose0_( pose ),
	alpha_beta_pattern_( alpha_beta_pattern ),
	dihedral_pattern_( dihedral_pattern )
{

	setup_minimization_graph( pose_, scorefxn_, min_map_ );

	// The alpha beta pattern and dihedral pattern
	// give a mapping between torsion IDs and the variable index that will control them
	// (hooray!)

	utility::vector1< char > uniqs;
	Size num;
	count_uniq_char( dihedral_pattern_, num, uniqs );
	utility::vector1< Size > starters( uniqs.size() );

	starters[ 1 ] = 1;
	for ( Size i = 2; i <= uniqs.size(); ++i ) {
		bool is_beta = false;
		for ( Size j = 1; j <= dihedral_pattern_.size(); ++j ) {
			if ( dihedral_pattern_[ j ] == uniqs[ i-1 ] ) {
				if ( alpha_beta_pattern_[ j ] == 'B' ) {
					is_beta = true;
				}
				j = dihedral_pattern_.size()+1;
			}
		}
		starters[ i ] = starters[ i-1 ] + ( is_beta ? 3 : 2 );
	}

	nvar_ = 0;
	for ( Size i = 1; i <= uniqs.size(); ++i ) {
		//TR << "What is the nature of the uniq " << uniqs[ i ] << "? " << std::endl;
		for ( Size j = 1; j <= dihedral_pattern_.size(); ++j ) {
			//TR << "Examinining dihedral pattern element j = " << j << " for comparison to uniqs[ " << i << " ] = " << uniqs[ i ] << std::endl;
			if ( dihedral_pattern_[ j ] == uniqs[ i ] ) {
				//TR << "Found an example thereof " << std::endl;
				if ( alpha_beta_pattern_[ j ] == 'A' || alpha_beta_pattern_[ j ] == 'P' ) {
					//TR << "It's an alpha!" << std::endl;
					nvar_ += 2;
				} else {
					nvar_ += 3;
				}
				j = dihedral_pattern_.size() + 1;
			}
		}
	}
	TR << "nvar_ = " << nvar_ << std::endl;

	for ( Size i = 1; i <= uniqs.size(); ++i ) {
		for ( Size j = 1; j <= dihedral_pattern_.size(); ++j ) {
			if ( dihedral_pattern_[ j ] == uniqs[ i ] ) {
				//TR << "Dihedral pattern at " << j << " equals uniqs[ " << i << " ] = " << uniqs[ i ] << std::endl;
				// we're looking at one of many residues represented by this uniq
				// therefore, identify this residue (resnum j!)
				// and push back its torsions
				id::TorsionID bb1( j, id::BB, 1 );
				vars_index_to_torsion_id_[ starters[ i ] ].push_back( bb1 );
				//TR << "So variable " << starters[ i ] << " controls residue " << j << " torsion 1 " << std::endl;
				id::TorsionID bb2( j, id::BB, 2 );
				vars_index_to_torsion_id_[ starters[ i ] + 1 ].push_back( bb2 );
				//TR << "So variable " << (starters[ i ]+1) << " controls residue " << j << " torsion 2 " << std::endl;
				if ( alpha_beta_pattern_[ j ] == 'B' ) {
					id::TorsionID bb3( j, id::BB, 3 );
					vars_index_to_torsion_id_[ starters[ i ] + 2 ].push_back( bb3 );
					//TR << "So variable " << (starters[ i ]+2) << " controls residue " << j << " torsion 3 " << std::endl;
				}
			}
		}
	}

	get_dofs_map();
	get_dofs_for_pose0();
}

Real
SecStructMinimizeMultiFunc::operator ()( Multivec const & vars ) const
{
	//TR << "Vars ";
	//for ( Size i = 1; i <= vars.size(); ++i ) {
	// TR << vars[i] << " ";
	//}
	//TR << std::endl;
	Multivec dofs( vars_to_dofs( vars ) );
	min_map_.copy_dofs_to_pose( pose_, dofs );

	return scorefxn_( pose_ );
}

//scoring::MinimizationGraphOP
void
SecStructMinimizeMultiFunc::setup_minimization_graph(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	MinimizerMap const & min_map
) const {
	scoring::MinimizationGraphOP mingraph( new scoring::MinimizationGraph( pose.total_residue() ) );

	scoring::EnergyMap dummy;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		sfxn.setup_for_minimizing_for_node( * mingraph->get_minimization_node( ii ),
			pose.residue( ii ), min_map, pose, false, dummy );
	}

	for ( graph::Graph::EdgeListIter
			eiter = mingraph->edge_list_begin(), eiter_end = mingraph->edge_list_end();
			eiter != eiter_end; ++eiter ) {
		Size const node1 = (*eiter)->get_first_node_ind();
		Size const node2 = (*eiter)->get_second_node_ind();

		scoring::MinimizationEdge & minedge( static_cast< scoring::MinimizationEdge & > (**eiter) );

		sfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
			pose.residue( node1 ), pose.residue( node2 ),
			minedge, min_map, pose, true, false,
			static_cast< scoring::EnergyEdge const * > ( 0 ), dummy );
	}

	/// Now initialize the long-range edges
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( scoring::ScoreFunction::LR_2B_MethodIterator
				iter = sfxn.long_range_energies_begin(),
				iter_end = sfxn.long_range_energies_end();
				iter != iter_end; ++iter ) {

			if ( (*iter)->minimize_in_whole_structure_context( pose ) ) continue;

			scoring::LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
			if ( !lrec || lrec->empty() ) continue;

			// Potentially O(N) operation...
			for ( scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_neighbor_iterator_begin( ii ), // traverse both upper and lower neighbors
					rniend = lrec->const_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const r1 = rni->lower_neighbor_id();
				Size const r2 = rni->upper_neighbor_id();
				Size const jj = ( r1 == ii ? r2 : r1 );
				bool const res_moving_wrt_eachother( true );

				if ( jj < ii ) continue; // only setup each edge once.

				conformation::Residue const & lower_res( r1 == ii ? pose.residue( ii ) : pose.residue( jj ) );
				conformation::Residue const & upper_res( r1 == ii ? pose.residue( jj ) : pose.residue( ii ) );
				sfxn.setup_for_lr2benmeth_minimization_for_respair(
					lower_res, upper_res, *iter, *mingraph, min_map, pose,
					res_moving_wrt_eachother, false, rni, dummy );
			}
		}
	}

	pose.energies().set_minimization_graph( mingraph );
	//return mingraph;
}


void
SecStructMinimizeMultiFunc::dfunc( Multivec const & vars, Multivec & dE_dvars ) const
{
	Multivec dE_ddofs;
	Multivec dofs = vars_to_dofs( vars );

	atom_tree_dfunc( pose_, min_map_, scorefxn_, dofs, dE_ddofs );
	dE_dvars = dEddofs_to_dEdvars( dE_ddofs );

	return;
}

void
SecStructMinimizeMultiFunc::get_dofs_for_pose0()
{
	dofs_for_pose0_.resize( min_map_.dof_nodes().size() );
	//TR << "Resized the pose0 dofs to hold " << min_map_.dof_nodes().size() << " vals " << std::endl;
	min_map_.copy_dofs_from_pose( pose0_, dofs_for_pose0_ );
}

Multivec
SecStructMinimizeMultiFunc::vars_to_dofs( Multivec const & vars ) const {

	Multivec dofs( dofs_for_pose0_ );

	for ( Size i_var = 1; i_var <= vars.size(); ++i_var ) {
		std::map< Size, utility::vector1< Size > >::const_iterator it = map_BB_to_DOF_.find( i_var );
		//TR << "Finding DOFs for " << i_var << std::endl;
		if ( it != map_BB_to_DOF_.end() ) {

			utility::vector1< Size > const i_dofs( it->second );
			//TR << "There are " << i_dofs.size() << " dofs! " << std::endl;
			// Overwrite values into dofs
			for ( Size iii = 1; iii <= i_dofs.size(); ++iii ) {
				//TR << "i_dofs[ iii ]  = i_dofs[ " << iii << " ] = " << i_dofs[ iii ] << std::endl;
				dofs[ i_dofs[ iii ] ] = vars[ i_var ];
				//TR << "Setting dofs[ " << i_dofs[ iii ] << " ] to " << vars[ i_var ] << std::endl;
			}
		}
	}

	return dofs;
} // end vars_to_dofs

Multivec
SecStructMinimizeMultiFunc::dofs_to_vars( Multivec const & dofs ) const
{

	// Next, iter over vars to convert dofs info to vars info
	Multivec vars( nvar_, 0.0 );

	for ( Size i_var = 1; i_var <= nvar_; ++i_var ) {
		std::map< Size, utility::vector1< Size > >::const_iterator it = map_BB_to_DOF_.find( i_var );
		if ( it != map_BB_to_DOF_.end() ) {
			utility::vector1< Size > const i_dofs( it->second );

			for ( Size iii = 1; iii <= i_dofs.size(); ++iii ) {
				vars[ i_var ] = dofs[ i_dofs[ iii ] ];
			}
		}
	}

	return vars;
} // end dof_to_vars

Multivec
SecStructMinimizeMultiFunc::dEddofs_to_dEdvars( Multivec const & dEddofs ) const
{
	Multivec dEdvars( nvar_, 0.0 );

	for ( Size i_var = 1; i_var <= nvar_; ++i_var ) {
		std::map< Size, utility::vector1< Size > >::const_iterator it = map_BB_to_DOF_.find( i_var );
		if ( it != map_BB_to_DOF_.end() ) {

			utility::vector1< Size > const i_dofs( it->second );

			// average!
			for ( Size iii = 1; iii <= i_dofs.size(); ++iii ) {
				dEdvars[ i_var ] += dEddofs[ i_dofs[ iii ] ];
			}
			dEdvars[ i_var ] /= i_dofs.size();
		}
	}

	return dEdvars;
}

void
SecStructMinimizeMultiFunc::get_dofs_map()
{
	std::list< DOF_NodeOP > dof_nodes( min_map_.dof_nodes() );

	// Map between var index and DOF
	map_BB_to_DOF_.clear();
	map_DOF_to_BB_.clear();

	Size imap = 1;

	for ( std::list< DOF_NodeOP >::const_iterator it=dof_nodes.begin(),
			it_end = dof_nodes.end(); it != it_end; ++it, ++imap ) {

		id::TorsionID const tor_id( (*it)->torsion_id() );

		for ( Size i = 1; i <= nvar_; ++i ) {
			//TR << "Initializing mapped DOFs for " << i << std::endl;

			utility::vector1< id::TorsionID > torsions = vars_index_to_torsion_id_[ i ];
			//TR << "Grabbed " << torsions.size() << " torsions " << std::endl;
			for ( Size j = 1; j <= torsions.size(); ++j ) {
				id::TorsionID id = torsions[ j ];

				if ( id == tor_id ) {
					//TR << "Passing " << i << " to get a dof will also get you dof number " << imap << std::endl;
					map_BB_to_DOF_[ i ].push_back( imap );

					//TR << "Passing " << imap << " to get a var number will get you var number " << i << std::endl;
					map_DOF_to_BB_[ imap ] = i;
					break;
				}
			}
		}

	}
}


/// @details Useful debugging code that can be re-enabled by changing the boolean
/// variables at the top of this function.
void
SecStructMinimizeMultiFunc::dump( Multivec const &/*vars*/, Multivec const &/*vars2*/ ) const {
	return;
}

} // namespace ncbb
} // namespace protocols

