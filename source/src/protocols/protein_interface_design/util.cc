// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/util.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit Headers
#include <protocols/protein_interface_design/util.hh>

// Package Headers


// Project Headers
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <boost/foreach.hpp>
#include <protocols/filters/Filter.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>


static thread_local basic::Tracer TR( "protocols.protein_interface_design.util" );

namespace protocols {
namespace protein_interface_design {

using namespace core::scoring;
using namespace protocols::moves;
using namespace core;
using namespace std;
using utility::vector1;

/// @brief generate a star fold tree (originating at residue 1)
core::kinematics::FoldTree
star_fold_tree( core::pose::Pose & pose )
{
	kinematics::FoldTree linear_ft;
	linear_ft.clear();
	core::Size const chain_num( pose.conformation().num_chains() );
	core::Size jump_num( 1 );
	for ( core::Size c=1; c<=chain_num; ++c ) {
		core::Size const chain_begin( pose.conformation().chain_begin( c ) );
		core::Size const chain_end( pose.conformation().chain_end( c ) );

		if ( !pose.residue( chain_end ).is_polymer() ) { //ligand
			linear_ft.add_edge( chain_end-1, chain_end, jump_num );
			++jump_num;
		} else {
			linear_ft.add_edge( chain_begin, chain_end, kinematics::Edge::PEPTIDE );
			if ( c > 1 ) {
				linear_ft.add_edge( 1,chain_begin, jump_num );
				++jump_num;
			}
		}
	}
	linear_ft.reorder( 1 );
	TR<<linear_ft<<std::endl;
	pose.fold_tree( linear_ft );
	pose.update_residue_neighbors();
	return( linear_ft );
}

// @details Connecting to the last carbon atom before the residue's action centre will allow                         // mimization to sample the residue's sidechain dofs without harming the interaction                                 // the stub with the target chain
//kdrew: HIS was defaulting to CB but CG makes more sense, also changing LEU from CG to CD2, and ILE, VAL, SER, THR, PRO, CYS not to be default CB
std::string
optimal_connection_point( std::string const & residue_type ){
	std::string connect_to( "CB" ); // to which atom to hook up the atom tree
	if ( residue_type == "GLN" || residue_type == "DGN"  || residue_type == "GLU"  || residue_type == "DGU" ) {
		connect_to = "CD";
	}
	if ( residue_type == "ARG" || residue_type == "DAR" ) {
		connect_to = "CZ";
	}
	if ( residue_type == "MET" || residue_type == "DME" ) {
		connect_to = "SD";
	}
	// AMW: cppcheck, thankfully, notices that DAS was included twice!
	if ( residue_type == "PHE" || residue_type == "DPH" || residue_type == "TRP" || residue_type == "DTR" || residue_type == "TYR" || residue_type == "DTY" || residue_type == "ASN" || residue_type == "DAN" || residue_type == "ASP"|| residue_type == "DAS" || residue_type == "HIS" || residue_type == "DHI" ) {
		connect_to = "CG";
	}
	if ( residue_type == "LEU" || residue_type == "DLE" ) {
		connect_to = "CD2";
	}
	if ( residue_type == "ILE" || residue_type == "DIL" ) {
		connect_to = "CD1";
	}
	if ( residue_type == "VAL" || residue_type == "DVA" ) {
		connect_to = "CG1";
	}
	if ( residue_type == "PRO" || residue_type == "DPR" ) {
		connect_to = "CD";
	}
	if ( residue_type == "SER" || residue_type == "DSE" ) {
		connect_to = "OG";
	}
	if ( residue_type == "THR" || residue_type == "DTH" ) {
		connect_to = "OG1";
	}
	if ( residue_type == "CYS" || residue_type == "DCS" ) {
		connect_to = "SG";
	}
	if ( residue_type == "LYS" || residue_type == "DLY"  ) {
		connect_to = "NZ";
	}
	if ( residue_type == "GLY" ) {
		connect_to = "CA";
	}

	TR.Debug << "residue_type: " << residue_type << " connect_to: " << connect_to << std::endl;

	return( connect_to );
}

/// @details generate a fold tree that connects each hotspot residue through its optimal connection point to
/// the nearest atom on chain 1
core::kinematics::FoldTree
make_hotspot_foldtree( core::pose::Pose const & pose )
{
	using namespace core::chemical;
	using namespace core::kinematics;
	FoldTree ft;
	ft.clear();

	utility::vector1< core::Size > connection_points;
	for ( core::Size jump( 1 ); jump<=pose.num_jump(); ++jump ) {
		core::Size const resi( pose.conformation().chain_begin( jump + 1) );
		std::string const residue_type( pose.residue( resi ).name3() );
		std::string const connect_to( optimal_connection_point( residue_type ) );
		core::Size const nearest_resi_on_target( find_nearest_residue( pose, 1/*target_chain*/, resi, connect_to ) );
		ft.add_edge( nearest_resi_on_target, resi, jump );
		ft.set_jump_atoms( jump, pose.residue( nearest_resi_on_target ).atom_type( pose.residue( nearest_resi_on_target ).nbr_atom() ).element(), connect_to );
		connection_points.push_back( nearest_resi_on_target );
	}
	std::sort( connection_points.begin(), connection_points.end() );
	utility::vector1< core::Size >::iterator last = std::unique( connection_points.begin(), connection_points.end() );
	connection_points.erase( last, connection_points.end() );
	core::Size upstream_position( 1 );
	BOOST_FOREACH ( core::Size const con, connection_points ) {
		ft.add_edge( upstream_position, con, Edge::PEPTIDE );
		upstream_position = con;
	}
	ft.add_edge( connection_points[ connection_points.size() ], pose.conformation().chain_end( 1 ), Edge::PEPTIDE );
	ft.delete_self_edges();
	ft.reorder( 1 );
	return( ft );
}

core::scoring::constraints::ConstraintCOPs
get_bbcsts( core::pose::Pose const & pose ) {
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( backbone_stub_constraint, 1.0 );
	core::pose::Pose  nonconst_pose = pose;
	// pre-score to get the active cst's
	(*scorefxn)(nonconst_pose);
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( nonconst_pose );

	// sort through pose's constraint set and pull out only active bbcst's
	ConstraintCOPs original_csts = nonconst_pose.constraint_set()->get_all_constraints();
	ConstraintCOPs new_csts;
	for ( ConstraintCOPs::const_iterator it = original_csts.begin(), end = original_csts.end(); it != end; ++it ) {
		ConstraintCOP cst( *it );
		if ( cst->type() == "AmbiguousConstraint" ) {
			AmbiguousConstraintCOP ambiguous_cst = AmbiguousConstraintCOP( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::AmbiguousConstraint const > ( cst ) ); //downcast to derived ambiguous constraint
			if ( ambiguous_cst ) { // safety check for downcasting
				if ( ambiguous_cst->active_constraint()->type() == "BackboneStub" ) {
					new_csts.push_back( ambiguous_cst->active_constraint() );
				}
			}
		}
	}
	return new_csts;
}

/// @details a utility function to evaluate backbone_stub_constraints for each residue in a chain and return a vector with the top n_return residue numbers by cst score
/// note that this function is NOT guaranteed to return n_return residues! It will return the best n<=n_return
utility::vector1< core::Size >
best_bbcst_residues( core::pose::Pose const & pose, core::Size const chain, core::Size const n_return )
{
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	core::pose::Pose  nonconst_pose = pose;
	utility::vector1< std::pair<core::Real, core::Size> > all_residues; // score, seqpos
	utility::vector1< core::Size > best_residues;

	// scorefxn containing only the constraint energy
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( backbone_stub_constraint, 1.0 );

	// score to make sure that the cst energies are current
	(*scorefxn)(nonconst_pose);
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( nonconst_pose );

	// get pairs of residue number and weighted cst energy
	for ( core::Size i=nonconst_pose.conformation().chain_begin( chain ); i<=nonconst_pose.conformation().chain_end( chain ); ++i ) {
		core::Real const score( nonconst_pose.energies().residue_total_energies( i )[ backbone_stub_constraint ] );
		core::Real const weight( (*scorefxn)[ backbone_stub_constraint ] );
		core::Real const curr_energy( weight * score );
		all_residues.push_back(std::make_pair( curr_energy, i ));
	}
	sort( all_residues.begin(), all_residues.end() );
	for ( core::Size i=1; i<=n_return; ++i ) {
		// only use it if the cst actually evaluates
		if ( all_residues[i].first < 0 ) best_residues.push_back( all_residues[i].second );
	}
	assert( best_residues.size() <= n_return );
	return best_residues;
}

void
find_lowest_constraint_energy_residue( core::pose::Pose const & pose, core::Size const chain, core::Size & resi, core::Real & lowest_energy )
{
	using namespace core::scoring;
	core::scoring::ScoreFunctionCOP scorefxn( get_score_function() );

	resi = 0;
	lowest_energy = 100000.0;
	for ( core::Size i=pose.conformation().chain_begin( chain ); i<=pose.conformation().chain_end( chain ); ++i ) {
		using namespace core::scoring;
		simple_filters::EnergyPerResidueFilter const eprf( i, scorefxn, backbone_stub_constraint, 10000.0/*dummy threshold*/ );
		core::Real const curr_energy( eprf.compute( pose ) );
		if ( curr_energy<=lowest_energy ) {
			lowest_energy = curr_energy;
			resi = i;
		}
	}
}

/// @details a utility function for removing ALL coordinate constraints from a pose.
/// returns the constraints that were removed
core::scoring::constraints::ConstraintCOPs
remove_coordinate_constraints_from_pose( core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;

	ConstraintCOPs original_csts = pose.constraint_set()->get_all_constraints() ;
	ConstraintCOPs crd_csts;
	for ( ConstraintCOPs::const_iterator it = original_csts.begin(), end = original_csts.end(); it != end; ++it ) {
		ConstraintCOP cst( *it );
		if ( cst->type() == "CoordinateConstraint" ) {
			ConstraintCOP crd_cst = utility::pointer::dynamic_pointer_cast< core::scoring::constraints::CoordinateConstraint const > ( cst ); //downcast to derived ambiguous constraint
			if ( crd_cst ) { // safety check for downcasting
				crd_csts.push_back( cst ); // add the entire ambiguous cst, since it contained a bbcst
			}
		}
	}
	pose.remove_constraints( crd_csts ); // remove all the ambigcsts that contain a bbcst
	return( crd_csts );
}


/// @details utility function for stub_based_atom_tree. tries to find an optimal cutpoint in a pose given two different boundaries. First looks for a 3-res loop stretch on the downstream partner and returns the middle residue. Then, does the same for the upstream chain. Then, becomes desperate and tries to find any loop residue on downstream chain, and then on upstream chain. Finally, if no success, returns 0 which means that no break was found
core::Size
best_cutpoint( core::pose::Pose & pose, core::Size const prev_u, core::Size const prev_d, core::Size const u, core::Size const d )
{
	// if the pose is all loops (mini default), then run dssp
	// this logic may cause a problem for miniproteins that are indeed all loop. But for now it's certainly OK
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		char const ss = pose.secstruct( i );
		if ( ss == 'H' || ss == 'S' ) break;
		if ( i == pose.total_residue() ) {
			core::scoring::dssp::Dssp dssp( pose );
			dssp.insert_ss_into_pose( pose );
		}
	}

	for ( core::Size res = prev_d; res <= d-1; ++res ) {
		if ( pose.secstruct( res ) == 'L' ) {
			if ( pose.secstruct( res + 1 ) == 'L' && pose.secstruct( res + 2 ) == 'L' ) return res;
		}
	}
	for ( core::Size res = prev_u; res <= u-1; ++res ) {
		if ( pose.secstruct( res ) == 'L' ) {
			if ( pose.secstruct( res + 1 ) == 'L' && pose.secstruct( res + 2 ) == 'L' ) return res;
		}
	}
	for ( core::Size res = prev_d; res <= d; ++res ) {
		if ( pose.secstruct( res ) == 'L' ) return res;
	}
	for ( core::Size res = prev_u; res <= u-1; ++res ) {
		if ( pose.secstruct( res ) == 'L' ) return res;
	}
	return 0; // sign of trouble
}

/// @details find the nearest residue on the target chain to res
core::Size
find_nearest_residue( core::pose::Pose const & pose, core::Size const target_chain, core::Size const res, std::string const & atom/*=="CA"*/ )
{
	core::Size nearest_resi( 0 );
	core::Real nearest_dist( 100000.0 );
	for ( core::Size resi( pose.conformation().chain_begin( target_chain ) ); resi<=pose.conformation().chain_end( target_chain ); ++resi ) {
		core::Real const distance( pose.residue(resi).xyz( pose.residue(resi).nbr_atom() ).distance( pose.residue( res ).xyz( atom ) ) );
		if ( distance<=nearest_dist ) {
			nearest_resi = resi;
			nearest_dist = distance;
		}
	}
	runtime_assert( nearest_resi );
	return( nearest_resi );
}


} //protein_interface_design
} //protocols
