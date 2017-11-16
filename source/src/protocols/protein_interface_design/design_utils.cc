// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/design_utils.cc
/// @brief various utilities for interface design.
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Project Headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/chemical/AA.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/symmetry/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/NonResidueTypeConstraint.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>

// Unit Headers
#include <protocols/protein_interface_design/design_utils.hh>
#include <protocols/simple_moves/ddG.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
#include <map>

#include <utility/vector0.hh>
#include <boost/lexical_cast.hpp>

// option key includes

namespace protocols {
namespace protein_interface_design {

using ObjexxFCL::format::A;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.protein_interface_design.design_utils" );

using namespace core;

typedef core::Real Real;
typedef core::Size Size;
typedef core::pose::Pose Pose;

// it is assumed that the pose is scored prior to calling this function
core::Real
sum_total_residue_energy( pose::Pose const & pose, core::Size const resid )
{
	using namespace core::scoring;

	typedef utility::vector1<ScoreType> ScoreTypeVec;

	ScoreTypeVec score_types;
	EnergyMap weights = pose.energies().weights();
	for ( core::Size i = 1; i <= n_score_types; ++i ) {
		ScoreType const st = ScoreType( i );
		if ( weights[ st ] != 0 ) score_types.push_back( st );
	}

	core::Real residue_total( 0.0 );
	for ( ScoreTypeVec::const_iterator it=score_types.begin(); it!=score_types.end(); ++it ) {
		residue_total+=(weights[*it] * pose.energies().residue_total_energies( resid )[ *it ]);
	}
	return( residue_total);
}

void
ReportSequenceDifferences::calculate( pose::Pose const & pose1_in, pose::Pose const & pose2_in )
{
	using namespace core::scoring;

	core::pose::Pose pose1( pose1_in );
	core::pose::Pose pose2( pose2_in );

	ScoreFunctionOP scorefxn1 = scorefxn_;
	ScoreFunctionOP scorefxn2 = scorefxn_;
	(*scorefxn1)(pose1);
	(*scorefxn2)(pose2);

	/// Now handled automatically.  scorefxn1->accumulate_residue_total_energies( pose1 );
	/// Now handled automatically.  scorefxn2->accumulate_residue_total_energies( pose2 );

	// core::scoring::EnergyMap weights1 = pose1.energies().weights(); // Unused variable causes a warning.
	// core::scoring::EnergyMap weights2 = pose2.energies().weights(); // Unused variable causes a warning.

	runtime_assert( pose1.size() == pose2.size() );
	for ( core::Size i = 1; i <= pose1.size(); ++i ) {
		if ( !pose1.residue(i).is_protein() ) continue;
		core::Size const restype1( pose1.residue(i).aa() );
		core::Size const restype2( pose2.residue(i).aa() );

		if ( restype1 != restype2 ) {
			res_energy1_.insert( std::make_pair( i, sum_total_residue_energy( pose1, i ) ));
			res_energy2_.insert( std::make_pair( i, sum_total_residue_energy( pose2, i ) ));

			res_name1_.insert( std::make_pair( i, pose1.residue(i).name3() ));
			res_name2_.insert( std::make_pair( i, pose2.residue(i).name3() ));
		}
	}
}

void
ReportSequenceDifferences::report( std::ostream & out ) const
{
	auto it_energy1=res_energy1_.begin();
	auto it_energy2=res_energy2_.begin();
	auto it_name1=res_name1_.begin();
	auto it_name2=res_name2_.begin();

	if ( it_energy1 == res_energy1_.end() ) {
		out<<"No changes\n";
		return;
	} else { out <<res_energy1_.size()<<" changes:\n"; }

	while ( it_energy1!=res_energy1_.end() ) {
		TR<<it_name1->second<<it_name1->first<<" "<<it_energy1->second<<'\t'<<
			it_name2->second<<it_name2->first<<" "<<it_energy2->second<<std::endl;
		++it_energy1; ++it_energy2; ++it_name1; ++it_name2;
	}
}

core::Real
ddG_cycles( pose::Pose const & pose, core::scoring::ScoreFunctionOP scorefxn, core::Size const cycles )
{
	pose::Pose temp_pose = pose;
	protocols::simple_moves::ddG ddg( scorefxn );
	core::Real ddG_val( 0.0 );
	for ( core::Size i = 1; i<=cycles; ++i ) {
		ddg.calculate( temp_pose );
		ddG_val += ddg.sum_ddG();
		temp_pose = pose; // reset the pose
	}
	ddG_val = ddG_val / cycles;
	return( ddG_val );
}

void
point_mutation( pose::Pose & pose, core::scoring::ScoreFunctionCOP scorefxn, core::Size const seqpos, core::Size const mutation )
{
	utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false );
	allowed_aas[ mutation ] = true;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( i == seqpos ) {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
		} else {
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	pack::pack_rotamers( pose, *scorefxn, task );
}

void
Revert::apply( pose::Pose & pose_wt, pose::Pose & pose_des ) const
{
	using namespace core::scoring;

	pose::Pose const saved_des = pose_des;

	ScoreFunctionOP scorefxn_wt = scorefxn_;
	ScoreFunctionOP scorefxn_des = scorefxn_;

	TR<<"Averaging all ddg calculations over "<<ddg_cycles_<<" iterations\n";
	core::Real const ddG_des_val( ddG_cycles( pose_des, scorefxn_des, ddg_cycles_));
	TR<<"average ddG for design: "<<ddG_des_val<<'\n';

	(*scorefxn_wt)(pose_wt);
	(*scorefxn_des)(pose_des);

	ReportSequenceDifferences seq_diff( scorefxn_ );
	seq_diff.calculate( pose_wt, pose_des );
	typedef std::map< core::Size, core::Real> EnMap;
	EnMap const *energy_map2( seq_diff.get_res_energy( 2 ));
	// now substitute each w/t sidechain on pose_des and measure ddg. if ddg doesn't change, apply it
	std::vector< int > revert_positions;
	std::vector< int > ala_positions;
	core::Size count_changes_orig( 0 ), count_changes_revert( 0 );
	std::string ignored_resids( "select resi " );
	std::string done_resids( "select resi " );
	bool first_pass_ignored( true ), first_pass_done( true );
	for ( auto const & it2 : *energy_map2 ) {
		using boost::lexical_cast;
		using std::string;
		core::Size const seqpos( it2.first );

		TR<<pose_des.residue(seqpos).name3()<<seqpos<<"->"<<pose_wt.residue(seqpos).name3()<<" ";
		pose_des.copy_segment( 1, pose_wt, seqpos, seqpos );
		core::Real const ddG_revert_val( ddG_cycles( pose_des, scorefxn_des, ddg_cycles_ ));
		TR<<"ddG change "<<ddG_revert_val-ddG_des_val<<". ";
		pose_des = saved_des;
		if ( ddG_revert_val <= ddG_des_val + ddg_tolerance_ ) {
			TR<<" Done.\n";
			if ( first_pass_done ) {
				first_pass_done = false;
			} else {
				done_resids += "+";
			}
			done_resids += lexical_cast< string >( seqpos );
			revert_positions.push_back( seqpos );
			++count_changes_orig;
		} else {
			TR<<" Ignored.\n";
			++count_changes_revert; ++count_changes_orig;
			if ( first_pass_ignored ) {
				first_pass_ignored = false;
			} else {
				ignored_resids += "+";
			}
			ignored_resids += lexical_cast< string >( seqpos );

			if ( it2.second > 0 ) {
				TR<<"but the total energy for "<<pose_des.residue(seqpos).name3()<<seqpos<<" is "<<it2.second<<" testing an Ala substitution\n";
				TR<<"mutation "<<pose_des.residue(seqpos).name3()<<seqpos<<"->ALA has ddG ";
				point_mutation( pose_des, scorefxn_des, seqpos, chemical::aa_ala );
				core::Real const ddG_ala( ddG_cycles( pose_des, scorefxn_des, ddg_cycles_ ));
				TR << ddG_ala<<" and will be ";
				if ( ddG_ala <= ddG_des_val ) {
					TR<<"kept\n";
					ala_positions.push_back( seqpos );
				} else { TR<<"ignored\n"; }
			}
		}
	}
	pose_des = saved_des;
	for ( std::vector<int>::const_iterator it_rev=revert_positions.begin(); it_rev!=revert_positions.end(); ++it_rev ) {
		point_mutation( pose_des, scorefxn_des, *it_rev, pose_wt.residue( *it_rev ).aa() );
	}
	for ( std::vector<int>::const_iterator it_ala=ala_positions.begin(); it_ala!=ala_positions.end(); ++it_ala ) {
		point_mutation( pose_des, scorefxn_des, *it_ala, chemical::aa_ala );
	}
	core::Real const ddG_all_changes( ddG_cycles( pose_des, scorefxn_des, ddg_cycles_ ));
	TR<<"Starting ddG "<<ddG_des_val<<" final ddG "<< ddG_all_changes<<"\n";
	TR<<"Sequences changes in design: "<<count_changes_orig<<". Sequence changes in reversion: "<<count_changes_revert<<'\n';
	TR<<"Ignored residues: \n"<<ignored_resids<<'\n';
	TR<<"Done residues: \n"<<done_resids<<'\n';
	TR.flush();
}


///////////////////////////////
FavorNativeResidue::FavorNativeResidue( core::pose::Pose & pose, core::Real const native_residue_bonus )
{
	core::Size const nres( pose.size() );
	for ( core::Size i = 1; i <= nres; ++i ) {
		native_residue_bonus_.push_back( native_residue_bonus );
	}
	add_residue_constraints( pose );
}

////////////////////////////
FavorNativeResidue::FavorNativeResidue( core::pose::Pose & pose, utility::vector1<core::Real> const & native_residue_bonus )
{
	core::Size const nres( pose.size() );
	for ( core::Size i = 1; i <= nres; ++i ) {
		native_residue_bonus_.push_back( native_residue_bonus[i] );
	}
	add_residue_constraints( pose );
}


///////////////////////////////////////////////////////////////////////////////////////////////
void
FavorNativeResidue::add_residue_constraints( core::pose::Pose & pose ) const {
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	core::Size const nres( pose.size() );
	for ( core::Size i=1; i<= nres;  ++i ) {
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new ResidueTypeConstraint( pose, i,  native_residue_bonus_[ i ]) ) ) );
	}

}

///////////////////////////////
FavorNonNativeResidue::FavorNonNativeResidue( Pose & pose, core::Real const non_native_residue_bonus )
{
	core::Size const nres( pose.size() );
	for ( core::Size i = 1; i <= nres; ++i ) {
		non_native_residue_bonus_.push_back( non_native_residue_bonus );
	}
	add_residue_constraints( pose );
}

////////////////////////////
FavorNonNativeResidue::FavorNonNativeResidue( Pose & pose, utility::vector1<core::Real> const & non_native_residue_bonus )
{
	core::Size const nres( pose.size() );
	for ( core::Size i = 1; i <= nres; ++i ) {
		non_native_residue_bonus_.push_back( non_native_residue_bonus[i] );
	}
	add_residue_constraints( pose );
}


///////////////////////////////////////////////////////////////////////////////////////////////
void
FavorNonNativeResidue::add_residue_constraints( pose::Pose & pose ) const {
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	core::Size const nres( pose.size() );
	for ( core::Size i=1; i<= nres;  ++i ) {
		pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new NonResidueTypeConstraint( pose, i,  non_native_residue_bonus_[ i ]) ) ) );
	}

}


/// @details minimize the interface between two partners. If target_residues is defined
/// the fold_tree for minimization is set up between the central residue in the target residues and the nearest residue on the partner.
/// if simultaneous minimization is true, then all dofs are minimized at once.
void
MinimizeInterface(
	pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	utility::vector1< bool > const & min_bb,
	utility::vector1< bool > const & min_sc,
	utility::vector1< bool > const & min_rb,
	bool const optimize_foldtree,
	utility::vector1< core::Size > const & target_residues,
	bool const simultaneous_minimization/* = false */ )
{
	using namespace optimization;

	runtime_assert( min_rb.size() == pose.num_jump() );

	core::Size const nres( pose.size() );
	runtime_assert( min_bb.size() == nres );
	runtime_assert( min_sc.size() == nres );
	runtime_assert( min_rb.size() == pose.num_jump() );

	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
	kinematics::FoldTree const saved_ft( pose.fold_tree() );

	if ( optimize_foldtree && target_residues.size() > 0 ) { //setup new fold_tree for better numerical behaviour between the residue at the centre of target_residues and the nearest residue on the partner
		core::Real min_mean_dist=10000;
		core::Size central_residue( *target_residues.begin() );
		for ( auto res_it1=target_residues.begin(); res_it1!=target_residues.end(); ++res_it1 ) {
			runtime_assert( *res_it1 <= pose.size() );

			core::conformation::Residue const res1( pose.residue(*res_it1) );
			core::Real mean_distance( 0.0 );
			for ( auto res_it2=res_it1+1; res_it2!=target_residues.end(); ++res_it2 ) {
				core::conformation::Residue const& res2( pose.residue(*res_it2) );

				mean_distance += res1.xyz( res1.nbr_atom() ).distance( res2.xyz( res2.nbr_atom() ) ) ;
			}
			if ( mean_distance<=min_mean_dist ) {
				central_residue = *res_it1;
				min_mean_dist = mean_distance;
			}
		}
		bool const central_res_in_chain1( central_residue < pose.conformation().chain_begin( 2 ) );
		core::Size const begin( central_res_in_chain1 ? pose.conformation().chain_begin( 2 ) : pose.conformation().chain_begin( 1 ) );
		core::Size const end  ( central_res_in_chain1 ? pose.conformation().chain_end( 2 )   : pose.conformation().chain_end( 1 ) );

		core::Real min_dist(10000);
		core::Size nearest_res( 0 );
		core::conformation::Residue const res_central( pose.residue( central_residue ) );
		for ( core::Size res=begin; res<=end; ++res ) {
			core::conformation::Residue const res2( pose.residue(res) );
			core::Real const distance( res_central.xyz( res_central.nbr_atom() ).distance( res2.xyz( res2.nbr_atom() ) ) );
			if ( distance<=min_dist ) {
				min_dist = distance;
				nearest_res = res;
			}
		}
		runtime_assert( nearest_res );

		kinematics::FoldTree new_ft;

		core::Size const rb_jump( 1 );
		core::Size const jump_pos1( central_res_in_chain1 ? central_residue : nearest_res );
		core::Size const jump_pos2( central_res_in_chain1 ? nearest_res : central_residue );
		new_ft.clear();
		new_ft.add_edge( jump_pos1, jump_pos2, rb_jump );
		new_ft.add_edge( 1, jump_pos1, kinematics::Edge::PEPTIDE );
		new_ft.add_edge( jump_pos1, pose.conformation().chain_end( 1 ), kinematics::Edge::PEPTIDE );
		new_ft.add_edge( pose.conformation().chain_begin( 2 ), jump_pos2, kinematics::Edge::PEPTIDE );
		new_ft.add_edge( jump_pos2, pose.size(), kinematics::Edge::PEPTIDE );
		new_ft.reorder( 1 );

		TR<<"setting fold_tree for minimization between "<<jump_pos1<<" and "<<jump_pos2<<"\n";
		pose.fold_tree( new_ft );
	} else { //setup new foldtree
		TR<<"Fold tree not optimized in minimize interface.\n"<<pose.fold_tree()<<std::endl;
	}

	core::kinematics::MoveMap mm;
	mm.set_bb( false );
	TR<<"minimizing sc of residues: ";
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		mm.set_chi( i, min_sc[ i ] );
		if ( min_sc[ i ] ) {
			TR<<i<<' ';
		}
	}
	TR<<'\n';
	if ( !simultaneous_minimization ) {
		AtomTreeMinimizer().run( pose, mm, *scorefxn,
			MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );
	}

	for ( core::Size rb_jump=1; rb_jump<=pose.num_jump(); ++rb_jump ) {
		mm.set_jump( rb_jump, min_rb[ rb_jump ] );
	}
	if ( !simultaneous_minimization && std::find( min_rb.begin(), min_rb.end(), true ) != min_rb.end() ) {
		TR<<"minimizing rigid body orientation\n";
		AtomTreeMinimizer().run( pose, mm, *scorefxn,
			MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );
	}

	if ( !simultaneous_minimization ) {
		TR<<"minimizing bb of residues (and sc of the previous subset): ";
	}
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		mm.set_bb( i, min_bb[ i ] );
		mm.set_chi( i, min_sc[ i ] );
		if ( min_bb[ i ] ) {
			TR<<i<<' ';
		}
	}
	TR<<"\nAnd now minimizing all dofs together\n";
	AtomTreeMinimizer().run( pose, mm, *scorefxn,
		MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );

	if ( 0/* target_residues.size() > 0*/ ) { //reset fold_tree
		pose.fold_tree( saved_ft );
	}
	TR.flush();
}

void
SymMinimizeInterface(
	pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	utility::vector1< bool > const & min_bb,
	utility::vector1< bool > const & min_sc,
	utility::vector1< bool > const & min_rb,
	//bool const optimize_foldtree,
	//utility::vector1< core::Size > const target_residues,
	bool const simultaneous_minimization/* = false */ )
{
	using namespace optimization;
	using namespace conformation::symmetry;

	runtime_assert( core::pose::symmetry::is_symmetric( pose ) );
	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`

	core::kinematics::MoveMap mm;
	mm.set_bb( false );
	TR<<"minimizing sc of residues: ";
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		mm.set_chi( i, min_sc[ i ] );
		if ( min_sc[ i ] ) {
			TR<<i<<' ';
		}
	}
	TR<<'\n';
	core::pose::symmetry::make_symmetric_movemap( pose, mm );

	if ( !simultaneous_minimization ) {
		optimization::symmetry::SymAtomTreeMinimizer().run( pose, mm, *scorefxn,
			MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );
	}

	mm.set_jump( true );
	core::pose::symmetry::make_symmetric_movemap( pose, mm );

	if ( min_rb.size() > 0 ) {
		TR<<"minimizing rigid body orientation\n";
		TR<<"By default all dofs in the symmetry input are used!!!Should change?\n";
		optimization::symmetry::SymAtomTreeMinimizer().run( pose, mm, *scorefxn,
			MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );
	}

	if ( !simultaneous_minimization ) {
		TR<<"minimizing bb of residues (and sc of the previous subset): ";
	}
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		mm.set_bb( i, min_bb[ i ] );
		mm.set_chi( i, min_sc[ i ] );
		if ( min_bb[ i ] ) {
			TR<<i<<' ';
		}
	}
	TR<<"\nAnd now minimizing all dofs together\n";
	core::pose::symmetry::make_symmetric_movemap( pose, mm );
	optimization::symmetry::SymAtomTreeMinimizer().run( pose, mm, *scorefxn,
		MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );
	TR.flush();
}

std::list< core::Size >
hbonded(
	Pose const & in_pose, core::Size const target_residue, std::set< core::Size > const & binders,
	bool const bb, bool const sc, core::Real const energy_thres, bool const bb_bb, core::scoring::ScoreFunctionOP sfxn )
{

	using namespace core::scoring::hbonds;

	std::list< core::Size > hbonded_list;
	core::scoring::ScoreFunctionOP scorefxn( sfxn );
	if ( !sfxn ) scorefxn=core::scoring::get_score_function();
	Pose pose( in_pose );
	(*scorefxn)(pose);

	HBondSet background_hbond_set;
	background_hbond_set.setup_for_residue_pair_energies( pose, false/*calculate_derivative*/, true/*backbone_only*/ );
	HBondDatabaseCOP hb_database( HBondDatabase::get_database( background_hbond_set.hbond_options().params_database_tag()));

	if ( bb_bb ) {
		TR << "decomposing bb hydrogen bond terms" << std::endl;
		core::scoring::methods::EnergyMethodOptionsOP energy_options( new core::scoring::methods::EnergyMethodOptions(scorefxn->energy_method_options()) );
		energy_options->hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options(*energy_options);
	}

	EnergyMap hbond_emap;
	core::conformation::Residue const resi( pose.residue( target_residue ));
	core::Real const distance_cutoff( 20.0 );
	for ( core::Size binder : binders ) {
		core::conformation::Residue const resj( pose.residue(binder) );

		core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
		if ( distance > distance_cutoff ) continue;

		HBondSet pair_hbond_set(2);
		identify_hbonds_1way(
			*hb_database,
			resi, resj, background_hbond_set.nbrs(resi.seqpos()), background_hbond_set.nbrs(resj.seqpos()),
			false /*calculate_derivative*/,
			!bb, !sc, !sc, !sc, pair_hbond_set);
		identify_hbonds_1way(
			*hb_database,
			resj, resi, background_hbond_set.nbrs(resj.seqpos()), background_hbond_set.nbrs(resi.seqpos()),
			false /*calculate_derivative*/,
			!bb, !sc, !sc, !sc, pair_hbond_set);

		hbond_emap.zero();
		get_hbond_energies( pair_hbond_set, hbond_emap );
		// The hbond_energies should be controlled by dotting hbond_emap
		// with the weights file used, but this cannot be done without
		// effecting hbond_energy_threshold_ which is calibrated for
		// unweighted hbond energies. Since STANDARD_WTS + SCORE12_PATCH
		// is hard coded, use this instead:
		//  hbond_emap[ hbond_sr_bb_sc ] = 0;
		//  hbond_emap[ hbond_lr_bb_sc ] = 0;

		core::Real total_hbond_energy( hbond_emap.sum() );


		// all of the bb / sc energies are lumped together
		// but it can be controlled whether bb - bb are included or not

		// counting the number of hbonds between the two residues
		if ( total_hbond_energy <= energy_thres ) {
			using namespace core::conformation;
			for ( core::Size i=1; i<=pair_hbond_set.nhbonds(); ++i ) {
				using namespace core::scoring::hbonds;
				HBond const & hb( pair_hbond_set.hbond( i ) );

				if ( !bb && ( hb.don_hatm_is_protein_backbone() && hb.acc_atm_is_protein_backbone() ) ) continue;
				if ( !sc && ( !hb.don_hatm_is_protein_backbone() || !hb.acc_atm_is_protein_backbone() ) ) continue;

				core::Size const don_res_i( hb.don_res() ), acc_res_i( hb.acc_res() );
				Residue const & don_rsd( pose.residue( don_res_i ) ),
					acc_rsd( pose.residue( acc_res_i ) );
				if ( (don_rsd.seqpos() == binder && acc_rsd.seqpos() == target_residue) ||
						( don_rsd.seqpos() == target_residue && acc_rsd.seqpos() == binder ) ) {
					hbonded_list.push_back( binder );
					core::Size const width( 10 );
					TR  << I( width, target_residue )
						<< I( width, binder )
						<< A( width, resi.name1() )
						<< A( width, resj.name1() )
						<< F( width, 3, hbond_emap[ hbond_sr_bb ] )
						<< F( width, 3, hbond_emap[ hbond_lr_bb ] )
						<< F( width, 3, hbond_emap[ hbond_sc ] )
						<< F( width, 3, hbond_emap[ hbond_bb_sc ] )
						<< F( width, 3, distance ) << "\n";
				}//correct residues
			}//hbond num
		}// if energy passes threshold
	} // residue j
	TR.flush();
	return( hbonded_list );
}

std::list< core::Size >
hbonded_atom(
	Pose const & in_pose, core::Size const target_residue, std::string const & target_atom, std::set< core::Size > const & binders,
	bool const bb, bool const sc, core::Real const energy_thres, bool const bb_bb, core::scoring::ScoreFunctionOP sfxn )
{

	using namespace core::scoring::hbonds;

	std::list< core::Size > hbonded_list;
	core::scoring::ScoreFunctionOP scorefxn( (sfxn ? sfxn : core::scoring::get_score_function()) );
	Pose pose( in_pose );
	(*scorefxn)(pose);

	HBondSet background_hbond_set;
	background_hbond_set.setup_for_residue_pair_energies( pose, false/*calculate_derivative*/, true/*backbone_only*/ );
	HBondDatabaseCOP hb_database( HBondDatabase::get_database( background_hbond_set.hbond_options().params_database_tag()));

	if ( bb_bb ) {
		TR << "decomposing bb hydrogen bond terms" << std::endl;
		core::scoring::methods::EnergyMethodOptionsOP energy_options( new core::scoring::methods::EnergyMethodOptions(scorefxn->energy_method_options()) );
		energy_options->hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options(*energy_options);
	}

	EnergyMap hbond_emap;
	core::conformation::Residue const resi( pose.residue( target_residue ));
	core::Real const distance_cutoff( 20.0 );
	for ( core::Size binder : binders ) {
		core::conformation::Residue const resj( pose.residue(binder) );

		core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
		if ( distance > distance_cutoff ) continue;

		HBondSet pair_hbond_set(2);
		identify_hbonds_1way(
			*hb_database,
			resi, resj, background_hbond_set.nbrs(resi.seqpos()), background_hbond_set.nbrs(resj.seqpos()),
			false /*calculate_derivative*/,
			!bb, !sc, !sc, !sc, pair_hbond_set);
		identify_hbonds_1way(
			*hb_database,
			resj, resi, background_hbond_set.nbrs(resj.seqpos()), background_hbond_set.nbrs(resi.seqpos()),
			false /*calculate_derivative*/,
			!bb, !sc, !sc, !sc, pair_hbond_set);

		hbond_emap.zero();
		get_hbond_energies( pair_hbond_set, hbond_emap );
		// The hbond_energies should be controlled by dotting hbond_emap
		// with the weights file used, but this cannot be done without
		// effecting hbond_energy_threshold_ which is calibrated for
		// unweighted hbond energies. Since STANDARD_WTS + SCORE12_PATCH
		// is hard coded, use this instead:
		//  hbond_emap[ hbond_sr_bb_sc ] = 0;
		//  hbond_emap[ hbond_lr_bb_sc ] = 0;

		core::Real total_hbond_energy( hbond_emap.sum() );

		// all of the bb / sc energies are lumped together
		// but it can be controlled whether bb - bb are included or not

		// counting the number of hbonds between the two residues
		if ( total_hbond_energy <= energy_thres ) {
			using namespace core::conformation;
			for ( core::Size i=1; i<=pair_hbond_set.nhbonds(); ++i ) {
				using namespace core::scoring::hbonds;
				HBond const & hb( pair_hbond_set.hbond( i ) );

				if ( !bb && ( hb.don_hatm_is_protein_backbone() && hb.acc_atm_is_protein_backbone() ) ) continue;
				if ( !sc && ( !hb.don_hatm_is_protein_backbone() || !hb.acc_atm_is_protein_backbone() ) ) continue;

				core::Size const don_res_i( hb.don_res() ), acc_res_i( hb.acc_res() );
				core::Size don_atom_i_base ( pose.residue( don_res_i ).atom_base( hb.don_hatm() ) );
				Residue const & don_rsd( pose.residue( don_res_i ) ),
					acc_rsd( pose.residue( acc_res_i ) );

				//access atom information
				core::Size const don_atom_i( hb.don_hatm() ), acc_atom_i( hb.acc_atm() );
				core::Size target_atom_id=resi.atom_index(target_atom);

				if ( ( (don_rsd.seqpos() == binder && acc_rsd.seqpos() == target_residue) && ( don_atom_i==target_atom_id || don_atom_i_base==target_atom_id ) ) ||
						( (don_rsd.seqpos() == binder && acc_rsd.seqpos() == target_residue) && acc_atom_i==target_atom_id ) ||
						( (don_rsd.seqpos() == target_residue && acc_rsd.seqpos() == binder) && ( don_atom_i==target_atom_id || don_atom_i_base==target_atom_id ) ) ||
						( (don_rsd.seqpos() == target_residue && acc_rsd.seqpos() == binder) && acc_atom_i==target_atom_id ) ) {
					hbonded_list.push_back( binder );
					core::Size const width( 10 );
					TR  << I( width, target_residue )
						<< I( width, target_atom )
						<< I( width, binder )
						<< A( width, resi.name1() )
						<< A( width, resj.name1() )
						<< F( width, 3, hbond_emap[ hbond_sr_bb ] )
						<< F( width, 3, hbond_emap[ hbond_lr_bb ] )
						<< F( width, 3, hbond_emap[ hbond_sc ] )
						<< F( width, 3, hbond_emap[ hbond_bb_sc ] )
						<< F( width, 3, distance ) << "\n";
				}//correct residues
			}//hbond num
		}// if energy passes threshold
	} // residue j
	TR.flush();
	return( hbonded_list );
}

} // namespace protocols
} // namespace protein_interface_design
