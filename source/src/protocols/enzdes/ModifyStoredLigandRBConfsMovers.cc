// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file .cc file for movers that mess around with additional ligand rigid body conformations
/// stored in the enzdes cacheable observer
/// @brief
/// @author Florian Richter, floric@u.washington.edu, oct 09


//unit headers
#include <protocols/enzdes/ModifyStoredLigandRBConfsMovers.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>

//project headers
//#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/pack/rtmin.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/Tracer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/MinMover.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

// C++ headers
#include <set>

#include <utility>
#include <utility/vector1.hh>


//debug headers
//#include <fstream>
//
//#include <utility/string_util.hh>

namespace protocols {
namespace enzdes {

static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.ModifyStoredLigandRBConfsMovers" );
//static core::Size applycalls_ = 0;

ModifyStoredRBConfs::ModifyStoredRBConfs( std::string const & name )
: parent( name ) {}

ModifyStoredRBConfs::~ModifyStoredRBConfs()= default;

std::string
ModifyStoredRBConfs::get_name() const {
	return "ModifyStoredRBConfs";
}

void
ModifyStoredRBConfs::swap_coordinates_in_pose(
	core::pose::Pose & pose,
	core::conformation::Residue & rescoords
) const
{
	Size seqpos( rescoords.seqpos() );
	runtime_assert( pose.residue( seqpos ).type().name3() == rescoords.type().name3() );
	runtime_assert( pose.residue( seqpos ).natoms() == rescoords.natoms() );
	core::Size res_atoms(  pose.residue( seqpos ).natoms() );

	for ( core::Size atm = 1; atm <= res_atoms; ++atm ) {
		core::PointPosition save_pos( pose.residue( seqpos ).xyz( atm ) );
		core::Size other_ind( rescoords.atom_index( pose.residue( seqpos ).atom_name( atm ) ) );
		pose.set_xyz( core::id::AtomID ( atm, seqpos ), rescoords.xyz( other_ind ) );
		rescoords.set_xyz( other_ind, save_pos );
	}
}

utility::vector1< ModifyStoredRBConfs::Real >
ModifyStoredRBConfs::closest_orient_atoms_msd(
	core::pose::Pose const & pose,
	utility::vector1< core::conformation::ResidueCOP > const & confs
) const
{
	utility::vector1< Real > to_return( confs.size(), 0.0 );
	if ( to_return.size() == 0 ) return to_return;
	Size center(0), nbr1(0), nbr2(0);
	core::conformation::Residue const & pose_res( pose.residue( confs[1]->seqpos() ) );
	runtime_assert( confs[1]->name3() == pose_res.name3() );
	confs[1]->select_orient_atoms( center, nbr1, nbr2 );
	runtime_assert( confs[1]->atom_name( center ) == pose_res.atom_name( center ) );
	runtime_assert( confs[1]->atom_name( nbr1 ) == pose_res.atom_name( nbr1 ) );
	runtime_assert( confs[1]->atom_name( nbr2 ) == pose_res.atom_name( nbr2 ) );

	for ( Size i = 1; i <= confs.size(); ++i ) {
		Real min_msd( pose_res.xyz( center ).distance_squared( confs[i]->xyz(center) ) + pose_res.xyz( nbr1 ).distance_squared( confs[i]->xyz(nbr1) ) + pose_res.xyz( nbr2 ).distance_squared( confs[i]->xyz(nbr2) ) );
		for ( Size j = 1; j < i; ++j ) {
			Real cur_msd( confs[j]->xyz( center ).distance_squared( confs[i]->xyz(center) ) + confs[j]->xyz( nbr1 ).distance_squared( confs[i]->xyz(nbr1) ) + confs[j]->xyz( nbr2 ).distance_squared( confs[i]->xyz(nbr2) ) );
			if ( cur_msd < min_msd ) min_msd = cur_msd;
		}
		to_return[i] = min_msd;
	}
	return to_return;
}


ModifyStoredRBConfs::RBConfLists
ModifyStoredRBConfs::get_rigid_body_confs( core::pose::Pose const & pose ) const
{
	RBConfLists to_return;
	toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) ); //toolbox::match_enzdes_util::get_enzdes_observer with const pose can return NULL
	if ( !enz_obs ) return to_return;
	std::map< core::Size, utility::vector1< core::conformation::ResidueCOP > > const & rb_confs( enz_obs->lig_rigid_body_confs() );
	for ( auto const & rb_conf : rb_confs ) {
		if ( pose.residue( rb_conf.first ).type().is_ligand() ) to_return.push_back( rb_conf.second );
	}
	return to_return;
}

void
ModifyStoredRBConfs::set_rigid_body_confs( RBConfLists rbs, core::pose::Pose & pose ) const
{
	for ( core::Size i = 1; i <= rbs.size(); ++i ) {
		if ( rbs[i].size() > 0 ) {
			set_rigid_body_confs_for_seqpos( rbs[i][1]->seqpos(), rbs[i], pose );
		}
	}
}

void
ModifyStoredRBConfs::set_rigid_body_confs_for_seqpos(
	core::Size seqpos,
	utility::vector1< core::conformation::ResidueCOP > & confs,
	core::pose::Pose & pose
) const
{
	toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_rigid_body_confs_for_lig( seqpos, confs );
}

GenerateStoredRBConfs::GenerateStoredRBConfs(
	Size num_total_rbconfs,
	bool include_metals
) :  parent( "GenerateStoredRBConfs" ),
	num_total_rbconfs_( num_total_rbconfs ), include_metals_(include_metals) {}

GenerateStoredRBConfs::~GenerateStoredRBConfs()= default;

std::string
GenerateStoredRBConfs::get_name() const {
	return "GenerateStoredRBConfs";
}

/// @brief two things happen:
/// 1. for ligands which already have multiple
/// conformations stored, their number will be increased to num_total_rbconfs_
/// in case there are already more rb_confs than num_total_rb_confs_, nothing happens
/// 2. for ligands that don't have multiple conformations stored, num_total_rbconfs_
/// will be generated
void
GenerateStoredRBConfs::apply(
	core::pose::Pose & pose )
{
	std::set< Size > present_confs_seqpos;
	DiversifyStoredRBConfs diversifier( 0.08 );

	RBConfLists rb_confs( get_rigid_body_confs( pose ) );
	//1.
	for ( Size i = 1; i <= rb_confs.size(); ++i ) {
		present_confs_seqpos.insert( rb_confs[i][1]->seqpos() );
		for ( Size j = rb_confs[i].size(); j <= num_total_rbconfs_; ++j ) {
			rb_confs[i].push_back( (rb_confs[i][1])->clone() );
		}
		diversifier.diversify_all_confs( pose, rb_confs[i] );
		set_rigid_body_confs_for_seqpos( rb_confs[i][1]->seqpos(), rb_confs[i], pose );
	}

	//2.
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue_type(i).is_ligand() && !( !include_metals_ && (pose.residue_type(i).natoms() <= 3) ) ) {
			if ( present_confs_seqpos.find( i ) != present_confs_seqpos.end() ) continue;

			tr << num_total_rbconfs_ << " random new rigid body conformations to be used in packing will be generated for residue " << pose.residue(i).name() << " " << i << "." << std::endl;
			utility::vector1< core::conformation::ResidueCOP > new_rb_confs;
			for ( Size j = 1; j <= num_total_rbconfs_; ++j ) {
				new_rb_confs.push_back( pose.residue( i ).clone() );
			}
			diversifier.diversify_all_confs( pose, new_rb_confs );
			set_rigid_body_confs_for_seqpos( i, new_rb_confs, pose );
		}
	}
}

ApplyRandomStoredRBConf::ApplyRandomStoredRBConf()
: parent( "ApplyRandomStoredRBConf" ) {}

ApplyRandomStoredRBConf::~ApplyRandomStoredRBConf()= default;


void
ApplyRandomStoredRBConf::apply(
	core::pose::Pose & pose )
{

	RBConfLists rb_confs( get_rigid_body_confs( pose ) );
	//pose.dump_pdb("random_apply_before.pdb");
	for ( core::Size i = 1; i <= rb_confs.size(); ++i ) {
		core::Size picked_conf( numeric::random::random_range( 1, rb_confs[i].size()) );
		tr << "Randomly appying conf " << picked_conf << " of " << rb_confs[i].size() <<" for respos " << rb_confs[i][picked_conf]->seqpos() << std::endl;
		core::conformation::ResidueOP exchange_res( rb_confs[i][picked_conf]->clone() );
		swap_coordinates_in_pose( pose, *exchange_res );
		rb_confs[i][picked_conf] = exchange_res;
	}
	//update the confs in the pose
	set_rigid_body_confs( rb_confs, pose );
	//pose.dump_pdb("random_apply_after.pdb");
}

std::string
ApplyRandomStoredRBConf::get_name() const {
	return "ApplyRandomStoredRBConf";
}


MinimizeStoredRBConfs::MinimizeStoredRBConfs( core::scoring::ScoreFunctionCOP sfxn )
: parent( "MinimizeStoredRBConfs" ), sfxn_(std::move(sfxn)), min_rms_(0.08) {}

MinimizeStoredRBConfs::~MinimizeStoredRBConfs()= default;


void
MinimizeStoredRBConfs::apply(
	core::pose::Pose & pose )
{
	RBConfLists rb_confs( get_rigid_body_confs( pose ) );
	for ( Size i = 1; i <= rb_confs.size(); ++i ) {
		if ( rb_confs[i].size() > 0 ) {
			rb_minimize_all_confs( pose, rb_confs[i] );
			set_rigid_body_confs_for_seqpos( rb_confs[i][1]->seqpos(), rb_confs[i], pose );
		}
	}
}

std::string
MinimizeStoredRBConfs::get_name() const {
	return "MinimizeStoredRBConfs";
}

void
MinimizeStoredRBConfs::rb_minimize_all_confs(
	core::pose::Pose const & pose,
	utility::vector1< core::conformation::ResidueCOP > & confs
) const
{
	//we need to make a copy of this pose that doesn't have the
	//rigid body confs
	//applycalls_++;
	core::pose::Pose mod_pose = pose;
	Size seqpos( confs[1]->seqpos() );
	Size natoms( confs[1]->natoms() );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	movemap->set_jump( pose.fold_tree().get_jump_that_builds_residue( seqpos ), true );
	protocols::simple_moves::MinMover minmover( movemap, sfxn_, "lbfgs_armijo_nonmonotone", 0.1, true );

	//I guess we need a task...
	core::pack::task::PackerTaskOP rtmin_task = core::pack::task::TaskFactory::create_packer_task( mod_pose );
	rtmin_task->initialize_from_command_line();
	for ( Size i = 1; i<= mod_pose.size(); ++i ) {
		if ( i == seqpos ) rtmin_task->nonconst_residue_task( i ).restrict_to_repacking();
		else rtmin_task->nonconst_residue_task( i ).prevent_repacking();
	}
	//if( applycalls_ == 1 ) mod_pose.dump_pdb("apply1_start.pdb");

	for ( Size conf = 1; conf <= confs.size(); ++conf ) {
		core::conformation::ResidueOP exchange_res( confs[conf]->clone() );
		swap_coordinates_in_pose( mod_pose, *exchange_res );
		//if( applycalls_ == 1 ) mod_pose.dump_pdb("before_rtmin_conf_"+utility::to_string( conf )+".pdb" );

		//note: ideally this would be rtmin, but i think it's too slow, especially
		//with a large number of rigid body confs/rotamers...
		//so let's do rotamer trials followed by a minimization
		//core::pack::RTMin rtmin( false, true );
		//rtmin.rtmin( mod_pose, *sfxn_, rtmin_task );
		core::pack::rotamer_trials( mod_pose, *sfxn_, rtmin_task );
		minmover.apply( mod_pose );
		//if( applycalls_ == 1 ) mod_pose.dump_pdb("after_rtmin_conf_"+utility::to_string( conf )+".pdb" );
		//now that the pose is rtmined, we need to get the new coordinates
		core::conformation::Residue const & newpos( mod_pose.residue( seqpos ) );
		for ( core::Size atm = 1; atm <= natoms; ++atm ) {
			exchange_res->set_xyz( atm, newpos.xyz( exchange_res->atom_name( atm ) ) );
		}
		confs[conf] = exchange_res;
	} //main rtmin loop over confs

	//debug
	/*
	if( applycalls_ == 1 ) {
	core::Size atcounter(0);
	std::ofstream out( "storeconfs_afterothermin.pdb" );
	for( Size conf = 1; conf <= confs.size(); ++conf ){
	out << "MODEL " << conf << "   " << std::endl;
	core::io::pdb::dump_pdb_residue( *confs[conf], atcounter, out );
	out << "ENDMDL " << std::endl;
	}
	out.close();
	}
	*/
	//debug over

	//now we need to diversify stuff, so not all rbconfs were minimized to the exact same local minima
	//diversify_rb_confs( pose, confs, min_rms_ );
	DiversifyStoredRBConfs diversifier( min_rms_ );
	diversifier.diversify_all_confs( pose, confs );

	//debug
	/*
	if( applycalls_ == 1 ) {
	core::Size atcounter(0);
	std::ofstream out( "storeconfs_afterdiversify.pdb" );
	for( Size conf = 1; conf <= confs.size(); ++conf ){
	out << "MODEL " << conf << "   " << std::endl;
	core::io::pdb::dump_pdb_residue( *confs[conf], atcounter, out );
	out << "ENDMDL " << std::endl;
	}
	out.close();
	}
	*/
	//debug over
}


DiversifyStoredRBConfs::DiversifyStoredRBConfs(
	Real min_rms )
: parent( "DiversifyStoredRBConfs" ), min_rms_(min_rms), max_trials_(10) {}

DiversifyStoredRBConfs::~DiversifyStoredRBConfs()= default;

void
DiversifyStoredRBConfs::apply(
	core::pose::Pose & pose )
{

	RBConfLists rb_confs( get_rigid_body_confs( pose ) );
	for ( Size i = 1; i <= rb_confs.size(); ++i ) {
		if ( rb_confs[i].size() > 0 ) {
			diversify_all_confs( pose, rb_confs[i] );
			set_rigid_body_confs_for_seqpos( rb_confs[i][1]->seqpos(), rb_confs[i], pose );
		}
	}
}


std::string
DiversifyStoredRBConfs::get_name() const {
	return "DiversifyStoredRBConfs";
}

void
DiversifyStoredRBConfs::diversify_all_confs(
	core::pose::Pose const & pose,
	utility::vector1< core::conformation::ResidueCOP > & confs
) const
{
	core::pose::Pose mod_pose = pose;
	Size natoms( confs[1]->natoms() );
	Size seqpos( confs[1]->seqpos() );
	protocols::rigid::RigidBodyPerturbMover simple_rigbod( mod_pose.fold_tree().get_jump_that_builds_residue( confs[1]->seqpos() ), numeric::conversions::degrees(0.05), min_rms_);

	for ( Size i = 1; i <= max_trials_; ++i ) {
		utility::vector1< Real > min_dist( closest_orient_atoms_msd( pose, confs ) );
		utility::vector1< Size > confs_to_change;
		//tr << "before trial " << i << ", confs have the following rms: " << std::endl;
		for ( Size j = 1; j <= min_dist.size(); ++j ) {
			//tr << j << " has rms of " << min_dist[j] << std::endl;
			if ( min_dist[j] <= min_rms_ ) confs_to_change.push_back( j );
		}
		if ( confs_to_change.size() == 0 ) break; // we're done

		for ( Size j = 1; j <= confs_to_change.size(); ++j ) {
			Size conf( confs_to_change[j] );
			core::conformation::ResidueOP exchange_res( confs[conf]->clone() );
			swap_coordinates_in_pose( mod_pose, *exchange_res );
			simple_rigbod.apply( mod_pose );
			core::conformation::Residue const & newpos( mod_pose.residue( seqpos ) );
			for ( core::Size atm = 1; atm <= natoms; ++atm ) {
				exchange_res->set_xyz( atm, newpos.xyz( exchange_res->atom_name( atm ) ) );
			}
			confs[conf] = exchange_res;
		}
	}
}

} // enzdes
} //protocols
