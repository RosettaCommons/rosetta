// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/RotamerSetOperations/AddGood2BPairEnergyRotamers.cc
/// @brief
/// @author Florian Richter, florian.richter.1@hu-berlin.de, june '14

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/AddGood2BPairEnergyRotamers.hh>

//Project headers

#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

// c++ headers
#include <fstream>
#include <map>

namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

static THREAD_LOCAL basic::Tracer tr( "protocols.toolbox.RotamerSetOperations.AddGood2BPairEnergyRotamers" );

AddGood2BPairEnergyRotamers::AddGood2BPairEnergyRotamers(
	Size seqpos,
	Size ex_level,
	utility::vector1< core::Size > const & target_seqpos,
	Real score_cutoff,
	bool drop_rots_above_cutoff
) : parent(),
	seqpos_(seqpos), ex_level_(ex_level),
	target_seqpos_(target_seqpos), score_cutoff_(score_cutoff),
	drop_rots_above_cutoff_(drop_rots_above_cutoff), debug_(false),
	disabled_(false)
{}

AddGood2BPairEnergyRotamers::AddGood2BPairEnergyRotamers( AddGood2BPairEnergyRotamers const & src )
: parent( src ),
	seqpos_(src.seqpos_), ex_level_(src.ex_level_),
	target_seqpos_(src.target_seqpos_), score_cutoff_(src.score_cutoff_),
	drop_rots_above_cutoff_(src.drop_rots_above_cutoff_), debug_(src.debug_), disabled_(false)
{}


AddGood2BPairEnergyRotamers::~AddGood2BPairEnergyRotamers(){}

core::pack::rotamer_set::RotamerSetOperationOP
AddGood2BPairEnergyRotamers::clone() const{
	return core::pack::rotamer_set::RotamerSetOperationOP( new AddGood2BPairEnergyRotamers( *this ) );
}


bool pair_first_sort_helper(
	std::pair< core::Real, core::conformation::ResidueCOP > const & pair1,
	std::pair< core::Real, core::conformation::ResidueCOP > const & pair2
)
{
	return (pair1.first < pair2.first );
}


/// @details
/// strategy: 1. note energy of existing rotamer towards target seqpos
///           2. build extra rotamers as specified by ex-flag
///           3. evaluate score towards target, keep rotamers that
///              have score below cutoff
///           4. little tricky: check against existing rotamer set
///              to make sure that no rotamers get added twice
///           5. in case of demand, output rotamers that were kept
///              to a file
void
AddGood2BPairEnergyRotamers::alter_rotamer_set(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::pack::task::PackerTask const & ptask,
	core::graph::GraphCOP packer_neighbor_graph,
	core::pack::rotamer_set::RotamerSet & rotamer_set
){

	if ( this->disabled_ ) return;
	core::Size this_pdbpos( pose.pdb_info()->number( seqpos_ ) );

	//1.
	Real existing_score(0.0);
	for ( Size i(1); i <= target_seqpos_.size(); ++i ) {
		existing_score += this->get_res_res_score( pose.residue( seqpos_), pose.residue( target_seqpos_[i]), pose, sfxn );
	}

	tr << "Beginning pdb position " <<  this_pdbpos << ", existing interactions is " << existing_score << "and existing rotamer set contains " << rotamer_set.num_rotamers() << " of " << rotamer_set.get_n_residue_groups() << " residue types." << std::endl;
	//2.
	core::pack::task::PackerTaskOP ptask_copy( ptask.clone() );
	core::pack::task::ExtraRotSample rot_explosion_sample( static_cast<core::pack::task::ExtraRotSample>(ex_level_) );
	ptask_copy->nonconst_residue_task( seqpos_).or_ex1( true ); ptask_copy->nonconst_residue_task( seqpos_).or_ex1_sample_level( rot_explosion_sample );
	ptask_copy->nonconst_residue_task( seqpos_).or_ex2( true ); ptask_copy->nonconst_residue_task( seqpos_).or_ex2_sample_level( rot_explosion_sample );
	ptask_copy->nonconst_residue_task( seqpos_).or_ex3( true ); ptask_copy->nonconst_residue_task( seqpos_).or_ex3_sample_level( rot_explosion_sample );
	ptask_copy->nonconst_residue_task( seqpos_).or_ex4( true ); ptask_copy->nonconst_residue_task( seqpos_).or_ex4_sample_level( rot_explosion_sample );

	this->disabled_ = true; //prevent infinite loop bc of packer task copy
	core::pack::rotamer_set::RotamerSetFactory rsf;
	core::pack::rotamer_set::RotamerSetOP extra_rotset( rsf.create_rotamer_set( pose.residue( seqpos_ ) ) ); //the residue isn't really used by the factory
	extra_rotset->set_resid( seqpos_ );
	extra_rotset->build_rotamers( pose, sfxn, *ptask_copy, packer_neighbor_graph );
	this->disabled_ = false;
	tr << "Finished building an extra rotset containing " << extra_rotset->num_rotamers() << " of " << extra_rotset->get_n_residue_groups() << " residue types." << std::endl;

	//3.
	//utlility::vector1< core::conformation::ResidueCOP > rots_to_keep;
	std::list< std::pair< Real, core::conformation::ResidueCOP > > rots_to_keep;
	Size number_passing_cutoff(0);
	for ( Size i(1); i <= extra_rotset->num_rotamers(); ++i ) {
		Real this_score(0.0);
		core::conformation::ResidueCOP candidate_rot( extra_rotset->rotamer( i ) );
		for ( Size j(1); j <= target_seqpos_.size(); ++j ) {
			this_score += this->get_res_res_score( *candidate_rot, pose.residue( target_seqpos_[j]), pose, sfxn );
		}
		if ( this_score <= score_cutoff_ ) {
			number_passing_cutoff++;
			//4. checking for rotamers existing twice
			if ( ! this->rotamer_set_contains_rotamer( rotamer_set, *candidate_rot ) ) {
				rots_to_keep.push_back( std::pair< Real, core::conformation::ResidueCOP > (this_score,candidate_rot) );
			}
		}
	}//loop over all extra rotamers

	for ( std::list< std::pair< Real, core::conformation::ResidueCOP > >::const_iterator list_it( rots_to_keep.begin() ), list_end( rots_to_keep.end() ); list_it != list_end; list_it++ ) {
		rotamer_set.add_rotamer( *(list_it->second) );
	}
	tr << "At pdb position " <<  this_pdbpos << ", " << rots_to_keep.size() << " additional rotamers were added. The standard set already contained " << number_passing_cutoff - rots_to_keep.size() << " rotamers passing the desired score cutoff." << std::endl;


	//5.
	if ( debug_ ) {

		std::string targetpdbpos_string;
		for ( Size i(1); i <= target_seqpos_.size(); ++i ) {
			targetpdbpos_string += utility::to_string(  pose.pdb_info()->number( target_seqpos_[i] ) ) + utility::to_string(  pose.pdb_info()->chain( target_seqpos_[i] ) );
		}

		rots_to_keep.sort( pair_first_sort_helper );
		std::map< std::string, std::list< std::pair< Real, core::conformation::ResidueCOP > > > restype_map;
		for ( std::list< std::pair< Real, core::conformation::ResidueCOP > >::const_iterator list_it( rots_to_keep.begin() ), list_end( rots_to_keep.end() ); list_it != list_end; list_it++ ) {

			std::string cur_name3( list_it->second->name3() );
			std::map< std::string, std::list< std::pair< Real, core::conformation::ResidueCOP > > >::iterator map_it( restype_map.find( cur_name3 ) );
			if ( map_it == restype_map.end() ) {
				restype_map.insert( std::pair< std::string, std::list< std::pair< Real, core::conformation::ResidueCOP > > >( cur_name3, std::list< std::pair< Real, core::conformation::ResidueCOP > > () ) );
				map_it = restype_map.find( cur_name3 );
			}
			map_it->second.push_back( *list_it );
		}

		for ( std::map< std::string, std::list< std::pair< Real, core::conformation::ResidueCOP > > >::const_iterator map_it( restype_map.begin() ), map_end( restype_map.end() ); map_it != map_end; ++map_it ) {
			std::string filename( "ExpRot_"+targetpdbpos_string+"_"+utility::to_string( this_pdbpos)+map_it->first+".pdb");
			std::ofstream file_out(filename.c_str() );
			Size atcounter(0);
			Size modelcounter(0);
			std::list< std::pair< Real, core::conformation::ResidueCOP > >::const_iterator list_it( map_it->second.begin() ), list_end( map_it->second.end() );
			for ( ; list_it != list_end; ++list_it ) {
				modelcounter++;
				core::conformation::Residue curres( *(list_it->second) );
				curres.seqpos( this_pdbpos );
				file_out << "MODEL " << modelcounter << " " << list_it->first << "\n";
				core::io::pdb::dump_pdb_residue( curres, atcounter, file_out );
				file_out << "ENDMDL\n";
			}
			file_out.close();
		}
	} //if debug_

} //alter rotamer set


/// @details two things need to be done
/// 1. easy: short range interactions
/// 2. more cumbersome: long range stuff,
/// need to duplicate some code from the scorefunction
AddGood2BPairEnergyRotamers::Real
AddGood2BPairEnergyRotamers::get_res_res_score(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn
) const
{
	using namespace core::scoring;

	EnergyMap tbemap;
	sfxn.eval_ci_2b( rsd1, rsd2, pose, tbemap );
	sfxn.eval_cd_2b( rsd1, rsd2, pose, tbemap );

	for ( ScoreFunction::CI_LR_2B_Methods::const_iterator iter = sfxn.ci_lr_2b_methods_begin(),
			iter_end = sfxn.ci_lr_2b_methods_end(); iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy( rsd1, rsd2, pose, sfxn, tbemap );
	}
	for ( ScoreFunction::CD_LR_2B_Methods::const_iterator iter = sfxn.cd_lr_2b_methods_begin(),
			iter_end = sfxn.cd_lr_2b_methods_end(); iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy( rsd1, rsd2, pose, sfxn, tbemap );
	}

	return tbemap.dot( sfxn.weights() );
}


bool
AddGood2BPairEnergyRotamers::rotamer_set_contains_rotamer(
	core::pack::rotamer_set::RotamerSet const & rotamer_set,
	core::conformation::Residue const & candidate_rot ) const
{
	Size num_resgroups( rotamer_set.get_n_residue_groups() );

	for ( Size resgroup(1); resgroup <= num_resgroups; ++resgroup ) {
		Size cur_rot_id( rotamer_set.get_residue_group_begin( resgroup ) );

		core::conformation::ResidueCOP cur_rot( rotamer_set.rotamer( cur_rot_id ) );
		if ( candidate_rot.name3() != cur_rot->name3() ) continue;

		Size resgroup_end( resgroup == num_resgroups ? rotamer_set.num_rotamers() : rotamer_set.get_residue_group_begin( resgroup + 1 ) - 1 );

		for ( ; cur_rot_id <= resgroup_end; ++cur_rot_id ) {
			cur_rot = rotamer_set.rotamer( cur_rot_id );
			bool all_chi_identical(true);
			for ( Size chi_id(1); chi_id <= candidate_rot.nchi(); ++chi_id ) {
				if ( candidate_rot.chi()[ chi_id ] != cur_rot->chi()[ chi_id ] ) {
					all_chi_identical = false;
					break;
				}
			}
			if ( all_chi_identical ) return true;
		} //loop over rotamers in resgroup
	} //loop over resgroups
	return false;
}


} //namespace rotamer_set_operations
} //namespace toolbox
} //namespace protocols
