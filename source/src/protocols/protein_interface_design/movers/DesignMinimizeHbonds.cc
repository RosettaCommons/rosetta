// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/DesignMinimizeHbonds.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/DesignMinimizeHbonds.hh>
#include <protocols/protein_interface_design/movers/DesignMinimizeHbondsCreator.hh>


#include <core/pack/pack_rotamers.hh>
#include <protocols/scoring/Interface.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/protein_interface_design/design_utils.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.DesignMinimizeHbonds" );

std::string
DesignMinimizeHbondsCreator::keyname() const
{
	return DesignMinimizeHbondsCreator::mover_name();
}

protocols::moves::MoverOP
DesignMinimizeHbondsCreator::create_mover() const {
	return new DesignMinimizeHbonds;
}

std::string
DesignMinimizeHbondsCreator::mover_name()
{
	return "DesignMinimizeHbonds";
}

DesignMinimizeHbonds::DesignMinimizeHbonds() :
	simple_moves::DesignRepackMover( DesignMinimizeHbondsCreator::mover_name() )
{
	min_rb_set_ = min_bb_set_ = min_sc_set_ = false;
	optimize_foldtree_ = true;
	automatic_repacking_definition_ = true;
}

DesignMinimizeHbonds::DesignMinimizeHbonds(
	ScoreFunctionCOP scorefxn_repack, ScoreFunctionCOP scorefxn_minimize,
	utility::vector1< core::Size > const target_residues,
	bool const donors,
	bool const acceptors,
	bool const bb_hbond,
	bool const sc_hbond,
	core::Real const hbond_energy_threshold,
	core::Real interface_distance_cutoff/*=8.0*/,
	bool const repack_partner1/*=true*/,
	bool const repack_partner2/*=false*/,
	bool const repack_non_ala/* = true*/
) :
	simple_moves::DesignRepackMover( DesignMinimizeHbondsCreator::mover_name() )
{
	scorefxn_repack_ = scorefxn_repack->clone();
	scorefxn_minimize_ = scorefxn_minimize->clone();
	target_residues_ = target_residues;
	repack_partner1_ = repack_partner1;
	repack_partner2_ = repack_partner2;
	donors_ = donors; acceptors_ = acceptors;
	bb_hbond_ = bb_hbond;
	sc_hbond_ = sc_hbond;
	hbond_energy_threshold_ = hbond_energy_threshold;
	interface_distance_cutoff_ = interface_distance_cutoff;
	repack_non_ala_ = repack_non_ala;
	// only requires donors_ || acceptors_ if doing sc hbonding
	runtime_assert( ( donors_ || acceptors_ ) || ( bb_hbond_ && !sc_hbond_ ) );
	runtime_assert( bb_hbond || sc_hbond );
	runtime_assert( hbond_energy_threshold_ <= 0 );
	runtime_assert (interface_distance_cutoff_ >= 0 );
}

DesignMinimizeHbonds::DesignMinimizeHbonds(
	ScoreFunctionOP scorefxn_repack,
	ScoreFunctionOP scorefxn_minimize,
	core::Size const target_residue,
	bool const donors,
	bool const acceptors,
	bool const bb_hbond,
	bool const sc_hbond,
	core::Real const hbond_energy_threshold,
	core::Real interface_distance_cutoff /*=8.0*/,
	bool const repack_partner1/*=true*/,
	bool const repack_partner2/*=false*/,
	bool const repack_non_ala/*=true*/
) :
	simple_moves::DesignRepackMover( DesignMinimizeHbondsCreator::mover_name() )
{
	  scorefxn_repack_ = scorefxn_repack;
	  scorefxn_minimize_ = scorefxn_minimize;
	  target_residues_.push_back( target_residue );
	  bb_hbond_ = bb_hbond;
	  sc_hbond_ = sc_hbond;
	  hbond_energy_threshold_ = hbond_energy_threshold;
	  repack_partner1_ = repack_partner1;
	  repack_partner2_ = repack_partner2;
	  interface_distance_cutoff_ = interface_distance_cutoff;
	  donors_ = donors; acceptors_ = acceptors;
		repack_non_ala_ = repack_non_ala;
}

DesignMinimizeHbonds::~DesignMinimizeHbonds() {}

protocols::moves::MoverOP
DesignMinimizeHbonds::clone() const {
	return( protocols::moves::MoverOP( new DesignMinimizeHbonds( *this ) ) );
}


/// @details Residues within a 10.0 Ang sphere around the set of target residues are designed with acceptors or donors
/// (or both) and the design is minimized. All residues that were changed but are not hbonded to the target residues
/// are reverted.  If only backbone hbonding is desired all residues except PRO, CYS, and GLY are used in design.
void
DesignMinimizeHbonds::apply( pose::Pose & pose )
{
	allowed_aas_.assign( chemical::num_canonical_aas, false );

	if ( bb_hbond_ && !sc_hbond_ ){
		allowed_aas_[ chemical::aa_ala ] = true;
		allowed_aas_[ chemical::aa_arg ] = true;
		allowed_aas_[ chemical::aa_asn ] = true;
		allowed_aas_[ chemical::aa_asp ] = true;
		allowed_aas_[ chemical::aa_gln ] = true;
		allowed_aas_[ chemical::aa_glu ] = true;
		allowed_aas_[ chemical::aa_his ] = true;
		allowed_aas_[ chemical::aa_ile ] = true;
		allowed_aas_[ chemical::aa_leu ] = true;
		allowed_aas_[ chemical::aa_lys ] = true;
		allowed_aas_[ chemical::aa_met ] = true;
		allowed_aas_[ chemical::aa_phe ] = true;
		allowed_aas_[ chemical::aa_ser ] = true;
		allowed_aas_[ chemical::aa_thr ] = true;
		allowed_aas_[ chemical::aa_trp ] = true;
		allowed_aas_[ chemical::aa_tyr ] = true;
		allowed_aas_[ chemical::aa_val ] = true;
	}
	else{
		if( donors_ ) {
			allowed_aas_[ chemical::aa_lys ] = true;
			allowed_aas_[ chemical::aa_asn ] = true;
			allowed_aas_[ chemical::aa_gln ] = true;
			allowed_aas_[ chemical::aa_ser ] = true;
			allowed_aas_[ chemical::aa_thr ] = true;
			allowed_aas_[ chemical::aa_trp ] = true;
			allowed_aas_[ chemical::aa_tyr ] = true;
			allowed_aas_[ chemical::aa_his ] = true;
		}
		if( acceptors_ ) {
			allowed_aas_[ chemical::aa_asp ] = true;
			allowed_aas_[ chemical::aa_glu ] = true;
			allowed_aas_[ chemical::aa_asn ] = true;
			allowed_aas_[ chemical::aa_gln ] = true;
			allowed_aas_[ chemical::aa_ser ] = true;
			allowed_aas_[ chemical::aa_thr ] = true;
		}
	}
	core::Size const rb_jump( 1 );

	runtime_assert( repack_partner1_ || repack_partner2_ );

	core::pose::Pose const saved_pose( pose );

	using ObjexxFCL::FArray1D_bool;
 	FArray1D_bool partner1( pose.total_residue() );
	pose.fold_tree().partition_by_jump( rb_jump, partner1 ); // partner1 is true for all residues in partner1; false o/w

  protocols::scoring::Interface interface_obj(rb_jump);
  interface_obj.distance( interface_distance_cutoff_ ); // to encourage longish residues
  pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
  interface_obj.calculate( pose );

	setup_packer_and_movemap( pose );
	// potential hbond partners will later on be reverted if they do not form hbonds
	std::set< core::Size > potential_hbond_partners;
	for( core::Size i = 1; i <= pose.total_residue(); ++i ){
		if( !pose.residue(i).is_protein() ) continue;
		core::Size const restype( pose.residue(i).aa() );
		if( (interface_obj.is_interface( i ) && // in interface
				(partner1( i ) && repack_partner1_ )) || ((!partner1(i) && repack_partner2_) && //designable
				( !( !repack_non_ala_ && (restype != chemical::aa_ala) ) || (restype == chemical::aa_pro) || (restype == chemical::aa_gly) || pose.residue(i).type().name() == "CYD" ))) { // design-allowed residues
	        core::conformation::Residue const resi( pose.residue( i ) );
	        for( utility::vector1< Size >::const_iterator target_it = target_residues_.begin();
	           target_it!=target_residues_.end(); ++target_it ) {
	          core::conformation::Residue const res_target( pose.residue( *target_it ) );

	          Real const distance( resi.xyz( resi.nbr_atom() ).distance( res_target.xyz( res_target.nbr_atom() ) ) );
	          Real const distance_cutoff( interface_distance_cutoff_ );
	          if ( distance <= distance_cutoff && automatic_repacking_definition_ )
	            potential_hbond_partners.insert( i );
					}
				}
		}

	pack::pack_rotamers( pose, *scorefxn_repack_, task_ );
	MinimizeInterface( pose, scorefxn_minimize_, curr_min_bb_, curr_min_sc_, curr_min_rb_, optimize_foldtree_, target_residues_ );
	pose.update_residue_neighbors();

	{ // replace any positions that were mutated but did not hbond with their previous identities
		std::set< core::Size > hbonded_residues;
		for( utility::vector1< core::Size >::const_iterator target_it = target_residues_.begin();
			 target_it!=target_residues_.end(); ++target_it ) {
			std::list< core::Size > new_list( hbonded( pose, *target_it, potential_hbond_partners, bb_hbond_, sc_hbond_,
													   hbond_energy_threshold_ ));
			hbonded_residues.insert( new_list.begin(), new_list.end() );
		}

		pack::task::PackerTaskOP to_Ala_task( pack::task::TaskFactory::create_packer_task( pose ));
		for( Size i=1; i<=pose.total_residue(); ++i ) {
			if( potential_hbond_partners.find( i ) == potential_hbond_partners.end() ) {
				to_Ala_task->nonconst_residue_task(i).prevent_repacking();
				continue;
			}
			if( hbonded_residues.find( i ) == hbonded_residues.end() ) { // revert
				TR<<"reverting "<< i <<'\n';
				utility::vector1< bool > revert_type( chemical::num_canonical_aas, false );
				revert_type[ saved_pose.residue( i ).aa() ] = true;
				to_Ala_task->nonconst_residue_task( i ).restrict_absent_canonical_aas( revert_type );
			}
			else
				to_Ala_task->nonconst_residue_task(i).prevent_repacking();
		}
		pack::pack_rotamers( pose, *scorefxn_repack_, to_Ala_task );
	} // end of replace non-hbonded residues scope
	(*scorefxn_minimize_)( pose );
	/// Now handled automatically.  scorefxn_minimize_->accumulate_residue_total_energies( pose );
	pose.update_residue_neighbors();
	TR.flush();
}

std::string
DesignMinimizeHbonds::get_name() const {
	return DesignMinimizeHbondsCreator::mover_name();
}

void
DesignMinimizeHbonds::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const & filters, Movers_map const & movers, core::pose::Pose const & pose )
{
	core::Real const hbond_weight( tag->getOption<core::Real>( "hbond_weight", 3.0 ) );
	TR<<"DesignMinimizeHbonds with the following parameters: "<<std::endl;
	donors_ = tag->getOption<bool>( "donors" );
	acceptors_ = tag->getOption<bool>( "acceptors" );
	bb_hbond_ = tag->getOption<bool>( "bb_hbond", 0 );
	sc_hbond_ = tag->getOption<bool>( "sc_hbond", 1 ) ;
	hbond_energy_threshold_ = tag->getOption<core::Real>( "hbond_energy", -0.5 );
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_cutoff_distance", 8.0 );

	simple_moves::DesignRepackMover::parse_my_tag( tag, data, filters, movers, pose );
	using namespace core::scoring;

  // change the weights on the hbonding terms
	scorefxn_repack_->set_weight( hbond_lr_bb, hbond_weight );
	scorefxn_repack_->set_weight( hbond_sr_bb, hbond_weight );
	scorefxn_repack_->set_weight( hbond_bb_sc, hbond_weight );
	scorefxn_repack_->set_weight( hbond_sc,   hbond_weight );

	scorefxn_minimize_->set_weight( hbond_lr_bb, hbond_weight );
	scorefxn_minimize_->set_weight( hbond_sr_bb, hbond_weight );
	scorefxn_minimize_->set_weight( hbond_bb_sc, hbond_weight );
	scorefxn_minimize_->set_weight( hbond_sc,   hbond_weight );

	runtime_assert( ( donors_ || acceptors_ ) || ( bb_hbond_ && !sc_hbond_ ) );
	runtime_assert( bb_hbond_ || sc_hbond_ );
	runtime_assert( hbond_energy_threshold_ <= 0 );
	runtime_assert(interface_distance_cutoff_ >= 0 );
	runtime_assert( target_residues_.size() );

	if( target_residues_.size() == 0 )
		TR<<"WARNING WARNING: no target residue defined for hbond design minimize"<<std::endl;
}

} //movers
} //protein_interface_design
} //protocols
