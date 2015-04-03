// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/TopologyBrokerMover.cc
/// @brief
/// @author Lei Shi (shilei@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/TopologyBrokerMover.hh>
#include <protocols/protein_interface_design/movers/TopologyBrokerMoverCreator.hh>

// Package headers
#include <protocols/rosetta_scripts/util.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>

//#include <protocols/topology_broker/TopologyBroker.hh>
//#include <protocols/topology_broker/util.hh>
#include <protocols/abinitio/BrokerMain.hh>
#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>

#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <core/pose/PDBInfo.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.TopologyBrokerMover" );

std::string
TopologyBrokerMoverCreator::keyname() const
{
	return TopologyBrokerMoverCreator::mover_name();
}

protocols::moves::MoverOP
TopologyBrokerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TopologyBrokerMover );
}

std::string
TopologyBrokerMoverCreator::mover_name()
{
	return "TopologyBrokerMover";
}

TopologyBrokerMover::TopologyBrokerMover() : simple_moves::DesignRepackMover( TopologyBrokerMoverCreator::mover_name() ) {}
TopologyBrokerMover::~TopologyBrokerMover() {}

void
TopologyBrokerMover::apply( pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace core::pack::task::operation;

	//separate the two chains
  core::pose::Pose pose1 = pose;
  core::pose::Pose pose2 = pose;
  Size chain1primet(pose.conformation().chain_begin( 1 ));
  Size chain1end(pose.conformation().chain_end( 1 ) );
  Size chain2primet(pose.conformation().chain_begin( 2 ));
  Size chain2end(pose.conformation().chain_end( 2 ) );
  pose1.conformation().delete_residue_range_slow( chain2primet, chain2end);
  pose2.conformation().delete_residue_range_slow( chain1primet, chain1end);
  core::pose::Pose ref_pose = pose2;

	protocols::protein_interface_design::movers::SaveAndRetrieveSidechainsOP srsc( new protocols::protein_interface_design::movers::SaveAndRetrieveSidechains(ref_pose, true, false, 0) );

  protocols::abinitio::AbrelaxApplication::register_options();
  protocols::abinitio::register_options_broker();
  protocols::abinitio::AbrelaxMoverOP abrelax( new protocols::abinitio::AbrelaxMover() );
	//abrelax->set_b_return_unrelaxed_fullatom(true);
	abrelax->apply( pose2 );
	//switch to fa output, relax is set in commond line
  //if ( !pose.is_fullatom() ) {
  //  core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
  //}
	srsc->apply(pose2);

	// superimpose
	core::id::AtomID_Map< id::AtomID > atom_map; 
	core::pose::initialize_atomid_map_AtomID ( atom_map, pose2, id::BOGUS_ATOM_ID ); 
	for ( core::Size i=start_; i<=end_; ++i ) { 
		core::id::AtomID const id1( pose2.residue(i-chain1end).atom_index("CA"), i-chain1end ); 
		core::id::AtomID const id2( ref_pose.residue(i-chain1end).atom_index("CA"), i-chain1end); 
		atom_map[ id1 ] = id2; 
	} 
  core::scoring::superimpose_pose( pose2, ref_pose, atom_map );
	Real rms=core::scoring::CA_rmsd(pose2,ref_pose,start_-chain1end,end_-chain1end);
	if (rms > 0.001) {
			utility_exit_with_message( "Something is wrong, the rigid part is not rigid, check your rigid_start, rigid_end and setup.tpb for broker !" );
	}
	//runtime_assert( rms < 0.001 );
	
	//coordinate constrained relax
	
	// append target
  pose1.append_pose_by_jump(pose2,chain1end);
	pose=pose1;
	//copy pdbinfo
  core::pose::PDBInfoOP pdbinfo( new core::pose::PDBInfo( pose) );
  pose.pdb_info(pdbinfo);
	
}

std::string
TopologyBrokerMover::get_name() const {
	return TopologyBrokerMoverCreator::mover_name();
}

void
TopologyBrokerMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	TR<<"Setup TopologyBrokerMover Mover " << std::endl;
	if (pose.conformation().num_chains()!=2) {
		utility_exit_with_message( "Expect your input pose -in:file:s to contain exactly two chains, A is target, B is design to be folded" );
	}

	align_ = tag->getOption< bool >( "realign", false );
	if ( align_ ) {
		if ( tag->hasOption("rigid_start") && tag->hasOption("rigid_end") ) {
	  	start_ = tag->getOption< core::Size >( "rigid_start", 1 );
	  	end_ = tag->getOption< core::Size >( "rigid_end", 1 );
			if (start_<pose.conformation().chain_begin( 2 ) || end_ < pose.conformation().chain_begin( 2 )
			 || start_>pose.conformation().chain_end( 2 ) || end_ > pose.conformation().chain_end( 2 ) 
			 || start_ >= end_ ) {
			TR << "rigid_start " << start_ << " and rigid_end " << end_ << " is not in range " << pose.conformation().chain_begin( 2 ) <<"-" << pose.conformation().chain_end( 2 ) << std::endl;
			utility_exit_with_message( "If you want to realign, please provide motif structure and target structure!" );
			}
		} else {
			utility_exit_with_message( "If you want to realign, please provide motif rigid_start and rigid_end in pose numbering!" );
		}
	} else  { 
		  utility_exit_with_message( "Expect to use realign option for futher design purpose" );
	}
}

protocols::moves::MoverOP
TopologyBrokerMover::clone() const {
    return( protocols::moves::MoverOP( new TopologyBrokerMover( *this ) ));
}

} //movers
} //protein_interface_design
} //protocols
