// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Bjorn Wallner
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/jumping/MembraneJump.hh>
#include <protocols/jumping/PairingLibrary.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <core/scoring/dssp/PairingsList.fwd.hh>

// Project Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedStubID.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/kinematics/FoldTree.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <cstdlib>
#include <string>


namespace protocols {
namespace jumping {

/// @details Auto-generated virtual destructor
MembraneJump::~MembraneJump() {}
using namespace core;
using namespace fragment;
static THREAD_LOCAL basic::Tracer tr( "protocols.jumping.MembraneJump" );

//default constructor
MembraneJump::MembraneJump()
{
	template_size_=0;
	pairings_size_=0;
}

// init given a template file and a pairings file
void
MembraneJump::init(std::string const& template_file,std::string const& pairings_file) {
	templates_.read_from_file_no_filters(template_file);
	read_pairing_list( pairings_file, pairings_);
	template_size_=templates_.size();
	pairings_size_=pairings_.size();

}

//this function will setup a fold tree to be used consisting of njumps using the info in templates_ and pairings_
void
MembraneJump::setup_fold_tree(core::pose::Pose & pose, core::Size njumps)
{
	using namespace ObjexxFCL;
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;

	if ( pairings_.size()==0 ) {
		return;
	}
	tr << "setting up fold_tree with " << njumps << " jump(s)\n";
	Size nres=pose.total_residue();
	core::kinematics::FoldTree f(nres);
	Size tries(0);
	core::scoring::MembraneTopology const & topology(*( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) )));
	FArray1D_int tmh(pose.total_residue());
	FArray1D_int tmh2(pose.total_residue(),0);
	Size total_tmhelix(topology.tmhelix());
	FArray1D_bool tmh_involved_in_jump(total_tmhelix,false);

	for ( Size j = 1; j <= pose.total_residue(); ++j ) {
		//bw change definition of membrane region to include jumps to non-tmh.
		if ( j<=topology.span_end(1) ) { //membrane_helix(1,2))
			tmh(j)=1;
		} else if ( j>topology.span_end(total_tmhelix) ) {
			tmh(j)=total_tmhelix;
		} else {
			for ( Size reg = 2; reg <= total_tmhelix; ++reg ) {
				if ( j>topology.span_end(reg-1) && j<=topology.span_end(reg) ) { //membrane_helix( reg-1, 2 ) && j<=membrane_helix(reg,2))
					tmh(j)=reg;
				}
			}
		}
	}
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		for ( Size j = 1; j <= total_tmhelix; ++j ) {
			if ( i>topology.span_begin(j) && i < topology.span_end(j) ) {
				tmh2(i)=j;
			}
		}
	}

	while ( selected_pairings_.size()<njumps && tries < 10 ) {
		Size index=static_cast< int >(numeric::random::rg().uniform()*pairings_.size()+1);
		std::cout << "Tries : " << tries << " " << index << ' ' << pairings_[index].Pos1()  << ' ' << tmh(pairings_[index].Pos1()) << ' ' << pairings_[index].Pos2() << ' ' << tmh(pairings_[index].Pos2()) <<std::endl;
		bool check_compatible=true;

		{
			if ( tmh_involved_in_jump(tmh(pairings_[index].Pos1())) ||
					tmh_involved_in_jump(tmh(pairings_[index].Pos2())) ) {
				check_compatible=false;
			}
		}
		for ( Size j = 1; j <= selected_pairings_.size(); ++j ) {
			if ( selected_pairings_[j].Pos1() == pairings_[index].Pos1() &&
					selected_pairings_[j].Pos2() == pairings_[index].Pos2() ) { // already in a jump
				check_compatible=false;
			}
		}

		if ( check_compatible ) {
			selected_pairings_.push_back(pairings_[index]);
			tmh_involved_in_jump(tmh(pairings_[index].Pos1()))=true;
			tmh_involved_in_jump(tmh(pairings_[index].Pos2()))=true;
		}
		++tries;
	}
	if ( selected_pairings_.size()<njumps ) {
		std::cout << "WARNING: Only picked " << selected_pairings_.size() << " given number was " << njumps << " only allow one jump between any two TMHs " << std::endl;
	}
	FArray2D_int jumps(2,selected_pairings_.size());
	for ( Size i=1; i<=selected_pairings_.size(); ++i ) {
		jumps(1,i)=selected_pairings_[i].Pos1();
		jumps(2,i)=selected_pairings_[i].Pos2();
	}
	FArray1D_float cut_bias(nres,0.0);
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( tmh2(i)==0 ) {
			cut_bias(i)=1;
		}
	}

	int num_jumps_in=selected_pairings_.size();
	f.random_tree_from_jump_points(nres,num_jumps_in,jumps,cut_bias);
	f.put_jump_stubs_intra_residue();

	std::cout <<  f;
	pose.fold_tree(f);
}
void
MembraneJump::rt_templates(core::pose::Pose& pose)
{

	for ( Size i=1; i<=selected_pairings_.size(); ++i ) {
		Size p1=selected_pairings_[i].Pos1();
		Size p2=selected_pairings_[i].Pos2();
		core::kinematics::FoldTree f(pose.fold_tree());
		core::kinematics::RT rt(templates_.get_random_tmh_jump(selected_pairings_[i].Orientation(),p1,p2));
		id::StubID up_stub(   core::pose::named_stub_id_to_stub_id( core::id::NamedStubID( "CA","N","CA","C", p1 ), pose ) );
		id::StubID down_stub( core::pose::named_stub_id_to_stub_id( core::id::NamedStubID( "CA","N","CA","C", p2 ), pose ) );
		pose.conformation().set_stub_transform( up_stub, down_stub, rt );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, f.cutpoint(i) );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, f.cutpoint(i)+1 );
	}
}

} //jumping
} //protocols
