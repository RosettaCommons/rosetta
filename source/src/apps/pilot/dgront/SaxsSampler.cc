// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 * SaxsSampler.cc
 *
 *  Created on: Jan 15, 2009
 *      Author: dgront
 */
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <protocols/Protocol.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <core/chemical/ChemicalManager.hh>


#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>

#include <utility/vector1.hh>

#include <cstdlib>

//Auto Headers
#include <core/pose/util.hh>


OPT_1GRP_KEY( IntegerVector, saxs, declare_domain )
OPT_1GRP_KEY( IntegerVector, saxs, declare_loop )
OPT_1GRP_KEY( Integer, saxs, declare_cutpoint )
OPT_1GRP_KEY( Boolean, saxs, randomize )
OPT_1GRP_KEY( Boolean, saxs, start_150 )

void register_options() {
	NEW_OPT( saxs::declare_cutpoint, "declare a cutpoint",0);
	NEW_OPT( saxs::declare_domain, "declare a protein fragment (possibly a domain) that will be treated as a rigid body",0);
	NEW_OPT( saxs::declare_loop, "declare a loose protein fragment (e.g. a loop) where fragment may be inserted",0);
	NEW_OPT( saxs::randomize, "Assign random conformations to loops. This option requires saxs::declare_loop to be used (user must say where the loops are)",0);
	NEW_OPT( saxs::start_150, "Assign -150,150 conformations to loops",0);
	OPT(out::file::silent);
	OPT(out::sf);
}

// #include "AbrelaxWithCutpoint.cc"

static THREAD_LOCAL basic::Tracer trSaxs( "SaxsSampler" );

int main(int argc, char * argv[]) {
    try {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using std::string;
	using utility::vector1;

//	protocols::abinitio::AbrelaxWithCutpoint::register_options();
	register_options();
	devel::init(argc, argv);

	//-------------- SET UP ABRELAX APPLICATION -------------
	protocols::abinitio::AbrelaxWithCutpoint *abrelax =
			new protocols::abinitio::AbrelaxWithCutpoint();

	// setup pose and abinitio
	protocols::Protocol* prot_ptr;
	pose::Pose init_pose;

	abrelax->setup();
	abrelax->setup_fold(init_pose, (protocols::ProtocolOP&) prot_ptr);
	init_pose.dump_pdb("init_pose.pdb");

	//-------------- Rotable parts for random restart
	utility::vector1<Size> rot_points;

	//-------------- SET UP DOMAINS -------------
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	movemap->set_bb(true);
	if (option[saxs::declare_domain].user()) {
		trSaxs.Info << "The fixed residues are:" << std::endl;
		utility::vector1<int> const& domainBoundaries(
				option[saxs::declare_domain]());
		for (Size i = 1; i <= domainBoundaries.size() / 2; i++) {
			for (Size iRes = domainBoundaries[i * 2 - 1]; iRes
					<= (Size) domainBoundaries[i * 2]; iRes++) {
				movemap->set_bb(iRes, false);
				trSaxs.Info << init_pose.residue(iRes).seqpos() << " "
						<< init_pose.residue(iRes).name() << std::endl;
			}
			for (Size iRes = domainBoundaries[i * 2 - 1] + 1; iRes
					< (Size) domainBoundaries[i * 2]; iRes++) {
				rot_points.push_back(iRes);
			}
		}
	}

	if (option[saxs::declare_loop].user()) {
		movemap->set_bb(false);
		trSaxs.Info << "The moveble residues are:" << std::endl;
		utility::vector1<int> const& domainBoundaries(
				option[saxs::declare_loop]());

		for (Size i = 1; i <= domainBoundaries.size() / 2; i++) {
			for (Size iRes = domainBoundaries[i * 2 - 1]; iRes
					<= (Size) domainBoundaries[i * 2]; iRes++) {
				movemap->set_bb(iRes, true);
				trSaxs.Info << init_pose.residue(iRes).seqpos() << " "
						<< init_pose.residue(iRes).name() << std::endl;
				//-------------- Random restart
				if (option[saxs::randomize].user()) {
					trSaxs.Debug << "Randomizing loops" << std::endl;
					Real phi = -((Real) (rand() % 50000)) - 120.0;
					Real psi = ((Real) (rand() % 50000)) + 120.0;
					init_pose.set_phi(iRes, phi);
					init_pose.set_psi(iRes, psi);
				}
				//-------------- Restart at 150
				if (option[saxs::start_150].user()) {
					trSaxs.Debug << "Expanding loops" << std::endl;
					Real phi = -150.0;
					Real psi = 150.0;
					init_pose.set_phi(iRes, phi);
					init_pose.set_psi(iRes, psi);
				}
			}
		}
	}

	//--------------- set jump ---------------
/*	if (option[saxs::declare_jump].user()) {
		Size res1, res2, cutpoint;
		utility::vector1<int> const& jump(option[saxs::declare_jump]());
		res1 = jump[1];
		res2 = jump[2];
		cutpoint = jump[3];
		tr.Info << "The jump positions are: " << res1 << " " << res2 << " "
				<< cutpoint << std::endl;
		chemical::ResidueType const& rt1(init_pose.residue_type(res1));
		chemical::ResidueType const& rt2(init_pose.residue_type(res2));
		id::AtomID a1(rt1.atom_index("N"), res1);
		id::AtomID a2(rt1.atom_index("CA"), res1);
		id::AtomID a3(rt1.atom_index("C"), res1);
		id::StubID down_stub(a1, a2, a3);

		id::AtomID b1(rt2.atom_index("N"), res2);
		id::AtomID b2(rt2.atom_index("CA"), res2);
		id::AtomID b3(rt2.atom_index("C"), res2);
		id::StubID up_stub(b1, b2, b3);
		kinematics::Stub up =
				init_pose.conformation().atom_tree().stub_from_id(up_stub);
		kinematics::Stub down =
				init_pose.conformation().atom_tree().stub_from_id(down_stub);
		kinematics::RT rt(up, down);
		kinematics::Jump nat_jump(rt);
		init_pose.set_jump(1,nat_jump);
		core::pose::add_variant_type_to_pose_residue(init_pose,
				chemical::CUTPOINT_LOWER, cutpoint);
		core::pose::add_variant_type_to_pose_residue(init_pose,
				chemical::CUTPOINT_UPPER, cutpoint + 1);
	}*/

	//-------------- RUN IT -------------
	trSaxs.Debug << "Runnig the sampler" << std::endl;
	protocols::abinitio::ClassicAbinitio* abinitio =
			static_cast<protocols::abinitio::ClassicAbinitio*> (prot_ptr);
	abinitio->bSkipStage3_ = true;
	abinitio->bSkipStage4_ = true;
	// plug in the new move map
	abinitio->set_movemap(movemap);
	// RUN
	abrelax->fold(init_pose, prot_ptr);
    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
       return 0;
}
