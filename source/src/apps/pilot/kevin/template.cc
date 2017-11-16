// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kevin/sandbox.cc
/// @brief
/// @details
/// @author Kevin Houlihan

#include <ctime>
#include <iostream>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
//#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>

//#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/scoring/sasa.hh>
#include <numeric/conversions.hh> //degrees-radians

//tracers
static basic::Tracer TR( "kevin.sandbox" );


//using namespace basic::options;
//using namespace core;
//using namespace core::conformation;
//using namespace chemical;
//using namespace utility;
//using namespace core::pose::metrics;
//using namespace core::scoring;

typedef numeric::xyzVector< core::Real > Vector;

//basic::options::RealOptionKey const sasa_cutoff("bunsat_sasa_cutoff");

/// @brief
class KHSandbox : public protocols::moves::Mover {
public:
	KHSandbox() {}

	virtual ~KHSandbox(){};

	virtual
	void
	apply(core::pose::Pose & pose)
	{
		utility::file::FileName filename(pose.pdb_info()->name());
		std::string pdbname_base = filename.base();

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scorefxn->score(pose);

		core::kinematics::AtomTree const & at = pose.atom_tree();
		core::kinematics::tree::AtomCOP rt = at.root();

		for ( core::Size i = 1; i < pose.size(); i++ ) {
			core::conformation::Residue const & rsd = pose.residue(i);
			TR << "Residue " << i << ": " << rsd.name3() << std::endl;
			for ( core::Size j = 1; j < rsd.natoms(); j++ ) {
				TR << "\t" << "Atom " << j << ": " << rsd.atom_name(j) << std::endl;
				//TR << "bonded_neighbors: " << rsd.bonded_neighbor(j) << std::endl;
				TR << "bonded_neighbors: " << rsd.atom_name(rsd.bonded_neighbor(j)[1]);
				for ( Size k = 2; k <= rsd.bonded_neighbor(j).size(); k ++ ) {
					TR << "," << rsd.atom_name(k);
				}
				TR << std::endl;
				TR << "atom_base: " << rsd.atom_name(rsd.atom_base(j)) << std::endl;
				TR << "atom_base of atom_base: " << rsd.atom_name(rsd.atom_base(rsd.atom_base(j))) << std::endl;
				if ( rsd.abase2(j) ) { TR << "abase2: " << rsd.atom_name(rsd.abase2(j)) << std::endl; }
				else { TR << "abase2: No abase2 for this atom" << std::endl; }
			}

			//for ( core::chemical::AtomIndices::const_iterator H_pos
			//  = restype->Hpos_polar().begin();
			//  H_pos != restype->Hpos_polar().end(); H_pos++) {

			//}
			for ( core::chemical::AtomIndices::const_iterator acc_pos
					= rsd.type().accpt_pos().begin();
					acc_pos != rsd.type().accpt_pos().end(); acc_pos++ ) {
				core::chemical::Hybridization const & hybrid( rsd.atom_type(*acc_pos).hybridization() );
				TR << "Atom # " << *acc_pos << ", Name: " << rsd.atom_name(*acc_pos) << std::endl;
				TR << "Hybridization: " << hybrid << std::endl;
			}
		}

		//rt->show();

		return;
	}

	virtual
	std::string
	get_name() const { return "KHSandbox"; }

private:

	// basic::MetricValue< id::AtomID_Map< Real > > atom_sasa_;

};

//typedef utility::pointer::owning_ptr<BuriedUnsatPolarsFinder> BuriedUnsatPolarsFinderOP;

int main(int argc, char* argv[])
{
	try {
		using basic::options::option;
		//option.add(sasa_cutoff, "SASA cutoff to be considered buried").def(0.01);

		TR << "devel::init()" << std::endl;
		std::time_t t = std::time(NULL);
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go(new KHSandbox);

		TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
