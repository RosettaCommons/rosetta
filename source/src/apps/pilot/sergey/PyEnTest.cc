// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
///
/// @brief
/// @author Sergey Lyskov

#include <core/scoring/EnergyMap.hh>

#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <devel/init.hh>

#include <core/import_pose/import_pose.hh>


#include <basic/Tracer.hh>

#include <iostream>

#include <core/scoring/methods/PyEnergyMethodRegisterer.hh>
#include <core/scoring/methods/EnergyMethodRegistrator.hh>

using core::pose::Pose;
using core::conformation::Residue;


basic::Tracer TR("PyEnTest");


using namespace core;
//using namespace protocols;


class MyCI1B_Creator : public core::scoring::methods::EnergyMethodCreator
{
public:
	virtual core::scoring::methods::EnergyMethodOP create_energy_method(core::scoring::methods::EnergyMethodOptions const &) const;// { return new MyCI1B; }

	virtual core::scoring::ScoreTypes score_types_for_method() const
	{
		std::cout << "MyCI1B_Creator::score_types_for_method()!" << std::endl;
		core::scoring::ScoreTypes sts;
		sts.push_back( core::scoring::PyRosettaEnergie_first );
		return sts;
	};

};


class MyCI1B : public core::scoring::methods::ContextIndependentOneBodyEnergy
{
	typedef ContextIndependentOneBodyEnergy parent;

public:
	MyCI1B() : parent(new MyCI1B_Creator) {
		std::cout << "MyCI1B::MyCI1B()!" << std::endl;
	}

	virtual void residue_energy(core::conformation::Residue const & /*rsd*/, core::pose::Pose const &, core::scoring::EnergyMap & emap) const
	{
		//std::cout << "MyCI1B residue_energy!" << std::endl;
		//emap[ core::scoring::PyRosettaEnergie_first ] = 1.0;
		emap.set(core::scoring::PyRosettaEnergie_first, 1.0);
	}

	core::scoring::methods::EnergyMethodOP clone() const
	{
		std::cout << "MyCI1B::clone()!" << std::endl;
		return new MyCI1B;
	}

	void indicate_required_context_graphs(utility::vector1< bool > &) const {}

	virtual core::Size version() const { return 1; }
};

core::scoring::methods::EnergyMethodOP MyCI1B_Creator::create_energy_method(core::scoring::methods::EnergyMethodOptions const &) const
{
	std::cout << "MyCI1B_Creator::create_energy_method()!" << std::endl;
	return new MyCI1B;
}


//static core::scoring::methods::EnergyMethodRegistrator< MyCI1B_Creator > ENC;

int main( int argc, char * argv [] )
{

	try {

	using namespace core;

	core::scoring::methods::PyEnergyMethodRegistrator ENC( new MyCI1B_Creator() );

	devel::init(argc, argv);

	core::pose::PoseOP pose = core::import_pose::pose_from_pdb("src/python/bindings/test/data/test_in.pdb");
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("python");

	//scorefxn->set_weight(core::scoring::python, 1.0);

	utility::vector1< double > v;  v.push_back(1.0);
	scorefxn->set_method_weights(core::scoring::python, v);

	std::cout << "Score:" << scorefxn->score(*pose) << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

