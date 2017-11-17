// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/jadolfbr/cluster_utilities/identify_cdr_clusters.cc
/// @brief This relaxes a CDR using cluster constraints.  North_AHO numbering required.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>

#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/constraints/util.hh>


#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

//Options
//#include <basic/options/option.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <core/types.hh>

using namespace protocols::antibody;
using namespace protocols::antibody::design;
using namespace protocols::antibody::constraints;

//using namespace basic::options;
//using namespace basic::options::OptionKeys;


///Documentation:
// This utility app relaxes all CDRs of a an input PDB or list of PDBs using cluster-based dihedral constraints.
// Use: JD2 enabled.  North_AHO numbering scheme required.  To get this: [insert IgClassify website]. No options other than input and output for now.

class RelaxCDRsMover : public protocols::moves::Mover {
public:
	RelaxCDRsMover(){};

	virtual ~RelaxCDRsMover(){};

	virtual
	std::string
	get_name() const {
		return "RelaxCDRs";
	}

	void
	apply(core::pose::Pose & pose){

		//Setup Instances
		AntibodyInfoOP ab_info( new AntibodyInfo(pose, AHO_Scheme, North) );
		ScoreFunctionOP scorefxn = core::scoring::get_score_function(true);
		AntibodyDesignModeler modeler = AntibodyDesignModeler(ab_info);

		//Setup Constraints
		scorefxn->set_weight(dihedral_constraint, 1.0);
		add_harmonic_cluster_constraints(ab_info, pose);

		//Run Relax.
		modeler.set_scorefunction(scorefxn);
		modeler.relax_cdrs(pose, false);
	}
};

int main(int argc, char* argv[]){
	try{
		//Make option for centroid or not.
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new RelaxCDRsMover ) );

		std::cout << "Done! -------------------------------\n";

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return(0);
}


//To be included:
// 1) Centroid-based relax vs. Minimization vs. relax.  For testing
// 2) Control of individual CDR's
