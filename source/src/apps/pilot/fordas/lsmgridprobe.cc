// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <string>
#include <fstream>
#include <vector>

#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>

#include <protocols/cluster/APCluster.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/xyz.functions.hh>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/MinimizingPatternSearch.hh>

#include <utility/excn/Exceptions.hh>


namespace protocols {
namespace hotspot_hashing {



bool tryParseLSMSpec(std::string lsmstring, core::Size & id, VectorPair & lsmspec)
{
	std::vector< std::string > vectors;
	boost::split( vectors, lsmstring, boost::is_any_of(":") );

	if(vectors.size() != 3){
		return false;
	}

	std::vector< std::string > position;
	std::vector< std::string > direction;

	boost::split( position, vectors[1], boost::is_any_of(",") );
	boost::split( direction, vectors[2], boost::is_any_of(",") );

	if(position.size() != 3 || direction.size() != 3){
		return false;
	}

	using boost::lexical_cast;

	id = lexical_cast<core::Size>(vectors[0]);

	lsmspec = VectorPair(
			Vector(
				lexical_cast<core::Real>(position[0]),
				lexical_cast<core::Real>(position[1]),
				lexical_cast<core::Real>(position[2])),
			Vector(
				lexical_cast<core::Real>(direction[0]),
				lexical_cast<core::Real>(direction[1]),
				lexical_cast<core::Real>(direction[2])));

	return true;
}

}
}

namespace basic{
	namespace options{
		namespace OptionKeys{
			basic::options::StringOptionKey const lsmdef("lsmdef");

			basic::options::RealOptionKey const angle_sampling("angle_sampling");
			basic::options::RealOptionKey const translocation_sampling("translocation_sampling");
			basic::options::RealOptionKey const max_radius("max_radius");
			basic::options::RealOptionKey const distance_sampling("distance_sampling");
			basic::options::RealOptionKey const max_distance("max_distance");

			basic::options::BooleanOptionKey const constrain_grid("constrain_grid");
		}
	}
}

int main( int argc, char * argv [] )
{
    try {
	using namespace protocols::hotspot_hashing;

	using basic::options::option;
  using namespace basic::options::OptionKeys;

  option.add( lsmdef, "lsmdef");
  option.add( angle_sampling, "angle_sampling");
  option.add( translocation_sampling, "translocation_sampling");
  option.add( max_radius, "max_radius");
  option.add( distance_sampling, "distance_sampling");
  option.add( max_distance, "max_distance");

  option.add( constrain_grid, "constrain_grid").def(false);

	devel::init( argc, argv );

	// Turn on extra rotamers
	option[ packing::ex1::ex1 ].value( true );
	option[ packing::ex2::ex2 ].value( true );
	option[ packing::extrachi_cutoff ].value( 0 );

	// Necessary to make sure NANs in hbonding don't cause an exit
	option[ in::file::fail_on_bad_hbond ].value( false );

	// Read target pose
	std::string targetFilename;
	core::pose::Pose targetPose;

	if ( option[hotspot::target].user() )
	{
		targetFilename = option[ hotspot::target ]();
		core::import_pose::pose_from_pdb( targetPose, targetFilename );
	}
	else
	{
		utility_exit_with_message("Must specify a target structure using -target <filename>");
	}

	// Read starting residue
	std::string resname;
	if (option[ hotspot::residue ].user() ) {
		resname = option[ hotspot::residue ]()[1];
	}
	else
	{
		utility_exit_with_message("You must specify probe residue identity using -residue <3Letter>");
	}

	// Initialize residue representation
	core::conformation::ResidueOP residue;
	core::chemical::ResidueTypeSet const & residue_set ( targetPose.residue(1).residue_type_set() );
	core::chemical::ResidueType const & restype( residue_set.name_map( resname ) );
	residue = core::conformation::ResidueFactory::create_residue( restype );

	// Read LSM definition
	std::string lsmstring;
	core::Size lsmid;
	VectorPair lsmspec;

	if (option[ lsmdef ].user() ) {
		lsmstring = option[ lsmdef ]();
	}
	else{
		utility_exit_with_message("You must specify lsm spec (id:location:normal) using -lsmdef <id>:<x>,<y>,<z>:<x>,<y>,<z>");
	}

	if (!tryParseLSMSpec(lsmstring, lsmid, lsmspec)){
		utility_exit_with_message("Invalid lsm spec, use format (id:location:normal)  <id>:<x>,<y>,<z>:<x>,<y>,<z>");
	}

	// Generate output filename
	utility::file::FileName output_file(targetPose.pdb_info()->name());
	output_file.ext(resname + "." + boost::lexical_cast<std::string>(lsmid));

	utility::file::PathName	output_directory( option[out::path::all]() );

	utility::file::FileName output_basename(output_file, output_directory);

	// Initialize score function
	core::scoring::ScoreFunctionOP scorefxn;
	if( option[ score::weights ].user() ) {
		scorefxn = core::scoring::get_score_function();
	}
	else {
		scorefxn = core::scoring::get_score_function_legacy( "score13" );
		scorefxn->set_weight( core::scoring::envsmooth, 0 );
	}

	core::scoring::methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.hbond_options().use_hb_env_dep( option[ hotspot::envhb]() );
	scorefxn->set_energy_method_options( options );

	LSMSearchPattern lsm(
			lsmspec,
			option[ angle_sampling ],
			option[ translocation_sampling ],
			option[ max_radius ],
			option[ distance_sampling ],
			option[ max_distance]
			);

	SearchPattern& pattern = lsm;

	MinimizingPatternSearch search(output_basename.name(), targetPose, residue, pattern, scorefxn, option[ constrain_grid ] );
	search.execute();
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;

}
