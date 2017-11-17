// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>


#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <utility/graph/graph.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>

static basic::Tracer tr( "apps.pilot.lsmgridprobe.cc" );

namespace basic {
namespace options {
namespace OptionKeys {
basic::options::FileOptionKey const stublist_filename("stublist");
}
}
}

std::ofstream stubIdentityFile;
std::ofstream stubstubEnergyFile;
std::ofstream stubtargetEnergyFile;

void ReportStubIdentity(std::string stubname, core::Size stubindex)
{
	stubIdentityFile << stubname << " " << stubindex << std::endl;
}

void ReportStubStubEnergy(core::Size stub_a, core::Size stub_b, core::Real energy)
{
	stubstubEnergyFile << stub_a << " " << stub_b << " " << energy << std::endl;
}

void ReportStubTargetEnergy(core::Size stub, core::Real energy)
{
	stubtargetEnergyFile << stub << " " << energy << std::endl;
}

int main( int argc, char * argv [] )
{
	try {
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		option.add( stublist_filename, "stublist");

		devel::init( argc, argv );

		// Read target pose
		std::string targetFilename;
		core::pose::Pose target_pose;
		if ( option[hotspot::target].user() ) {
			targetFilename = option[ hotspot::target ]();
			core::import_pose::pose_from_file( target_pose, targetFilename , core::import_pose::PDB_file);
		} else {
			utility_exit_with_message("Must specify a target structure using -target <filename>");
		}

		// Read input placement filename
		std::ifstream stubinput;
		std::vector<core::conformation::ResidueCOP> stublist;

		if ( option[ stublist_filename ].user() ) {
			stubinput.open( option[ stublist_filename ]().name().c_str(), std::ifstream::in);
		} else {
			utility_exit_with_message("Must specify a target stublist using -target <filename>");
		}

		// Open report files
		utility::file::FileName stubidentity_filename(option[ stublist_filename ]());
		stubidentity_filename.ext(".stubidentity");
		stubIdentityFile.open(stubidentity_filename.name().c_str());

		utility::file::FileName stubstubenergy_filename(option[ stublist_filename ]());
		stubstubenergy_filename.ext(".stubstubenergy");
		stubstubEnergyFile.open(stubstubenergy_filename.name().c_str());

		utility::file::FileName stubtargetenergy_filename(option[ stublist_filename ]());
		stubtargetenergy_filename.ext(".stubtargetenergy");
		stubtargetEnergyFile.open(stubtargetenergy_filename.name().c_str());

		// Read input placements
		while ( !stubinput.eof() )
				{
			std::string filename;
			std::getline(stubinput, filename);

			core::pose::Pose stubpose;
			core::import_pose::pose_from_file( stubpose, filename, core::import_pose::PDB_file);

			if ( !stubpose.size() == 1 ) {
				utility_exit_with_message( "Stub pose must contain only one residue:" + filename );
			}

			core::conformation::ResidueCOP stub_residue( stubpose.residue( 1 ));
			stublist.push_back( stub_residue );

			ReportStubIdentity(filename, stublist.size());
		}

		// Place stub structures onto target and store residue indicies
		core::Size maxTargetResidue = target_pose.size();

		std::vector<core::Size> stub_residue_indicies;
		std::map<core::Size, core::Size> residue_indicies_to_stub;

		for ( core::Size i = 0; i < stublist.size(); i++ ) {
			target_pose.append_residue_by_jump(*stublist[i], target_pose.size(), "", stublist[i]->atom_name(stublist[i]->nbr_atom()), true);

			stub_residue_indicies.push_back(target_pose.size());
			residue_indicies_to_stub[target_pose.size()] = i;
		}

		// Score the pose
		core::scoring::ScoreFunctionOP scorefxn;
		scorefxn = core::scoring::get_score_function();

		(*scorefxn)(target_pose);

		core::scoring::EnergyMap const graphWeights = target_pose.energies().weights();
		core::scoring::EnergyGraph const & targetEnergyGraph = target_pose.energies().energy_graph();

		for ( core::Size i = 0; i < stub_residue_indicies.size(); i++ ) {
			core::Size stub_residue = stub_residue_indicies[i];
			core::Real targetEnergy = 0;

			for ( utility::graph::EdgeListConstIterator egraph_it = targetEnergyGraph.get_node( stub_residue )->const_edge_list_begin(); egraph_it != targetEnergyGraph.get_node( stub_residue )->const_edge_list_end(); ++egraph_it ) {
				core::scoring::EnergyEdge const * eedge = static_cast< core::scoring::EnergyEdge const * > (*egraph_it);

				core::Real eedge_energy = eedge->dot( graphWeights );
				core::Size eedge_other_residue = eedge->get_other_ind(stub_residue);

				if ( eedge_other_residue <= maxTargetResidue ) {
					targetEnergy += eedge_energy;
				} else if ( eedge_other_residue > stub_residue ) {
					// Only emit edge when passing to greater indexed node so edges are only emitted once.
					ReportStubStubEnergy(i, residue_indicies_to_stub[eedge_other_residue], eedge_energy);
				}
			}

			ReportStubTargetEnergy(i, targetEnergy);
		}

		stubIdentityFile.close();
		stubstubEnergyFile.close();
		stubtargetEnergyFile.close();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
