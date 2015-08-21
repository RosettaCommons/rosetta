// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/xml_parsing.bench.cc
///
/// @brief  Performance benchmark for parsing XML for RosettaScripts etc.
/// @brief  parse a representative set of rosetta scripts documents from the integration tests
/// @author Matthew O'Meara

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>
#include <sstream>

using namespace core;

class XMLParseBenchmark : public PerformanceBenchmark
{
public:
	std::stringstream rosetta_script_;

	XMLParseBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		rosetta_script_.str( std::string() );

		rosetta_script_
			<< "This protocol will simply do low-resolution followed by high-resolution docking.\n"
			<< "It will also report the binding energy (ddg) and buried-surface area (sasa) in the score file.\n"
			<< "<ROSETTASCRIPTS>\n"
			<< "\t<SCOREFXNS>\n"
			<< "\t\t<ligand_soft_rep weights=ligand_soft_rep>\n"
			<< "\t\t\t<Reweight scoretype=fa_elec weight=0.42/>\n"
			<< "\t\t\t<Reweight scoretype=hbond_bb_sc weight=1.3/>\n"
			<< "\t\t\t<Reweight scoretype=hbond_sc weight=1.3/>\n"
			<< "\t\t\t<Reweight scoretype=rama weight=0.2/>\n"
			<< "\t\t</ligand_soft_rep>\n"
			<< "\t\t<hard_rep weights=ligand>\n"
			<< "\t\t\t<Reweight scoretype=fa_intra_rep weight=0.004/>\n"
			<< "\t\t\t<Reweight scoretype=fa_elec weight=0.42/>\n"
			<< "\t\t\t<Reweight scoretype=hbond_bb_sc weight=1.3/>\n"
			<< "\t\t\t<Reweight scoretype=hbond_sc weight=1.3/>\n"
			<< "\t\t\t<Reweight scoretype=rama weight=0.2/>\n"
			<< "\t\t</hard_rep>\n"
			<< "\t</SCOREFXNS>\n"
			<< "\t<SCORINGGRIDS>\n"
			<< "\t\t<atr grid_type=\"AtrGrid\" weight=\"1.0\"/>\n"
			<< "\t\t<rep grid_type=\"RepGrid\" weight=\"1.0\"/>\n"
			<< "\t\t<hba grid_type=\"HbaGrid\" weight=\"1.0\"/>\n"
			<< "\t\t<hbd grid_type=\"HbdGrid\" weight=\"1.0\"/>\n"
			<< "\t\t<vdw grid_type=\"VdwGrid\" weight=\"1.0\"/>\n"
			<< "\t\t<classic grid_type=\"ClassicGrid\" weight=\"1.0\"/>\n"
			<< "\t\t<charge grid_type=\"ChargeGrid\" weight=\"1.0\"/>\n"
			<< "\t\t\n"
			<< "\t</SCORINGGRIDS>\n"
			<< "\t<LIGAND_AREAS>\n"
			<< "\t\t<docking_sidechain chain=X cutoff=6.0 add_nbr_radius=true all_atom_mode=true minimize_ligand=10/>\n"
			<< "\t\t<final_sidechain chain=X cutoff=6.0 add_nbr_radius=true all_atom_mode=true/>\n"
			<< "\t\t<final_backbone chain=X cutoff=7.0 add_nbr_radius=false all_atom_mode=true Calpha_restraints=0.3/>\n"
			<< "\t</LIGAND_AREAS>\n"
			<< "\t<INTERFACE_BUILDERS>\n"
			<< "\t\t<side_chain_for_docking ligand_areas=docking_sidechain/>\n"
			<< "\t\t<side_chain_for_final ligand_areas=final_sidechain/>\n"
			<< "\t\t<backbone ligand_areas=final_backbone extension_window=3/>\n"
			<< "\t</INTERFACE_BUILDERS>\n"
			<< "\t<MOVEMAP_BUILDERS>\n"
			<< "\t\t<docking sc_interface=side_chain_for_docking minimize_water=true/>\n"
			<< "\t\t<final sc_interface=side_chain_for_final bb_interface=backbone minimize_water=true/>\n"
			<< "\t</MOVEMAP_BUILDERS>\n"
			<< "\t<MOVERS>\n"
			<< "\tsingle movers\t\t\n"
			<< "\t\t<StartFrom name=start_from chain=X>\n"
			<< "\t\t\t<Coordinates x=-1.731 y=32.589 z=-5.039/>\n"
			<< "\t\t</StartFrom>\n"
			<< "\t\t<SlideTogether name=slide_together chains=X/>\n"
			<< "\t\t<Transform name=\"transform\" chain=\"X\" box_size=\"5.0\" move_distance=\"1.0\" angle=\"45\" cycles=\"5000\" temperature=\"100\"/>\n"
			<< "\t\t<InterfaceScoreCalculator name=add_scores chains=X scorefxn=hard_rep native=\"inputs/7cpa_7cpa_native.pdb\"/>\n"
			<< "\tcompound movers\n"
			<< "\t\tA stride of 5 is used to cut down on integration test file size. In production use a stride of 1 or 2\n"
			<< "\t\t<RenderGridsToKinemage name=\"kineAtr\" file_name=\"output.kin\" grid_name=\"atr\" color=\"1.0,0.0,0.0\" stride=\"5\"/>\n"
			<< "\t\t<RenderGridsToKinemage name=\"kineRep\" file_name=\"output.kin\" grid_name=\"rep\" color=\"0.0,1.0,0.0\" stride=\"5\"/>\n"
			<< "\t\t<RenderGridsToKinemage name=\"kineHba\" file_name=\"output.kin\" grid_name=\"hba\" low_color=\"1.0,1.0,1.0\" high_color=\"1.0,0.0,0.0\" stride=\"5\"/>\n"
			<< "\t\t<RenderGridsToKinemage name=\"kineHbd\" file_name=\"output.kin\" grid_name=\"hbd\" low_color=\"1.0,1.0,1.0\" high_color=\"1.0,1.0,0.0\" stride=\"5\"/>\n"
			<< "\t\t<RenderGridsToKinemage name=\"kineVdw\" file_name=\"output.kin\" grid_name=\"vdw\" low_color=\"1.0,1.0,1.0\" high_color= \"1.0,0.0,1.0\" stride=\"5\"/>\n"
			<< "\t\t<RenderGridsToKinemage name=\"kineClassic\" file_name=\"output.kin\" grid_name=\"classic\" low_color=\"0.0,1.0,0.0\" high_color = \"0.0,0.0,1.0\" stride=\"5\"/>\n"
			<< "\t\t<ParsedProtocol name=low_res_dock>\n"
			<< "\t\t\t<Add mover_name=start_from/>\n"
			<< "\t\t\t<Add mover_name=transform/>\n"
			<< "\t\t</ParsedProtocol>\n"
			<< "\t\t<ParsedProtocol name=output_grids>\n"
			<< "\t\t\t<Add mover_name=\"kineAtr\"/>\n"
			<< "\t\t\t<Add mover_name=\"kineRep\"/>\n"
			<< "\t\t\t<Add mover_name=\"kineHba\"/>\n"
			<< "\t\t\t<Add mover_name=\"kineHbd\"/>\n"
			<< "\t\t\t<Add mover_name=\"kineVdw\"/>\n"
			<< "\t\t\t<Add mover_name=\"kineClassic\"/>\n"
			<< "\t\t</ParsedProtocol>\n"
			<< "\t</MOVERS>\n"
			<< "\t<PROTOCOLS>\n"
			<< "\t\t<Add mover_name=low_res_dock/>\n"
			<< "\t\t<Add mover_name=output_grids/>\n"
			<< "\t</PROTOCOLS>\n"
			<< "</ROSETTASCRIPTS>\n";

	}

	virtual void run(core::Real scaleFactor) {
		core::Size reps( (core::Size)(100*scaleFactor) ); // amw 10000 to 1000

		if ( reps == 0 ) { reps = 1; } // do at least one rep, regardless of scaling factor
		for ( core::Size i=0; i<reps; i++ ) {
			utility::tag::Tag::create(rosetta_script_);
			rosetta_script_.clear();
			rosetta_script_.seekg(std::ios_base::beg);
		}
	}

	virtual void tearDown() {}
};
