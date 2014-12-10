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
			 << "	<SCOREFXNS>\n"
			 << "		<ligand_soft_rep weights=ligand_soft_rep>\n"
			 << "			<Reweight scoretype=fa_elec weight=0.42/>\n"
			 << "			<Reweight scoretype=hbond_bb_sc weight=1.3/>\n"
			 << "			<Reweight scoretype=hbond_sc weight=1.3/>\n"
			 << "			<Reweight scoretype=rama weight=0.2/>\n"
			 << "		</ligand_soft_rep>\n"
			 << "		<hard_rep weights=ligand>\n"
			 << "			<Reweight scoretype=fa_intra_rep weight=0.004/>\n"
			 << "			<Reweight scoretype=fa_elec weight=0.42/>\n"
			 << "			<Reweight scoretype=hbond_bb_sc weight=1.3/>\n"
			 << "			<Reweight scoretype=hbond_sc weight=1.3/>\n"
			 << "			<Reweight scoretype=rama weight=0.2/>\n"
			 << "		</hard_rep>\n"
			 << "	</SCOREFXNS>\n"
			 << "	<SCORINGGRIDS>\n"
			 << "		<atr grid_type=\"AtrGrid\" weight=\"1.0\"/>\n"
			 << "		<rep grid_type=\"RepGrid\" weight=\"1.0\"/>\n"
			 << "		<hba grid_type=\"HbaGrid\" weight=\"1.0\"/>\n"
			 << "		<hbd grid_type=\"HbdGrid\" weight=\"1.0\"/>\n"
			 << "		<vdw grid_type=\"VdwGrid\" weight=\"1.0\"/>\n"
			 << "		<classic grid_type=\"ClassicGrid\" weight=\"1.0\"/>\n"
			 << "		<charge grid_type=\"ChargeGrid\" weight=\"1.0\"/>\n"
			 << "		\n"
			 << "	</SCORINGGRIDS>\n"
			 << "	<LIGAND_AREAS>\n"
			 << "		<docking_sidechain chain=X cutoff=6.0 add_nbr_radius=true all_atom_mode=true minimize_ligand=10/>\n"
			 << "		<final_sidechain chain=X cutoff=6.0 add_nbr_radius=true all_atom_mode=true/>\n"
			 << "		<final_backbone chain=X cutoff=7.0 add_nbr_radius=false all_atom_mode=true Calpha_restraints=0.3/>\n"
			 << "	</LIGAND_AREAS>\n"
			 << "	<INTERFACE_BUILDERS>\n"
			 << "		<side_chain_for_docking ligand_areas=docking_sidechain/>\n"
			 << "		<side_chain_for_final ligand_areas=final_sidechain/>\n"
			 << "		<backbone ligand_areas=final_backbone extension_window=3/>\n"
			 << "	</INTERFACE_BUILDERS>\n"
			 << "	<MOVEMAP_BUILDERS>\n"
			 << "		<docking sc_interface=side_chain_for_docking minimize_water=true/>\n"
			 << "		<final sc_interface=side_chain_for_final bb_interface=backbone minimize_water=true/>\n"
			 << "	</MOVEMAP_BUILDERS>\n"
			 << "	<MOVERS>\n"
			 << "	single movers		\n"
			 << "		<StartFrom name=start_from chain=X>\n"
			 << "			<Coordinates x=-1.731 y=32.589 z=-5.039/>\n"
			 << "		</StartFrom>\n"
			 << "		<SlideTogether name=slide_together chains=X/>\n"
			 << "		<Transform name=\"transform\" chain=\"X\" box_size=\"5.0\" move_distance=\"1.0\" angle=\"45\" cycles=\"5000\" temperature=\"100\"/>\n"
			 << "		<InterfaceScoreCalculator name=add_scores chains=X scorefxn=hard_rep native=\"inputs/7cpa_7cpa_native.pdb\"/>\n"
			 << "	compound movers\n"
			 << "		A stride of 5 is used to cut down on integration test file size. In production use a stride of 1 or 2\n"
			 << "		<RenderGridsToKinemage name=\"kineAtr\" file_name=\"output.kin\" grid_name=\"atr\" color=\"1.0,0.0,0.0\" stride=\"5\"/>\n"
			 << "		<RenderGridsToKinemage name=\"kineRep\" file_name=\"output.kin\" grid_name=\"rep\" color=\"0.0,1.0,0.0\" stride=\"5\"/>\n"
			 << "		<RenderGridsToKinemage name=\"kineHba\" file_name=\"output.kin\" grid_name=\"hba\" low_color=\"1.0,1.0,1.0\" high_color=\"1.0,0.0,0.0\" stride=\"5\"/>\n"
			 << "		<RenderGridsToKinemage name=\"kineHbd\" file_name=\"output.kin\" grid_name=\"hbd\" low_color=\"1.0,1.0,1.0\" high_color=\"1.0,1.0,0.0\" stride=\"5\"/>\n"
			 << "		<RenderGridsToKinemage name=\"kineVdw\" file_name=\"output.kin\" grid_name=\"vdw\" low_color=\"1.0,1.0,1.0\" high_color= \"1.0,0.0,1.0\" stride=\"5\"/>\n"
			 << "		<RenderGridsToKinemage name=\"kineClassic\" file_name=\"output.kin\" grid_name=\"classic\" low_color=\"0.0,1.0,0.0\" high_color = \"0.0,0.0,1.0\" stride=\"5\"/>\n"
			 << "		<ParsedProtocol name=low_res_dock>\n"
			 << "			<Add mover_name=start_from/>\n"
			 << "			<Add mover_name=transform/>\n"
			 << "		</ParsedProtocol>\n"
			 << "		<ParsedProtocol name=output_grids>\n"
			 << "			<Add mover_name=\"kineAtr\"/>\n"
			 << "			<Add mover_name=\"kineRep\"/>\n"
			 << "			<Add mover_name=\"kineHba\"/>\n"
			 << "			<Add mover_name=\"kineHbd\"/>\n"
			 << "			<Add mover_name=\"kineVdw\"/>\n"
			 << "			<Add mover_name=\"kineClassic\"/>\n"
			 << "		</ParsedProtocol>\n"
			 << "	</MOVERS>\n"
			 << "	<PROTOCOLS>\n"
			 << "		<Add mover_name=low_res_dock/>\n"
			 << "		<Add mover_name=output_grids/>\n"
			 << "	</PROTOCOLS>\n"
			 << "</ROSETTASCRIPTS>\n";

	}

	virtual void run(core::Real scaleFactor) {
		core::Size reps( (core::Size)(10000*scaleFactor) );
		if( reps == 0 ) { reps = 1; } // do at least one rep, regardless of scaling factor
		for(core::Size i=0; i<reps; i++) {
			utility::tag::Tag::create(rosetta_script_);
			rosetta_script_.clear();
			rosetta_script_.seekg(std::ios_base::beg);
		}
	}

	virtual void tearDown() {}
};
