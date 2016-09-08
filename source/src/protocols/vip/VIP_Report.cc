// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/vip/VIP_Report.cc
/// @brief Prints information about changes implemented in vip mover

#include <protocols/vip/VIP_Utils.hh>
#include <protocols/vip/VIP_Report.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>

#include <core/pose/PDBInfo.hh>

#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/ScoreMover.hh>

#include <fstream>

namespace protocols {
namespace vip {

VIP_Report::VIP_Report() = default;
VIP_Report::~VIP_Report()= default;

void
VIP_Report::get_GOE_repack_report(
	core::pose::Pose & goe_native,
	utility::vector1<core::Real> & goe_repack_e,
	utility::vector1<core::conformation::ResidueOP> & goe_repack_res,
	utility::vector1<core::Size> & goe_repack_pos,
	core::Size it,
	bool use_stored,
	core::Real stored_e
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::ofstream output;
	std::string fname = option[cp::vipReportFile];
	output.open( fname.c_str(), std::ios::out | std::ios::app );

	output << "Iteration " << it << " :  Found candidate mutations:" << std::endl;

	core::scoring::ScoreFunctionOP sf2 = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::pack_sfxn] );
	protocols::simple_moves::ScoreMoverOP score_em( new protocols::simple_moves::ScoreMover(sf2) );
	score_em->apply( goe_native );

	core::Real check_E( use_stored ? stored_e : goe_native.energies().total_energy() );

	core::Size num_accepted( 0 );
	for ( core::Size i = 1; i <= goe_repack_e.size(); i++ ) {
		if ( goe_repack_e[i] < check_E ) {
			if ( goe_repack_res[i]->name() != goe_native.residue(goe_repack_pos[i]).name() ) {
				num_accepted++;
				output << "Position: " << goe_repack_pos[i] << " Native AA: " << goe_native.residue(goe_repack_pos[i]).name() << "  Mutant AA: " << goe_repack_res[i]->name() << "  ddEgoe: " << goe_repack_e[i] - check_E << std::endl;
			}
		}
	}

	if ( num_accepted == 0 ) {
		output << "Iteration  :  No candidate mutations found!" << std::endl;
	}

	output.close();
}


void
VIP_Report::get_GOE_relaxed_report(
	core::pose::Pose & goe_native,
	utility::vector1<core::Real> & goe_relax_e,
	utility::vector1<core::conformation::ResidueOP> & goe_relax_res,
	utility::vector1<core::Size> & goe_relax_pos,
	core::Size it,
	bool use_stored,
	core::Real stored_e
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::fstream output;
	std::string fname = option[cp::vipReportFile];
	output.open( fname.c_str(), std::ios::out | std::ios::app );

	output << "Iteration " << it << " :  The following mutations were accomodated after relaxation:" << std::endl;

	core::scoring::ScoreFunctionOP sf2 = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::relax_sfxn] );
	protocols::simple_moves::ScoreMoverOP score_em( new protocols::simple_moves::ScoreMover(sf2) );
	score_em->apply( goe_native );

	core::Real check_E( use_stored ? stored_e : goe_native.energies().total_energy() );

	core::Size num_accepted( 0 );
	core::Size best_index( 0 );
	core::Real best_score( 9999.0 );
	for ( core::Size i = 1; i <= goe_relax_e.size(); i++ ) {
		//    output << "Check All Score Position: " << goe_native.pdb_info()->number( goe_relax_pos[i] ) << " chain:  " << goe_native.pdb_info()->chain( goe_relax_pos[i] )  << " Native AA: " << goe_native.residue(goe_relax_pos[i]).name() << "  Mutant AA: " << goe_relax_res[i]->name() << "  relax score " << goe_relax_e[i] << " compare with stored " << check_E << std::endl;
		if ( goe_relax_e[i] < check_E ) {
			if ( goe_relax_res[i]->name() != goe_native.residue(goe_relax_pos[i]).name() ) {
				num_accepted++;
				core::Real Ediff( goe_relax_e[i] - check_E );
				if ( Ediff < best_score ) {
					best_index = i;
					best_score = Ediff;
				}
				output << "Position: " << goe_native.pdb_info()->number( goe_relax_pos[i] ) << " chain:  " << goe_native.pdb_info()->chain( goe_relax_pos[i] )  << " Native AA: " << goe_native.residue(goe_relax_pos[i]).name() << "  Mutant AA: " << goe_relax_res[i]->name() << "  ddEgoe: " << Ediff << std::endl;
			}
		}
	}


	if ( num_accepted > 0 ) {
		output << "Accepted mutation from " << goe_native.residue(goe_relax_pos[best_index]).name() << " to " << goe_relax_res[best_index]->name() << " at position " << goe_native.pdb_info()->number( goe_relax_pos[best_index] ) <<  "  chain:  " << goe_native.pdb_info()->chain( goe_relax_pos[best_index] ) << std::endl;
	} else {
		output << "Iteration  :  No mutations were accommodated after relaxation!" << std::endl;
	}

	output.close();
}

void
VIP_Report::get_GOE_packstat_report(
	core::pose::Pose & goe_native,
	utility::vector1<core::pose::Pose> & goe_relax
){
	std::fstream output;
	using namespace core::scoring::packstat;
	//std::string filename;
	//std::string filename0 = "Intermediate";

	output << " GOE packstat report: " << std::endl;

	for ( core::Size i = 1; i <= goe_relax.size(); i++ ) {
		for ( core::Size j = 1; j <= goe_relax[i].size(); j++ ) {
			if ( goe_relax[i].residue(j).name() != goe_native.residue(j).name() ) {
				// filename = filename0 + ".designed" + utility::to_string<int>(i) + ".pdb";
				//        goe_relax[i].dump_pdb( filename );

				output << " Position: " << j << "  Native AA: " << goe_native.residue(j).name() << "  Mutant AA: " << goe_relax[i].residue(j).name() << "  Packing Score: " << output_packstat(goe_relax[i]) << std::endl;
			}
		}
	}
}

}}
