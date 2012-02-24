// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/vip/VIP_Report.cc
/// @brief Prints information about changes implemented in vip mover


#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/init.hh>
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/packing_score_params.hh>
#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <core/scoring/packstat/SimplePDB.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <core/conformation/Residue.hh>
#include <numeric/all.hh>
#include <utility/all.hh>
#include <core/pose/Pose.hh>
//#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>
#include <utility/options/keys/OptionKey.hh>
#include <fstream>
#include <iostream>
#include <protocols/vip/VIP_Utils.hh>
#include <protocols/vip/VIP_Report.hh>

namespace protocols {
namespace vip {

                VIP_Report::VIP_Report() {}
                VIP_Report::~VIP_Report(){}

void
VIP_Report::get_GOE_repack_report(
	core::pose::Pose & goe_native,
	utility::vector1<core::Real> & goe_repack_e,
	utility::vector1<core::conformation::ResidueOP> & goe_repack_res,
	utility::vector1<core::Size> & goe_repack_pos
){
        using namespace basic::options;
        using namespace basic::options::OptionKeys;
	std::ofstream output;
	std::string fname = option[cp::vipReportFile];
        output.open( fname.c_str(), std::ios::out | std::ios::app );

	output << " GOE point mutant report: " << std::endl;

        core::scoring::ScoreFunctionOP sf2 = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::pack_sfxn] );
	protocols::simple_moves::ScoreMoverOP score_em = new protocols::simple_moves::ScoreMover(sf2);
	score_em->apply( goe_native );

	for( core::Size i = 1; i <= goe_repack_e.size(); i++ ){
	    if( goe_repack_e[i] < goe_native.energies().total_energy() ){
		    if( goe_repack_res[i]->name() != goe_native.residue(goe_repack_pos[i]).name() ){
			output << "Position: " << goe_repack_pos[i] << " Native AA: " << goe_native.residue(goe_repack_pos[i]).name() << "  Mutant AA: " << goe_repack_res[i]->name() << "  ddEgoe: " << goe_native.energies().total_energy() - goe_repack_e[i] << std::endl;}}}
output.close();
}


void
VIP_Report::get_GOE_relaxed_report(
	core::pose::Pose & goe_native,
	utility::vector1<core::Real> & goe_relax_e,
	utility::vector1<core::conformation::ResidueOP> & goe_relax_res,
	utility::vector1<core::Size> & goe_relax_pos
  ){
        using namespace basic::options;
        using namespace basic::options::OptionKeys;
	std::fstream output;
	std::string fname = option[cp::vipReportFile];
        output.open( fname.c_str(), std::ios::out | std::ios::app );

        output << " GOE after relax report: " << std::endl;

        core::scoring::ScoreFunctionOP sf2 = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::relax_sfxn] );
        protocols::simple_moves::ScoreMoverOP score_em = new protocols::simple_moves::ScoreMover(sf2);
        score_em->apply( goe_native );

	for( core::Size i = 1; i <= goe_relax_e.size(); i++ ){
	   if( goe_relax_e[i] < goe_native.energies().total_energy() ){
		    if( goe_relax_res[i]->name() != goe_native.residue(goe_relax_pos[i]).name() ){
			output << "Position: " << goe_relax_pos[i] << " Native AA: " << goe_native.residue(goe_relax_pos[i]).name() << "  Mutant AA: " << goe_relax_res[i]->name() << "  ddEgoe: " << goe_native.energies().total_energy() - goe_relax_e[i] << std::endl;}}}
	output.close();
}

void
VIP_Report::get_GOE_packstat_report(
	core::pose::Pose & goe_native,
	utility::vector1<core::pose::Pose> & goe_relax
  ){
	std::fstream output;
        using namespace core::scoring::packstat;
	std::string filename;
	std::string filename0 = "Intermediate";

        output << " GOE packstat report: " << std::endl;

	for( core::Size i = 1; i <= goe_relax.size(); i++ ){
	        for( core::Size j = 1; j <= goe_relax[i].total_residue(); j++ ){
			if( goe_relax[i].residue(j).name() != goe_native.residue(j).name() ){
//	filename = filename0 + ".designed" + utility::to_string<int>(i) + ".pdb";
//        goe_relax[i].dump_pdb( filename );

	output << " Position: " << j << "  Native AA: " << goe_native.residue(j).name() << "  Mutant AA: " << goe_relax[i].residue(j).name() << "  Packing Score: " << output_packstat(goe_relax[i]) << std::endl;
}}}
}

}}
