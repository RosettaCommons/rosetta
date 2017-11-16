// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file rotamerize_disulfides.cc
/// @brief Reads in a pdb to rosetta
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @created October 2008
/// @usage rotamerize_disulfides -database db -extra_res_fa DSF.params
///    -s input.pdb -prefix output_prefix
///    -target_res 1 [-nstructs 9999]


//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>

//Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/blivens.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>


#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/exit.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/protein_interface_design/movers/TryRotamers.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/ScoreType.hh>

#include <algorithm>

using namespace core;
using namespace std;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring;

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
static basic::Tracer TR( "pilot_apps.blivens.rotamerize_disulfides" );

int
usage(string msg)
{
	TR.Error
		<< "usage: rotamerize_disulfides -database db -extra_res_fa DSF.params \\" << ends
		<< "         -s input.pdb -prefix output_prefix                        \\" << ends
		<< "         -target_res 1                                             \\" << ends
		<< "         [-score_type total_score # filter score. default none     \\" <<ends
		<< "          -threshold -1.0] # default -1.0                          \\"
		<< "         [-nstructs 9999] # default unlimited                      " << ends
		<< msg << endl;
	utility_exit();
	exit(1);
}

int main( int argc, char * argv [] )
{
	try {
		//init options system
		devel::init(argc, argv);

		if ( ! option[ in::file::s ].user() ) {
			return usage("No in file given: Use -s to designate pdb files to search for disulfides");
		}
		string pdb = basic::options::start_file();

		if ( ! option[ out::prefix ].user() ) {
			return usage("No out file given: Use -prefix to designate an output file prefix");
		}
		string outprefix = option[ out::prefix ]();

		if ( ! option[ hotspot::target_res ].user() ) {
			return usage("No target_res specified");
		}
		Size target = option[ hotspot::target_res ]();

		bool use_filter( option[ blivens::score_type ].user() );
		ScoreType score_type( score_type_from_name(option[ blivens::score_type ]) );

		Real energy_threshold( option[ hotspot::threshold ]() );

		pose::Pose pose;

		core::import_pose::pose_from_file( pose, pdb , core::import_pose::PDB_file);
		pose.update_residue_neighbors();

		// string res_name = pose.residue(target).name3();

		// ResidueTypeSetCAP restype_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		// ResidueType const& restype(restype_set->name_map(res_name));
		//
		// //mutate the two residues
		// pose.replace_residue(target, dsf, true /*orient backbone*/);

		ResidueType const& restype( pose.residue_type(target));

		//Iterate through all possible rotamers
		scoring::ScoreFunctionOP sfxn = scoring::get_score_function();

		protocols::simple_filters::EnergyPerResidueFilter energy_filter(
			target, sfxn, score_type, energy_threshold );
		// protocols::protein_interface_design::movers::TryRotamers rot_mover(
		//  target, *sfxn, 1/*explosion*/, 0 /*no jumps*/ );
		// rot_mover.set_resnum(target);
		// rot_mover.set_scorefxn(sfxn);

		pack::rotamers::SingleResidueRotamerLibraryCAP rot_lib( RotamerLibrary::get_instance()->get_rsd_library( restype ));
		// pack::rotamers::SingleLigandRotamerLibrary const& lig_rot_lib(dynamic_cast< pack::rotamers::SingleLigandRotamerLibrary const& > (*rot_lib) );

		utility::vector1< conformation::ResidueOP> rotamers;
		utility::vector1< utility::vector1< Real > > extra_chi_steps( restype.nchi() );
		utility::graph::GraphOP graph( new utility::graph::Graph( pose.energies().energy_graph() ));
		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );

		rot_lib->fill_rotamer_vector(pose, *sfxn, *task, graph, ResidueTypeCOP(&restype),
			pose.residue(target),extra_chi_steps, false, //guess its not buried
			rotamers);
		Size n(rotamers.size());

		TR << "Found "<< n << " rotamers." << endl;

		Size max_count(0);
		if ( option[ OptionKeys::out::nstruct ].user() ) {
			max_count = option[ OptionKeys::out::nstruct ]();
		} else {
			max_count = n;
		}
		TR << "Outputing at most " << max_count << " structures." << endl;

		Size count(1);
		for ( vector1< conformation::ResidueOP >::const_iterator rot = rotamers.begin(),
				end_rot = rotamers.end();
				rot != end_rot && count <= max_count;
				++rot ) {
			pose.replace_residue(target, * *rot, true /*orient backbone*/);
			if ( use_filter && !energy_filter.apply(pose) ) {
				continue; // ignore clashing rotamers
			}
			ostringstream outfile;
			outfile << outprefix << '_';
			outfile << right << setw(4) << setfill('0') << count;
			outfile << ".pdb";
			pose.dump_pdb(outfile.str());
			++count;
		}
		// core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
		// task->initialize_from_command_line();
		// task->restrict_to_repacking();
		// vector1<bool> repackable(pose.size() );
		// fill(repackable.begin(),repackable.end(), false);
		// repackable[target] = true;
		// task->restrict_to_residues(repackable);
		//
		// Size const packer_runs( option[ OptionKeys::packing::ndruns ]() );
		// // containers in which to store the results of packing
		// utility::vector1< std::pair< Real, std::string > > results;
		// utility::vector1< pose::PoseOP > pose_list;
		//
		// pack::pack_rotamers_loop( pose, *sfxn, task, packer_runs, results, pose_list );
		//
		// TR.Info << "Found " << pose_list.size() << " rotamers." << endl;
		//
		// //Output pdbs
		// pose.dump_pdb(outprefix + "_0000.pdb");
		// Size count = 1;
		//
		// Size max_count;
		// if( option[ OptionKeys::out::nstruct ].user() ) {
		//  max_count = option[ OptionKeys::out::nstruct ]();
		// }
		// else {
		//  max_count = pose_list.size() + 1;
		// }
		// TR << "Outputing at most " << max_count << " structures." << endl;
		//
		// for(vector1< pose::PoseOP >::const_iterator new_pose = pose_list.begin(),
		//   end_new_pose = pose_list.end();
		//   new_pose != end_new_pose && count < max_count;
		//   ++new_pose )
		// {
		//  ostringstream outfile;
		//  outfile << outprefix << '_';
		//  outfile << right << setw(4) << setfill('0') << count;
		//  outfile << ".pdb";
		//  (*new_pose)->dump_pdb(outfile.str());
		//  ++count;
		// }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // end main
