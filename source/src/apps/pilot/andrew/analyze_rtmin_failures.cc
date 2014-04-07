#include <devel/init.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/ResidueKinWriter.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/id/AtomID.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/packer_neighbors.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <core/kinematics/MoveMap.hh>

// C++ headers
#include <fstream>

OPT_1GRP_KEY( Integer, rtminfail, residue_of_interest )

bool
is_nat_rot(
	core::pose::Pose const & ref_pose,
	core::pose::Pose const & min_pose,
	core::Size const resid
)
{
	using core::pack::dunbrack::subtract_chi_angles;

	core::Real const natcut( 40 );
	core::conformation::Residue const & refres( ref_pose.residue( resid ) );
	core::conformation::Residue const & minres( min_pose.residue( resid ) );

	assert( refres.aa() == minres.aa() );

	for ( core::Size ii = 1; ii <= refres.nchi(); ++ii ) {
		if ( refres.type().is_proton_chi( ii ) ) continue;
		if ( std::abs(subtract_chi_angles( refres.chi(ii), minres.chi(ii), refres.aa(), ii )) > natcut ) {
			return false;
		}
	}
	return true;
}

int main( int argc, char * argv [] )
{
  try {

	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::io::pdb;
	using namespace core::graph;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace basic::options;

	NEW_OPT( rtminfail::residue_of_interest, "Residue index (pose numbering) for which rotamers should be examined", 1 );

	devel::init( argc, argv );

	Size const residue_of_interest = option[ rtminfail::residue_of_interest ];

	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

	if ( input_jobs.size() != 1 ) {
		utility_exit_with_message( "Expected exactly one pdb to be specified from the -s or -l flags" );
	}

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, input_jobs[ 1 ]->input_tag() );
	ScoreFunctionOP sfxn = getScoreFunction();
	(*sfxn)(pose);

	PackerTaskOP task = TaskFactory::create_packer_task( pose );
	task->initialize_from_command_line();
	task->set_bump_check( false );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ii != residue_of_interest ) {
			task->nonconst_residue_task( ii ).prevent_repacking();
		} else {
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}
	}

	sfxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

	graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn, task );

	rotamer_set::RotamerSetsOP rotsets( new rotamer_set::RotamerSets() );
	rotsets->set_task( task );
	rotsets->build_rotamers( pose, *sfxn, packer_neighbor_graph );
	rotsets->prepare_sets_for_packing( pose, *sfxn );

	rotamer_set::RotamerSetCOP roi_rotset = rotsets->rotamer_set_for_residue( residue_of_interest );

	kinematics::MoveMapOP mm = new kinematics::MoveMap; mm->set_chi( residue_of_interest, true );
	protocols::simple_moves::MinMover minmover( mm, sfxn, "dfpmin", 0.1, true, false, false );

	//utility::vector1< ResidueOP > minimized_rotamers( roi_rotset->num_rotamers() );
	//utility::vector1< EnergyMap > minimized_scores( roi_rotset->num_rotamers() );
	//utility::vector1< Real >      minimized_tots( roi_rotset->num_rotamers() );

	Real best_natrot_energy = 0.0;
	EnergyMap best_natrot_energies;
	ResidueOP best_natrot = 0;

	Real best_nonnatrot_energy = 0.0;
	EnergyMap best_nonnatrot_energies;
	ResidueOP best_nonnatrot = 0;

	std::cout << "Minimizing " << roi_rotset->num_rotamers() << " rotamers" << std::endl;
	for ( Size ii = 1; ii <= roi_rotset->num_rotamers(); ++ii ) {
		Pose minpose( pose );
		minpose.replace_residue( residue_of_interest, *roi_rotset->rotamer( ii ), false );
		minmover.apply( minpose );

		Real totalE = minpose.energies().total_energy();
		if ( is_nat_rot( pose, minpose, residue_of_interest )) {
			if ( best_natrot() == 0 || totalE < best_natrot_energy ) {
				best_natrot = minpose.residue( residue_of_interest ).clone();
				best_natrot_energy = totalE;
				best_natrot_energies = minpose.energies().total_energies();
			}
		} else {
			if ( best_nonnatrot() == 0 || totalE < best_nonnatrot_energy ) {
				best_nonnatrot = minpose.residue( residue_of_interest ).clone();
				best_nonnatrot_energy = totalE;
				best_nonnatrot_energies = minpose.energies().total_energies();
			}
		}
	}

	bool compare = true;
	if ( best_natrot() == 0 ) {
		compare = false;
		std::cout << "No rotamer minimized to within the cutoff of the native rotamer!" << std::endl;
	}
	if ( best_nonnatrot() == 0 ) {
		compare = false;
		std::cout << "No rotamer minimized away from within the cutoff of the native rotamer!" << std::endl;
	}

	if ( compare ) {
		if ( best_nonnatrot_energy < best_natrot_energy ) {
			std::cout << "Native rotamer is worse than the non-native rotamer" << std::endl;
			EnergyMap diff = best_natrot_energies;
			diff -= best_nonnatrot_energies;
			std::cout << "Differences in energies between the best native and the best non-native" << std::endl;
			diff.show_weighted( std::cout, sfxn->weights() );
			std::cout << std::endl;
		} else {
			std::cout << "The native rotamer is the best rotamer" << std::endl;
		}
		pose.replace_residue( residue_of_interest, *best_natrot, false );
		(*sfxn)(pose);
		pose.dump_pdb( "rtmin_failure_best_native_rotamer.pdb" );

		pose.replace_residue( residue_of_interest, *best_nonnatrot, false );
		(*sfxn)(pose);
		pose.dump_pdb( "rtmin_failure_best_nonnative_rotamer.pdb" );
	}
  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
}

