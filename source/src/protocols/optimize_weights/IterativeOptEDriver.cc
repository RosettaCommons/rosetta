// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/optimize_weights/IterativeOptEDriver.cc
/// @brief  Implementation of iterative weight fitting protocol
/// @author Andrew Leaver-Fay -- emulating a protocol by Jim Havranek and Brian Kuhlman.

// Unit headers
#include <protocols/optimize_weights/IterativeOptEDriver.hh>

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/UnfoldedStatePotential.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/scoring/rms_util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/ReplicateTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>


#include <protocols/optimize_weights/OptEData.hh>
#include <protocols/optimize_weights/NestedEnergyTermOptEData.hh>
#include <protocols/optimize_weights/DGBindOptEData.hh>
#include <protocols/optimize_weights/DDGBindOptEData.hh>
#include <protocols/optimize_weights/PNatLigPoseOptEData.hh>

#include <utility/graph/Graph.hh>

#include <core/kinematics/FoldTree.hh>


#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <protocols/optimize_weights/OptEMultifunc.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/util.hh>
#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/statistics/functions.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>

//silent file stuff
#include <core/io/silent/SilentFileData.hh>

#ifdef USEMPI
/// MPI
#include <mpi.h>
#endif

// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>

// option key includes

#include <basic/options/option.hh>
#include <basic/options/keys/optE.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

namespace protocols {
namespace optimize_weights {

using namespace core;
using namespace scoring;
using namespace optimization;

using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace utility;

using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "protocols.optimize_weights.IterativeOptEDriver" );
static THREAD_LOCAL basic::Tracer TR_VERBOSE( "protocols.optimize_weights.IterativeOptEDriver.verbose" );


void attach_debugger();

class ScaleAnnealerTemperatureOperation : public core::pack::task::operation::TaskOperation {

public:
	typedef TaskOperation parent;

	ScaleAnnealerTemperatureOperation(
		ScaleAnnealerTemperatureOperation const & other
	) :
		parent(),
		scale_factor_( other.scale_factor_ )
	{}
	ScaleAnnealerTemperatureOperation( core::Real scale ) : scale_factor_( scale ) {}

	~ScaleAnnealerTemperatureOperation() override = default;


	core::pack::task::operation::TaskOperationOP
	clone() const override {
		return core::pack::task::operation::TaskOperationOP( new ScaleAnnealerTemperatureOperation( *this ) );
	}


	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & task ) const override {
		task.low_temp(    0.3 * scale_factor_ );
		task.high_temp( 100.0 * scale_factor_ );
	}

	void scale_factor( core::Real setting ) {
		scale_factor_ = setting;
	}

private:
	core::Real scale_factor_;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

///
/// @brief
/// Main constructor for the IterativeOptEDriver class.  Note that mpi_rank and mpi_nprocs get set even if
/// USEMPI is not defined.  These values are then used to set MPI_rank_ and MPI_nprocs_.
/// Also calls the initialize_free_and_fixed_terms method.
///
IterativeOptEDriver::IterativeOptEDriver() :
	ligand_repack_pdbs_(),
	ligand_repack_native_poses_(),
	decoy_discrim_data_( /* 0 */ ),
	ligand_discrim_data_( /* NULL */ ),
	dG_binding_data_( /* NULL */ ),
	ddG_bind_optE_data_( /* NULL */ ),
	include_count_( 0 ),
	fixed_count_( 0 ),
	free_count_( 0 ),
	component_weights_( n_optE_data_types, 1.0 ),
#ifdef USEMPI
	tag_( 1 ),
#endif
	outer_loop_counter_( 0 ),
	inner_loop_counter_( 1 ),
	total_positions_( 0 ), // new params by APL for sequence entropy optimization
	count_recovered_( 0 ), // entropy
	aa_obs_( core::chemical::num_canonical_aas, 0 ), // entropy
	aa_exp_( core::chemical::num_canonical_aas, 0 ), // entropy
	aa_freq_obs_( core::chemical::num_canonical_aas, 0.0 ), // entropy
	aa_freq_exp_( core::chemical::num_canonical_aas, 0.0 ), // entropy
	mixing_factor_( 0.0 ),
	outer_loop_last_sequence_recovery_rate_( 0.0 ),
	outer_loop_seq_profile_cross_entropy_( 0.0 ), // entropy
	inner_loop_sequence_recovery_rate_( 0.0 ),
	using_unfolded_energy_term_( false )
{
	// default task factory, generates 'vanilla' PackerTasks
	task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );

	// load custom TaskOperations according to an xml-like utility::tag file
	if ( option[ optE::parse_tagfile ].user() ) {
		using namespace core::pack::task::operation;
		std::string tagfile_name( option[ optE::parse_tagfile ]() );
		read_tagfile_to_taskfactory( tagfile_name, task_factory_ );
		// else use default TaskOperation(s)
	} else if ( ! option[ optE::design_with_minpack ] ) {
		using core::pack::task::operation::TaskOperationCOP;
		task_factory_->push_back( TaskOperationCOP( new pack::task::operation::InitializeFromCommandline ) );
	}

	int mpi_rank( 0 ), mpi_nprocs( 1 );
#ifdef USEMPI
	MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &mpi_nprocs);/* get number of processes */
#endif

	MPI_rank_ = mpi_rank;
	MPI_nprocs_ = mpi_nprocs;

	// only init the refE vector1's if we're using reference energies
	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		before_minimization_reference_energies_.resize( chemical::num_canonical_aas, 0.0 );
		after_minimization_reference_energies_.resize( chemical::num_canonical_aas, 0.0 );
		reference_energies_inner_loop_.resize( chemical::num_canonical_aas, 0.0 );
	}

	intialize_free_and_fixed_energy_terms();
	load_component_weights( component_weights_ );

	// set the using unfolded boolean, by iterating over the score types in both the fixed and free lists and
	// checking for "unfolded". checking if the emaps have non-zero weights associated with "unfolded" would also work.
	for ( auto & score_type_iter : free_score_list_ ) {
		if ( name_from_score_type( score_type_iter ) == "unfolded" ) {
			TR << "IterativeOptEDriver(): setting 'using_unfolded_energy_term_' to true." << std::endl;
			using_unfolded_energy_term_ = true;
		}
	}
	for ( auto & score_type_iter : fixed_score_list_ ) {
		if ( name_from_score_type( score_type_iter ) == "unfolded" ) {
			TR << "IterativeOptEDriver(): setting 'using_unfolded_energy_term_' to true." << std::endl;
			using_unfolded_energy_term_ = true;
		}
	}

}

///
IterativeOptEDriver::~IterativeOptEDriver() = default;

void
IterativeOptEDriver::task_factory( core::pack::task::TaskFactoryCOP tf )
{
	task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory( *tf ) );
}

///
/// @brief
/// Reads in an XML formatted task operation and puts builds a task factory from it.
///
void
IterativeOptEDriver::read_tagfile_to_taskfactory(std::string tagfile_name,
	core::pack::task::TaskFactoryOP task_factory){
	using namespace core::pack::task::operation;
	TaskOperationFactory::TaskOperationOPs tops;
	basic::datacache::DataMap datamap;
	TaskOperationFactory::get_instance()->newTaskOperations( tops, datamap, tagfile_name );
	for ( auto & top : tops ) {
		task_factory->push_back( top );
	}
}

///
/// @brief
/// loads structure into pose - decides between silent or pdb


// PTC - this is a quick and dirty function to intercept the file name intended for pose_from_pdb and retrieve it from a silent file instead
// it dramatically speeds up decoy discrimination (more than 70% of the time is spent on loading pdbs!)
// it uses the path of requested pdb to find silent file, each PDB needs to have all of its structures in its own folder (ie: 1agy/pdb_set.silent)
// it looks within each folder for the filename passed to optE::load_from_silent option
// only used in optimize_decoy_discrimination and use of optE::load_from_silent option is not exhaustively tested!

void
IterativeOptEDriver::load_pose( pose::Pose & pose, std::string const & filename, bool ignore_centroid_input_flag=false )
{
	if ( option[ optE::load_from_silent ].user() ) {
		/// APL -- refactor this.  Static data is unacceptible here.
		static std::string prev_path = "";
		static core::io::silent::SilentFileData * sfd;

		Size slash_index = filename.find_last_of("/\\");
		std::string path = filename.substr(0, slash_index);
		std::string tag = filename.substr(slash_index+1);
		std::string filename = option[ optE::load_from_silent ];
		TR_VERBOSE << "loading: " << tag << "from " << path << "/" << filename << std::endl;

		if ( prev_path != path ) {
			prev_path = path;
			delete sfd;
			sfd = new core::io::silent::SilentFileData();
			sfd->read_file( path + "/" + filename );
		}
		if ( ! sfd ) { utility_exit_with_message("Error loading SilentFileData in IterativeOptEDriver"); }
		(*sfd)[ tag ]->fill_pose( pose );
	} else {
		if ( option[ in::file::centroid_input ] && !ignore_centroid_input_flag ) {
			core::import_pose::centroid_pose_from_pdb( pose, filename , core::import_pose::PDB_file);
		} else {
			core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);

			if ( option[ corrections::score::cenrot ] ) {
				using namespace protocols::simple_moves;
				SwitchResidueTypeSetMoverOP to_cenrot( new SwitchResidueTypeSetMover("centroid_rot") );
				to_cenrot->apply(pose);
			}
		}
	}
}

///
/// @brief
/// The head node has to send out to all the work nodes the list of pdb files they have to do their thing on.
/// It itself doesn't do any of the calculations, right?
/// Work nodes get their list of pdb to work on.
///
void
IterativeOptEDriver::divide_up_pdbs()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	if ( MPI_rank_ == 0 ) {
		utility::vector1< std::string > all_filenames = get_native_pdb_names();

		Size const num_pdbs_per_cpu = all_filenames.size() / MPI_nprocs_;
		Size const nextra = all_filenames.size() - num_pdbs_per_cpu * MPI_nprocs_;

		Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_pdbs_per_cpu;
		for ( Size ii = 1; ii <= my_njobs; ++ii ) {
			native_pdbs_.push_back(all_filenames[ ii ] );
			//#ifndef USEMPI
			next_iteration_pdbs_.push_back(all_filenames[ ii ] );
			//#else
			//next_iteration_pdbs_.push_back( "workdir_" + to_string( MPI_rank_ ) + "/" + all_filenames[ ii ] );
			//#endif
			//TR << "divide_up_pdbs(): PROC #" << MPI_rank_ << " has native pdb " << native_pdbs_[ ii ] << std::endl;
		}

#ifdef USEMPI
		//std::cout << " number of nodes " << MPI_nprocs_ << std::endl;
		Size ii_offset = my_njobs;
		for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
			Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_pdbs_per_cpu;
			MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
			for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
				send_string_to_node( ii, all_filenames[ jj ] );
			}
			ii_offset += ii_njobs;
		}
	} else {
		Size my_njobs;
		MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
		native_pdbs_.reserve( my_njobs );
		for ( Size ii = 1; ii <= my_njobs; ++ii ) {
			native_pdbs_.push_back( receive_string_from_node( 0 ) );
			next_iteration_pdbs_.push_back( native_pdbs_[ ii ] );
			//next_iteration_pdbs_.push_back( "workdir_" + to_string( MPI_rank_ ) + "/" + all_filenames[ my_offset + ii ] );
			//TR << "divide_up_pdbs(): PROC #" << MPI_rank_ << " has native pdb " << native_pdbs_[ ii ] << std::endl;
		}
#endif
	}

#ifdef USEMPI
	for ( Size ii = 1; ii <= native_pdbs_.size(); ++ii ) {
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		//printf("Structure %s assigned to %s (rank = %d)\n", native_pdbs_[ ii ].c_str(), hostname, (int) MPI_rank_);
		//fflush( stdout );
		TR_VERBOSE << "divide_up_pdbs(): structure '" << native_pdbs_[ii] << "' assigned to " << hostname << " (rank = " << MPI_rank_ << ")" << std::endl;
	}
#endif


	// Decoy discrimination option processing...
	if ( option[ optE::optimize_decoy_discrimination ].user() ) {
		if ( MPI_rank_ == 0 ) {
			utility::vector1< std::pair< std::string, std::string > > file_lists;
			utility::vector1< std::string > crystal_native_list;
			std::ifstream native_and_decoy_lists( option[ optE::optimize_decoy_discrimination ]()().c_str() );
			while ( native_and_decoy_lists ) {
				std::string native_files;
				std::string decoy_files;
				std::string crystal_native_file;
				native_and_decoy_lists >> native_files >> decoy_files >> crystal_native_file;
				if ( native_files != "" && decoy_files != "" && crystal_native_file != "" ) {
					file_lists.push_back( std::make_pair( native_files, decoy_files ));
					crystal_native_list.push_back( crystal_native_file );
				}
			}

			Size const num_pdbs_per_cpu = file_lists.size() / MPI_nprocs_;
			Size const nextra = file_lists.size() - num_pdbs_per_cpu * MPI_nprocs_;

			Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_pdbs_per_cpu;
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				decdisc_native_decoy_pairs_.push_back(file_lists[ ii ] );
				decdisc_crystal_natives_.push_back( crystal_native_list[ ii ] );
				//TR << " PROC #" << MPI_rank_ << " "<<  ii << " decdiscrim:  "
				// << decdisc_native_decoy_pairs_[ ii ].first << " "
				// << decdisc_native_decoy_pairs_[ ii ].second << " "
				// << decdisc_crystal_natives_[ ii ] << std::endl;
			}

#ifdef USEMPI
			Size ii_offset = my_njobs;
			for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
				Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_pdbs_per_cpu;
				MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
					send_string_to_node( ii, file_lists[ jj ].first );
					send_string_to_node( ii, file_lists[ jj ].second );
					send_string_to_node( ii, crystal_native_list[ jj ] );
				}
				ii_offset += ii_njobs;
			}
		} else {
			Size my_njobs( 0 );
			MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
			decdisc_native_decoy_pairs_.reserve( my_njobs );
			decdisc_crystal_natives_.reserve( my_njobs );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				std::string native_pdb_list_name = receive_string_from_node( 0 );
				std::string decoy_pdb_list_name = receive_string_from_node( 0 );
				decdisc_native_decoy_pairs_.push_back( std::make_pair( native_pdb_list_name, decoy_pdb_list_name )  );
				std::string crystal_native = receive_string_from_node( 0 );
				decdisc_crystal_natives_.push_back( crystal_native );

				//TR << " PROC #" << MPI_rank_ << " "<<  ii << " decdiscrim:  "
				//	<< native_decoy_pairs_[ ii ].first << " "
				//	<< native_decoy_pairs_[ ii ].second << " "
				//	<< crystal_natives_[ ii ] << std::endl;
			}
#endif
		}

	}

	if ( option[ optE::optimize_ligand_discrimination ].user() ) {
		if ( MPI_rank_ == 0 ) {
			utility::vector1< std::pair< std::string, std::string > > file_lists;
			utility::vector1< std::string > crystal_native_list;
			std::ifstream native_and_decoy_lists( option[ optE::optimize_ligand_discrimination ]()().c_str() );
			while ( native_and_decoy_lists ) {
				std::string native_files;
				std::string decoy_files;
				std::string crystal_native_file;
				native_and_decoy_lists >> native_files >> decoy_files >> crystal_native_file;
				if ( native_files != "" && decoy_files != "" && crystal_native_file != "" ) {
					file_lists.push_back( std::make_pair( native_files, decoy_files ));
					crystal_native_list.push_back( crystal_native_file );
				}
			}

			Size const num_pdbs_per_cpu = file_lists.size() / MPI_nprocs_;
			Size const nextra = file_lists.size() - num_pdbs_per_cpu * MPI_nprocs_;

			Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_pdbs_per_cpu;
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				ligand_native_decoy_pairs_.push_back(file_lists[ ii ] );
				ligand_crystal_natives_.push_back( crystal_native_list[ ii ] );
				//TR << " PROC #" << MPI_rank_ << " "<<  ii << " lig discrim:  "
				// << ligand_native_decoy_pairs_[ ii ].first << " "
				// << ligand_native_decoy_pairs_[ ii ].second << " "
				// << ligand_crystal_natives_[ ii ] << std::endl;
			}

#ifdef USEMPI
			Size ii_offset = my_njobs;
			for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
				Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_pdbs_per_cpu;
				MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
					send_string_to_node( ii, file_lists[ jj ].first );
					send_string_to_node( ii, file_lists[ jj ].second );
					send_string_to_node( ii, crystal_native_list[ jj ] );
				}
				ii_offset += ii_njobs;
			}
		} else {
			Size my_njobs( 0 );
			MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
			ligand_native_decoy_pairs_.reserve( my_njobs );
			ligand_crystal_natives_.reserve( my_njobs );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				std::string native_pdb_list_name = receive_string_from_node( 0 );
				std::string decoy_pdb_list_name = receive_string_from_node( 0 );
				ligand_native_decoy_pairs_.push_back( std::make_pair( native_pdb_list_name, decoy_pdb_list_name )  );
				std::string crystal_native = receive_string_from_node( 0 );
				ligand_crystal_natives_.push_back( crystal_native );

				//TR << " PROC #" << MPI_rank_ << " "<<  ii << " lig discrim:  "
				//	<< ligand_native_decoy_pairs_[ ii ].first << " "
				//	<< ligand_native_decoy_pairs_[ ii ].second << " "
				//	<< ligand_crystal_natives_[ ii ] << std::endl;
			}
#endif
		}

	}

	if ( option[ optE::optimize_ligand_rot ].user() ) {
		if ( MPI_rank_ == 0 ) {
			utility::vector1< std::string > file_list;
			std::ifstream pdb_list_file( option[ optE::optimize_ligand_rot ]()().c_str() );
			while ( pdb_list_file ) {
				std::string pdb_file;
				pdb_list_file >> pdb_file;
				if ( pdb_file != "" ) {
					file_list.push_back( pdb_file );
				}
			}

			Size const num_pdbs_per_cpu = file_list.size() / MPI_nprocs_;
			Size const nextra = file_list.size() - num_pdbs_per_cpu * MPI_nprocs_;

			Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_pdbs_per_cpu;
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				ligand_repack_pdbs_.push_back( file_list[ ii ] );
			}

#ifdef USEMPI
			Size ii_offset = my_njobs;
			for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
				Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_pdbs_per_cpu;
				MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
					send_string_to_node( ii, file_list[ jj ] );
				}
				ii_offset += ii_njobs;
			}
		} else {
			Size my_njobs( 0 );
			MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
			ligand_repack_pdbs_.reserve( my_njobs );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				std::string pdb_name = receive_string_from_node( 0 );
				ligand_repack_pdbs_.push_back( pdb_name );
			}
#endif
		}

	}

	if ( option[ optE::optimize_dGbinding ].user() ) {
		//TR << "divide_up_pdbs(): node " << MPI_rank_ << " reading optimize_dGbinding input file..." << std::endl;
		if ( MPI_rank_ == 0 ) {
			utility::vector1< std::pair< std::string, std::string > > dG_pdb_files;
			utility::vector1< Real > dgs;
			utility::io::izstream dg_data( option[ optE::optimize_dGbinding ]() );
			while ( dg_data ) {
				std::string bound_pdb, unbound_pdb;
				Real dg_experimental;
				dg_data >> bound_pdb;
				if ( bound_pdb == "" ) break;
				dg_data >> unbound_pdb;
				dg_data >> dg_experimental;
				dG_pdb_files.push_back( std::make_pair( bound_pdb, unbound_pdb ) );
				dgs.push_back( dg_experimental );
			}

			Size const num_dgs_per_cpu = dgs.size() / MPI_nprocs_;
			Size const nextra = dgs.size() - num_dgs_per_cpu * MPI_nprocs_;

			//TR << "divide_up_pdbs(): node " << MPI_rank_ << " read " << dgs.size() << " dG pairs: sending " << num_dgs_per_cpu << " to slaves" << std::endl;

			Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_dgs_per_cpu;
			dG_bound_unbound_pairs_.reserve( my_njobs ); dG_binding_.reserve( my_njobs );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				dG_bound_unbound_pairs_.push_back( dG_pdb_files[ ii ] );
				dG_binding_.push_back( dgs[ ii ] );

				//TR << "divide_up_pdbs(): node " << MPI_rank_ << " dG_bind:  "
				// << dG_bound_unbound_pairs_[ ii ].first << " " << dG_bound_unbound_pairs_[ ii ].second << " " << dG_binding_[ ii ] << std::endl;
			}

#ifdef USEMPI
			Size ii_offset = my_njobs;
			for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
				Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_dgs_per_cpu;
				MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
					send_string_to_node( ii, dG_pdb_files[ jj ].first );
					send_string_to_node( ii, dG_pdb_files[ jj ].second );
					MPI_Send( & dgs[ jj ], 1, MPI_DOUBLE, ii, tag_, MPI_COMM_WORLD );
				}
				ii_offset += ii_njobs;
			}
		} else {
			Size my_njobs( 0 );
			MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
			dG_bound_unbound_pairs_.reserve( my_njobs );
			dG_binding_.resize( my_njobs, 0.0 );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				std::string bound_pdb = receive_string_from_node( 0 );
				std::string unbound_pdb = receive_string_from_node( 0 );
				dG_bound_unbound_pairs_.push_back( std::make_pair( bound_pdb, unbound_pdb ));
				MPI_Recv( &dG_binding_[ ii ], 1, MPI_DOUBLE, 0, tag_, MPI_COMM_WORLD, & stat_ );

				//TR << "divide_up_pdbs(): node " << MPI_rank_ << " dG_bind:  "
				//	<< dG_bound_unbound_pairs_[ ii ].first << " " << dG_bound_unbound_pairs_[ ii ].second << " " << dG_binding_[ ii ] << std::endl;
			}
#endif

		}
		//TR << "Exiting divide_up_pdbs for dG bind." << std::endl;
	}

	if ( option[ optE::optimize_ddGmutation ].user() ) {
		if ( MPI_rank_ == 0 ) {
			utility::vector1< std::pair< std::string, std::string > > ddG_mut_files;
			utility::vector1< Real > ddgs;
			utility::io::izstream ddg_data( option[ optE::optimize_ddGmutation ]() );
			while ( ddg_data ) {
				std::string wt_file, mut_file;
				Real ddg_experimental;
				ddg_data >> wt_file;
				if ( wt_file == "" ) break;
				ddg_data >> mut_file;
				ddg_data >> ddg_experimental;
				ddG_mut_files.push_back( std::make_pair( wt_file, mut_file ) );
				ddgs.push_back( ddg_experimental );
			}

			Size const num_ddgs_per_cpu = ddgs.size() / MPI_nprocs_;
			Size const nextra = ddgs.size() - num_ddgs_per_cpu * MPI_nprocs_;

			//TR << "Node 0 read " << ddgs.size() << " ddG pairs: sending " << num_ddgs_per_cpu << " to slaves" << std::endl;

			Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_ddgs_per_cpu;
			ddg_mut_wt_pairs_.reserve( my_njobs ); ddGs_.reserve( my_njobs );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				ddg_mut_wt_pairs_.push_back( ddG_mut_files[ ii ] );
				ddGs_.push_back( ddgs[ ii ] );
				//TR << " PROC #" << MPI_rank_ << " "<<  ii << " ddGmut:  "
				// << ddg_mut_wt_pairs_[ ii ].first << " "
				// << ddg_mut_wt_pairs_[ ii ].second << " "
				// << ddGs_[ ii ] << std::endl;
			}

#ifdef USEMPI
			Size ii_offset = my_njobs;
			for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
				Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_ddgs_per_cpu;
				MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
					send_string_to_node( ii, ddG_mut_files[ jj ].first );
					send_string_to_node( ii, ddG_mut_files[ jj ].second );
					MPI_Send( & ddgs[ jj ], 1, MPI_DOUBLE, ii, tag_, MPI_COMM_WORLD );
				}
				ii_offset += ii_njobs;
			}
		} else {
			Size my_njobs( 0 );
			MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
			ddg_mut_wt_pairs_.reserve( my_njobs );
			ddGs_.resize( my_njobs, 0.0 );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				std::string wt_list = receive_string_from_node( 0 );
				std::string mut_list = receive_string_from_node( 0 );
				ddg_mut_wt_pairs_.push_back( std::make_pair( wt_list, mut_list ));
				MPI_Recv( &ddGs_[ ii ], 1, MPI_DOUBLE, 0, tag_, MPI_COMM_WORLD, & stat_ );
				//TR << " PROC #" << MPI_rank_ << " "<<  ii << " ddGmut:  "
				//	<< ddg_mut_wt_pairs_[ ii ].first << " "
				//	<< ddg_mut_wt_pairs_[ ii ].second << " "
				//	<< ddGs_[ ii ] << std::endl;
				//TR_VERBOSE << "divide_up_pdbs(): node " << MPI_rank_ << " has ddG wt file: '" << ddg_mut_wt_pairs_[ ii ].first << "', ddG mutant file '"
				//	<< ddg_mut_wt_pairs_[ ii ].second << "' and experimental ddG: " << ddGs_[ ii ] << std::endl;
			}

#endif

		}
	}

	if ( option[ optE::optimize_ddG_bind_correlation ].user() ) {
		if ( MPI_rank_ == 0 ) {
			utility::vector1< utility::vector1< std::string > > ddG_bind_files;
			utility::vector1< Real > ddGs_binding;
			utility::io::izstream ddG_bind_data( option[ optE::optimize_ddG_bind_correlation ]() );
			while ( ddG_bind_data ) {
				std::string wt_complexes_file, mut_complexes_file, wt_unbounds_file, mut_unbounds_file;
				Real ddG_experimental;
				ddG_bind_data >> wt_complexes_file; if ( wt_complexes_file == "" ) break;
				ddG_bind_data >> mut_complexes_file;
				ddG_bind_data >> wt_unbounds_file;
				ddG_bind_data >> mut_unbounds_file;
				ddG_bind_data >> ddG_experimental;
				// not sure of a better way to make a vector out of these files
				utility::vector1< std::string > files;
				files.push_back( wt_complexes_file ); files.push_back( mut_complexes_file );
				files.push_back( wt_unbounds_file ); files.push_back( mut_unbounds_file );
				ddG_bind_files.push_back( files );
				ddGs_binding.push_back( ddG_experimental );
			}

			Size const num_ddGs_per_cpu = ddGs_binding.size() / MPI_nprocs_;
			Size const nextra = ddGs_binding.size() - num_ddGs_per_cpu * MPI_nprocs_;

			Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_ddGs_per_cpu;
			// resize the class member variables to the sizes we read in
			ddG_bind_files_.reserve( my_njobs ); ddGs_binding_.reserve( my_njobs );
			// the local variables have the same name as the class member variables, minus a trailing underscore
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				// ddG_bind_files_ is a vector of vectors
				ddG_bind_files_.push_back( ddG_bind_files[ ii ] );
				ddGs_binding_.push_back( ddGs_binding[ ii ] );
			}

#ifdef USEMPI
			Size ii_offset = my_njobs;
			for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
				Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_ddGs_per_cpu;
				MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
					send_string_to_node( ii, ddG_bind_files[ jj ][ DDGBindOptEData::WT_COMPLEXES_LIST_FILE ] );
					send_string_to_node( ii, ddG_bind_files[ jj ][ DDGBindOptEData::MUT_COMPLEXES_LIST_FILE ] );
					send_string_to_node( ii, ddG_bind_files[ jj ][ DDGBindOptEData::WT_UNBOUNDS_LIST_FILE ] );
					send_string_to_node( ii, ddG_bind_files[ jj ][ DDGBindOptEData::MUT_UNBOUNDS_LIST_FILE ] );
					MPI_Send( & ddGs_binding[ jj ], 1, MPI_DOUBLE, ii, tag_, MPI_COMM_WORLD );
				}
				ii_offset += ii_njobs;
			}
		} else {
			// this code is what the slave nodes will execute; basically, receive the work unit
			Size my_njobs( 0 );
			MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
			ddG_bind_files_.reserve( my_njobs );
			ddGs_binding_.resize( my_njobs, 0.0 );
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				std::string wt_complexes_file = receive_string_from_node( 0 );
				std::string mut_complexes_file = receive_string_from_node( 0 );
				std::string wt_unbounds_file = receive_string_from_node( 0 );
				std::string mut_unbounds_file = receive_string_from_node( 0 );

				utility::vector1< std::string > files;
				files.push_back( wt_complexes_file ); files.push_back( mut_complexes_file );
				files.push_back( wt_unbounds_file ); files.push_back( mut_unbounds_file );

				ddG_bind_files_.push_back( files );
				MPI_Recv( & ddGs_binding_[ ii ], 1, MPI_DOUBLE, 0, tag_, MPI_COMM_WORLD, & stat_ );
				//TR << " PROC #" << MPI_rank_ << " "<<  ii << " ddGmut:  "
				//	<< ddg_mut_wt_pairs_[ ii ].first << " "
				//	<< ddg_mut_wt_pairs_[ ii ].second << " "
				//	<< ddGs_[ ii ] << std::endl;
				TR_VERBOSE << "divide_up_pdbs(): node " << MPI_rank_ << " has "
					<< "ddG bind wt complexes file: '" << ddG_bind_files_[ ii ][ DDGBindOptEData::WT_COMPLEXES_LIST_FILE ]
					<< "', ddG bind mut complexes file '"  << ddG_bind_files_[ ii ][ DDGBindOptEData::MUT_COMPLEXES_LIST_FILE ]
					<< "', ddG bind wt unbounds file '"  << ddG_bind_files_[ ii ][ DDGBindOptEData::WT_UNBOUNDS_LIST_FILE ]
					<< "', ddG bind mut unbounds file '"  << ddG_bind_files_[ ii ][ DDGBindOptEData::MUT_UNBOUNDS_LIST_FILE ]
					<< "' and experimental ddG bind: " << ddGs_binding_[ ii ] << std::endl;
			}
#endif
		}
	}


	if ( option[ optE::rescore::context_round ].user() ) {
		Size context_round = option[ optE::rescore::context_round ]();
		if ( context_round != 0 ) {
			setup_pdbnames_next_round( context_round, next_iteration_pdbs_, native_pdbs_ );
		}
	}

}

///
/// @brief
///
///
void
IterativeOptEDriver::collect_rotamer_energies()
{
	using namespace pack::rotamer_set;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	/// Do this once per iteration through the outer loop
	/// Make sure this happens before we abort rotamer-energy collection, if we
	/// are using the -design_first flag.
	pdbs_this_round_ = next_iteration_pdbs_;
	setup_pdbnames_next_round( outer_loop_counter_, next_iteration_pdbs_, native_pdbs_ );


	/// Don't bother collecting rotamer optE data if we're not going to optimize
	/// weights.  Just skip straight to the design step.
	//if ( outer_loop_counter_ == 1 && option[ optE::design_first ].user() ) {
	// TR_VERBOSE << "collect_rotamer_energies(): design_first flag in use. leaving method" << std::endl;
	// return;
	//}
	// taking this out because otherwise the rescore.log isn't made correctly (-ronj)

	compute_rotamer_energies_for_assigned_pdbs();
	compute_rotamers_around_ligands();
	collect_decoy_discrimination_data();
	collect_ligand_discrimination_data();
	collect_dG_of_binding_data();
	collect_ddG_of_mutation_data();
	collect_ddG_of_binding_data();

	if ( ! option[ optE::mpi_weight_minimization ] ) {
		if ( MPI_rank_ == 0 ) {
			collect_rotamer_energies_from_slave_cpus();
		} else {
			send_rotamer_energies_to_master_cpu();
		}
	} // else, keep the data on the original cpus for now; send the data after minimization completes.

	/// If we're simply rescoring the optE data for a particular weight set,
	/// quit as soon as rotamer data has been collected and written to a file.
	if ( option[ optE::rescore::weights ].user() ) {
		score_position_data();

		// extra option for testing sequence recovery with a weight set
		if ( option[ optE::rescore::measure_sequence_recovery ].user() ) {
			// optimize_weights();  // leaving this line to be extra clear; definitely don't optimize the weights!

			// score_position_data() leaves all the weights in the vars and fixed_terms vectors. need to move the values
			// out of those containers and into the after_minimization containers because that's what the write_new()
			// function expects the values to be in.
			free_weights_after_minimization_ = free_parameters_;  // fixed params are set in score_position_data()

			write_new_scorefile();
			test_sequence_recovery();
		}

		exit_gracefully();
	}

}

void
IterativeOptEDriver::exit_gracefully()
{
	barrier();
#ifdef USEMPI
		MPI_Finalize();
#endif
	exit( 0 );
}

///
/// @brief
/// include_terms_ is an EnergyMap, as well.  I think this function sets up the free and fixed score lists which are just
/// a vector1 of ScoreType objects.  include_, fixed_ and free_count_ are just (Size) member variables.
///
void
IterativeOptEDriver::setup_derived_free_and_fixed_data()
{
	include_terms_ = free_parameters_;
	include_terms_ += fixed_parameters_;

	include_count_ = fixed_count_ = free_count_ = 0;

	for ( int i=1 ; i <= n_score_types ; ++i ) {
		if ( include_terms_[ ScoreType(i) ] != 0.0 ) {
			if ( fixed_parameters_[ ScoreType(i) ] == 0.0 ) {
				free_score_list_.push_back( ScoreType(i) );
				++free_count_;
			}
			++include_count_;
		}
	}

	for ( int i = 1; i <= n_score_types; ++i ) {
		if ( fixed_parameters_[ ScoreType(i) ] != 0.0 ) {
			fixed_score_list_.push_back( ScoreType(i) );
			++fixed_count_;
		}
	}

}

///
/// @brief
/// Computes the rotamer energies for all positions for all pdbs given (in the call to get_nat_aa_opte_data()).
/// Also, optionally, does the same for native rotamer recovery data.  Creates an unweighted score function using
/// the the class method 'create_unweighed_scorefunction'.  The scorefunction used to get the interaction energies
/// of the rotamers is the one that's created here.
///
/// Note: Surface is a PackerTask option.  No repacking/design is done here. Just scoring.  So the surface
/// energies would not be included here in the optEdata.  If surface was made into a ScoreType, then it could be
/// added. -ronj
/// Update: Surface score is now it's own EnergyMethod. If 'surface' is detected in the ScoreFunction, then the
/// EnergyMethod calculates the surface energy of the pose. You don't get per-residue surface energies this way, but
/// you do get the total pose surface energy. -ronj
///
void
IterativeOptEDriver::compute_rotamer_energies_for_assigned_pdbs()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::rotamer_set;

	if ( MPI_rank_ == 0 ) TR << "compute_rotamer_energies_for_assigned_pdbs(): entered method" << std::endl;

	optE_data_ = OptEDataOP( new OptEData ); // get rid of old optEdata...

	if ( MPI_rank_ == 0 && option[ optE::constrain_weights ].user() ) {
		ConstraintedOptimizationWeightFuncOP cst( new ConstraintedOptimizationWeightFunc( free_score_list_ ) );
		std::string cstfilename = option[ optE::constrain_weights ]();
		std::ifstream input( cstfilename.c_str() );
		cst->initialize_constraints_from_file( input );
		optE_data_->add_position_data( cst );
	}

	ScoreFunctionOP scorefxn = create_unweighted_scorefunction();
	if ( MPI_rank_ == 0 ) {
		TR << "compute_rotamer_energies_for_assigned_pdbs(): created scorefxn for calculating rotamer energies" << std::endl;
		scorefxn->show( std::cout );
	}

	// for when unfolded term is in use... -ronj
	// the thing to be careful about here is to create a scorefunction which has a weight for the unfolded term, but doesn't
	// set the unfolded term's method weights. if you do that, then the unfolded energy method will return actual values.
	// that's only meant to be used during the design steps of optE, though, not the rotamer energy collection steps.

	// do this loop for every pdb we have in the native_pdbs_ list
	for ( Size n=1; n<= native_pdbs_.size(); ++n ) {
		//std::string const & filename( pdbs_this_round_[n] );
		std::string const & native_filename( native_pdbs_[n] );

		if ( option[ optE::optimize_pssm ] ) {
			load_pssm_data( native_filename, n );
		}

		core::pose::Pose pose, native_pose;

		if ( outer_loop_counter_ == 1 ) {
			if ( option[ in::file::centroid_input ] ) {
				core::import_pose::centroid_pose_from_pdb( native_pose, native_filename , core::import_pose::PDB_file);
			} else {
				core::import_pose::pose_from_file( native_pose, native_filename , core::import_pose::PDB_file);

				if ( option[ corrections::score::cenrot ] ) {
					using namespace protocols::simple_moves;
					SwitchResidueTypeSetMoverOP to_cenrot( new SwitchResidueTypeSetMover("centroid_rot") );
					to_cenrot->apply(native_pose);
				}

			}
			pose = native_pose;

			/// these are stored regardless of whether or not no_design is on the command line...
			/// useful if rotamer recovery alone is being measured
			native_poses_.push_back( native_pose );
			TR_VERBOSE << "compute_rotamer_energies_for_assigned_pdbs(): pushing " << native_filename << " onto native_poses_ vector." << std::endl;
			context_poses_.push_back( native_pose );
			if ( option[ optE::recover_nat_rot ] ) rotamer_recovery_context_poses_.push_back( native_pose );

		} else {
			pose = context_poses_[ n ];
			native_pose = native_poses_[ n ];
		}

		TR_VERBOSE << "compute_rotamer_energies_for_assigned_pdbs(): " << node_name( MPI_rank_ ) << " calling get_opte_data for " << native_filename << std::endl;

		(*scorefxn)( pose );
		(*scorefxn)( native_pose );

		utility::file::FileName natname( native_filename );

		if ( option[ optE::optimize_nat_aa ] || option[ optE::optimize_pssm ] ) {
			get_nat_aa_opte_data(
				natname.base(), pose, native_pose,
				*scorefxn, free_score_list_, fixed_score_list_,
				*optE_data_ );
		}

		if ( option[ optE::optimize_nat_rot ]() ) {
			core::pose::Pose context_pose;

			/// Should we use the previously repacked pose, or the native pose to gather data from?
			if ( option[ optE::recover_nat_rot ] ) {
				context_pose = rotamer_recovery_context_poses_[ n ];
			} else {
				context_pose = native_pose;
			}

			utility::vector1<bool> include_rsd( context_pose.size(), true );

			for ( Size j(1); j <= context_pose.size(); ++j ) {
				include_rsd[j] = pose.residue_type(j).is_protein();
			}

			get_nat_rot_opte_data(
				natname.base(), context_pose, native_pose, include_rsd, *scorefxn,
				free_score_list_, fixed_score_list_, *optE_data_ );
		}

	} // now do it all over again for another pdb in the list native_pdbs_

}

///
void
IterativeOptEDriver::load_pssm_data(
	std::string const & native_filename,
	Size const which_protein // which of the several proteins that this node is responsible for redesigning
)
{
	if ( outer_loop_counter_ == 1 ) {

		std::string native_substr = native_filename.substr( 0, native_filename.size() - 4 );
		std::string pssm_file_name = native_substr + ".fasta.probs";
		//std::cerr << "Openning PSSM File " << pssm_file_name << std::endl;
		std::ifstream pssm_file( pssm_file_name.c_str() );

		std::list< std::pair< chemical::AA, utility::vector1< Real > > > pssm_data;
		utility::vector1< Real > pssm_prob_dist( chemical::num_canonical_aas, 0.0 );
		Size linenum( 0 );
		while ( pssm_file ) {
			++linenum;
			char line_aa;
			pssm_file >> line_aa;
			chemical::AA aa( chemical::aa_from_oneletter_code( line_aa ));
			Real sum( 0.0 );
			for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
				pssm_file >> pssm_prob_dist[ ii ];
				sum += pssm_prob_dist[ ii ];
			}
			if ( std::abs( sum - 1 ) > 0.001 ) {
				TR << "Warning: pssm probability distribution does not sum to 1.0: " << sum << std::endl;
				TR << "Problem on line " << linenum << " of " << pssm_file_name << std::endl;
			}
			pssm_data.push_back( std::make_pair( aa, pssm_prob_dist ));
		}
		pssm_data_.clear();
		pssm_data_.resize( pssm_data.size() );

		// copy to vector
		std::copy(  pssm_data.begin(), pssm_data.end(), pssm_data_.begin() );
		all_pssm_data_.push_back( pssm_data_ );

		if ( pssm_data_.size() == 0 ) { std::cerr << "Did not read file -- possibly not found" << std::endl; }
	} else {
		pssm_data_ = all_pssm_data_[ which_protein ];
	}
}

#ifdef USEMPI
///
///
/// @brief
/// Takes a std::string and a destination and constructs the MPI_Send call.
///
void
IterativeOptEDriver::send_string_to_node( int destination, std::string const & string_to_send )
{
	int tag( 1 );
	int len( string_to_send.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	MPI_Send( const_cast< char * > (string_to_send.c_str()), len, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
}

///
/// @brief
/// Receive a string from the master node.  First find out how long the message is, then allocate space for it and
/// actually receive the message. Returns to the calling function the string that was received.
///
std::string
IterativeOptEDriver::receive_string_from_node( int source )
{
	int len( 0 );
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &len, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
	char * str = new char[ len + 1 ];
	str[ len ] = '\0'; // ? do I need null terminated strings?
	MPI_Recv( str, len, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
	std::string return_string( str, len );
	delete [] str;
	return return_string;

}

#endif

///
/// @brief
/// Used by all slave nodes; sends the rotamer energies (according to the protocol we've set up here) to the master node.
////
void
IterativeOptEDriver::send_rotamer_energies_to_master_cpu()
{
#ifdef USEMPI

	using namespace core::pack::rotamer_set;
	//std::cout << " PROC #" << MPI_rank_ << " send_rotamer_energies_to_master_cpu" << std::endl;

	/// 1. Sanity: send a "boolean" vector of the free and fixed energy terms
	int * free_energy_terms = new int[ core::scoring::n_score_types ];
	int * fixed_energy_terms = new int[ core::scoring::n_score_types ];
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		free_energy_terms[  ii - 1 ] = (int) (free_parameters_[  (ScoreType) ii ] != 0.0);
		fixed_energy_terms[ ii - 1 ] = (int) (fixed_parameters_[ (ScoreType) ii ] != 0.0);
	}
	MPI_Send( free_energy_terms, n_score_types, MPI_INT, 0, tag_, MPI_COMM_WORLD );
	MPI_Send( fixed_energy_terms, n_score_types, MPI_INT, 0, tag_, MPI_COMM_WORLD );

	delete [] free_energy_terms;   free_energy_terms = 0;
	delete [] fixed_energy_terms; fixed_energy_terms = 0;

	/// 2. Number of positions on which OptE data has been gathered
	Size n_pos = optE_data_->num_positions();
	MPI_Send( & n_pos, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );

	for ( OptEPositionDataOPs::const_iterator
			iter = optE_data_->position_data_begin(),
			iter_end = optE_data_->position_data_end();
			iter != iter_end; ++iter ) {
		int position_data_type = (*iter)->type();

		MPI_Send( & position_data_type, 1, MPI_INT, 0, tag_, MPI_COMM_WORLD );
		(*iter)->send_to_node( 0, tag_ );
	}
#endif
}

///
/// @brief
/// Helper method for collecting energies. Calls collect_rotamer_energies_from_slave_cpu for all CPU's being used.
////
void IterativeOptEDriver::collect_rotamer_energies_from_slave_cpus()
{
	using namespace core::pack::rotamer_set;

	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		collect_rotamer_energies_from_slave_cpu( ii );
	}
	TR << "collect_rotamer_energies_from_slave_cpus(): master node with " << optE_data_->num_positions() << " positions" << std::endl;

	Size total( 0 );
	for ( auto
			iter = optE_data_->position_data_begin(),
			iter_end = optE_data_->position_data_end();
			iter != iter_end; ++iter ) {
		total += (*iter)->memory_use();
	}
	TR << "collect_rotamer_energies_from_slave_cpus(): master node using " << total << " bytes for "
		<< optE_data_->num_positions() << " positions." << std::endl;

}

///
/// @brief
/// Collect optE data for decoy discrimination. Similar to the get_nat_aa_opte_data method.
////
void
IterativeOptEDriver::collect_decoy_discrimination_data()
{
	using namespace core::io::pdb;
	using namespace core::pack::rotamer_set;
	using namespace core::pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;

	utility::vector1< Real > free_data( free_score_list_.size() );
	utility::vector1< Real > fixed_data( fixed_score_list_.size() );

	bool const calc_decoy_score_sd( option[ optE::normalize_decoy_score_spread ].user() );


	if ( decoy_discrim_data_ == nullptr ) {

		if ( option[ optE::ramp_nativeness ] ) {
			if ( option[ optE::min_decoy_rms_to_native ].user()  ) {
				PNatStructureOptEData::set_nativeness_high( option[ optE::min_decoy_rms_to_native ] );
			}
			if ( option[ optE::max_rms_from_native ].user() ) {
				PNatStructureOptEData::set_nativeness_low( option[ optE::max_rms_from_native ] );
			}
		}

		decoy_discrim_data_ = OptEDataOP( new OptEData );

		ScoreFunctionOP scorefxn = create_unweighted_scorefunction();
		ScoreFunctionOP weighted_sfxn = create_weighted_scorefunction();

		/// Collect decoy energies to compute standard deviation
		ScoreFunctionOP decoy_spread_reference_sfxn;

		if ( option[ optE::normalize_decoy_score_spread ].user() ) {
			/// special case if you say "SCORE12" as the weights file, then create the score function as the standard.wts + score12.wts_patch
			if ( option[ optE::normalize_decoy_score_spread ]() == "SCORE12" ) {
				decoy_spread_reference_sfxn = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, SCORE12_PATCH );
			} else {
				decoy_spread_reference_sfxn = ScoreFunctionFactory::create_score_function( option[ optE::normalize_decoy_score_spread ] );
			}
		}

		if ( option[ optE::repack_and_minimize_decoys ] ) {
			decdisc_native_poses_.resize( decdisc_native_decoy_pairs_.size() );
			decdisc_decoy_poses_.resize( decdisc_native_decoy_pairs_.size() );
			decdisc_xtal_natives_.resize( decdisc_native_decoy_pairs_.size() );
		}

		for ( Size ii = 1; ii <= decdisc_native_decoy_pairs_.size(); ++ii ) {


			PNatStructureOptEDataOP structure_data( new PNatStructureOptEData );
			if ( option[ optE::n_top_natives_to_optimize ].user() ) {
				structure_data->n_top_natives_to_score( option[ optE::n_top_natives_to_optimize ] );
			}

			{//scope
				utility::file::FileName cryst_fname( decdisc_crystal_natives_[ ii ] );
				structure_data->tag( cryst_fname.base() ); // trim path and extension data -- beautiful utility
			}

			utility::vector1< std::string > native_pdb_names, decoy_pdb_names;
			//std::cout << " PROC #" << MPI_rank_ << " reading pdb lists " ;
			//std::cout << decdisc_native_decoy_pairs_[ ii ].first << " and ";
			//std::cout << decdisc_native_decoy_pairs_[ ii ].second << std::endl;

			std::ifstream native_pdblist( decdisc_native_decoy_pairs_[ ii ].first.c_str() );
			while ( native_pdblist ) {
				std::string native_pdb;
				native_pdblist >> native_pdb;
				if ( native_pdb != "" ) native_pdb_names.push_back( native_pdb );
			}
			if ( native_pdb_names.size() == 0 ) {
				std::cerr << "ERROR: no native structures specified in " << decdisc_native_decoy_pairs_[ ii ].first << " on node " << MPI_rank_ << std::endl;
			}

			std::ifstream decoy_pdblist( decdisc_native_decoy_pairs_[ ii ].second.c_str() );
			while ( decoy_pdblist ) {
				std::string decoy_pdb;
				decoy_pdblist >> decoy_pdb;
				if ( decoy_pdb != "" ) decoy_pdb_names.push_back( decoy_pdb );
			}
			if ( decoy_pdb_names.size() == 0 ) {
				std::cerr << "ERROR: no native structures specified in " << decdisc_native_decoy_pairs_[ ii ].second << " on node " << MPI_rank_ << std::endl;
			}

			TR_VERBOSE << "collect_decoy_discrimination_data(): scoring natives and decoys of " << structure_data->tag() << std::endl;

			/// Collect decoy energies to compute standard deviation
			utility::vector1< Real > decoy_energies;
			if ( option[ optE::normalize_decoy_score_spread ].user() ) {
				decoy_energies.reserve( decoy_pdb_names.size() + native_pdb_names.size() );
			}


			core::pose::Pose crystal_native;
			load_pose( crystal_native, decdisc_crystal_natives_[ ii ], false );
			if ( option[ optE::repack_and_minimize_decoys ] ) {
				decdisc_xtal_natives_[ ii ] = crystal_native;
			}

			Size first_total_residue( 0 );
			for ( Size jj = 1; jj <= native_pdb_names.size(); ++jj ) {
				//std::cout << " PROC #" << MPI_rank_ << " reading pdb: #" << jj << " " << native_pdb_names[ jj ] << std::endl;
				/// read the pdb into a pose
				core::pose::Pose pose;
				load_pose( pose, native_pdb_names[ jj ], false );

				if ( option[ optE::repack_and_minimize_decoys ] ) {
					decdisc_native_poses_[ ii ].push_back( pose );
				}

				if ( jj == 1 ) {
					structure_data->set_total_residue( pose.size() );
					first_total_residue = pose.size();
				} else if ( first_total_residue != pose.size() ) {
					std::cerr << "Warning: total_residue for " << native_pdb_names[ jj ];
					std::cerr << "not equal to native #1 total_residue: " << first_total_residue << " vs ";
					std::cerr << pose.size() << std::endl;
					std::cerr << "Excluding structure!" << std::endl;
					continue;
				}

				SingleStructureDataOP ssd = single_structure_data_for_pose(
					scorefxn, pose, crystal_native, free_data, fixed_data, native_pdb_names[ jj ] );

				AddStatus added_native = add_structure_based_on_rms( ssd, structure_data, true /* intended native */ );
				if ( calc_decoy_score_sd && added_native == ADDED_STRUCTURE_OPPOSITE_AS_INTENDED ) {
					decoy_energies.push_back( (*decoy_spread_reference_sfxn)( pose ) );
				}

				if ( option[ optE::repack_and_minimize_input_structures ] ) {
					repack_and_minimize_pose( pose, weighted_sfxn );
					SingleStructureDataOP ssd = single_structure_data_for_pose(
						scorefxn, pose, crystal_native, free_data, fixed_data, "rpmin_0_"+native_pdb_names[ jj ]  );

					add_structure_based_on_rms( ssd, structure_data, true /* intended native */ );
				}
			}

			for ( Size jj = 1; jj <= decoy_pdb_names.size(); ++jj ) {

				/// read the pdb into a pose
				core::pose::Pose pose;
				//std::cout << " PROC #" << MPI_rank_ << " reading pdb: #" << jj << " " << decoy_pdb_names[ jj ] << std::endl;
				load_pose( pose, decoy_pdb_names[ jj ], false );

				if ( first_total_residue != pose.size() ) {
					std::cerr << "Warning: total_residue for " << decoy_pdb_names[ jj ];
					std::cerr << "not equal to native #1 total_residue: " << first_total_residue << " vs ";
					std::cerr << pose.size() << std::endl;
					std::cerr << "Excluding structure!" << std::endl;
					continue;
				}

				if ( option[ optE::repack_and_minimize_decoys ] ) {
					decdisc_decoy_poses_[ ii ].push_back( pose );
				}

				SingleStructureDataOP ssd = single_structure_data_for_pose(
					scorefxn, pose, crystal_native, free_data, fixed_data, decoy_pdb_names[ jj ] );

				AddStatus added_decoy = add_structure_based_on_rms( ssd, structure_data, false /* intended native */ );
				if ( calc_decoy_score_sd && added_decoy == ADDED_STRUCTURE_AS_INTENDED ) {
					decoy_energies.push_back( (*decoy_spread_reference_sfxn)( pose ) );
				}


				if ( option[ optE::repack_and_minimize_input_structures ] ) {
					repack_and_minimize_pose( pose, weighted_sfxn );
					SingleStructureDataOP ssd2 = single_structure_data_for_pose(
						scorefxn, pose, crystal_native, free_data, fixed_data, "rpmin_0_" + decoy_pdb_names[ jj ]  );

					add_structure_based_on_rms( ssd2, structure_data, false /* intended native */ );
				}

			}
			if ( calc_decoy_score_sd ) {
				Real decoy_score_sd = numeric::statistics::std_dev( decoy_energies.begin(), decoy_energies.end(), Real( 0.0 ) );
				if ( decoy_score_sd != 0 ) {
					structure_data->set_normalize_decoy_stddev( true );
					structure_data->set_initial_decoy_stddev( decoy_score_sd );
				}
			}
			decoy_discrim_data_->add_position_data( structure_data );
		}
	}

	if ( option[ optE::repack_and_minimize_decoys ] && outer_loop_counter_ != 1 ) {
		Size count = 0;
		ScoreFunctionOP weighted_sfxn = create_weighted_scorefunction();
		ScoreFunctionOP unweighted_sfxn = create_unweighted_scorefunction();

		for ( auto
				iter = decoy_discrim_data_->position_data_begin(),
				iter_end = decoy_discrim_data_->position_data_end();
				iter != iter_end; ++iter ) {

			utility::vector1< Pose > new_nats, new_decs;
			utility::vector1< Real > new_nats_scores, new_decs_scores;

			++count;
			runtime_assert( dynamic_cast< PNatStructureOptEData * > ( (*iter).get() ) );
			PNatStructureOptEDataOP structure_data(
				utility::pointer::static_pointer_cast< PNatStructureOptEData > ( (*iter) ) );
			/// Create new natives
			for ( Size ii = 1, iie = decdisc_native_poses_[ count ].size(); ii <= iie; ++ii ) {
				core::pose::Pose pose = decdisc_native_poses_[ count ][ ii ];
				repack_and_minimize_pose( pose, weighted_sfxn );
				new_nats.push_back( pose );
				new_nats_scores.push_back( ( *weighted_sfxn )( pose ) );
				SingleStructureDataOP ssd = single_structure_data_for_pose(
					unweighted_sfxn, pose, decdisc_xtal_natives_[ count ], free_data, fixed_data,
					"rpmin_nat_" + utility::to_string( outer_loop_counter_ ) );
				add_structure_based_on_rms( ssd, structure_data, true /* intended native */ );
			}
			/// Create new decoys
			for ( Size ii = 1, iie = decdisc_decoy_poses_[ count ].size(); ii <= iie; ++ii ) {
				core::pose::Pose pose = decdisc_decoy_poses_[ count ][ ii ];
				repack_and_minimize_pose( pose, weighted_sfxn );
				new_decs.push_back( pose );
				new_decs_scores.push_back( ( *weighted_sfxn )( pose ) );
				SingleStructureDataOP ssd = single_structure_data_for_pose(
					unweighted_sfxn, pose, decdisc_xtal_natives_[ count ], free_data, fixed_data ,
					"rpmin_dec_" + utility::to_string( outer_loop_counter_ ) );
				add_structure_based_on_rms( ssd, structure_data, false /* intended native */ );
			}

			if ( option[ optE::output_top_n_new_decoys ].user() ) {
				Size n_to_output = option[ optE::output_top_n_new_decoys ];
				utility::vector1< Size > top_decoy_inds( n_to_output, 0 );
				utility::arg_least_several( new_decs_scores, top_decoy_inds );
				for ( Size ii = 1; ii <= top_decoy_inds.size(); ++ii ) {
					new_decs[ top_decoy_inds[ ii ] ].dump_pdb( "workdir_" + to_string( MPI_rank_ ) +
						"/" + structure_data->tag() + "_" + to_string( outer_loop_counter_ ) + "_"
						+ to_string( ii ) + ".pdb" );
				}
			}
		}
	}

	for ( auto
			iter = decoy_discrim_data_->position_data_begin(),
			iter_end = decoy_discrim_data_->position_data_end();
			iter != iter_end; ++iter ) {
		optE_data_->add_position_data( *iter );
	}

}

///
SingleStructureDataOP
IterativeOptEDriver::single_structure_data_for_pose(
	core::scoring::ScoreFunctionOP scorefxn,
	core::pose::Pose & pose,
	core::pose::Pose const & crystal_native,
	utility::vector1< Real > & free_data, // scratch space; avoids new
	utility::vector1< Real > & fixed_data, // scratch space; avoids new
	std::string const & structure_tag
) const
{
	(*scorefxn)( pose );

	for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
		free_data[ kk ] = pose.energies().total_energies()[ free_score_list_[ kk ] ];
	}
	for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
		fixed_data[ kk ] = pose.energies().total_energies()[ fixed_score_list_[ kk ] ];
	}
	SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
	core::Real native_rms = core::scoring::CA_rmsd( crystal_native, pose );
	ssd->rms( native_rms );
	ssd->tag( structure_tag );
	return ssd;
}

///
/// @details Returns 1 if added as intended, 0 if not added at all, and -1 if added in the opposite of its intention.
///
AddStatus
IterativeOptEDriver::add_structure_based_on_rms(
	SingleStructureDataOP ssd,
	PNatStructureOptEDataOP structure_data,
	bool intended_native
) const
{
	bool const ramp_nativeness( basic::options::option[ basic::options::OptionKeys::optE::ramp_nativeness ] );


	if ( intended_native ) {
		if ( ramp_nativeness ) {
			if ( ssd->rms() > PNatStructureOptEData::nativeness_high() ) {
				structure_data->add_decoy( ssd );
				return ADDED_STRUCTURE_OPPOSITE_AS_INTENDED;
			}
			/// else
			structure_data->add_native( ssd );
			return ADDED_STRUCTURE_AS_INTENDED;
		}

		if ( option[ optE::max_rms_from_native ].user() ) {
			if ( ssd->rms() > option[ optE::max_rms_from_native ]() ) {
				if ( option[ optE::min_decoy_rms_to_native ]() ) {
					if ( ssd->rms() > option[ optE::min_decoy_rms_to_native ]() ) {
						structure_data->add_decoy( ssd );
						return ADDED_STRUCTURE_OPPOSITE_AS_INTENDED;
						//TR << "Excluding decoy " << decoy_pdb_names[ jj ] << " with rms: " << decoy_rms << std::endl;
					}
				}
				return DID_NOT_ADD_STRUCTURE; // Do not count this structure as a native.
			}
		}
		structure_data->add_native( ssd );
	} else {

		if ( ramp_nativeness ) {
			if ( ssd->rms() > PNatStructureOptEData::nativeness_high() ) {
				structure_data->add_decoy( ssd );
				return ADDED_STRUCTURE_AS_INTENDED;
			}
			/// else
			structure_data->add_native( ssd );
			return ADDED_STRUCTURE_OPPOSITE_AS_INTENDED;
		}


		if ( option[ optE::min_decoy_rms_to_native ].user() ) {
			if ( ssd->rms() <  option[ optE::min_decoy_rms_to_native ]() ) {
				if ( option[ optE::max_rms_from_native ].user() ) {
					if ( ssd->rms() < option[ optE::max_rms_from_native ]() ) {
						structure_data->add_native( ssd );
						return ADDED_STRUCTURE_OPPOSITE_AS_INTENDED;
					}
				}
				//TR << "Excluding decoy " << decoy_pdb_names[ jj ] << " with rms: " << decoy_rms << std::endl;
				return DID_NOT_ADD_STRUCTURE; // Do not treat this structure as a decoy
			}
		}
		//std::cout << "decoy rms: " << structure_data->tag() << " " << decoy_rms << std::endl;
		structure_data->add_decoy( ssd );
	}

	return ADDED_STRUCTURE_AS_INTENDED;
}

///
/// @brief
/// Ligand stuff.
////
void
IterativeOptEDriver::compute_rotamers_around_ligands()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::rotamer_set;

	ScoreFunctionOP scorefxn = create_unweighted_scorefunction();

	for ( Size n=1; n <= ligand_repack_pdbs_.size(); ++n ) {
		std::string const & native_filename( ligand_repack_pdbs_[n] );

		core::pose::Pose native_pose;//, pose;
		if ( outer_loop_counter_ == 1 ) {
			core::import_pose::pose_from_file( native_pose, native_filename , core::import_pose::PDB_file);
			//pose = native_pose;
			ligand_repack_native_poses_.push_back( native_pose );
			//ligand_repack_context_poses_.push_back( native_pose );
		} else {
			//pose = ligand_repack_context_poses_[ n ];
			native_pose = ligand_repack_native_poses_[ n ];
		}

		//(*scorefxn)( pose );
		(*scorefxn)( native_pose );

		utility::file::FileName natname( native_filename );
		core::pose::Pose context_pose;
		/// Should we use the previously repacked pose, or the native pose to gather data from?
		//if ( option[ optE::recover_nat_rot ] ) {
		// context_pose = rotamer_recovery_context_poses_[ n ];
		//} else {
		context_pose = native_pose;
		//}

		// Only include protein residues within 6A of touching the ligand
		utility::vector1<bool> include_rsd( context_pose.size(), false );
		int const jump_id = context_pose.num_jump(); // assume ligand is last jump
		ObjexxFCL::FArray1D_bool is_upstream ( context_pose.size(), false );
		context_pose.fold_tree().partition_by_jump( jump_id, is_upstream );
		for ( core::Size i = 1, i_end = context_pose.size(); i <= i_end; ++i ) {
			// Nothing on ligand side can move
			if ( ! is_upstream(i) ) continue;
			// on protein side, have to do distance check
			core::conformation::Residue const & prot_rsd = context_pose.residue(i);
			if ( ! prot_rsd.is_protein() ) continue;
			for ( core::Size j = 1, j_end = context_pose.size(); j <= j_end; ++j ) {
				if ( is_upstream(j) ) continue; // compare against only ligand residues
				core::conformation::Residue const & lig_rsd = context_pose.residue(j);
				for ( core::Size k = 1, k_end = lig_rsd.nheavyatoms(); k <= k_end; ++k ) {
					double dist2 = lig_rsd.xyz(k).distance_squared( prot_rsd.xyz(prot_rsd.nbr_atom()) );
					double cutoff = prot_rsd.nbr_radius() + 6.0;
					if ( dist2 <= cutoff * cutoff ) {
						include_rsd[i] = true;
						goto END_LIGRES_LOOP; // C++ lacks multi-level break  :(
					}
				}
			}
			END_LIGRES_LOOP: ; // compiler needs ; as a no-op before end of loop
		}

		get_nat_rot_opte_data(
			natname.base(), context_pose, native_pose, include_rsd, *scorefxn,
			free_score_list_, fixed_score_list_, *optE_data_ );
	}
}

///
/// @brief
/// Ligand stuff.
////
void
IterativeOptEDriver::collect_ligand_discrimination_data()
{
	using namespace core::pack::rotamer_set;

	if ( ligand_discrim_data_ == nullptr ) {

		ligand_discrim_data_ = OptEDataOP( new OptEData );

		/// Refactor this
		ScoreFunctionOP scorefxn = create_unweighted_scorefunction();

		//std::string scorelog_name( "workdir_" + to_string( MPI_rank_ ) + "/decdisc_scores.dat" );
		//std::ofstream scorelog( scorelog_name.c_str() );

		for ( Size ii = 1; ii <= ligand_native_decoy_pairs_.size(); ++ii ) {
			PNatLigPoseOptEDataOP structure_data( new PNatLigPoseOptEData );

			{//scope
				utility::file::FileName cryst_fname( ligand_crystal_natives_[ ii ] );
				//structure_data->tag( cryst_fname.base() ); // trim path and extension data -- beautiful utility
				structure_data->tag( ligand_crystal_natives_[ ii ] ); // I need full path name for debugging
			}

			utility::vector1< std::string > native_pdb_names, decoy_pdb_names;
			//std::cout << " PROC #" << MPI_rank_ << " reading pdb lists " ;
			//std::cout << native_decoy_pairs_[ ii ].first << " and ";
			//std::cout << native_decoy_pairs_[ ii ].second << std::endl;

			std::ifstream native_pdblist( ligand_native_decoy_pairs_[ ii ].first.c_str() );
			if ( native_pdblist.bad() ) utility_exit_with_message("Cannot open file "+ligand_native_decoy_pairs_[ ii ].first);
			while ( native_pdblist ) {
				std::string native_pdb;
				native_pdblist >> native_pdb;
				if ( native_pdb != "" ) native_pdb_names.push_back( native_pdb );
			}

			std::ifstream decoy_pdblist( ligand_native_decoy_pairs_[ ii ].second.c_str() );
			if ( decoy_pdblist.bad() ) utility_exit_with_message("Cannot open file "+ligand_native_decoy_pairs_[ ii ].second);
			while ( decoy_pdblist ) {
				std::string decoy_pdb;
				decoy_pdblist >> decoy_pdb;
				if ( decoy_pdb != "" ) decoy_pdb_names.push_back( decoy_pdb );
			}
			//std::cout << native_pdb_names.size() << " natives and " << decoy_pdb_names.size() << " decoys" << std::endl;

			if ( native_pdb_names.size() == 0 ) {
				TR << "[rank " << MPI_rank_ << "] No native entries in " << ligand_native_decoy_pairs_[ ii ].first << "; skipping ligand decoy discrimination for this target." << std::endl;
				continue;
			}
			if ( decoy_pdb_names.size() == 0 ) {
				TR << "[rank " << MPI_rank_ << "] No decoy entries in " << ligand_native_decoy_pairs_[ ii ].second << "; skipping ligand decoy discrimination for this target." << std::endl;
				continue;
			}

			core::pose::Pose crystal_native;
			//if ( option[ in::file::centroid_input ] ) {
			// core::import_pose::centroid_pose_from_pdb( crystal_native, crystal_natives_[ ii ] , core::import_pose::PDB_file);
			//} else {
			core::import_pose::pose_from_file( crystal_native, ligand_crystal_natives_[ ii ] , core::import_pose::PDB_file);
			//}

			utility::vector1< Real > free_data( free_score_list_.size() );
			utility::vector1< Real > fixed_data( fixed_score_list_.size() );
			Size first_total_residue( 0 );
			for ( Size jj = 1; jj <= native_pdb_names.size(); ++jj ) {
				//std::cout << " PROC #" << MPI_rank_ << " reading pdb: #" << jj << " " << native_pdb_names[ jj ] << std::endl;
				/// read the pdb into a pose
				core::pose::Pose pose;
				//if ( option[ in::file::centroid_input ] ) {
				// core::import_pose::centroid_pose_from_pdb( pose, native_pdb_names[ jj ] , core::import_pose::PDB_file);
				//} else {
				core::import_pose::pose_from_file( pose, native_pdb_names[ jj ] , core::import_pose::PDB_file);
				//}

				if ( jj == 1 ) {
					structure_data->set_total_residue( pose.size() );
					first_total_residue = pose.size();
				} else if ( first_total_residue != pose.size() ) {
					std::cerr << "Warning [node " << MPI_rank_ << "]: total_residue for " << native_pdb_names[ jj ];
					std::cerr << " not equal to native #1 total_residue: " << first_total_residue << " vs ";
					std::cerr << pose.size() << std::endl;
					continue;
				}

				///*Real score = */(*scorefxn)( pose );
				//scorelog << "Decoy Discrimination NATIVE " << native_pdb_names[ jj ] << " " << score << "\n";

				EnergyMap emap = score_ligand_interface(*scorefxn, pose);
				for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
					free_data[ kk ] = emap[ free_score_list_[ kk ] ];
				}
				for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
					fixed_data[ kk ] = emap[ fixed_score_list_[ kk ] ];
				}
				SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
				structure_data->add_native( ssd );
				//std::cout << "Adding native, size = " << structure_data->size() << std::endl;
			}

			for ( Size jj = 1; jj <= decoy_pdb_names.size(); ++jj ) {

				/// read the pdb into a pose
				core::pose::Pose pose;
				//std::cout << " PROC #" << MPI_rank_ << " reading pdb: #" << jj << " " << decoy_pdb_names[ jj ] << std::endl;
				//if ( option[ in::file::centroid_input ] ) {
				// core::import_pose::centroid_pose_from_pdb( pose, decoy_pdb_names[ jj ] , core::import_pose::PDB_file);
				//} else {
				core::import_pose::pose_from_file( pose, decoy_pdb_names[ jj ] , core::import_pose::PDB_file);
				//}

				if ( first_total_residue != pose.size() ) {
					std::cerr << "Warning [node " << MPI_rank_ << "]: total_residue for " << decoy_pdb_names[ jj ];
					std::cerr << " not equal to native #1 total_residue: " << first_total_residue << " vs ";
					std::cerr << pose.size() << std::endl;
					continue;
				}

				///*Real score = */ (*scorefxn)( pose );
				//scorelog << "Decoy Discrimination DECOY " << decoy_pdb_names[ jj ] << " " << score << "\n";

				//if ( option[ optE::min_decoy_rms_to_native ].user() ) {
				// Real decoy_rms = core::scoring::CA_rmsd( crystal_native, pose );
				// if ( decoy_rms <  option[ optE::min_decoy_rms_to_native ]() ) {
				//  //TR << "Excluding decoy " << decoy_pdb_names[ jj ] << " with rms: " << decoy_rms << std::endl;
				//  continue;
				// }
				//}

				EnergyMap emap = score_ligand_interface(*scorefxn, pose);
				for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
					free_data[ kk ] = emap[ free_score_list_[ kk ] ];
				}
				for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
					fixed_data[ kk ] = emap[ fixed_score_list_[ kk ] ];
				}
				SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
				if ( structure_data->size() != 0 ) structure_data->add_decoy( ssd );
				//std::cout << "Adding decoy, size = " << structure_data->size() << std::endl;
			}
			ligand_discrim_data_->add_position_data( structure_data );
		}
	}
	for ( auto
			iter = ligand_discrim_data_->position_data_begin(),
			iter_end = ligand_discrim_data_->position_data_end();
			iter != iter_end; ++iter ) {
		optE_data_->add_position_data( *iter );
	}
}

///
/// @brief
/// Ligand stuff.
////
core::scoring::EnergyMap
IterativeOptEDriver::score_ligand_interface( core::scoring::ScoreFunction const & scorefxn, core::pose::Pose & pose )
{
	// For the plain old total score:
	//scorefxn(pose);
	//return pose.energies().total_energies();

	// For the interface score:
	int const jump_id = pose.num_jump();
	core::pose::Pose split_pose( pose ); // make a copy
	protocols::rigid::RigidBodyTransMover trans_mover( split_pose, jump_id );
	// Default direction is to move centroids apart
	trans_mover.step_size(500); // make sure they're fully separated!
	trans_mover.apply( split_pose );
	scorefxn(pose);
	EnergyMap emap( pose.energies().total_energies() ); // make a copy
	scorefxn(split_pose);
	EnergyMap const & smap = split_pose.energies().total_energies();
	//std::cout << "delta_scores";
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		ScoreType i = (ScoreType) ii;
		emap[i] -= smap[i];
		//if( emap[i] != 0 ) std::cout << " " << name_from_score_type(i) << " " << emap[i];
	}
	//std::cout << "\n";
	return emap;
}

///
/// @brief
/// For a calling master node, collects the rotamer energies that were calculated on a slave CPU.
////
void
IterativeOptEDriver::collect_rotamer_energies_from_slave_cpu
(
#ifdef USEMPI
	Size const which_cpu
#else
	Size const
#endif
)
{
#ifdef USEMPI
	using namespace core::pack::rotamer_set;

	//std::cout << "Node 0 preparing to receive from " << which_cpu << std::endl;

	/// 1. Sanity: check that "boolean" vector of the free and fixed energy terms
	/// from the source cpu matches the boolean vectors we are expecting
	int * free_energy_terms = new int[ core::scoring::n_score_types ];
	int * fixed_energy_terms = new int[ core::scoring::n_score_types ];

	//std::cout << " PROC #" << MPI_rank_ << "collect_rotamer_energies_from_slave_cpu: " << which_cpu << std::endl;

	MPI_Recv( free_energy_terms, n_score_types, MPI_INT, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );
	MPI_Recv( fixed_energy_terms, n_score_types, MPI_INT, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );

	//std::cout << " PROC #" << MPI_rank_ << "received free and fixed from " << which_cpu << std::endl;

	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( free_energy_terms[  ii - 1 ] != (int) ( free_parameters_[  (ScoreType) ii ] != 0.0 )) {
			std::cerr << "Free energy term mismatch! " << ScoreType( ii ) << " " << free_energy_terms[ ii - 1 ] << " & " << free_parameters_[  (ScoreType) ii ]  << std::endl;
			utility_exit_with_message( "Free energy term on Node 0 does not match free energy term remotely");
		}
		if ( fixed_energy_terms[ ii - 1 ] != (int) ( fixed_parameters_[ (ScoreType) ii ] != 0.0 )) {
			std::cerr << "Fixed energy term mismatch! " << ScoreType( ii ) << " " << fixed_energy_terms[ ii - 1 ] << " & " << fixed_parameters_[  (ScoreType) ii ]  << std::endl;
			utility_exit_with_message( "Free energy term on Node 0 does not match fixed energy term remotely");
		}
	}

	delete [] free_energy_terms;   free_energy_terms = 0;
	delete [] fixed_energy_terms; fixed_energy_terms = 0;

	/// 2. Number of positions on which OptE data has been gathered remotely
	Size n_pos;
	MPI_Recv( & n_pos, 1, MPI_UNSIGNED_LONG, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );
	//std::cout << " PROC #" << MPI_rank_ << "npos from " << which_cpu << " " << n_pos << std::endl;

	for ( Size ii = 1; ii <= n_pos; ++ii ) {
		//std::cout << "Waiting to receive position data type" << std::endl;
		int position_data_type;
		MPI_Recv( & position_data_type, 1, MPI_INT, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );

		OptEPositionDataOP ii_data = OptEPositionDataFactory::create_position_data( (OptEPositionDataType) position_data_type );
		//std::cout << "Node 0 about to receive from node " << which_cpu << std::endl;
		ii_data->receive_from_node( which_cpu, tag_ );

		optE_data_->add_position_data( ii_data );
		//std::cout << "Node 0 added position data object " << std::endl;
	}
#endif

}


Size IterativeOptEDriver::num_outer_iterations() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return option[ optE::n_design_cycles ];
}

Size IterativeOptEDriver::num_inner_iterations() const { return 6; }

///
/// @brief
/// Calls the method initialize_free_and_fixed() which reads in the files free and fixed and sets the EnergyMap vectors
/// free_parameters_ and fixed_parameters_.  Also here the reference energies array gets init'd.  Finally,
/// setup_derived_free_and_fixed_data gets called which
///
void IterativeOptEDriver::intialize_free_and_fixed_energy_terms() {

	initialize_free_and_fixed( free_parameters_, fixed_parameters_ );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		if ( MPI_rank_ == 0 ) {
			//for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			// if ( free_parameters_[ (ScoreType) ii ] != 0.0 ) {
			//  free_parameters_[ (ScoreType) ii ] = numeric::random::rg().uniform() + 0.01; // random non-zero starting point
			// }
			//}
			Real const rpp_refs[20] = {
				0.16, 1.7, -0.67, -0.81, 0.63, -0.17, 0.56, 0.24, -0.65, -0.1,
				-0.34, -0.89, 0.02, -0.97, -0.98, -0.37, -0.27, 0.29, 0.91, 0.51
				};

			for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
				before_minimization_reference_energies_[ ii ] = rpp_refs[ ii - 1 ];
			}
			std::cout << "INITIALIZED before_minimization_reference_energies_ REFERENCE ENERGIES" << std::endl;
		}
	}
	setup_derived_free_and_fixed_data();
}

///
/// @brief
/// Optimizes weights using either a standard minimizer or a ParticleSwarmMinimizer (which is significantly better!).
///
/// Assuming the ParticleSwarmMinimizer is used...
/// Each particle traverses weight space and comes up with a set of values, then it evaluates the fitness function
/// (which basically calls all of the appropriate underlying get_score() methods) to see how good the new weights are.
/// The particles then update their direction and velocity and come up with a new set of weights to evaluate. This is
/// done with a set number of particles and a set number of cycles.  Currently using on the order of 100 particles and
/// 20 cycles. (No tests have been done to determine the optimum number of particles and/or cycles, although ronj likes
/// to use more particles whenever possible.)
///
void IterativeOptEDriver::optimize_weights()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::optimization;

	barrier();

	bool const optimize_in_parallel( option[ optE::mpi_weight_minimization ] );
	if ( MPI_rank_ != 0 && !optimize_in_parallel ) {
		return; // do nothing if we're not the master cpu AND we're not distributing weight optimization over multiple cpus.
	}

	if ( MPI_rank_ == 0 ) TR << "optimize_weights(): entered method" << std::endl;

	//attach_debugger();

	free_weights_before_minimization_ = free_parameters_;
	if ( option[ optE::wrap_dof_optimization ].user() ) {
		minimizer_dofs_before_minimization_ = minimizer_dofs_mixed_;
	}

	/// Skip weight optimization in the context of the native protein if the
	/// "design first" flag is on the command line, but only after
	/// free_weights_before_mininmization_ has been updated.
	if ( outer_loop_counter_ == 1 && option[ optE::design_first ].user() ) {
		TR << "optimize_weights(): design_first flag in use. node " << MPI_rank_ << " leaving method." << std::endl;
		return;
	}

	if ( option[ optE::starting_refEs ].user() && outer_loop_counter_ == 1 ) {
		std::cout << "READING REFERENCE ENERGIES FROM FILE" << std::endl;
		std::cout << "FREE COUNT: " << free_count_ << std::endl;
		after_minimization_reference_energies_ = read_reference_energies_from_file( option[ optE::starting_refEs ] );
		before_minimization_reference_energies_ = after_minimization_reference_energies_;
		//std::cout << "after_minimization_reference_energies_: ";
		//for ( Size ii = 1; ii <= after_minimization_reference_energies_.size(); ++ii ) {
		// std::cout << after_minimization_reference_energies_[ ii ] << " ";
		//}
		//std::cout << std::endl;
		if ( free_count_ == 0 ) {
			TR << "optimize weights early exit" << std::endl;
			return; // don't bother trying to optimize the weights
		}
	}

	// Create an "optE data minimizer" object
	// Do the actual weight optimization
	OptEMultifuncOP opt_min_ptr( new OptEMultifunc(
		*optE_data_, fixed_parameters_,
		(int) free_count_,
		free_score_list_,
		fixed_score_list_,
		before_minimization_reference_energies_,
		component_weights_ ) );
	OptEMultifunc & opt_min( * opt_min_ptr );


	if ( option[ optE::fit_reference_energies_to_aa_profile_recovery ] && option[ optE::dont_use_reference_energies ] ) {
		utility_exit_with_message("optimize_weights(): can't do fitting of reference energies to profile recovery when 'dont_use_reference_energies' flag in use.");
	}

	if ( option[ optE::fit_reference_energies_to_aa_profile_recovery ] && !(outer_loop_counter_ == 1 && inner_loop_counter_ == 1) ) {
		opt_min.fix_reference_energies( true );
	} else if ( option[ optE::starting_refEs ].user() && outer_loop_counter_ == 1 ) {
		opt_min.fix_reference_energies( true );
	}

	if ( MPI_rank_ == 0 ) {

		clock_t starttime = clock();

		/// Set include terms according to current weight set
		include_terms_ = free_parameters_;
		include_terms_ += fixed_parameters_;

		//TR << " about to start getting dofs " <<std::endl;
		optimization::Multivec start_dofs  = opt_min.get_dofs_from_energy_map( include_terms_ );
		Real start_fitness = opt_min( start_dofs );
		Real end_fitness = 0.0;

		TR << "optimize_weights(): objective function start: " << start_fitness << " dofs: [ ";
		for ( Size ii = 1; ii <= start_dofs.size(); ++ii ) { TR << F( 8,4,start_dofs[ii] ) << ", "; }
		TR << "]" << std::endl;

		MultifuncOP opt_min2;

		Size ndofs = start_dofs.size();
		if ( option[ optE::wrap_dof_optimization ].user() ) {
			if ( outer_loop_counter_ == 1 ) {
				wrapped_opt_min_ = WrapperOptEMultifuncOP( new WrapperOptEMultifunc() );
				wrapped_opt_min_->init(
					free_score_list_, (int) free_count_,
					fixed_score_list_, fixed_parameters_,
					opt_min_ptr );
				minimizer_dofs_before_minimization_.resize( wrapped_opt_min_->n_real_dofs(), 0.0 );
				minimizer_dofs_after_minimization_.resize( wrapped_opt_min_->n_real_dofs(), 0.0 );
				minimizer_dofs_mixed_.resize( wrapped_opt_min_->n_real_dofs(), 0.0 );
			} else {
				wrapped_opt_min_->set_multifunc( opt_min_ptr );
			}
			opt_min2 = wrapped_opt_min_;
			ndofs = wrapped_opt_min_->n_real_dofs();
		} else {
			opt_min2 = opt_min_ptr;
		}

		if ( option[ optE::optimize_starting_free_weights ] ) {

			using namespace core::optimization;
			/// High tolerance -- don't over minimize before pswarm gets to explore
			optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone_atol", 1, true, false, false );
			optimization::Minimizer minimizer( *opt_min2, options );

			/// Low tolerance -- drill down!
			optimization::MinimizerOptions options2( "lbfgs_armijo_nonmonotone", 1e-4, true, false, false );
			optimization::Minimizer minimizer2( *opt_min2, options2 );

			/// Lowest tolerance -- drill down!
			optimization::MinimizerOptions options3( "lbfgs_armijo_nonmonotone", 1e-9, true, false, false );
			optimization::Minimizer minimizer3( *opt_min2, options3 );

			// create two Multivec's (vector1 Real) that will hold the minimum and maximum weights a Particle can achieve
			// sometimes the particles leave the weight on fa_intra_rep as zero since the score12 weight is so low. this messes
			// things up down the line when sending things via MPI because only non-zero values are sent to speed up the run.
			// if the particles leave the weight at zero, reset it to a very small value so it's still present.
			// or just enforce a minimum weight of 0.0001 on the particles. -ronj
			Multivec min( ndofs, 0.001), max( ndofs, 5.0);

			ParticleSwarmMinimizer psm(min, max);

			// run 20 cycles of the PS algorithm, using opt_min2 as the fitness function and some number of particles
			// the option has a default value of 100, if the user doesn't specify it.
			TR_VERBOSE << "optimize_weights(): node " << MPI_rank_ << " beginning round 1 of PSO." << std::endl;
			ParticleOPs particles = psm.run( option[ optE::number_of_swarm_cycles ], *opt_min2, option[ optE::number_of_swarm_particles ] );

			//   psm.print_particles( particles, "optimize_weights(): round 1" );

			for ( core::Size i = 1; i <= particles.size(); ++i ) {
				Real start( (*opt_min2)( particles[i]->p_ ) );
				minimizer.run( particles[i]->p_ );
				TR.Trace << "first round minimization: " << i << " start: " << start <<  ", end: " << (*opt_min2)( particles[ i ]->p_ ) << std::endl;
			}
			psm.run( 0, *opt_min2, particles );
			//   psm.print_particles( particles, "optimize_weights(): round 1 post min" );

			TR_VERBOSE << "optimize_weights(): node " << MPI_rank_ << " beginning round 2 of PSO using [minimized] round 1 particles." << std::endl;

			psm.run( option[ optE::number_of_swarm_cycles ], *opt_min2, particles );
			//   psm.print_particles( particles, "optimize_weights(): round 2" );

			clock_t min_starttime = clock();
			//for(core::Size ii = 1; ii <= particles.size(); ++ii) {
			for ( core::Size ii = 1; ii <= 15; ++ii ) { // minimizing is really slow so only do some particles!
				Real start( (*opt_min2)( particles[ii]->p_ ) );
				TR.Trace << "starting minimization of particle " << ii << std::endl;
				minimizer2.run( particles[ii]->p_ );
				TR.Trace << "second round minimization: " << ii << " start: " << F( 9,5,start ) << ", end: " << F( 9,5,(*opt_min2)( particles[ ii ]->p_ )) << std::endl;
				std::cout.flush();
			}
			clock_t min_stoptime = clock();
			TR_VERBOSE << "optimize_weights(): particle minimization took " << ((double)min_stoptime-min_starttime) / CLOCKS_PER_SEC << " seconds." << std::endl;

			// This will re-sort the particles and update pbest_
			psm.run(0, *opt_min2, particles);
			//   psm.print_particles( particles, "optimize_weights(): round 2 post min" );

			ParticleOP p = particles[1];
			TR_VERBOSE << "optimize_weights(): best particle fitness: " << F( 9,5,-1.0 * p->fitness_pbest() ) << " dofs: [";
			for ( core::Size j=1; j <= p->pbest().size(); ++j ) { TR_VERBOSE << F(8,4,p->pbest()[j]) << ", "; }
			TR_VERBOSE << " ]" << std::endl;

			start_dofs = particles[1]->pbest();
			//TR << "Final round of gradient-based minimization, score before: " << -1 * particles[1]->fitness_pbest() << std::endl;
			//minimizer3.run( start_dofs ); //ronj don't use minimization when using the new unfolded state energy
			//TR << "Final round of gradient-based minimization, score after: " << (*opt_min2)( start_dofs ) << std::endl;

			end_fitness = (*opt_min2)( start_dofs );

			if ( ( end_fitness > start_fitness ) && option[ optE::repeat_swarm_optimization_until_fitness_improves ] ) {
				// try one more round with a new set of particles that have starting values closer to the start dofs

				// run 40 cycles of the PS algorithm, using opt_min2 as the fitness function and some number of particles
				// the option has a default value of 100, if the user doesn't specify it.
				ParticleOPs particles = psm.run( 2 * option[ optE::number_of_swarm_cycles ], *opt_min2, option[ optE::number_of_swarm_particles ],
					opt_min.get_dofs_from_energy_map( free_weights_before_minimization_ ) );

				//    psm.print_particles( particles, "optimize_weights(): round extra innings" );

				clock_t min_starttime = clock();
				for ( core::Size ii = 1; ii <= 10; ++ii ) { // minimizing is really slow so only do some of the particles
					Real start( (*opt_min2)( particles[ii]->p_ ) );
					TR_VERBOSE << "starting minimization of particle " << ii << std::endl;
					//minimizer2.run( particles[ii]->p_ );
					minimizer.run( particles[ii]->p_ );
					TR_VERBOSE << "second round minimization: " << ii << " start: " << F( 9,5,start ) << ", end: " << F( 9,5,(*opt_min2)( particles[ ii ]->p_ )) << std::endl;
					std::cout.flush();
				}
				clock_t min_stoptime = clock();
				TR_VERBOSE << "optimize_weights(): particle minimization took " << ((double)min_stoptime-min_starttime) / CLOCKS_PER_SEC << " seconds." << std::endl;

				// re-score and re-sort the particles to get the best scoring one
				psm.run( 0, *opt_min2, particles );
				//    psm.print_particles( particles, "optimize_weights(): round extra innings post min" );

				ParticleOP p = particles[1];
				TR_VERBOSE << "optimize_weights(): best particle fitness: " << F( 9,5,-1.0 * p->fitness_pbest() ) << " dofs: [";
				for ( core::Size j=1; j <= p->pbest().size(); ++j ) { TR_VERBOSE << F(8,4,p->pbest()[j]) << ", "; }
				TR_VERBOSE << " ]" << std::endl;

				start_dofs = particles[1]->pbest();
			}


			if ( option[ optE::wrap_dof_optimization ].user() ) {
				minimizer_dofs_after_minimization_ = start_dofs;
				/// From this point forward start_dofs needs to be the size the OptEMultifunc expects;
				start_dofs = wrapped_opt_min_->derived_dofs( minimizer_dofs_after_minimization_ );
			}

		} else {
			// somebody is crazy and not using the Particle Swarm to do weight space exploration...
			if ( using_unfolded_energy_term_ ) {
				TR << "optimize_weights(): minimization not recommended when using unfolded state energy" << std::endl;
			}
			optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1e-9, true, false, false );
			optimization::Minimizer minimizer( opt_min, options );
			minimizer.run( start_dofs );
		}

		TR << "optimize_weights(): end: " << ( *opt_min2 )( start_dofs ) << ", dofs: [ ";
		for ( Size ii = 1; ii <= start_dofs.size(); ++ii ) { TR << F(8,4,start_dofs[ii]) << ", "; }
		TR << " ]" << std::endl;

		if ( option[ optE::wrap_dof_optimization ].user() ) {
			TR << "Wrapped weights after minimization" << std::endl;
			wrapped_opt_min_->print_dofs( minimizer_dofs_after_minimization_, TR );
			TR << std::endl;
		}

		test_weight_sensitivity(opt_min, start_dofs);

		// set the after_min refE vector
		if ( ! option[ optE::dont_use_reference_energies ].user() ) {
			//for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {   // 20 should not be hardcoded here!
			// after_minimization_reference_energies_[ ii ] = start_dofs[ free_count_+ii ]; // save them in non-negated form
			//}
			// apparently, now, the opt_min object can be queried to get the reference energy values instead of what the
			// particles get out
			after_minimization_reference_energies_ = opt_min.get_reference_energies_from_dofs( start_dofs );
		}

		// create an EnergyMap from the most recent set of DOFs, but make sure to reset the energies for the fixed terms to 0.
		free_weights_after_minimization_ = opt_min.get_energy_map_from_dofs( start_dofs );
		for ( Size ii = 1 ; ii <= fixed_score_list_.size(); ++ii ) {
			free_weights_after_minimization_[ fixed_score_list_[ ii ]] = 0; // reset fixed parameters
		}

		// initialize some vectors that will later be used by the PositionData print/get_score methods
		optimization::Multivec vars( free_count_ + after_minimization_reference_energies_.size(), 0.0 );
		optimization::Multivec dE_dvars( free_count_ + after_minimization_reference_energies_.size(), 0.0 );

		Size num_energy_dofs( free_count_ );
		int num_ref_dofs( after_minimization_reference_energies_.size() );
		int num_total_dofs( num_energy_dofs + num_ref_dofs );
		scoring::EnergyMap fixed_terms = fixed_parameters_;
		scoring::ScoreTypes score_list( free_score_list_ );
		scoring::ScoreTypes fixed_score_list( fixed_score_list_ );

		// set the vars Mulitvec to contain the free weights and reference weights
		for ( Size ii = 1; ii <= free_score_list_.size(); ++ii ) {
			vars[ ii ] = free_weights_after_minimization_[ free_score_list_[ ii ] ] ;
			TR_VERBOSE << "optimize_weights(): free weights before/after minimization_: [ " << name_from_score_type( free_score_list_[ii] ) << " ]: "
				<< F(8,4,free_weights_before_minimization_[ free_score_list_[ ii ] ]) << " -> "
				<< F(8,4,free_weights_after_minimization_[ free_score_list_[ ii ] ]) << std::endl;
		}

		if ( ! option[ optE::dont_use_reference_energies ].user() ) {
			for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
				vars[ ii + free_count_ ] = after_minimization_reference_energies_[ ii ];
			}
		}

		/// NOW write all the score data to the log.  This means first grabbing the data from distant nodes,
		/// if that data is not on the master node already
		/// There's no reason, of course, why the optE data should have to live on a single CPU all at once...
		clock_t stoptime = clock();

		TR << "optimize_weights(): optimization took " << ((double) stoptime-starttime)/CLOCKS_PER_SEC << " seconds." << std::endl;

		// release hold of the multifunc that the wrapper wraps, but
		/// hold on to the wrapper itself for later reuse.  No more calls to func() after this point.
		if ( option[ optE::wrap_dof_optimization ].user() ) {
			wrapped_opt_min_->set_multifunc( nullptr );
		}

		/// Tell remote nodes that function evaulation is over.
		if ( optimize_in_parallel ) {
			opt_min.declare_minimization_over(); // tell remote nodes to stop waiting for func/dfuncs to evaluate
			collect_rotamer_energies_from_slave_cpus();
		}

		///////////////////// Log the results after minimization is complete ////////////////////////////////

		/// 2. open log file.
		std::string logname = "logdir/minimization_dat_" + utility::to_string( outer_loop_counter_ ) + ".dat";
		std::ofstream outlog( logname.c_str() );

		/// 3. score all the position data
		utility::vector1< Real > cumulative_score_list(n_optE_data_types,0.0);
		//int ii = 1;
		for ( auto itr = optE_data_->position_data_begin(),
				e_itr = optE_data_->position_data_end(); itr != e_itr ; ++itr  ) {
			(*itr)->print_score( outlog, component_weights_, vars, dE_dvars,
				num_energy_dofs, num_ref_dofs, num_total_dofs,
				fixed_terms, score_list, fixed_score_list );
			cumulative_score_list[(*itr)->type()] +=  (*itr)->get_score(
				component_weights_, vars, dE_dvars,
				num_energy_dofs, num_ref_dofs, num_total_dofs,
				fixed_terms, score_list, fixed_score_list );
			//++ii;
			//if ( ii % 10 == 0 ) {
			// //TR << ".";
			//}
		}
		for ( Size ii=1; ii <= n_optE_data_types; ii++ ) {
			TR << "optimize_weights(): energy component: " <<  OptEPositionDataFactory::optE_type_name( OptEPositionDataType( ii ))
				<< " " << cumulative_score_list[ii] << std::endl;
		}


		/// 4. prints some extra information about the range of energies observed to the minimization data file

		// db //
		utility::vector1< EnergyMap > rawE_min( n_optE_data_types );
		utility::vector1< EnergyMap > rawE_max( n_optE_data_types );
		Real const faux_max( -1234 ); Real const faux_min( 1234 );
		for ( Size ii = 1; ii <= n_optE_data_types; ++ii ) {
			//for ( EnergyMap::iterator iter = rawE_min[ ii ].begin(); iter != rawE_min[ ii ].end(); ++iter ) {
			// *iter = faux_min;
			//}
			//for ( EnergyMap::iterator iter = rawE_max[ ii ].begin(); iter != rawE_max[ ii ].end(); ++iter ) {
			// *iter = faux_max;
			//}
			for ( Size jj = 1; jj <= n_score_types; ++jj ) {
				rawE_min[ ii ][ (ScoreType) jj ] = faux_min;
				rawE_max[ ii ][ (ScoreType) jj ] = faux_max;
			}

		}

		for ( auto itr = optE_data_->position_data_begin(),
				e_itr = optE_data_->position_data_end(); itr != e_itr ; ++itr  ) {
			(*itr)->range( score_list, fixed_score_list,
				rawE_min[(*itr)->type()], rawE_max[(*itr)->type()] );
		}

		ScoreTypes free_and_fixed( score_list );
		for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
			free_and_fixed.push_back( fixed_score_list[ ii ] );
		}
		std::sort(free_and_fixed.begin(), free_and_fixed.end() );

		for ( Size ii = 1; ii <= n_optE_data_types; ++ii ) {
			outlog << "DATA RANGE: ";
			outlog << OptEPositionDataFactory::optE_type_name( OptEPositionDataType( ii ));
			outlog << " ";
			for ( Size jj = 1; jj <= free_and_fixed.size(); ++jj ) {
				if ( rawE_min[ ii ][ free_and_fixed[ jj ] ] > rawE_max[ ii ][ free_and_fixed[ jj ] ] ) continue;
				outlog << "( " << name_from_score_type( free_and_fixed[ jj ]) << ", ";
				outlog << rawE_min[ ii ][ free_and_fixed[ jj ] ] << ", ";
				outlog << rawE_max[ ii ][ free_and_fixed[ jj ] ] << " ) ";
			}
			outlog << "\n";
		}
	} else {
		opt_min.wait_for_remote_vars();
		/// once minimization has completed, send this data to the master cpu so it can output it to the log.
		/// There's no reason, of course, why the optE data should have to live on a single CPU all at once...
		send_rotamer_energies_to_master_cpu();
	}
	TR << "leaving optimize_weights()" << std::endl;
}


void
IterativeOptEDriver::test_weight_sensitivity(
	OptEMultifunc const & func,
	core::optimization::Multivec const & dofs
) const
{
	std::string logname = "weightdir/sensitivity_" + utility::to_string( outer_loop_counter_ ) + ".dat";
	std::ofstream out( logname.c_str() );
	Real const minval = func(dofs);
	out << "Minimum function value " << minval << std::endl;
	// Later DOFs include the AA ref energies, so we don't want ALL of them:
	//for(Size dof_idx = 1; dof_idx <= dofs.size(); ++dof_idx) {
	for ( Size dof_idx = 1; dof_idx <= free_score_list_.size(); ++dof_idx ) {
		Real maxval = minval;
		out << "term_" << dof_idx << " " << free_score_list_[dof_idx];
		for ( Real val = 0.0; val <= 2.0; val += 0.1 ) {
			//for(Real scale = 0.0; scale <= 2.0; scale += 0.5) {
			//if(scale == 1.0) continue;
			core::optimization::Multivec dofs_copy = dofs;
			//dofs_copy[dof_idx] *= scale;
			dofs_copy[dof_idx] = val;
			Real const newval = func(dofs_copy);
			maxval = std::max(maxval, newval);
			out << "  " << (newval - minval);
		}
		out << "  maxDelta " << (maxval - minval) << std::endl;
	}
	out.close();
}

utility::vector1< Real >
IterativeOptEDriver::read_reference_energies_from_file( std::string const & fname ) const
{
	utility::vector1< Real > reference_energies( chemical::num_canonical_aas, 0.0 );

	utility::io::izstream weight_file( fname.c_str() );
	bool read_refEs = false;
	while ( weight_file ) {
		std::string tag;
		weight_file >> tag;
		if ( tag == "METHOD_WEIGHTS" ) {
			weight_file >> tag; // "ref"
			if ( tag != "ref" ) continue;
			for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
				Real aa_refE;
				weight_file >> aa_refE;
				if ( !weight_file ) break;
				reference_energies[ ii ] = aa_refE;
			}
			read_refEs = true;
		}
	}
	if ( ! read_refEs ) {
		utility_exit_with_message( "Failed to read METHOD_WEIGHTS from file " + fname + ".  IterativeOptEDriver::read_reference_energies_from_file() " );
	}
	return reference_energies;

}

void
IterativeOptEDriver::score_position_data()
{
	using namespace core::pack::rotamer_set;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	runtime_assert( option[ optE::rescore::weights ].user() );

	bool const optimize_in_parallel( option[ optE::mpi_weight_minimization ] );

	if ( MPI_rank_ != 0 ) {
		//opt_min.wait_for_remote_vars();
		/// once minimization has completed, send this data to the master cpu so it can output it to the log.
		/// There's no reason, of course, why the optE data should have to live on a single CPU all at once...
		if ( optimize_in_parallel ) {
			send_rotamer_energies_to_master_cpu();
		} else {
			return;
		}

	} else {

		/// Tell remote nodes that function evaulation is over.
		if ( optimize_in_parallel ) {
			//opt_min.declare_minimization_over(); // tell remote nodes to stop waiting for func/dfuncs to evaluate
			collect_rotamer_energies_from_slave_cpus();
		}

		/// 1. read in the weight set and reference energies
		/// and initialize the mutliVects needed to call OptEPositionData::print_score()
		utility::vector1< Real > reference_energies( chemical::num_canonical_aas, 0.0 );
		EnergyMap weight_set;

		std::ifstream weight_file( option[ optE::rescore::weights ]()().c_str() );
		while ( weight_file ) {
			std::string tag;
			weight_file >> tag;
			if ( tag == "METHOD_WEIGHTS" ) {
				weight_file >> tag; // "ref"
				for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
					Real aa_refE;
					weight_file >> aa_refE;
					reference_energies[ ii ] = aa_refE;
				}
			} else if ( ScoreTypeManager::is_score_type( tag ) ) {
				ScoreType st = ScoreTypeManager::score_type_from_name( tag );
				Real weight;
				weight_file >> weight;
				TR << "Reading weight: " << st << " " << weight << std::endl;
				weight_set[ st ] = weight;
			} else {
				TR << "Warning:: ignoring tag " << tag << " from " << option[ optE::rescore::weights ]()() << std::endl;
			}

		}

		/// The following variables are input parameters to the print_score function in optEPositionData
		optimization::Multivec vars( free_count_ + before_minimization_reference_energies_.size(), 0.0 );
		optimization::Multivec dE_dvars( free_count_ + before_minimization_reference_energies_.size(), 0.0 );

		Size num_energy_dofs( free_count_ );
		int num_ref_dofs( before_minimization_reference_energies_.size() );
		int num_total_dofs( num_energy_dofs + num_ref_dofs );

		scoring::EnergyMap fixed_terms;
		scoring::ScoreTypes score_list( free_score_list_ );
		scoring::ScoreTypes fixed_score_list( fixed_score_list_ );

		for ( Size ii = 1; ii <= score_list.size(); ++ii ) {
			TR_VERBOSE << "setting vars " << ii << " to " << weight_set[ free_score_list_[ ii ] ] << std::endl;
			vars[ ii ] = weight_set[ free_score_list_[ ii ] ];
			// also set the free_parameters Map in case the user isn't optimizing weights but just rescoring a weight set
			free_parameters_[ free_score_list_[ ii ] ] = weight_set[ free_score_list_[ ii ] ];
		}

		if ( ! option[ optE::dont_use_reference_energies ].user() ) {
			TR_VERBOSE << "setting reference energies ";
			for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
				TR_VERBOSE << ii << " " << reference_energies[ ii ] << ", ";
				vars[ ii + free_count_ ] = reference_energies[ ii ];
				after_minimization_reference_energies_[ ii ] = reference_energies[ ii ];
			}
			TR_VERBOSE << std::endl;
		}

		for ( Size ii = 1; ii <= fixed_score_list_.size(); ++ii ) {
			TR_VERBOSE << "setting fixed term " << ii << " " << fixed_score_list_[ ii ] << " " << weight_set[ fixed_score_list_[ ii ] ] << std::endl;
			fixed_terms[ fixed_score_list_[ ii ] ] = weight_set[ fixed_score_list_[ ii ] ];
			// also set the fixed_parameters Map in case the user isn't optimizing weights but just rescoring a weight set
			fixed_parameters_[ fixed_score_list_[ ii ] ] = weight_set[ fixed_score_list_[ ii ] ];
		}

		/// 2. open log file.
		std::ofstream outlog( option[ optE::rescore::outlog ]()().c_str() );
		TR_VERBOSE << "log file opened, scoring position data" << std::endl;

		/// 3. score all the position data
		utility::vector1< Real > cumulative_score_list( n_optE_data_types, 0.0 );
		//int ii = 1;
		for ( auto itr = optE_data_->position_data_begin(), e_itr = optE_data_->position_data_end() ; itr != e_itr ; ++itr  ) {
			(*itr)->print_score( outlog, component_weights_, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, score_list, fixed_score_list );
			cumulative_score_list[(*itr)->type()] +=  (*itr)->get_score( component_weights_, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, score_list, fixed_score_list );
			//++ii;
			//if ( ii % 10 == 0 ) {
			// //TR << ".";
			//}
		}
		for ( Size ii=1; ii <= n_optE_data_types; ii++ ) {
			if ( cumulative_score_list[ii] != 0 ) {
				TR_VERBOSE << "optimization function energy component: " << OptEPositionDataFactory::optE_type_name( OptEPositionDataType( ii ))
					<< " " << cumulative_score_list[ii] << std::endl;
			}
		}

		/// 4. close log file.
		/// noop -- happens automatically when the ofstream goes out of scope.
	}
}

///
/// @brief
/// send new score file via MPI instead of writing to disk.
///
/// @remarks
/// The reference energy term is automatically added to the end of the scorefile with a weight of 1, so don't put it as
/// a free or fixed param or it will be counted twice.  The actual reference energies per AA are varied during the
/// protocol and also written out here.
///
void IterativeOptEDriver::write_new_scorefile()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	optE_data_.reset(); /// clear memory!

	if ( MPI_rank_ == 0 ) {

		mixing_factor_ = opte_weight_mixing_factor( outer_loop_counter_, inner_loop_counter_ );
		Real alpha = 1 - mixing_factor_;

		// here we begin a significant branching of write_new_scorefile
		// we have four possibilities: writing with and without consideration of sequence entropy and writing
		// with and without use of the wrapped dof minimizer

		if ( option[ optE::fit_reference_energies_to_aa_profile_recovery ] ) {
			using namespace core::chemical;

			// the below is too confusing and should be against mini coding guidelines (-ronj)
			//Real const fa_atr_weight(
			// fixed_parameters_[ fa_atr ] == 0.0 ?
			// (free_weights_inner_loop_[ fa_atr ] == 0.0 ?  1.0 : free_weights_inner_loop_[ fa_atr ] )
			// : fixed_parameters_[ fa_atr ] );

			Real fa_atr_weight;
			if ( fixed_parameters_[ fa_atr ] == 0.0 ) {
				if ( free_weights_inner_loop_[ fa_atr ] == 0.0 ) {
					fa_atr_weight = 1.0;
				} else {
					fa_atr_weight = free_weights_inner_loop_[ fa_atr ];
				}
			} else {
				fa_atr_weight = fixed_parameters_[ fa_atr ];
			}

			//// Maybe I'll add flags for these so that I can tweak them from
			//// the command line.
			Real const pretty_close( 1.10 );
			Real const way_off( 1.5 );
			Real const way_way_off( 3.0 );
			Real const nowhere_close( 10.0 );
			Real const biggest_step( 0.8 * fa_atr_weight );
			Real const bigger_step( 0.4 * fa_atr_weight );
			Real const big_step( 0.1 * fa_atr_weight );
			Real const small_step( 0.05 * fa_atr_weight );
			Real const tiny_step( 0.01 * fa_atr_weight );

			if ( inner_loop_counter_ == 1 && ( outer_loop_counter_ == 1 ||  free_count_ != 0 ) ) {
				/// all the way, baby!
				///free_weights_inner_loop_ = free_weights_after_minimization_;

				if ( option[ optE::wrap_dof_optimization ].user() ) {
					/// Interpolate the DOFs, reconstruct the weight sets from the dofs.
					if ( ! wrapped_opt_min_ ) {
						utility_exit_with_message( "ERROR in IterativeOptEDriver::write_new_scorefile(); wrapped_opt_min_ is NULL");
					}
					for ( Size kk = 1; kk <= minimizer_dofs_before_minimization_.size(); ++kk ) {
						minimizer_dofs_mixed_[ kk ] =
							alpha * minimizer_dofs_before_minimization_[ kk ] +
							mixing_factor_ * minimizer_dofs_after_minimization_[ kk ];
					}
					utility::vector1< Real > vars = wrapped_opt_min_->derived_dofs( minimizer_dofs_mixed_ );
					free_weights_and_refEs_from_vars( vars, free_weights_inner_loop_, reference_energies_inner_loop_ );
					if ( ! option[ optE::no_design ] ) {
						TR << "Wrapped weights round: " << outer_loop_counter_ << ", " << inner_loop_counter_ << std::endl;
						wrapped_opt_min_->print_dofs( minimizer_dofs_mixed_, TR );
						TR << std::endl;
					}

					reference_energies_inner_loop_ = after_minimization_reference_energies_;
					// moving this call up to here, b/c if wrap dof not in use, the reference energies get set correctly

				} else {

					if ( outer_loop_counter_ == 1 && option[ optE::design_first ].user() ) { mixing_factor_ = 0.0; alpha = 1.0; }
					if ( option[ optE::rescore::measure_sequence_recovery ].user() ) { mixing_factor_ = 1.0; alpha = 0.0; } // the right wts are in after_min vector
					if ( option[ optE::optimize_ddGmutation ].user() && ! ( option[ optE::optimize_nat_aa ].user() ) ) { mixing_factor_ = 1.0; alpha = 0.0; }

					for ( Size kk = 1; kk <= n_score_types; ++kk ) {
						free_weights_inner_loop_[ (ScoreType) kk ] = alpha * free_weights_before_minimization_[ (ScoreType) kk ];
						free_weights_inner_loop_[ (ScoreType) kk ] += mixing_factor_ * free_weights_after_minimization_[ (ScoreType) kk ];
					}

					// output to the terminal the weights before and after mixing
					for ( Size ii = 1; ii <= n_score_types; ++ii ) {
						if ( free_weights_inner_loop_[ (ScoreType) ii ] != 0.0 ) {
							TR_VERBOSE << "write_new_scorefile(): free weights before/after mixing: [ " << name_from_score_type( (ScoreType) ii ) << " ]: "
								<< F(8,4,free_weights_after_minimization_[ (ScoreType) ii ]) << " -> "
								<< F(8,4,free_weights_inner_loop_[ (ScoreType) ii ]) << std::endl;
						}
					}

					for ( Size kk = 1; kk <= before_minimization_reference_energies_.size(); ++kk ) {
						reference_energies_inner_loop_[ kk ] = ( 1.0 - mixing_factor_ ) * before_minimization_reference_energies_[ kk ];
						reference_energies_inner_loop_[ kk ] += mixing_factor_ * after_minimization_reference_energies_[ kk ];
					}

					TR_VERBOSE << "write_new_scorefile(): reference energies after mixing: ";
					for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
						TR_VERBOSE << core::chemical::oneletter_code_from_aa( core::chemical::AA( ii ) ) << ":" << F(5,2,reference_energies_inner_loop_[ ii ]) << ", ";
					}
					TR_VERBOSE << std::endl;

					// special section for the unfolded state energy term; this block of code could be removed if the wrapped multifunc
					// optimizer was being used, but I added this prior to APL's release of the wrapped multifunc optimizer and I'm used
					// to running things this way. (-ronj)
					output_weighted_unfolded_energies();
				}

			} else { // inner loop counter is not 1 or we have no other free terms besides refEs
				TR << "Tuning reference energies using amino acid profile recovery data: round " <<
					outer_loop_counter_ << " " << inner_loop_counter_ << std::endl;
				TR << "write_new_scorefile(): reference energies before entropy: " << inner_loop_counter_ << ": ";
				for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
					TR << core::chemical::oneletter_code_from_aa( core::chemical::AA( ii ) ) << ":" << F(5,2,reference_energies_inner_loop_[ ii ]) << ", ";
				}
				TR << std::endl;

				/// Only modify the reference energies so that the designed sequence profile
				/// matches the experimentally observed sequence profile.
				for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
					if ( aa_freq_exp_[ ii ] != 0 ) {
						if ( aa_freq_obs_[ ii ] > aa_freq_exp_[ ii ] ) {
							// if the designed freq is greater than what's observed in nature, INCREASE the values on the reference energies
							if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] > nowhere_close ) {
								reference_energies_inner_loop_[ ii ] += biggest_step;
							} else if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] > way_way_off ) {
								reference_energies_inner_loop_[ ii ] += bigger_step;
							} else if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] > way_off ) {
								reference_energies_inner_loop_[ ii ] += big_step;
							} else if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] > pretty_close ) {
								reference_energies_inner_loop_[ ii ] += small_step;
							} else {
								reference_energies_inner_loop_[ ii ] += tiny_step; // we're < 1.1% away, take a tiny step
							}
						} else {
							// if the designed freq is less than what's observed in nature, DECREASE the values on the reference energies
							if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] < ( 1.0 / nowhere_close ) ) {
								reference_energies_inner_loop_[ ii ] -= biggest_step;
							} else if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] < ( 1.0 / way_way_off ) ) {
								reference_energies_inner_loop_[ ii ] -= bigger_step;
							} else if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] < ( 1.0 / way_off ) ) {
								reference_energies_inner_loop_[ ii ] -= big_step;
							} else if ( aa_freq_obs_[ ii ] / aa_freq_exp_[ ii ] < ( 1.0 / pretty_close ) ) {
								reference_energies_inner_loop_[ ii ] -= small_step;
							} else {
								reference_energies_inner_loop_[ ii ] -= tiny_step;
							}
						}
					}
				}
			}

			TR << "write_new_scorefile(): reference energies after entropy:  " << inner_loop_counter_ << ": ";
			for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
				TR << core::chemical::oneletter_code_from_aa( core::chemical::AA( ii ) ) << ":" << F(5,2,reference_energies_inner_loop_[ ii ]) << ", ";
			}
			TR << std::endl;

			Real total_refE( 0.0 );
			for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) total_refE += reference_energies_inner_loop_[ ii ];
			for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) reference_energies_inner_loop_[ ii ] -= total_refE / num_canonical_aas;

			TR << "write_new_scorefile(): ref energies after normalization:  " << inner_loop_counter_ << ": ";
			for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) { TR << F(5,2,reference_energies_inner_loop_[ ii ]) << ", "; }
			TR << std::endl;

		} else {  // not fitting reference energies to a profile

			if ( option[ optE::wrap_dof_optimization ].user() ) {
				/// Interpolate the DOFs, reconstruct the weight sets from the dofs.
				if ( ! wrapped_opt_min_ ) {
					utility_exit_with_message( "ERROR in IterativeOptEDriver::write_new_scorefile(); wrapped_opt_min_ is NULL");
				}
				for ( Size kk = 1; kk <= minimizer_dofs_before_minimization_.size(); ++kk ) {
					minimizer_dofs_mixed_[ kk ] =
						alpha * minimizer_dofs_before_minimization_[ kk ] +
						mixing_factor_ * minimizer_dofs_after_minimization_[ kk ];
				}
				utility::vector1< Real > vars = wrapped_opt_min_->derived_dofs( minimizer_dofs_mixed_ );
				free_weights_and_refEs_from_vars( vars, free_weights_inner_loop_, reference_energies_inner_loop_ );
				if ( ! option[ optE::no_design ] ) {
					TR << "Wrapped weights round: " << outer_loop_counter_ << ", " << inner_loop_counter_ << std::endl;
					wrapped_opt_min_->print_dofs( minimizer_dofs_mixed_, TR );
					TR << std::endl;
				}

			} else { // don't use the wrapped multifunc optimizer

				if ( outer_loop_counter_ == 1 && option[ optE::design_first ].user() ) { mixing_factor_ = 0.0; alpha = 1.0; }
				if ( option[ optE::rescore::measure_sequence_recovery ].user() ) { mixing_factor_ = 1.0; alpha = 0.0; } // the right wts are in after_min vector
				if ( option[ optE::optimize_ddGmutation ].user() && ! ( option[ optE::optimize_nat_aa ].user() ) ) { mixing_factor_ = 1.0; alpha = 0.0; }

				for ( Size kk = 1; kk <= n_score_types; ++kk ) {
					free_weights_inner_loop_[ (ScoreType) kk ] = alpha * free_weights_before_minimization_[ (ScoreType) kk ];
					free_weights_inner_loop_[ (ScoreType) kk ] += mixing_factor_ * free_weights_after_minimization_[ (ScoreType) kk ];
				}

				// output to the terminal the weights before and after mixing
				for ( Size ii = 1; ii <= n_score_types; ++ii ) {
					if ( free_weights_inner_loop_[ (ScoreType) ii ] != 0.0 ) {
						TR_VERBOSE << "write_new_scorefile(): free weights before/after mixing: [ " << name_from_score_type( (ScoreType) ii ) << " ]: "
							<< F(8,4,free_weights_after_minimization_[ (ScoreType) ii ]) << " -> "
							<< F(8,4,free_weights_inner_loop_[ (ScoreType) ii ]) << std::endl;
					}
				}

				if ( ! option[ optE::dont_use_reference_energies ].user() ) {
					for ( Size kk = 1; kk <= before_minimization_reference_energies_.size(); ++kk ) {
						reference_energies_inner_loop_[ kk ] = ( 1.0 - mixing_factor_ ) * before_minimization_reference_energies_[ kk ];
						reference_energies_inner_loop_[ kk ] += mixing_factor_ * after_minimization_reference_energies_[ kk ];
					}
				}

				// special section for the unfolded state energy term; this block of code could be removed if the wrapped multifunc
				// optimizer was being used, but I added this prior to APL's release of the wrapped multifunc optimizer and I'm used
				// to running things this way. (-ronj)
				output_weighted_unfolded_energies();

			}
		}

#ifdef USEMPI
		/// Send the weights we've just computed to the other nodes
		Real * free_wts = new Real[ n_score_types ];
		for ( Size ii = 0; ii < n_score_types; ++ii ) {
			free_wts[ ii ] = free_weights_inner_loop_[ ScoreType( ii + 1 ) ];
		}

		// if refE's not in use, this code will create arrays of size 0 - that's fine
		Size n_ref_Es = reference_energies_inner_loop_.size();
		Real * ref_Es = new Real[ n_ref_Es ];
		for ( Size ii = 0; ii < n_ref_Es; ++ii ) {
			ref_Es[ ii ] = reference_energies_inner_loop_[ ii + 1 ];
		}

		for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
			MPI_Send( free_wts, n_score_types, MPI_DOUBLE, ii, tag_, MPI_COMM_WORLD );
			if ( ! option[ optE::dont_use_reference_energies ].user() ) {
				MPI_Send( & n_ref_Es, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				MPI_Send( ref_Es, n_ref_Es, MPI_DOUBLE, ii, tag_, MPI_COMM_WORLD );
			}
		}

		delete [] ref_Es; ref_Es = 0;
		delete [] free_wts; free_wts = 0;

	} else {
		/// MPI_rank_ != 0; receive energies from the master node.
		Real * free_wts = new Real[ n_score_types ];
		MPI_Recv( free_wts, n_score_types, MPI_DOUBLE, 0, tag_, MPI_COMM_WORLD, &stat_ );

		for ( Size ii = 0; ii < n_score_types; ++ii ) {
			free_weights_inner_loop_[ ScoreType( ii + 1 ) ] = free_wts[ ii ];
		}
		delete [] free_wts; free_wts = 0;

		// don't bother with the reference energies if the user doesn't want them
		if ( ! option[ optE::dont_use_reference_energies ].user() ) {

			Size n_ref_Es( 0 );
			MPI_Recv( & n_ref_Es, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, &stat_ );
			Real * ref_Es = new Real[ n_ref_Es ];
			MPI_Recv( ref_Es, n_ref_Es, MPI_DOUBLE, 0, tag_, MPI_COMM_WORLD, &stat_ );

			reference_energies_inner_loop_.resize( n_ref_Es );
			for ( Size ii = 0; ii < n_ref_Es; ++ii ) {
				reference_energies_inner_loop_[ ii + 1 ] = ref_Es[ ii ];
			}

			delete[] ref_Es; ref_Es = 0;
		}
#endif // USEMPI
	}

	if ( MPI_rank_ == 0  ) {
		/// for posterity.
		EnergyMap combined_weights( fixed_parameters_ );
		combined_weights += free_weights_inner_loop_;

		std::string scorefile_name = get_scorefile_name();
		// Make sure directory exists:
		utility::file::create_directory_recursive( utility::file::PathName(scorefile_name).parent() );
		std::ofstream fout( scorefile_name.c_str() );

		/// Ensure the score file includes soft rep if its requested.
		if ( option[ optE::optE_soft_rep ].user() ) {
			fout << "ETABLE FA_STANDARD_SOFT" << std::endl;
		}

		if ( ! option[ optE::dont_use_reference_energies ].user() ) {
			/// Weight file output:
			fout << "METHOD_WEIGHTS ref ";
			for ( Size ii = 1; ii <= reference_energies_inner_loop_.size(); ++ii ) {
				fout << reference_energies_inner_loop_[ ii ] << " ";
			}
			fout << "\n";
		}
		if ( option[ optE::no_hb_env_dependence ] ) {
			fout << "NO_HB_ENV_DEP\n";
		}

		for ( Size ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
			if ( combined_weights[ ScoreType( ii ) ] != 0 ) {
				fout << name_from_score_type( ScoreType( ii ) ) << " " << combined_weights[ ScoreType( ii ) ] << "\n";
			}
		}

		if ( ! option[ optE::dont_use_reference_energies ].user() ) {
			fout << "ref 1\n"; // DONT FORGET TO USE THE REFERENCE ENERGIES YOU JUST CALCULATED!
		}

		fout.close();
	}

}

///
/// @brief Multiply out the unweighted unfolded energies with the current free and fixed term weights
///
/// @details
/// Dots the unweighted, unfolded energies with the current free and fixed term weights and prints out the weighted
/// unfolded energies to stdout.
///
void
IterativeOptEDriver::output_weighted_unfolded_energies() {

	utility::vector1< Real > wtd_unfE( chemical::num_canonical_aas, 0.0 );

	if ( using_unfolded_energy_term_ ) {

		UnfoldedStatePotential const & unfE_potential = ScoringManager::get_instance()->get_UnfoldedStatePotential( scoring::UNFOLDED_SCORE12 );
		utility::vector1< EnergyMap > unweighted_unfolded_emap( chemical::num_canonical_aas );
		for ( Size aa=1; aa <= chemical::num_canonical_aas; ++aa ) {
			unweighted_unfolded_emap[ aa ].zero();
			unfE_potential.raw_unfolded_state_energymap( chemical::name_from_aa( (chemical::AA) aa ), unweighted_unfolded_emap[ aa ] );
		}

		for ( Size aa = 1; aa <= chemical::num_canonical_aas; ++aa ) {

			Real unfolded_energy_for_one_aa = 0.0;
			Real weighted_unfolded_energy_for_one_aa = 0.0;

			// free weights first
			for ( Size ii = 1; ii <= n_score_types; ++ii ) {
				if ( free_weights_inner_loop_[ (ScoreType) ii ] != 0.0 ) {
					TR_VERBOSE << "output_weighted_unfolded_energies(): adding unfolded energy for aa '" << chemical::name_from_aa( (chemical::AA) aa )
						<< "' unweighted free '" << name_from_score_type( (ScoreType)ii ) << "' energy: "
						<< unweighted_unfolded_emap[ aa ][ (ScoreType)ii ] << " * '" << name_from_score_type( (ScoreType)ii )
						<< "' weight: " << free_weights_inner_loop_[ (ScoreType)ii ]
						<< " = " << unweighted_unfolded_emap[ aa ][ (ScoreType)ii ] * free_weights_inner_loop_[ (ScoreType)ii ]
						<< std::endl;
					unfolded_energy_for_one_aa += (unweighted_unfolded_emap[ aa ][ (ScoreType)ii ] * free_weights_inner_loop_[ (ScoreType)ii ]);
				}
			}

			// then fixed weights
			for ( Size ii = 1; ii <= fixed_score_list_.size(); ++ii ) {
				if ( fixed_parameters_[ fixed_score_list_[ ii ] ] != 0.0 ) {
					TR_VERBOSE << "output_weighted_unfolded_energies(): adding unfolded energy for aa '" << chemical::name_from_aa( (chemical::AA) aa )
						<< "' unweighted fixed '" << name_from_score_type( fixed_score_list_[ ii ] ) << "' energy: "
						<< unweighted_unfolded_emap[ aa ][ fixed_score_list_[ ii ]] << " * '"
						<<  name_from_score_type( fixed_score_list_[ ii ] )
						<< "' weight: " << fixed_parameters_[ fixed_score_list_[ ii ] ]
						<< " = " << unweighted_unfolded_emap[ aa ][ fixed_score_list_[ ii ]] * fixed_parameters_[ fixed_score_list_[ ii ] ]
						<< std::endl;
					unfolded_energy_for_one_aa += (unweighted_unfolded_emap[ aa ][ fixed_score_list_[ ii ] ] * fixed_parameters_[ fixed_score_list_[ ii ] ]);
				}
			}

			if ( free_weights_inner_loop_[ unfolded ] == 0.0 && fixed_parameters_[ unfolded ] == 0.0 ) {
				TR << "output_weighted_unfolded_energies(): unfolded term has no weight! using 1.0 to avoid errors." << std::endl;
				weighted_unfolded_energy_for_one_aa = unfolded_energy_for_one_aa * 1.0;
			} else if ( free_weights_inner_loop_[ unfolded ] != 0.0 ) {
				// unfolded term weight is variable
				TR << "output_weighted_unfolded_energies(): weighting unfolded energy '" << unfolded_energy_for_one_aa << "' by unfolded term weight: '"
					<< free_weights_inner_loop_[ unfolded ] << "' gives weighted unfolded energy for one aa of "
					<< F(4,2,unfolded_energy_for_one_aa * free_weights_inner_loop_[ unfolded ]) << std::endl;
				weighted_unfolded_energy_for_one_aa = unfolded_energy_for_one_aa * free_weights_inner_loop_[ unfolded ];
			} else if ( fixed_parameters_[ unfolded ] != 0.0 ) {
				TR << "output_weighted_unfolded_energies(): weighting unfolded energy '" << unfolded_energy_for_one_aa << "' by unfolded term weight: '"
					<< fixed_parameters_[ unfolded ] << "' gives weighted unfolded energy for one aa of "
					<< F(4,2,unfolded_energy_for_one_aa * fixed_parameters_[ unfolded ]) << std::endl;
				weighted_unfolded_energy_for_one_aa = unfolded_energy_for_one_aa * fixed_parameters_[ unfolded ];
			} else {
				TR << "output_weighted_unfolded_energies(): error with checking the weight of the unfolded term. using 1.0 to avoid errors." << std::endl;
				weighted_unfolded_energy_for_one_aa = unfolded_energy_for_one_aa * 1.0;
			}

			wtd_unfE[ aa ] = weighted_unfolded_energy_for_one_aa;
		}

		TR << "output_weighted_unfolded_energies(): weighted unfoldedE by aa: [ ";
		for ( Size aa = 1; aa <= chemical::num_canonical_aas; ++aa ) {
			TR << wtd_unfE[ aa ] << " ";
		}
		TR << std::endl;
	}
}

///
/// @brief Expand free variables and combine with fixed to make an Energy Map
///
/// @details This dofs Multivec is the list of weights that the OptEMultifunc
/// sees; do not confuse this set of dofs with the set of dofs that the
/// Minimizer and the WrappedOptEMultifunc use.
///
core::scoring::EnergyMap
IterativeOptEDriver::free_terms_energy_map_from_dofs(
	core::optimization::Multivec const & dofs
) const
{
	EnergyMap return_map;

	// This covers the variable weights
	Size dof_index( 1 );
	for ( auto itr : free_score_list_ ) {
		return_map[ itr ] = dofs[ dof_index++ ];
	}

	return return_map;
}

///
void
IterativeOptEDriver::free_weights_and_refEs_from_vars(
	utility::vector1< Real > const & vars,
	core::scoring::EnergyMap & weights,
	utility::vector1< Real > & reference_energies
) const
{
	// conditional needed because if reference_energies is accessed when refE's are not being optimized, errors will occur
	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			reference_energies[ ii ] = vars[ free_count_ + ii ]; // save them in non-negated form
		}
	}
	weights.zero();
	weights = free_terms_energy_map_from_dofs( vars );
}

///
/// @brief
/// Sets functional forms (e.g. soft-rep) but doesn't set any weights.
/// If the option -optE::optE_soft_rep is specified, then an empty scorefunction with the FA_STANDARD_SOFT etable
/// is returned.  Another option is -optE::optE_no_protein_fa_elec.  This excludes protein_protein_fa_elec in the
/// the scorefunction energy method options.
///
ScoreFunctionOP
IterativeOptEDriver::configure_new_scorefunction() const
{
	ScoreFunctionOP scorefxn( new ScoreFunction );
	if ( option[ optE::optE_soft_rep ].user() ) {
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.etable_type( FA_STANDARD_SOFT );
		scorefxn->set_energy_method_options( options );
	}
	if ( option[ optE::optE_no_protein_fa_elec ]() ) {
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.exclude_protein_protein_fa_elec( true );
		scorefxn->set_energy_method_options( options );
	}
	if ( option[ optE::no_hb_env_dependence ] ) {
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.hbond_options().use_hb_env_dep( false );
		scorefxn->set_energy_method_options( options );
	}
	if ( option[ optE::no_hb_env_dependence_DNA ] ) {
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.hbond_options().use_hb_env_dep_DNA( false );
		scorefxn->set_energy_method_options( options );
	}
	if ( using_unfolded_energy_term_ ) {
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.unfolded_energies_type( scoring::UNFOLDED_SCORE12 );
		scorefxn->set_energy_method_options( options );
	}
	if ( option[ optE::centroid_rot ] ) {
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.atom_vdw_atom_type_set_name( "centroid_rot" );
		scorefxn->set_energy_method_options( options );
	}
	if ( option[ optE::centroid_rot_min ] ) {
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.atom_vdw_atom_type_set_name( "centroid_rot:min" );
		scorefxn->set_energy_method_options( options );
	}
	return scorefxn;
}

///
/// @brief
/// Takes a std::string and a destination and constructs the MPI_Send call.  Does this include the reference energy
/// term somehow?  I don't believe it does.
///
ScoreFunctionOP
IterativeOptEDriver::create_unweighted_scorefunction() const
{
	ScoreFunctionOP scorefxn = configure_new_scorefunction();
	for ( int i=1 ; i <= n_score_types ; ++i ) {
		if ( include_terms_[ ScoreType(i) ] != 0.0 ) {
			scorefxn->set_weight( ScoreType(i), include_terms_[ ScoreType(i) ] );
		}
	}
	return scorefxn;
}

/// @details Construct a score function: set etable type, set weights, set reference weights.
core::scoring::ScoreFunctionOP
IterativeOptEDriver::create_weighted_scorefunction() const
{
	ScoreFunctionOP sfxn = configure_new_scorefunction();
	for ( int i=1 ; i <= n_score_types ; ++i ) {
		if ( free_weights_inner_loop_[ ScoreType(i) ] != 0.0 ) {
			//std::cout << " PROC #" << MPI_rank_ << " include term: " << ScoreType(i) << std::endl;
			sfxn->set_weight( ScoreType(i), free_weights_inner_loop_[ ScoreType(i) ] );
		} else if ( fixed_parameters_[ ScoreType(i) ] != 0.0 ) {
			sfxn->set_weight( ScoreType(i), fixed_parameters_[ ScoreType(i) ] );
		}
	}
	// sfxn->energy_method_options().set_method_weights( ref, reference_energies_inner_loop_ );
	sfxn->set_method_weights( ref, reference_energies_inner_loop_ );
	sfxn->set_weight( ref, 1.0 );

	return sfxn;
}

///
/// @remarks
/// IMPORTANT IMPORTANT IMPORTANT: requires weightdir having been created before launching the program.
///
std::string
IterativeOptEDriver::get_scorefile_name()
{
	std::stringstream instream;
	instream << outer_loop_counter_;
	std::string scorefile_name = "weightdir/optE_scorefile_" + instream.str() + ".wts";
	return scorefile_name;
}

///
/// @brief
/// Calls run_design on all pdbs and collects the results from slave cpus if MPI is in use.
///
void IterativeOptEDriver::test_sequence_recovery()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ optE::no_design ]() ) {
		exit_gracefully();
	}

	if ( MPI_rank_ == 0 ) {
		run_design_on_assigned_pdbs(); std::cout.flush();
		repack_assigned_pdbs(); std::cout.flush();
		collect_sequence_recovery_data_from_slave_cpus(); std::cout.flush();
		collect_rotamer_recovery_data_from_slave_cpus();
	} else {
		run_design_on_assigned_pdbs(); std::cout.flush();
		repack_assigned_pdbs(); std::cout.flush();
		send_sequence_recovery_data_to_master_cpu(); std::cout.flush();
		send_rotamer_recovery_data_to_master_cpu();
	}
}

///
/// @brief
/// Helper method for master node.  Calls collect_recovery_data_from_slave_cpu on all slave CPUs.
///
void IterativeOptEDriver::collect_sequence_recovery_data_from_slave_cpus() {

	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		collect_sequence_recovery_data_from_slave_cpu( ii );
	}

	TR << "collect_sequence_recovery_data_from_slave_cpus(): collected data from " << total_positions_ << " residues: " << count_recovered_
		<< " native amino acids recovered" << std::endl;
	inner_loop_sequence_recovery_rate_ = (( Real ) count_recovered_ ) / total_positions_;

	TR  << "collect_sequence_recovery_data_from_slave_cpus(): overall sequence recovery rate: " << inner_loop_sequence_recovery_rate_ << std::endl;

	/// This is now a good time to compute the designed frequency and the experimental frequency of the various aa's.
	for ( Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
		aa_freq_obs_[ ii ] = ((Real) aa_obs_[ ii ] ) / total_positions_;
		aa_freq_exp_[ ii ] = ((Real) aa_exp_[ ii ] ) / total_positions_;
	}

	TR << "amino acid counts: observed: ";
	for ( Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
		TR << core::chemical::oneletter_code_from_aa( core::chemical::AA( ii ) ) << ": " << aa_obs_[ ii ] << ", ";
		//if ( ii % 5 == 0 ) TR << std::endl;
	}
	TR << std::endl;

	TR << "amino acid counts: expected: ";
	for ( Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
		TR << core::chemical::oneletter_code_from_aa( core::chemical::AA( ii ) ) << ": " << aa_exp_[ ii ] << ", ";
	}
	TR << std::endl;


	Real cross_entropy( 0.0 );
	TR << "amino acid frequency: obs (exp):";
	for ( Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
		TR << " " << core::chemical::oneletter_code_from_aa( core::chemical::AA( ii ) )
			<< ": " << ((Real) aa_obs_[ ii ] ) / total_positions_ << " (" << ((Real) aa_exp_[ ii ] ) / total_positions_<< ")";
		//if ( ii % 5 == 0 ) TR << std::endl;
		cross_entropy -= (((Real) aa_exp_[ ii ] ) / total_positions_) * std::log(  ((Real) aa_obs_[ ii ] ) / total_positions_ + 1e-5  );
	}
	TR << std::endl;

	//TR << "Cross Entropy: " << cross_entropy << std::endl;
}

///
void IterativeOptEDriver::collect_sequence_recovery_data_from_slave_cpu(
#ifdef USEMPI
	Size const which_cpu
#else
	Size const
#endif
)
{
#ifdef USEMPI
	Size cpu_positions;
	Size cpu_recovered;
	MPI_Recv( & cpu_positions, 1, MPI_UNSIGNED_LONG, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );
	MPI_Recv( & cpu_recovered, 1, MPI_UNSIGNED_LONG, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );

	total_positions_ += cpu_positions;
	count_recovered_ += cpu_recovered;

	Size aa_counts[ core::chemical::num_canonical_aas ];

	/// 1. Send counts of amino acid types coming out of design (observed)
	MPI_Recv( aa_counts, core::chemical::num_canonical_aas, MPI_UNSIGNED_LONG, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );
	for ( Size ii = 1, iim1 = 0; ii <= core::chemical::num_canonical_aas; ++ii, ++iim1 ) aa_obs_[ ii ] += aa_counts[ iim1 ];

	/// 2. Send counts of amino acid types in the input data (expected)
	MPI_Recv( aa_counts, core::chemical::num_canonical_aas, MPI_UNSIGNED_LONG, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );
	for ( Size ii = 1, iim1 = 0; ii <= core::chemical::num_canonical_aas; ++ii, ++iim1 ) aa_exp_[ ii ] += aa_counts[ iim1 ];

#endif
}

///
/// @brief
/// Helper method for master node.  Calls collect_rotamer_recovery_data_from_slave_cpu on all slave CPUs.
///
void IterativeOptEDriver::collect_rotamer_recovery_data_from_slave_cpus()
{
	if ( basic::options::option[ basic::options::OptionKeys::in::file::centroid_input ] ) return;
	if ( ! basic::options::option[ basic::options::OptionKeys::optE::recover_nat_rot ] ) return;

	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		collect_rotamer_recovery_data_from_slave_cpu( ii );
	}
	TR << "collect_rotamer_recovery_data_from_slave_cpus(): collected rotamer recovery data from " << total_rotamer_positions_
		<< " residues: " << count_rotamers_recovered_ << " native amino acids recovered" << std::endl;
	inner_loop_rotamer_recovery_rate_ = (( Real ) count_rotamers_recovered_ ) / total_rotamer_positions_;
	TR  << "collect_rotamer_recovery_data_from_slave_cpus(): rotamer recovery rate: " << inner_loop_rotamer_recovery_rate_ << std::endl;

}

///
void IterativeOptEDriver::collect_rotamer_recovery_data_from_slave_cpu(
#ifdef USEMPI
	Size const which_cpu
#else
	Size const
#endif
)
{
#ifdef USEMPI
	Size cpu_rotamer_positions;
	Size cpu_rotamer_recovered;
	MPI_Recv( & cpu_rotamer_positions, 1, MPI_UNSIGNED_LONG, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );
	MPI_Recv( & cpu_rotamer_recovered, 1, MPI_UNSIGNED_LONG, which_cpu, tag_, MPI_COMM_WORLD, &stat_ );

	total_rotamer_positions_ += cpu_rotamer_positions;
	count_rotamers_recovered_ += cpu_rotamer_recovered;

#endif
}

///
/// @brief
/// Runs design on the pdbs assigned to this node/cpu.
///
void IterativeOptEDriver::run_design_on_assigned_pdbs()
{
	/// NO MORE: don't read from disk, instead, create a score function
	/// based on the *_inner_loop weight/reference energy arrays
	///ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( get_scorefile_name() );

	/// Construct a score function: set etable type, set weights, set reference weights.
	//ScoreFunctionOP sfxn = create_weighted_scorefunction();

	/// Veto'ing APL's code above because users may not want refE's; need to rewrite this section to use the new
	/// configure_new_scorefunction method but then sets the refE weight to zero.

	EnergyMap wts_map;
	ScoreFunctionOP sfxn = configure_new_scorefunction();

	for ( int i=1 ; i <= n_score_types ; ++i ) {
		if ( free_weights_inner_loop_[ ScoreType(i) ] != 0.0 ) {
			sfxn->set_weight( ScoreType(i), free_weights_inner_loop_[ ScoreType(i) ] );
			wts_map[ ScoreType(i) ] = free_weights_inner_loop_[ ScoreType(i) ];
		} else if ( fixed_parameters_[ ScoreType(i) ] != 0.0 ) {
			sfxn->set_weight( ScoreType(i), fixed_parameters_[ ScoreType(i) ] );
			wts_map[ ScoreType(i) ] = fixed_parameters_[ ScoreType(i) ];
		}
	}

	// Adding a special check here for the unfolded state energy term. Like the reference energy term, the unfolded
	// term needs to have extra weights set.  These can be set by passing in a vector1 of Reals of size n_score_types
	// which is the weights for every score term desired in the unfolded state energy.  So what needs to be done is
	// an EnergyMap needs to be created to hold all of the free and fixed weights - and then the emap needs to be
	// converted into a vector1 of Reals.

	// the vector of weights and the set_method_weights function will only get called if 'unfolded' is being used,
	// i.e. if the term has a nonzero weight during the run

	if ( free_weights_inner_loop_[ scoring::unfolded ] != 0.0 || fixed_parameters_[ scoring::unfolded ] != 0.0 ) {
		// obtaining a vector of Reals from an EnergyMap should probably be an EnergyMap method, but whatever.
		utility::vector1< Real > wts_vector( scoring::n_score_types, 0.0 );
		for ( int ii=1; ii < n_score_types; ++ii ) {
			wts_vector[ ii ] = wts_map[ ScoreType(ii) ];
		}
		// set method weights should cause the 'unfolded' EnergyMethod object living inside the score function to be
		// recreated with the weights in the wts_vector.  this should allow the unfolded state term to actually be
		// used during design (unlike in collect_rotamer_energies() where it just returns zero).
		sfxn->set_method_weights( scoring::unfolded, wts_vector );
	}

	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		sfxn->set_method_weights( ref, reference_energies_inner_loop_ );
		sfxn->set_weight( ref, 1.0 );
	}

	if ( MPI_rank_ == 0 ) {
		TR_VERBOSE << "run_design_on_assigned_pdbs(): created scorefxn for running design" << std::endl;
		//sfxn->show( std::cout );
	}

	/// prep arrays so that we can recover data from
	zero_aa_counts();

	measure_sequence_recovery(
		native_pdbs_,
		next_iteration_pdbs_,
		sfxn,
		total_positions_,
		count_recovered_
	);

}

///
void IterativeOptEDriver::repack_assigned_pdbs()
{
	if ( basic::options::option[ basic::options::OptionKeys::in::file::centroid_input ] ) return;
	if ( ! basic::options::option[ basic::options::OptionKeys::optE::recover_nat_rot ] ) return;

	total_rotamer_positions_ = count_rotamers_recovered_ = 0;

	/// NO MORE: don't read from disk, instead, create a score function
	/// based on the *_inner_loop weight/reference energy arrays
	///ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( get_scorefile_name() );


	/// Construct a score function: set etable type, set weights, set reference weights.
	ScoreFunctionOP sfxn = configure_new_scorefunction();
	for ( int i=1 ; i <= n_score_types ; ++i ) {
		if ( free_weights_inner_loop_[ ScoreType(i) ] != 0.0 ) {
			//std::cout << " PROC #" << MPI_rank_ << " include term: " << ScoreType(i) << std::endl;
			sfxn->set_weight( ScoreType(i), free_weights_inner_loop_[ ScoreType(i) ] );
		} else if ( fixed_parameters_[ ScoreType(i) ] != 0.0 ) {
			sfxn->set_weight( ScoreType(i), fixed_parameters_[ ScoreType(i) ] );
		}
	}

	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		// sfxn->energy_method_options().set_method_weights( ref, reference_energies_inner_loop_ );
		sfxn->set_method_weights( ref, reference_energies_inner_loop_ );
		sfxn->set_weight( ref, 1.0 );
	}

	measure_rotamer_recovery(
		native_pdbs_,
		next_iteration_pdbs_,
		sfxn,
		total_rotamer_positions_,
		count_rotamers_recovered_
	);

}

///
void IterativeOptEDriver::send_sequence_recovery_data_to_master_cpu()
{
#ifdef USEMPI
	MPI_Send( & total_positions_, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );
	MPI_Send( & count_recovered_, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );

	Size aa_counts[ core::chemical::num_canonical_aas ];
	/// 1. Send counts of amino acid types coming out of design (observed)
	for ( Size ii = 1, iim1 = 0; ii <= core::chemical::num_canonical_aas; ++ii, ++iim1 ) aa_counts[ iim1 ] = aa_obs_[ ii ];
	MPI_Send( aa_counts, core::chemical::num_canonical_aas, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );

	/// 2. Send counts of amino acid types in the input data (expected)
	for ( Size ii = 1, iim1 = 0; ii <= core::chemical::num_canonical_aas; ++ii, ++iim1 ) aa_counts[ iim1 ] = aa_exp_[ ii ];
	MPI_Send( aa_counts, core::chemical::num_canonical_aas, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );

#endif

}

///
void IterativeOptEDriver::send_rotamer_recovery_data_to_master_cpu()
{
	if ( basic::options::option[ basic::options::OptionKeys::in::file::centroid_input ] ) return;
	if ( ! basic::options::option[ basic::options::OptionKeys::optE::recover_nat_rot ] ) return;
#ifdef USEMPI
	MPI_Send( & total_rotamer_positions_, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );
	MPI_Send( & count_rotamers_recovered_, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );
#endif
}

///
/// @brief
/// The final function call of the go() method.  After all the pdbs have been designed and repacked, and the recovery
/// data collected, decide if this new set of weights improved the sequence recovery.
///
bool IterativeOptEDriver::decide_if_sequence_recovery_improved()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	int accept_new_weight_set( 0 ); // 0 == weight set rejected; 1 == weight set accepted

	if ( MPI_rank_ != 0 ) {
#ifdef USEMPI
		MPI_Recv( & accept_new_weight_set, 1, MPI_INT, 0, tag_, MPI_COMM_WORLD, & stat_ );
#endif
	} else {

		if ( option[ optE::fit_reference_energies_to_aa_profile_recovery ] ) {

			Real const Wcross_ent = -0.1;
			// typically x-entropy is ~3. typical seq. rec is ~0.30. make x-entropy 2x as important as seq. recovery

			if ( outer_loop_counter_ <= 2 && inner_loop_counter_ == 4 ) { /// perform 4 rounds the first AND SECOND time thru the outer loop
				accept_new_weight_set = 1;
			}

			// but after this, accept a weight set if the sequence recovery and the entropy improve - don't force a certain number of
			// inner iterations before accepting a weight set
			// the weight on the entropy term below decided how important entropy is relative to overall sequence recovery
			// the weight is negative so that we can MAXIMIZE the cross entropy (when typically you minimize cross entropy)
			// there's two ways to make the decision whether or not to accept a weight set
			// 1) add seq recovery and wtd cross entropy and accept if the sum is better
			// 2) accept only if seq recovery is better and wtd cross entropy is better
			//
			// Reverting my change to something new: force 3 runs of the outer loop before we accept a weight set. the entropy just
			// doesn't kick in until the later rounds.

			Real cross_entropy( 0.0 );
			for ( Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
				cross_entropy -= aa_freq_exp_[ ii ] * std::log( aa_freq_obs_[ ii ] + 1e-5 );
			}
			Real weighted_cross_entropy = cross_entropy * Wcross_ent;

			TR << "decide_if_sequence_recovery_improved(): inner loop recovery rate: " << inner_loop_sequence_recovery_rate_
				<< ", cross entropy: " << cross_entropy << ", weighted cross entropy: " << weighted_cross_entropy
				<< ", inner loop count: " << inner_loop_counter_ << std::endl;

			if ( outer_loop_counter_ > 2 &&
					/* inner_loop_counter_ > 3 && */
					accept_new_weight_set == 0 &&
					weighted_cross_entropy + inner_loop_sequence_recovery_rate_ > ( outer_loop_seq_profile_cross_entropy_ * Wcross_ent ) + outer_loop_last_sequence_recovery_rate_ ) {
				TR << "decide_if_sequence_recovery_improved(): accepting new weight set: "
					<< "inner loop recovery rate: " << inner_loop_sequence_recovery_rate_
					<< ", outer loop recovery rate: " << outer_loop_last_sequence_recovery_rate_
					<< ", inner loop weighted cross entropy: " << weighted_cross_entropy
					<< ", outer loop weighted cross entropy: " << outer_loop_seq_profile_cross_entropy_ * Wcross_ent << std::endl;
				accept_new_weight_set = 1;
			}

			//if ( inner_loop_counter_ == num_inner_iterations() ) {
			// accept_new_weight_set = 1;
			//}
			// ronj just because we hit the max number of inner iterations, don't accept it it's worse. the inner loop will kill itself and
			// ronj we'll just reswarm and mix with the last good weight set.

			if ( accept_new_weight_set != 0 ) {
				outer_loop_seq_profile_cross_entropy_ = cross_entropy;
			}
		} else {
			if ( outer_loop_counter_ == 1 || inner_loop_sequence_recovery_rate_ > outer_loop_last_sequence_recovery_rate_ ) {
				TR << "decide_if_sequence_recovery_improved(): accepting new weight set: "
					<< "inner loop recovery rate: " << inner_loop_sequence_recovery_rate_
					<< ", outer loop recovery rate: " << outer_loop_last_sequence_recovery_rate_
					<< ", mixing factor: " << mixing_factor_ << std::endl;
				accept_new_weight_set = 1;
			} else {
				accept_new_weight_set = 0;
			}
		}

#ifdef USEMPI
		for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
			MPI_Send( & accept_new_weight_set, 1, MPI_INT, ii, tag_, MPI_COMM_WORLD );
		}
#endif
	}

	/// if we're going to accept this weight set, prepare for the next round of sequence optimization
	if ( accept_new_weight_set != 0 && MPI_rank_ == 0 ) {

		if ( ! option[ optE::dont_use_reference_energies ].user() ) {
			before_minimization_reference_energies_ = reference_energies_inner_loop_;
		}

		TR << "decide_if_sequence_recovery_improved(): free_parameters_ before: ";
		free_parameters_.show_nonzero( TR );
		TR << std::endl;

		TR << "decide_if_sequence_recovery_improved(): free_parameters_ after: ";
		free_weights_inner_loop_.show_nonzero( TR  );
		TR << std::endl;

		outer_loop_last_sequence_recovery_rate_ = inner_loop_sequence_recovery_rate_;
		free_weights_before_minimization_ = free_weights_inner_loop_;
		free_parameters_ = free_weights_inner_loop_;

	} else if ( MPI_rank_ == 0 ) {
		TR << "decide_if_sequence_recovery_improved(): rejected weight set:\n";
		free_weights_inner_loop_.show_nonzero( TR  );
		TR << std::endl;

		// accept_new_weight_set is 0
		//ronj accept the refEs BUT NOT THE OTHER FREE TERMS on the last iteration to avoid redoing the same entropy work
		if ( inner_loop_counter_ == num_inner_iterations() && option[ optE::fit_reference_energies_to_aa_profile_recovery ] ) {
			if ( ! option[ optE::dont_use_reference_energies ].user() ) {
				before_minimization_reference_energies_ = reference_energies_inner_loop_;
			}
		}
	}
	return accept_new_weight_set != 0;

}

///
/// @brief
/// Main loop for the optE protocol.  This is function the apps call to do optE.
///
void
IterativeOptEDriver::go()
{

	barrier();
	//intialize_free_and_fixed_energy_terms();
	TR << "go(): " << node_name( MPI_rank_ ) << std::endl;
	divide_up_pdbs();

	for ( outer_loop_counter_ = 1; outer_loop_counter_ <= num_outer_iterations(); ++outer_loop_counter_ ) {
		TR.Debug << "Node " << MPI_rank_ << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ") collect_rotamer_energies ..." << std::endl;
		collect_rotamer_energies();
		TR.Debug << "Node " << MPI_rank_ << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ") optimize_weights ..." << std::endl;
		optimize_weights();
		for ( inner_loop_counter_ = 1; inner_loop_counter_ <= num_inner_iterations(); ++inner_loop_counter_ ) {
			TR.Debug << "Node " << MPI_rank_ << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ") write_new_scorefile ..." << std::endl;
			if ( MPI_rank_ == 0 ) {
				TR << "go(): " << node_name( MPI_rank_ ) << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ")" << std::endl;
			}
			write_new_scorefile();
			TR.Debug << "Node " << MPI_rank_ << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ") barrier [1] ..." << std::endl;
			barrier();
			if ( ! ( option[ optE::optimize_nat_aa ] || option[ optE::optimize_pssm ] ) ) { break; }
			TR.Debug << "Node " << MPI_rank_ << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ") test_sequence_recovery ..." << std::endl;
			test_sequence_recovery();
			TR.Debug << "Node " << MPI_rank_ << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ") barrier [2] ..." << std::endl;
			barrier();
			TR.Debug << "Node " << MPI_rank_ << " (" << outer_loop_counter_ << "," << inner_loop_counter_ << ") decide_if_sequence_recovery_improved ..." << std::endl;
			if ( decide_if_sequence_recovery_improved() ) break;
		}
	}

	if ( MPI_rank_ == 0 ) { TR << "go(): DONE with weight optimization." << std::endl; }
}

void
IterativeOptEDriver::barrier()
{
#ifdef USEMPI
	MPI_Barrier( MPI_COMM_WORLD );
#endif
	std::cout.flush();
}

///
/// @brief
/// Reads in the list of pdb file names to use for the optimization. Uses basic::options::start_file which returns a vector1 of strings.
/// Unfortunately, this method reads data out of the return value of start_file() which is why you have to use a listfile with the
/// -s option for things to work correctly.
///
utility::vector1< std::string >
IterativeOptEDriver::get_native_pdb_names()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// read list file of pdbs
	utility::vector1< std::string > filenames;


	std::string const listfile( start_file() );
	std::ifstream data( listfile.c_str() );
	std::string line;
	while ( getline( data,line ) ) {
		filenames.push_back( line );
	}
	data.close();

	return filenames;
}

///
/// @brief
/// This function is the heart of the optE protocol.
/// For each position of the protein, we're going to create a new PackerTask that only tries all of the various amino
/// acids at that position.  (It also uses any packer task related flags on the command line.)
/// opte_data is a container class for optE data objects.  So each position is going to have a PNatAAOptEPositionData
/// object.  For each rotamer that's built for each position, the energies for the fixed and free energy terms are
/// stored in a PNatAAOptERotamerData object which gets added to the PNatAAOptEPositionData object.  So at the end
/// of this function we have one optEData object which has position info for all position and rotamer data for all
/// rotamers at each position.
///
void
IterativeOptEDriver::get_nat_aa_opte_data(
	std::string const & pdb_name,
	pose::Pose & pose,
	pose::Pose & native_pose,
	ScoreFunction const & scorefxn,
	ScoreTypes & score_list,
	ScoreTypes & fixed_score_vec,
	OptEData & opte_data
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace pose;
	using namespace pack::task;
	using namespace pack::rotamer_set;

	TR.Trace << "Getting native aa opte data for " << pdb_name << ":";

	PackerTaskOP design_task;
	// If desired lock the task to the starting pose to prevent changes as we go on.
	// Prevent the task from changing which residues are designed over the course of the runs.
	if ( option[ optE::constant_logic_taskops_file ].user() ) {
		design_task = copy_native_packertask_logic( native_pose,
			pose,
			task_factory_);
	} else if ( option[ optE::design_with_minpack ] ) {
		// Don't use extra rotamers during design, but do use extra rotamers when coming up with the
		// rotamers used in the weight fitting step.  Basically, with design_with_minpack on the command
		// line, command line flags that effect the packer no longer effect the design phase.
		design_task = task_factory_->create_task_and_apply_taskoperations( pose );
		pack::task::operation::InitializeFromCommandline ifcl_operation;
		ifcl_operation.apply( pose, *design_task );
	} else {
		design_task = task_factory_->create_task_and_apply_taskoperations( pose ) ;
	}
	design_task->set_bump_check( true );
	//design_task->or_include_current( true ); // WHOA we definately don't want to include current.

	scorefxn.setup_for_packing( pose, design_task->repacking_residues(), design_task->designing_residues() );

	// create some unfolded state potential objects here before we go into the pose residue loop to speed things up
	UnfoldedStatePotential const & unfE_potential = ScoringManager::get_instance()->get_UnfoldedStatePotential( scoring::UNFOLDED_SCORE12 );
	utility::vector1< EnergyMap > e;
	if ( using_unfolded_energy_term_ ) {
		e.resize( chemical::num_canonical_aas );
		for ( Size aa=1; aa <= chemical::num_canonical_aas; ++aa ) {
			e[ aa ].zero();
			unfE_potential.raw_unfolded_state_energymap( chemical::name_from_aa( (chemical::AA) aa ), e[ aa ] );
		}
	}

	// used to restrict design to one position at a time
	utility::vector1< bool > task_mask( pose.size(), false );
	Size num_diffs_between_native_and_input( 0 );

	for ( Size resi = 1; resi <= pose.size(); ++resi ) {

		if ( ! pose.residue(resi).is_protein() ) continue;
		// do not consider residues that are not designable
		if ( ! design_task->residue_task( resi ).being_designed() ) continue;

		TR.Trace << " ";
		if ( pose.pdb_info() ) TR.Trace << pose.pdb_info()->number(resi);
		else TR.Trace << resi;
		TR.Trace << "." << pose.residue_type(resi).name3();

		// use new naive PackerTask for getting single-residue data
		PackerTaskOP single_residue_task = TaskFactory::create_packer_task( pose );
		single_residue_task->initialize_from_command_line();
		task_mask[ resi ] = true;
		// the following turns off packing of all other residues
		single_residue_task->restrict_to_residues( task_mask );
		task_mask[ resi ] = false;

		PNatAAOptEPositionDataOP this_pos_data;

		if ( option[ optE::optimize_pssm ] && resi <= pssm_data_.size() ) {
			if ( pssm_data_[ resi ].first == native_pose.residue( resi ).aa() ) {
				PSSMOptEPositionDataOP data( new PSSMOptEPositionData );
				data->set_pssm_probabilities( pssm_data_[ resi ].second );
				this_pos_data = data;
			} else {
				std::cerr << "Warning position " << resi << " in " << pdb_name << " pssm data does not match native amino acid: ";
				std::cerr << pssm_data_[ resi ].first << " vs " << native_pose.residue( resi ).aa()  << std::endl;
				std::cerr << "Falling back on PNatAAOptEPositionData" << std::endl;
				this_pos_data = PNatAAOptEPositionDataOP( new PNatAAOptEPositionData );
			}
		} else {
			if ( option[ optE::optimize_pssm ]() ) {
				TR << "Warning: " << pdb_name << ".fasta.probs is shorter than PDB file!\n";
				TR << "Falling back on PNatAAOptEPositionData for residue " << resi << std::endl;
			}

			// Create a position data object of the type that has special processing for unfolded state energy calculations.
			// If not for this special check, then we'll always be creating the standard position data objects and won't ever
			// get to the code that deals with the unfolded state energy.
			// Note: This special position data class is not compatible with PSSM optimization.
			if ( using_unfolded_energy_term_ ) {
				this_pos_data = PNatAAOptEPositionDataOP( new NestedEnergyTermPNatAAOptEPositionData );
				(utility::pointer::dynamic_pointer_cast< protocols::optimize_weights::NestedEnergyTermPNatAAOptEPositionData > ( this_pos_data ))->set_unfolded_energy_emap_vector( e );

			} else {
				this_pos_data = PNatAAOptEPositionDataOP( new PNatAAOptEPositionData );
			}
		}

		this_pos_data->tag( pdb_name );
		this_pos_data->set_position( resi );
		this_pos_data->set_native_aa( native_pose.residue( resi ).aa() );
		if ( native_pose.residue( resi ).aa() != pose.residue( resi ).aa() ) {
			//std::cout << "native_residue # " << resi << " of " << native_pose.residue( resi ).aa() << " differs with pose residue " << pose.residue( resi ).aa() << std::endl;
			++num_diffs_between_native_and_input;
		}
		this_pos_data->set_neighbor_count(
			pose.energies().tenA_neighbor_graph().get_node( resi )->num_neighbors_counting_self() );

		graph::GraphCOP packer_neighbor_graph(  pack::create_packer_graph( pose, scorefxn, single_residue_task ) );

		RotamerSetFactory rsf;
		RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue( resi ) );
		//  RotamerSetOP rotset = RotamerSetFactory::create_rotamer_set( pose.residue( resi ) );

		rotset->set_resid( resi );
		rotset->build_rotamers( pose, scorefxn, *single_residue_task, packer_neighbor_graph );
		scorefxn.prepare_rotamers_for_packing( pose, *rotset );

		// First, need a vector of energy maps
		utility::vector1< EnergyMap > emap_vector( rotset->num_rotamers() );

		// Call the new energy map fn
		rotset->compute_one_body_energy_maps( pose, scorefxn, *single_residue_task, packer_neighbor_graph, emap_vector );

		for ( Size jj = 1; jj <= rotset->num_rotamers(); ++jj ) {
			EnergyMap & emap_total( emap_vector[jj] );

			// Hacky limit for fa_rep
			if ( emap_total[ fa_rep ] > 10.0 ) {
				emap_total[ fa_rep ] = 10.0;
			}

			utility::vector1< Real > energy_info;
			utility::vector1< Real > fixed_energy_info;

			// put all the energies for the free energy terms into the 'energy_info' vector
			for ( auto & score_type_iter : score_list ) {
				energy_info.push_back( emap_total[ score_type_iter ] );
			}

			// put all the energies for the fixed energy terms into the 'fixed_energy_info' vector
			for ( auto & score_type_iter : fixed_score_vec ) {
				fixed_energy_info.push_back( emap_total[ score_type_iter ] );
			}

			// the data held inside a PNatAAOptERotamerData object is the two vectors we just set above and the amino acid type
			PNatAAOptERotamerDataOP new_rot_line( new PNatAAOptERotamerData( (*rotset->rotamer( jj )).aa(), jj, energy_info, fixed_energy_info ) );

			// this rotamer line information gets added to the PNatAAOptEPositionData object created above
			this_pos_data->add_rotamer_line_data( new_rot_line );
		}

		// Now that we have rotamer data for all rotamers in the position data object, store the position data object in the
		// optE data container (which holds position data objects for all positions).
		opte_data.add_position_data( this_pos_data );
	}
	TR.Trace << std::endl;

	//TR << "get_nat_aa_opte_data(): num_diffs_between_native_and_input: " << num_diffs_between_native_and_input << std::endl;

}

///
/// @brief
/// Similar to get_nat_aa_opte_data.  See comments there for more info.
///
void
IterativeOptEDriver::get_nat_rot_opte_data(
	std::string const & pdb_name,
	pose::Pose & pose,
	pose::Pose & native_pose,
	utility::vector1<bool> include_rsd,
	ScoreFunction const & scorefxn,
	ScoreTypes & score_list,
	ScoreTypes & fixed_score_vec,
	OptEData & opte_data
)
{
	using namespace pose;
	using namespace pack::task;
	using namespace pack::rotamer_set;
	using namespace pack::dunbrack;
	using namespace pack::rotamers;
	using namespace chemical;

	TR.Trace << "Getting native rotamer opte data for " << pdb_name << ":";
	PackerTaskOP packing_task;
	// If desired lock the task to the starting pose to prevent changes as we go on.
	if ( option[ optE::constant_logic_taskops_file ].user() ) {
		packing_task = copy_native_packertask_logic( native_pose, pose, task_factory_);
	} else if ( option[ optE::design_with_minpack ] ) {
		// Don't use extra rotamers during design, but do use extra rotamers when coming up with the
		// rotamers used in the weight fitting step.  Basically, with design_with_minpack on the command
		// line, command line flags that effect the packer no longer effect the design phase.
		packing_task = task_factory_->create_task_and_apply_taskoperations( pose );
		pack::task::operation::InitializeFromCommandline ifcl_operation;
		ifcl_operation.apply( pose, *packing_task );
	} else {
		packing_task = task_factory_->create_task_and_apply_taskoperations( pose ) ;
	}

	packing_task->set_bump_check( false );
	packing_task->restrict_to_repacking();

	scorefxn( native_pose );
	scorefxn( pose );

	scorefxn.setup_for_packing( pose, packing_task->repacking_residues(), packing_task->designing_residues() );

	utility::vector1< bool > task_mask( pose.size(), false );
	for ( Size resi = 1; resi <= pose.size(); ++resi ) {
		// only consider residues that the master task considers packable (via task_factory_)
		if ( ! packing_task->residue_task( resi ).being_packed() ) continue;

		//if ( residue_has_unacceptably_bad_dunbrack_energy( native_pose, resi )) continue;
		if ( residue_has_bad_bfactor( native_pose, resi ) ) continue;

		//if ( ! pose.residue(resi).is_protein() ) continue;
		if ( ! include_rsd[ resi ] ) continue;

		TR.Trace << " ";
		if ( pose.pdb_info() ) TR.Trace << pose.pdb_info()->number(resi);
		else TR.Trace << resi;
		TR.Trace << "." << pose.residue_type(resi).name3();

		// use new naive PackerTask to get data for one residue at a time
		PackerTaskOP task = TaskFactory::create_packer_task( pose );
		task_mask[ resi ] = true;
		task->restrict_to_residues( task_mask );
		task->restrict_to_repacking();
		task->initialize_from_command_line();
		task->set_bump_check( false );
		//task->or_include_current( true ); // apl TO DO -- do we want to include the native rotamer?
		task_mask[ resi ] = false;

		utility::vector1< Size > rot_wells;

		SingleResidueDunbrackLibrary::n_rotamer_bins_for_aa( pose.residue_type( resi ).aa(), rot_wells );

		if ( rot_wells.size()  == 0 ) continue;

		SingleResidueRotamerLibraryCOP srlib(
			core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue_type( resi ) ) );

		SingleResidueDunbrackLibraryCOP srdlib( utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const > ( srlib ));
		runtime_assert( srdlib != nullptr );

		PNatRotOptEPositionDataOP this_pos_data( new PNatRotOptEPositionData );
		this_pos_data->aa() = pose.residue( resi ).aa();
		std::string tag_to_assign =
			pdb_name + " " +
			chemical::oneletter_code_from_aa( pose.residue( resi ).aa() ) + " " +
			utility::to_string( resi );
		this_pos_data->tag( tag_to_assign );
		this_pos_data->phi() = pose.phi( resi );
		this_pos_data->psi() = pose.psi( resi );
		this_pos_data->set_rotamer_well_counts( rot_wells );
		set_aa_periodicity( this_pos_data, native_pose.residue(resi).aa() );
		this_pos_data->set_native_rotamer_chi( native_pose.residue(resi).chi() );

		runtime_assert( native_pose.residue( resi ).aa() == pose.residue( resi ).aa() );

		utility::vector1< Size > nat_rot_indices;
		srdlib->get_rotamer_from_chi(
			native_pose.residue( resi ).chi(),
			nat_rot_indices );
		nat_rot_indices.resize( rot_wells.size() );
		this_pos_data->set_native_rotamer_index( nat_rot_indices );

		graph::GraphCOP packer_neighbor_graph(  pack::create_packer_graph( pose, scorefxn, task ) );

		RotamerSetFactory rsf;
		RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue( resi ) );
		//  RotamerSetOP rotset = RotamerSetFactory::create_rotamer_set( pose.residue( resi ) );

		rotset->set_resid( resi );
		rotset->build_rotamers( pose, scorefxn, *task, packer_neighbor_graph );
		scorefxn.prepare_rotamers_for_packing( pose, *rotset );

		// First, need a vector of energy maps
		utility::vector1< EnergyMap > emap_vector( rotset->num_rotamers() );

		// Call the new energy map fn
		rotset->compute_one_body_energy_maps( pose, scorefxn, *task, packer_neighbor_graph, emap_vector );

		//std::cout << "Nrotamers for " << pdb_name << " " << resi << " " << rotset->num_rotamers() << " " << pose.residue_type( resi ).name() << std::endl;

		for ( Size jj = 1; jj <= rotset->num_rotamers(); ++jj ) {
			EnergyMap & emap_total( emap_vector[jj] );

			// Hacky limit for fa_rep
			if ( emap_total[ fa_rep ] > 10.0 ) emap_total[ fa_rep ] = 10.0;

			utility::vector1< Real > free_energy_info;
			utility::vector1< Real > fixed_energy_info;

			for ( auto & score_type_iter : score_list ) {
				free_energy_info.push_back( emap_total[ score_type_iter ] );
			}

			for ( auto & score_type_iter : fixed_score_vec ) {
				fixed_energy_info.push_back( emap_total[ score_type_iter ] );
			}

			utility::vector1< Size > rot_index_vector;
			srdlib->get_rotamer_from_chi(
				rotset->rotamer( jj )->chi(),
				rot_index_vector );

			rot_index_vector.resize( rot_wells.size() );

			PNatRotOptERotamerDataOP new_rot_line( new PNatRotOptERotamerData(
				rot_index_vector,
				rotset->rotamer( jj )->chi(),
				free_energy_info,
				fixed_energy_info ) );

			this_pos_data->add_rotamer_line_data( new_rot_line );
		}
		// Done with rotamers for this position, store this position data object
		opte_data.add_position_data( this_pos_data );

		//TR << "Added rotamer position data for residue " << resi << " in pose " << pdb_name << std::endl;

	}
	TR.Trace << std::endl;

}

///
/// @remarks Andrew?
///
void
IterativeOptEDriver::set_aa_periodicity( PNatRotOptEPositionDataOP pos_data, core::chemical::AA aa ) const
{
	using namespace core::chemical;
	switch ( aa ) {
	case aa_ala: case aa_gly : break;
	case aa_cys: case aa_ser: case aa_thr: case aa_val : {
		utility::vector1< Real > asym1(1, 360);
		pos_data->set_native_chi_periodicity( asym1 );
	}
		break;
	case aa_asp: case aa_phe: case aa_tyr : {
		utility::vector1< Real > sym2( 2 ); sym2[ 1 ] = 360; sym2[ 2 ] = 180;
		pos_data->set_native_chi_periodicity( sym2 );
	}
		break;
	case aa_his: case aa_ile: case aa_leu: case aa_asn: case aa_trp : {
		utility::vector1< Real > asym2(2, 360);
		pos_data->set_native_chi_periodicity( asym2 );
	}
		break;
	case aa_glu : {
		utility::vector1< Real > sym3( 3, 360 ); sym3[ 3 ] = 180;
		pos_data->set_native_chi_periodicity( sym3 );
	}
		break;

	case aa_met: case aa_gln: case aa_pro : {
		utility::vector1< Real > asym3( 3, 360 );
		pos_data->set_native_chi_periodicity( asym3 );
	}
		break;
	case aa_arg: case aa_lys : {
		utility::vector1< Real > asym4( 4, 360 );
		pos_data->set_native_chi_periodicity( asym4 );
	}
	default :
		break;

	}
}

///
/// @details Precondition: pose must have been scored
///
bool
IterativeOptEDriver::residue_has_unacceptably_bad_dunbrack_energy( core::pose::Pose const & pose, Size const resid ) const
{
	using namespace core::chemical;
	using namespace core::scoring;
	switch ( pose.residue_type( resid ).aa() ) {
	case aa_ala: case aa_gly : return false; break;
	case aa_cys: case aa_ser: case aa_thr: case aa_val: case aa_pro :
		if ( pose.energies().residue_total_energies( resid )[ fa_dun ] > 10 ) return true;
		break;
	case aa_asp: case aa_phe: case aa_his: case aa_ile: case aa_leu: case aa_asn: case aa_trp: case aa_tyr :
		if ( pose.energies().residue_total_energies( resid )[ fa_dun ] > 15 ) return true;
		break;
	case aa_glu: case aa_met: case aa_gln :
		if ( pose.energies().residue_total_energies( resid )[ fa_dun ] > 18 ) return true;
		break;
	case aa_arg: case aa_lys :
		if ( pose.energies().residue_total_energies( resid )[ fa_dun ] > 22 ) return true;
		break;
	default :
		break;
	}
	return false;
}

///
/// @details pose must have been read from a pdb.
///
bool
IterativeOptEDriver::residue_has_bad_bfactor( core::pose::Pose const & pose, Size const resid ) const
{
	using namespace core::pose;
	PDBInfoCOP info = pose.pdb_info();
	if ( !info ) return false;

	for ( Size ii = 1; ii <= info->natoms( resid ); ++ii ) {
		//std::cout << "Temperature on " << resid << " " << ii << " " << info->temperature( resid, ii ) << std::endl;
		if ( info->temperature( resid, ii ) > 40 ) {
			return true;
		}

	}
	return false;
}

///
/// @brief
/// Helper function to reduce code duplication.
////
SingleStructureDataOP
IterativeOptEDriver::make_simple_ssd_from_pdb( std::string const & pdb_filename,
	core::scoring::ScoreFunctionOP sfxn, bool pretend_no_fa_rep ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::pose::Pose structure;
	if ( option[ in::file::centroid_input ] ) {
		core::import_pose::centroid_pose_from_pdb( structure, pdb_filename , core::import_pose::PDB_file);
	} else {
		core::import_pose::pose_from_file( structure, pdb_filename , core::import_pose::PDB_file);
	}

	/// score this pose, create SingleStructureData.
	(*sfxn)( structure );
	utility::vector1< Real > free_data( free_score_list_.size(), 0.0 );
	utility::vector1< Real > fixed_data( fixed_score_list_.size(), 0.0 );

	for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
		if ( !pretend_no_fa_rep || free_score_list_[ kk ] != fa_rep ) {
			free_data[ kk ] = structure.energies().total_energies()[ free_score_list_[ kk ] ];
		}
	}
	for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
		if ( !pretend_no_fa_rep || fixed_score_list_[ kk ] != fa_rep ) {
			fixed_data[ kk ] = structure.energies().total_energies()[ fixed_score_list_[ kk ] ];
		}
	}

	SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
	return ssd;
}

///
/// @brief
/// dG optimization optE data collection.
////
void
IterativeOptEDriver::collect_dG_of_binding_data()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( dG_binding_data_ == nullptr ) {
		dG_binding_data_ = OptEDataOP( new OptEData() );
		ScoreFunctionOP sfxn = create_unweighted_scorefunction();

		bool const no_fa_rep = option[ optE::pretend_no_ddG_repulsion ]();
		for ( Size ii = 1; ii <= dG_bound_unbound_pairs_.size(); ++ii ) {
			DGBindOptEDataOP dg_data( new DGBindOptEData() );
			dg_data->deltaG_bind( dG_binding_[ ii ] );
			dg_data->bound_struct(   make_simple_ssd_from_pdb( dG_bound_unbound_pairs_[ ii ].first,  sfxn, no_fa_rep ) );
			dg_data->unbound_struct( make_simple_ssd_from_pdb( dG_bound_unbound_pairs_[ ii ].second, sfxn, no_fa_rep ) );
			dG_binding_data_->add_position_data( dg_data );
		}
	}

	for ( auto
			iter = dG_binding_data_->position_data_begin(),
			iter_end = dG_binding_data_->position_data_end();
			iter != iter_end; ++iter ) {
		optE_data_->add_position_data( *iter );
	}
}

///
/// @brief
/// ddG optimization optE data collection.
////
void
IterativeOptEDriver::collect_ddG_of_mutation_data()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( ddG_mutation_data_ == nullptr ) {
		ddG_mutation_data_ = OptEDataOP( new OptEData );

		// rescore wt's less often
		std::map< std::string, std::pair< SingleStructureDataOP, std::string > > structure_map;
		std::map< std::string, EnergyMap > structure_energy_map;
		std::map< std::string, Energy > structure_bestfarep_map;


		ScoreFunctionOP sfxn = create_unweighted_scorefunction();

		if ( MPI_rank_ == 0 ) {
			TR_VERBOSE << "collect_ddG_of_mutation_data(): created sfxn for calculating ddGs" << std::endl;
			sfxn->show( std::cout );
		}

		// for each line of the input ddG data file...
		// the first column is a text file that contains the names of wt pdbs
		// the second column is another text file that contains the names of mutant pdbs
		// the third column gets added as the experimental ddG for this mutation

		// this input system allows one to have multiple structures of the same exact sequence. this way, if a mutation is causing a slight clash
		// or maybe creating a hydrogen bond in one structure but not in another, the better structure will be used for weight optimization.
		// this can make a big difference in energies depending on what the weight on the fa_rep or hbond term, etc, is.
		for ( Size ii = 1; ii <= ddg_mut_wt_pairs_.size(); ++ii ) {

			DDGMutationOptEDataOP ddg_data;
			UnfoldedStatePotential const & unfE_potential = ScoringManager::get_instance()->get_UnfoldedStatePotential( scoring::UNFOLDED_SCORE12 );

			// Create a position data object of the type that has special processing for unfolded state energy calculations.
			// If not for this special check, then we'll always be creating the standard position data objects and won't ever
			// get to the code that deals with the unfolded state energy.
			if ( using_unfolded_energy_term_ ) {
				ddg_data = DDGMutationOptEDataOP( new NestedEnergyTermDDGMutationOptEData );
			} else {
				ddg_data = DDGMutationOptEDataOP( new DDGMutationOptEData );
			}

			// save the experimental ddg for this wt/mut list-of-files pair
			ddg_data->set_experimental_ddg( ddGs_[ ii ] );

			utility::file::FileName wts(ddg_mut_wt_pairs_[ ii ].first);
			std::string file_extension = wts.ext();

			utility::vector1< std::string > wt_pdb_names, mut_pdb_names;

			bool read_silent( false );
			core::io::silent::SilentFileData sfd_wt;
			core::io::silent::SilentFileData sfd_mut;

			if ( file_extension.compare("list") == 0 ) {

				/// Read names of wt pdbs; i.e. open up the filename given in the 1st column of the ddG data file
				/// and read out the strings of pdbs listed there?
				TR << "collect_ddG_of_mutation_data(): reading file '" << wts() << "' to get list of wt pdb names." << std::endl;
				std::ifstream wt_pdblist( wts().c_str() );
				while ( wt_pdblist ) {
					std::string wt_pdb;
					wt_pdblist >> wt_pdb;
					if ( wt_pdb != "" ) wt_pdb_names.push_back( wt_pdb );
				}
			} else if ( file_extension.compare("out") == 0 ) { //add in silent file capabilities
				read_silent=true;
				sfd_wt.set_filename(wts()); //for now assume binary
				if ( !sfd_wt.read_file(wts()) ) {
					std::cout << "[ERROR ERROR ERROR] did not read in silent file properly! " << wts() << std::endl;
				}
				wt_pdb_names = sfd_wt.tags();
			} else {
				//file extension not recognized
				std::cerr << "ERROR! file " << wts() << " has un-recognized extension " << file_extension  << std::endl;
				utility_exit();
			}

			/// Read names of mut pdbs
			bool no_tag_yet_assigned( true );
			utility::file::FileName muts(ddg_mut_wt_pairs_[ii].second);
			file_extension = muts.ext();

			if ( file_extension.compare("list") == 0 ) {
				TR << "collect_ddG_of_mutation_data(): reading file '" << muts() << "' to get list of mutant pdb names." << std::endl;
				std::ifstream mut_pdblist( muts().c_str() );
				while ( mut_pdblist ) {
					std::string mut_pdb;
					mut_pdblist >> mut_pdb;
					if ( mut_pdb != "" ) mut_pdb_names.push_back( mut_pdb );
					if ( no_tag_yet_assigned ) {
						utility::file::FileName mut1( mut_pdb );
						ddg_data->tag( mut1.base() );
						no_tag_yet_assigned = false;
					}
				}
			} else if ( file_extension.compare("out") == 0 ) {
				read_silent=true;
				sfd_mut.set_filename(muts());
				if ( !sfd_mut.read_file(muts()) ) {
					std::cout << "[ERROR ERROR ERROR] did not read in silent file properly! " << muts() << std::endl;
				}
				mut_pdb_names = sfd_mut.tags();
				if ( no_tag_yet_assigned ) {
					ddg_data->tag(mut_pdb_names[1]);
					no_tag_yet_assigned = false;
				}
			} else {
				std::cerr << "ERROR! file " << muts() << " has un-recognized extension " << file_extension << std::endl;
				utility_exit();
				//file extension not recognized
			}
			std::string wt_seq, mut_seq; // wt and mutant sequences; must differ at exactly one position
			Real best_wt_rep( 12345678 ), best_mut_rep( 12345678 );
			bool collect_best_rep( option[ optE::exclude_badrep_ddGs ].user() );

			for ( Size jj = 1; jj <= wt_pdb_names.size(); ++jj ) {


				// since a given wt protein might have a few hundred characterized mutants, there's no point in scoring the
				// wild type structure for each of those hundred mutants. we can score the wt once and save that score for
				// all of the mutants of that structure.  ingenious time-saver thanks to APL.  -ronj
				if ( structure_map.find( wts()+wt_pdb_names[ jj ] ) == structure_map.end() ) { //wt_pdb_name not already in structure_map

					core::pose::Pose wt_structure;
					if ( !read_silent ) {
						if ( option[ in::file::centroid_input ] ) {
							core::import_pose::centroid_pose_from_pdb( wt_structure, wt_pdb_names[ jj ] , core::import_pose::PDB_file);
						} else {
							core::import_pose::pose_from_file( wt_structure, wt_pdb_names[ jj ]  , core::import_pose::PDB_file);
						}
					} else {
						core::io::silent::SilentStructOP ss = sfd_wt[wt_pdb_names[ jj ]];
						core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
						ss->fill_pose(wt_structure,*rsd_set);
					}

					/// make sure sequences match across poses...
					if ( jj == 1 ) {
						wt_seq = wt_structure.sequence();
						//wt_structure.dump_pdb("wt_structure.pdb");
					} else {
						if ( wt_seq != wt_structure.sequence() ) {
							std::cerr << "Wild type sequence inconsistent across wts " << std::endl;
							std::cerr << wt_seq << std::endl;
							std::cerr << wt_structure.sequence() << std::endl;
							std::cerr << "Node " << MPI_rank_ << " " << wt_pdb_names[ jj ] << " and " << wt_pdb_names[ 1 ] << std::endl;
							wt_structure.dump_pdb("offending_wt_sequence.pdb");
							utility_exit();
						}
					}

					/// score this pose, create SingleStructureData.
					utility::vector1< Real > free_data( free_score_list_.size(), 0.0 );
					utility::vector1< Real > fixed_data( fixed_score_list_.size(), 0.0 );

					(*sfxn)( wt_structure );

					// Adding a special check for whether or not surface scoring is in use. If yes, then we need to zero out the value for
					// the surface energy.  It appears to be hurting the correlations for ddG rather than helping. The EnergyMethod returns
					// 0.0 when the residue_energy() method is called, but at the end of scoring in finalize_total_energy(), the total
					// surface score for the pose is calculated and placed into the pose.energies() object.
					wt_structure.energies().total_energies()[ scoring::surface ] = 0.0;

					if ( collect_best_rep ) {
						if ( wt_structure.energies().total_energies()[ fa_rep ] < best_wt_rep ) {
							best_wt_rep = wt_structure.energies().total_energies()[ fa_rep ];
						}
						structure_bestfarep_map[ wts()+wt_pdb_names[ jj ] ] = best_wt_rep;
					}

					for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
						if ( !option[ optE::pretend_no_ddG_repulsion ] || free_score_list_[ kk ] != fa_rep ) {
							free_data[ kk ] = wt_structure.energies().total_energies()[ free_score_list_[ kk ] ];
						}
					}
					for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
						if ( !option[ optE::pretend_no_ddG_repulsion ] || fixed_score_list_[ kk ] != fa_rep ) {
							fixed_data[ kk ] = wt_structure.energies().total_energies()[ fixed_score_list_[ kk ] ];
						}
					}
					SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
					ddg_data->add_wt( ssd );
					structure_map[ wts()+wt_pdb_names[ jj ] ] = std::make_pair( ssd, wt_seq );

					// The 'unfolded' energy term, if in use, also needs special handling. However, at this point, we don't have a
					// set of weights. So we can only store the unweighted energy term energies that go into the unfolded energy
					// here and not actually calculate a final unfolded energy that can be placed in the pose emap. This case will be
					// handled in the same way as PNatAA handles it, with a class that extends DDGOptEData.  But, we need to set that
					// extra variable here.  Regardless of whether the wild-type structures have already been encountered or not,
					// set the emap here, on this first iteration through the list of wild-type structures.
					if ( using_unfolded_energy_term_ ) {
						if ( jj == 1 ) {
							EnergyMap e;
							unfE_potential.pose_raw_unfolded_state_energymap( wt_structure, e );
							(utility::pointer::dynamic_pointer_cast< protocols::optimize_weights::NestedEnergyTermDDGMutationOptEData > ( ddg_data ))->set_wt_unfolded_energies_emap( e );
							structure_energy_map[ wts()+wt_pdb_names[ jj ] ] = e;
						}
					}

				} else {
					// else, this wild-type structure (or list of wild-type structures) has already been encountered previously. So don't waste
					// time scoring it again.  Just store the results from the previous one.  But make sure to also store the best rep and
					// unfolded energy.
					if ( jj == 1 ) {
						wt_seq = structure_map[ wts()+wt_pdb_names[ jj ] ].second;
					}

					if ( collect_best_rep ) {
						best_wt_rep = structure_bestfarep_map[ wts()+wt_pdb_names[ jj ] ];
					}

					if ( using_unfolded_energy_term_ ) {
						if ( jj == 1 ) {
							(utility::pointer::dynamic_pointer_cast< protocols::optimize_weights::NestedEnergyTermDDGMutationOptEData > ( ddg_data ))->set_wt_unfolded_energies_emap( structure_energy_map[ wts()+wt_pdb_names[ jj ] ] );
						}
					}

					ddg_data->add_wt( structure_map[ wts()+wt_pdb_names[ jj ] ].first );

				}

			}

			for ( Size jj = 1; jj <= mut_pdb_names.size(); ++jj ) {
				core::pose::Pose mut_structure;
				if ( !read_silent ) {
					if ( option[ in::file::centroid_input ] ) {
						core::import_pose::centroid_pose_from_pdb( mut_structure, mut_pdb_names[ jj ] , core::import_pose::PDB_file);
					} else {
						core::import_pose::pose_from_file( mut_structure, mut_pdb_names[ jj ]  , core::import_pose::PDB_file);
					}
				} else {
					core::io::silent::SilentStructOP ss = sfd_mut[mut_pdb_names[ jj ]];
					ss->fill_pose(mut_structure, *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)) );
				}

				/// make sure sequences match across poses...
				if ( jj == 1 ) {
					mut_seq = mut_structure.sequence();
				} else {
					if ( mut_seq != mut_structure.sequence() ) {
						std::cerr << "Mutant sequence inconsistent across muts " << std::endl;
						std::cerr << mut_seq << std::endl;
						std::cerr << mut_structure.sequence() << std::endl;
						std::cerr << "Node " << MPI_rank_ << " " << mut_pdb_names[ jj ] << " and " << mut_pdb_names[ 1 ] << std::endl;
						mut_structure.dump_pdb("offending_mutant_structure_seq.pdb");
						utility_exit();
					}
				}

				/// score this pose, create SingleStructureData.
				utility::vector1< Real > free_data( free_score_list_.size(), 0.0 );
				utility::vector1< Real > fixed_data( fixed_score_list_.size(), 0.0 );

				(*sfxn)( mut_structure );

				// Taking out the calculation of the surface score for ddG optimization. It appears to be hurting the correlations rather
				// than helping.
				mut_structure.energies().total_energies()[ scoring::surface ] = 0.0;

				if ( collect_best_rep ) {
					if ( mut_structure.energies().total_energies()[ fa_rep ] < best_mut_rep ) {
						best_mut_rep = mut_structure.energies().total_energies()[ fa_rep ];
					}
				}


				for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
					if ( !option[ optE::pretend_no_ddG_repulsion ] || free_score_list_[ kk ] != fa_rep ) {
						free_data[ kk ] = mut_structure.energies().total_energies()[ free_score_list_[ kk ] ];
					}
				}
				for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
					if ( !option[ optE::pretend_no_ddG_repulsion ] || fixed_score_list_[ kk ] != fa_rep ) {
						fixed_data[ kk ] = mut_structure.energies().total_energies()[ fixed_score_list_[ kk ] ];
					}
				}
				SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
				ddg_data->add_mutant( ssd );

				// See note above for explanation of what's needed to correctly handle the unfolded state energy
				if ( using_unfolded_energy_term_ ) {
					if ( jj == 1 ) {
						EnergyMap e;
						unfE_potential.pose_raw_unfolded_state_energymap( mut_structure, e );
						(utility::pointer::dynamic_pointer_cast< protocols::optimize_weights::NestedEnergyTermDDGMutationOptEData > ( ddg_data ))->set_mut_unfolded_energies_emap( e );
					}
				}

			}

			/// Find the discrepant position;
			//TR << "collect_ddG_of_mutation_data(): looking for discrepant position between wt and mut structures." << std::endl;
			int discrepant_position( -1 );
			if ( wt_seq.size() != mut_seq.size() ) {
				std::cerr << "ERROR: Wild type and mutant sequences of different length: " << std::endl;
				std::cerr << wt_seq << std::endl;
				std::cerr << mut_seq << std::endl;
				std::cerr << "Rank " << MPI_rank_ << " " << wt_pdb_names[ 1 ] << " vs " << mut_pdb_names[ 1 ] << std::endl;
				utility_exit();
			} else {
				for ( Size jj = 0; jj < wt_seq.size(); ++jj ) {
					if ( discrepant_position == -1 ) {
						if ( wt_seq[ jj ] != mut_seq[ jj ] ) {
							discrepant_position = jj;
						}
					} else if ( wt_seq[ jj ] != mut_seq[ jj ] ) {
						std::cerr << "Error: Wild type and mutant sequences differ at more than one position: "
							<< discrepant_position << " " << jj << std::endl;
						std::cerr << wt_seq << std::endl;
						std::cerr << mut_seq << std::endl;
						std::cerr << "Rank " << MPI_rank_ << " " << wt_pdb_names[ 1 ] << " vs " << mut_pdb_names[ 1 ] << std::endl;
						utility_exit();
					}
				}
			}
			if ( discrepant_position == -1 ) {
				std::cerr << "ERROR: mutant and wild type sequences must differ by at least one position";
				std::cerr << "Rank " << MPI_rank_ << " " << wt_pdb_names[ 1 ] << " vs " << mut_pdb_names[ 1 ] << std::endl;
				utility_exit();
			}
			ddg_data->set_wt_aa( chemical::aa_from_oneletter_code( wt_seq[ discrepant_position ] ));
			ddg_data->set_mut_aa( chemical::aa_from_oneletter_code( mut_seq[ discrepant_position ] ));

			if ( ! collect_best_rep  || best_mut_rep - best_wt_rep < option[ optE::exclude_badrep_ddGs ]() ) {
				ddG_mutation_data_->add_position_data( ddg_data );
			} else { /// else, discard this position...
				TR << "Rank " << MPI_rank_ << " Excluding ddG data from " << ddg_mut_wt_pairs_[ ii ].second
					<< " with mt-wt rep delta: "<< best_mut_rep << " - " << best_wt_rep
					<< " = " << best_mut_rep - best_wt_rep << std::endl;
			}
		}
	} // if ddg_mutation_data_ == 0

	for ( auto
			iter = ddG_mutation_data_->position_data_begin(),
			iter_end = ddG_mutation_data_->position_data_end();
			iter != iter_end; ++iter ) {
		optE_data_->add_position_data( *iter );
	}
}

///
/// @brief
/// The calculations for ddG of binding for interfaces doesn't really fit with the ddG stability and dG binding
/// optE modes. This functions loads all of the necessary structures for optimization of ddG of binding into the
/// optE framework.
////
void
IterativeOptEDriver::collect_ddG_of_binding_data()
{
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( ddG_bind_optE_data_ == nullptr ) {
		ddG_bind_optE_data_ = OptEDataOP( new OptEData );

		// rescore wt's less often
		std::map< std::string, std::pair< SingleStructureDataOP, std::string > > structure_map;
		std::map< std::string, Energy > structure_bestfarep_map;

		ScoreFunctionOP sfxn = create_unweighted_scorefunction();

		if ( MPI_rank_ == 0 ) {
			TR_VERBOSE << "collect_ddG_of_binding_data(): created sfxn for calculating ddGs of binding" << std::endl;
			sfxn->show( std::cout );
		}

		// for each line of the input ddG data file...
		// the first and second columns are text files containing the names of wt and mutant complex pdbs
		// the third and fourth columns are text files containing the names of wt and mutant unbounded pdbs
		// the fifth column gets added as the experimental ddG of binding for this mutation

		// this input system allows one to have multiple structures of the same exact sequence. this way, if a mutation is causing a slight clash
		// or maybe creating a hydrogen bond in one structure but not in another, the better structure will be used for weight optimization.
		// this can make a big difference in energies depending on what the weight on the fa_rep or hbond term, etc, is.
		for ( Size ii = 1; ii <= ddG_bind_files_.size(); ++ii ) {

			DDGBindOptEDataOP ddg_bind_position_data( new DDGBindOptEData );

			// save the experimental ddg for this wt/mut list-of-files pair
			ddg_bind_position_data->set_experimental_ddg_bind( ddGs_binding_[ ii ] );


			utility::file::FileName wt_complexes_list_file( ddG_bind_files_[ ii ][ DDGBindOptEData::WT_COMPLEXES_LIST_FILE ] );
			utility::file::FileName mut_complexes_list_file( ddG_bind_files_[ ii ][ DDGBindOptEData::MUT_COMPLEXES_LIST_FILE ] );
			utility::file::FileName wt_unbounds_list_file( ddG_bind_files_[ ii ][ DDGBindOptEData::WT_UNBOUNDS_LIST_FILE ] );
			utility::file::FileName mut_unbounds_list_file( ddG_bind_files_[ ii ][ DDGBindOptEData::MUT_UNBOUNDS_LIST_FILE ] );

			utility::vector1< std::string > wt_complex_pdb_names, mut_complex_pdb_names, wt_unbounded_pdb_names, mut_unbounded_pdb_names;

			/// Read names of wt complexes pdbs
			TR << "collect_ddG_of_binding_data(): reading file '" << wt_complexes_list_file() << "' to get list of wt complex pdb names." << std::endl;
			std::ifstream wt_complex_pdblist( wt_complexes_list_file().c_str() );
			while ( wt_complex_pdblist ) {
				std::string wt_complex_pdb;
				wt_complex_pdblist >> wt_complex_pdb;
				if ( wt_complex_pdb != "" ) wt_complex_pdb_names.push_back( wt_complex_pdb );
			}

			/// Read names of mut complexes pdbs
			bool no_tag_yet_assigned( true );

			TR << "collect_ddG_of_binding_data(): reading file '" << mut_complexes_list_file() << "' to get list of mutant complex pdb names." << std::endl;
			std::ifstream mut_complex_pdblist( mut_complexes_list_file().c_str() );
			while ( mut_complex_pdblist ) {
				std::string mut_complex_pdb;
				mut_complex_pdblist >> mut_complex_pdb;
				if ( mut_complex_pdb != "" ) mut_complex_pdb_names.push_back( mut_complex_pdb );
				if ( no_tag_yet_assigned ) {
					utility::file::FileName mut1( mut_complex_pdb );
					ddg_bind_position_data->tag( mut1.base() );
					TR << "collect_ddG_of_binding_data(): assigned tag: '" << mut1.base() << "' to this set of ddG bind files." << std::endl;
					no_tag_yet_assigned = false;
				}
			}

			/// Read names of wt unbounded pdbs
			TR << "collect_ddG_of_binding_data(): reading file '" << wt_unbounds_list_file() << "' to get list of wt unbounded pdb names." << std::endl;
			std::ifstream wt_unbounded_pdblist( wt_unbounds_list_file().c_str() );
			while ( wt_unbounded_pdblist ) {
				std::string wt_unbounded_pdb;
				wt_unbounded_pdblist >> wt_unbounded_pdb;
				if ( wt_unbounded_pdb != "" ) wt_unbounded_pdb_names.push_back( wt_unbounded_pdb );
			}

			/// Read names of mut unbounded pdbs
			TR << "collect_ddG_of_binding_data(): reading file '" << mut_unbounds_list_file() << "' to get list of mutant unbounded pdb names." << std::endl;
			std::ifstream mut_unbounded_pdblist( mut_unbounds_list_file().c_str() );
			while ( mut_unbounded_pdblist ) {
				std::string mut_unbounded_pdb;
				mut_unbounded_pdblist >> mut_unbounded_pdb;
				if ( mut_unbounded_pdb != "" ) mut_unbounded_pdb_names.push_back( mut_unbounded_pdb );
			}


			// make sure the wt complexes and wt unbounds have the same sequence; same for the mutant
			// however, don't require that the wt and mutant complex (or unbounds) differ in only one position as for ddG stability
			std::string wt_complex_seq, mut_complex_seq;
			Real best_wt_complex_rep( 12345678 ), best_mut_complex_rep( 12345678 );
			bool filtering_bad_ddGs( option[ optE::exclude_badrep_ddGs ].user() );

			for ( Size jj = 1; jj <= wt_complex_pdb_names.size(); ++jj ) {

				// since a given wt protein might have a few hundred characterized mutants, there's no point in scoring the
				// wild type structure for each of those hundred mutants. we can score the wt once and save that score for
				// all of the mutants of that structure.  ingenious time-saver thanks to APL.  -ronj
				if ( structure_map.find( wt_complexes_list_file()+wt_complex_pdb_names[ jj ] ) == structure_map.end() ) { //wt_pdb_name not already in structure_map

					pose::Pose wt_complex;
					core::import_pose::pose_from_file( wt_complex, wt_complex_pdb_names[ jj ]  , core::import_pose::PDB_file);

					/// make sure sequences match across the list of poses...
					if ( jj == 1 ) {
						wt_complex_seq = wt_complex.sequence();
					} else {
						if ( wt_complex_seq != wt_complex.sequence() ) {
							std::cerr << "wild type complex sequence inconsistent across wts " << std::endl;
							std::cerr << wt_complex_seq << std::endl;
							std::cerr << wt_complex.sequence() << std::endl;
							std::cerr << "Node " << MPI_rank_ << " " << wt_complex_pdb_names[ jj ] << " and " << wt_complex_pdb_names[ 1 ] << std::endl;
							wt_complex.dump_pdb("offending_wt_complex_sequence.pdb");
							utility_exit();
						}
					}

					/// score this pose, create SingleStructureData.
					utility::vector1< Real > free_data( free_score_list_.size(), 0.0 );
					utility::vector1< Real > fixed_data( fixed_score_list_.size(), 0.0 );

					(*sfxn)( wt_complex );

					for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
						free_data[ kk ] = wt_complex.energies().total_energies()[ free_score_list_[ kk ] ];
					}
					for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
						fixed_data[ kk ] = wt_complex.energies().total_energies()[ fixed_score_list_[ kk ] ];
					}

					if ( filtering_bad_ddGs ) {
						if ( wt_complex.energies().total_energies()[ fa_rep ] < best_wt_complex_rep ) {
							best_wt_complex_rep = wt_complex.energies().total_energies()[ fa_rep ];
						}
						structure_bestfarep_map[ wt_complexes_list_file()+wt_complex_pdb_names[ jj ] ] = best_wt_complex_rep;
					}

					SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
					ddg_bind_position_data->add_wt_complex( ssd );
					structure_map[ wt_complexes_list_file()+wt_complex_pdb_names[ jj ] ] = std::make_pair( ssd, wt_complex_seq );

				} else {

					// else, this wild-type structure (or list of wild-type structures) has already been encountered previously. So don't waste
					// time scoring it again.  Just store the results from the previous one.  But make sure to also store the best rep and
					// unfolded energy.
					if ( jj == 1 ) {
						wt_complex_seq = structure_map[ wt_complexes_list_file()+wt_complex_pdb_names[ jj ] ].second;
					}

					if ( filtering_bad_ddGs ) {
						best_wt_complex_rep = structure_bestfarep_map[ wt_complexes_list_file()+wt_complex_pdb_names[ jj ] ];
					}

					ddg_bind_position_data->add_wt_complex( structure_map[ wt_complexes_list_file()+wt_complex_pdb_names[ jj ] ].first );
				}
			}

			for ( Size jj = 1; jj <= mut_complex_pdb_names.size(); ++jj ) {

				pose::Pose mut_complex;
				core::import_pose::pose_from_file( mut_complex, mut_complex_pdb_names[ jj ]  , core::import_pose::PDB_file);

				/// make sure sequences match across poses...
				if ( jj == 1 ) {
					mut_complex_seq = mut_complex.sequence();
				} else {
					if ( mut_complex_seq != mut_complex.sequence() ) {
						std::cerr << "mutant complex sequence inconsistent across muts " << std::endl;
						std::cerr << mut_complex_seq << std::endl;
						std::cerr << mut_complex.sequence() << std::endl;
						std::cerr << "Node " << MPI_rank_ << " " << mut_complex_pdb_names[ jj ] << " and " << mut_complex_pdb_names[ 1 ] << std::endl;
						mut_complex.dump_pdb("offending_mutant_complex_structure_seq.pdb");
						utility_exit();
					}
				}

				/// score this pose, create SingleStructureData.
				utility::vector1< Real > free_data( free_score_list_.size(), 0.0 );
				utility::vector1< Real > fixed_data( fixed_score_list_.size(), 0.0 );

				(*sfxn)( mut_complex );

				for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
					free_data[ kk ] = mut_complex.energies().total_energies()[ free_score_list_[ kk ] ];
				}
				for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
					fixed_data[ kk ] = mut_complex.energies().total_energies()[ fixed_score_list_[ kk ] ];
				}

				if ( filtering_bad_ddGs ) {
					if ( mut_complex.energies().total_energies()[ fa_rep ] < best_mut_complex_rep ) {
						best_mut_complex_rep = mut_complex.energies().total_energies()[ fa_rep ];
					}
				}

				SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
				ddg_bind_position_data->add_mutant_complex( ssd );
			}

			for ( Size jj = 1; jj <= wt_unbounded_pdb_names.size(); ++jj ) {

				if ( structure_map.find( wt_unbounds_list_file()+wt_unbounded_pdb_names[ jj ] ) == structure_map.end() ) {

					pose::Pose wt_unbounded;
					core::import_pose::pose_from_file( wt_unbounded, wt_unbounded_pdb_names[ jj ]  , core::import_pose::PDB_file);

					/// make sure these unbounded structure sequences match the complex sequences...
					/// this will also check to make sure all of the unbounded structure match in sequence
					//if ( wt_unbounded.sequence() != wt_complex_seq ) {
					// std::cerr << "wild type unbounded sequence inconsistent with complex structure " << std::endl;
					// std::cerr << wt_unbounded.sequence() << std::endl;
					// std::cerr << wt_complex_seq << std::endl;
					// std::cerr << "Node " << MPI_rank_ << " " << wt_unbounded_pdb_names[ jj ] << " and " << wt_complex_pdb_names[ 1 ] << std::endl;
					// wt_unbounded.dump_pdb("offending_wt_unbounded_sequence.pdb");
					// utility_exit();
					//}

					/// score this pose, create SingleStructureData.
					utility::vector1< Real > free_data( free_score_list_.size(), 0.0 );
					utility::vector1< Real > fixed_data( fixed_score_list_.size(), 0.0 );

					(*sfxn)( wt_unbounded );

					for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
						free_data[ kk ] = wt_unbounded.energies().total_energies()[ free_score_list_[ kk ] ];
					}
					for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
						fixed_data[ kk ] = wt_unbounded.energies().total_energies()[ fixed_score_list_[ kk ] ];
					}

					SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
					ddg_bind_position_data->add_wt_unbounds( ssd );
					structure_map[ wt_unbounds_list_file()+wt_unbounded_pdb_names[ jj ] ] = std::make_pair( ssd, wt_complex_seq ); // use the wt_complex_seq here

				} else {
					ddg_bind_position_data->add_wt_unbounds( structure_map[ wt_unbounds_list_file()+wt_unbounded_pdb_names[ jj ] ].first );
				}
			}

			for ( Size jj = 1; jj <= mut_unbounded_pdb_names.size(); ++jj ) {

				pose::Pose mut_unbounded;
				core::import_pose::pose_from_file( mut_unbounded, mut_unbounded_pdb_names[ jj ]  , core::import_pose::PDB_file);

				/// make sure these unbounded structure sequences match the complex sequences...
				/// this will also check to make sure all of the unbounded structure match in sequence
				//if ( mut_unbounded.sequence() != mut_complex_seq ) {
				// std::cerr << "mutant unbounded sequence inconsistent with complex structure " << std::endl;
				// std::cerr << mut_unbounded.sequence() << std::endl;
				// std::cerr << mut_complex_seq << std::endl;
				// std::cerr << "Node " << MPI_rank_ << " " << mut_unbounded_pdb_names[ jj ] << " and " << mut_complex_pdb_names[ 1 ] << std::endl;
				// mut_unbounded.dump_pdb("offending_mutant_unbounded_structure_seq.pdb");
				// utility_exit();
				//}

				/// score this pose, create SingleStructureData.
				utility::vector1< Real > free_data( free_score_list_.size(), 0.0 );
				utility::vector1< Real > fixed_data( fixed_score_list_.size(), 0.0 );

				(*sfxn)( mut_unbounded );

				for ( Size kk = 1; kk <= free_score_list_.size(); ++kk ) {
					free_data[ kk ] = mut_unbounded.energies().total_energies()[ free_score_list_[ kk ] ];
				}
				for ( Size kk = 1; kk <= fixed_score_list_.size(); ++kk ) {
					fixed_data[ kk ] = mut_unbounded.energies().total_energies()[ fixed_score_list_[ kk ] ];
				}

				SingleStructureDataOP ssd( new SingleStructureData( free_data, fixed_data ) );
				ddg_bind_position_data->add_mutant_unbounds( ssd );
			}


			/// Report the discrepant positions
			TR << "collect_ddG_of_binding_data(): looking for discrepant positions between wt and mut complex structures: ";
			bool one_discrepancy_found( false );
			if ( wt_complex_seq.size() != mut_complex_seq.size() ) {
				std::cerr << "ERROR: Wild type and mutant complex sequences of different length: " << std::endl;
				std::cerr << wt_complex_seq << std::endl;
				std::cerr << mut_complex_seq << std::endl;
				utility_exit();
			} else {
				for ( Size jj = 0; jj < wt_complex_seq.size(); ++jj ) {
					if ( wt_complex_seq[ jj ] != mut_complex_seq[ jj ] ) {
						// we're going to have to store all of the mutations if the get_score method is to work correctly
						// but we want to keep the wt and mut and position together
						// vector of pairs, with first being the position, and second being a pair of wt and mut aa's
						ddg_bind_position_data->add_mutation(
							std::make_pair( jj, std::make_pair( chemical::aa_from_oneletter_code( wt_complex_seq[ jj ] ), chemical::aa_from_oneletter_code( mut_complex_seq[ jj ] ) ) ) );
						TR << wt_complex_seq[ jj ] << jj << mut_complex_seq[ jj ] << ", ";
						one_discrepancy_found = true;
					}
				}
				TR << std::endl;
			}
			if ( one_discrepancy_found == false ) {
				std::cerr << "ERROR: mutant and wild type complex sequences must differ by at least one position";
				utility_exit();
			}

			if ( filtering_bad_ddGs && ( best_mut_complex_rep - best_wt_complex_rep > option[ optE::exclude_badrep_ddGs ]() ) ) {
				TR << "Rank " << MPI_rank_ << " Excluding ddG bind data from " << ddG_bind_files_[ ii ][ DDGBindOptEData::MUT_COMPLEXES_LIST_FILE ]
					<< " with mut-wt rep delta: "<< best_mut_complex_rep << " - " << best_wt_complex_rep
					<< " = " << best_mut_complex_rep - best_wt_complex_rep << std::endl;
			} else {
				ddG_bind_optE_data_->add_position_data( ddg_bind_position_data );
			}
		}
	} // if ddG_bind_optE_data_ == 0

	for ( auto iter = ddG_bind_optE_data_->position_data_begin(),
			iter_end = ddG_bind_optE_data_->position_data_end(); iter != iter_end; ++iter ) {
		optE_data_->add_position_data( *iter );
	}
}

///
/// @brief
/// Set the counts for the amino acid frequencies (observed and expected) to zero.
///
void
IterativeOptEDriver::zero_aa_counts() {
	std::fill( aa_obs_.begin(), aa_obs_.end(), 0 );
	std::fill( aa_exp_.begin(), aa_exp_.end(), 0 );
	std::fill( aa_freq_obs_.begin(), aa_freq_obs_.end(), 0.0 );
	std::fill( aa_freq_exp_.begin(), aa_freq_exp_.end(), 0.0 );
}

///
/// @detail iterate across all the native pdbs,
///
Real
IterativeOptEDriver::measure_sequence_recovery(
	utility::vector1< std::string > const & native_pdb_names,
	utility::vector1< std::string > const & names_for_output_pdbs,
	core::scoring::ScoreFunctionOP sfxn,
	//std::list< core::pack::task::operation::TaskOperationOP > operation_list,
	Size & nresidues_designed,
	Size & nresidues_recovered
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::scoring;

	//sleep( MPI_rank_ );
	//std::cout << "NODE: " << MPI_rank_ << " with reference weight " << sfxn->weights()[ ref ] << " and refE's: ";
	//for ( Size ii = 1; ii <= sfxn->energy_method_options().method_weights( ref ).size(); ++ii ) {
	// std::cout << sfxn->energy_method_options().method_weights( ref )[ ii ] << " ";
	//}
	//std::cout << std::endl;

	nresidues_designed = 0;
	nresidues_recovered = 0;

	//ScoreFunctionOP sfxn2 = ScoreFunctionFactory::create_score_function( get_scorefile_name() );
	using core::pack::task::operation::TaskOperationCOP;
	TaskFactoryOP task_factory_for_design( new TaskFactory( *task_factory_ ) );
	task_factory_for_design->push_back( TaskOperationCOP( new ScaleAnnealerTemperatureOperation( sfxn->weights()[ fa_atr ] / 0.8 ) ) );

	if ( MPI_rank_ == 0 ) {
		if ( sfxn->get_weight( scoring::surface ) != 0.0 ) {
			TR << "measure_sequence_recovery(): designing with surface term, with weight: "
				<< F(8,4, sfxn->get_weight( scoring::surface )) << std::endl;
		}
	}

	protocols::moves::MoverOP design_mover;

	if ( option[ optE::design_with_minpack ] ) {
		protocols::simple_moves::MinPackMoverOP minpack_mover( new protocols::simple_moves::MinPackMover );
		minpack_mover->task_factory( task_factory_for_design );
		minpack_mover->score_function( sfxn );
		design_mover = minpack_mover;
	} else {
		protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
		pack_mover->task_factory( task_factory_for_design );
		pack_mover->score_function( sfxn );
		design_mover = pack_mover;
	}

	for ( Size poseindex = 1; poseindex <= native_pdb_names.size(); ++poseindex ) {
		TR << "begin measure_sequence_recovery(): PDB: " << native_pdb_names[ poseindex ] << std::endl;
		/// read the pdb into a pose
		core::pose::Pose pose;

		// the native_poses_ vector1 gets set in compute_energies,
		// but that doesn't get called if design_first is on the command line
		// so either set the native_poses_ array or load the pdb correctly
		if ( option[ optE::design_first ].user() && outer_loop_counter_ == 1 ) {
			TR << "measure_sequence_recovery(): design_first in use! pushing "
				<< native_pdb_names[poseindex] << " onto native_poses_ vector." << std::endl;
			core::pose::Pose native_pose;
			core::import_pose::pose_from_file( native_pose, native_pdb_names[ poseindex ] , core::import_pose::PDB_file);
			native_poses_.push_back( native_pose );
			context_poses_.push_back( native_pose );

			if ( option[ optE::recover_nat_rot ] ) {
				rotamer_recovery_context_poses_.push_back( native_pose );
			}
		}

		// read in the native pose to do design on
		pose = native_poses_[ poseindex ];

		//if ( option[ in::file::centroid_input ] ) {
		// core::import_pose::centroid_pose_from_pdb( pose, native_pdb_names[ poseindex ] , core::import_pose::PDB_file);
		//} else {
		// core::import_pose::pose_from_file( pose, native_pdb_names[ poseindex ] , core::import_pose::PDB_file);
		//}

		//Real score2 = (*sfxn2)( pose );
		//Real score1 = (*sfxn)( pose );

		//if ( std::abs( score1 - score2 ) > 1e-8 ) {
		// std::cerr << "Score discrepancy for " << native_pdb_names[ poseindex ]
		// << "score1: " << score1 << " score2: " << score2 << std::endl;
		//}

		Size const nresidues_pose( pose.size() );

		/// record original sequence (for all residues in pose)
		utility::vector1< chemical::AA > full_input_sequence( nresidues_pose );
		for ( Size resi = 1; resi <= nresidues_pose; ++resi ) {
			full_input_sequence[ resi ] = pose.residue(resi).aa();
		}

		/// redesign the native pose
		design_mover->apply( pose );

		context_poses_[ poseindex ] = pose; // save this for the next round of optimization

		// use a 'dummy' PackerTask to determine which residues were designed in this pose
		// this should reflect the same PackerTask that PackRotamersMover used for design
		core::pack::task::PackerTaskOP ptask;
		// If desired make sure the task logic is aligned to the native.
		if ( option[ optE::constant_logic_taskops_file ].user() ) {
			ptask = copy_native_packertask_logic( native_poses_[ poseindex ],
				pose,
				task_factory_for_design);
		} else {
			ptask = task_factory_for_design->create_task_and_apply_taskoperations( pose ) ;
		}
		//PackerTaskOP task_for_design = task_factory_for_design->create_task_and_apply_taskoperations( pose );
		//pack_mover->task(task_for_design)


		/// measure seq recov
		for ( Size resi = 1; resi <= nresidues_pose; ++resi ) {
			// do not compile statistics for residues that were not designed
			if ( ! ptask->being_designed(resi) ) continue;
			if ( ! pose.residue(resi).is_protein() ) continue;
			++nresidues_designed;
			++aa_exp_[ full_input_sequence[ resi ]];
			++aa_obs_[ pose.residue(resi).aa() ];
			if ( full_input_sequence[ resi ] == pose.residue(resi).aa() ) ++nresidues_recovered;
		}

		/// write out new pdb for posterity
		if ( ! option[ optE::no_design_pdb_output ] ) {
			pose.dump_scored_pdb( names_for_output_pdbs[ poseindex ], *sfxn );
		}

		// print out score information for the redesign
		(*sfxn)( pose );
		print_energies( pose, sfxn, TR.Trace );

	}

	Real recovery(0.0);
	if ( nresidues_designed != 0 ) {
		recovery = ( static_cast< Real > (nresidues_recovered) ) / nresidues_designed;
		TR_VERBOSE << "measure_sequence_recovery(): recovery: " << recovery << std::endl;
	}

	return recovery;
}

///
Real
IterativeOptEDriver::measure_rotamer_recovery(
	utility::vector1< std::string > const & native_pdb_names,
	utility::vector1< std::string > const & , // names_for_output_pdbs,
	core::scoring::ScoreFunctionOP sfxn,
	//std::list< core::pack::task::operation::TaskOperationOP > operation_list,
	Size & nresidues_repacked,
	Size & nrotamers_recovered
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::dunbrack;

	//sleep( MPI_rank_ );
	//std::cout << "NODE: " << MPI_rank_ << " with reference weight " << sfxn->weights()[ ref ] << " and refE's: ";
	//for ( Size ii = 1; ii <= sfxn->energy_method_options().method_weights( ref ).size(); ++ii ) {
	// std::cout << sfxn->energy_method_options().method_weights(ref )[ ii ] << " ";
	//}
	//std::cout << std::endl;

	nresidues_repacked = 0;
	nrotamers_recovered = 0;

	//ScoreFunctionOP sfxn2 = ScoreFunctionFactory::create_score_function( get_scorefile_name() );

	using namespace core::pack::task;
	using core::pack::task::operation::TaskOperationCOP;
	TaskFactoryOP task_factory_for_repacking( new TaskFactory( *task_factory_ ) );
	task_factory_for_repacking->push_back( TaskOperationCOP( new operation::RestrictToRepacking ) );

	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
	pack_mover->task_factory( task_factory_for_repacking );
	pack_mover->score_function( sfxn );

	for ( Size poseindex = 1; poseindex <= native_pdb_names.size(); ++poseindex ) {
		/// read the pdb into a pose
		core::pose::Pose pose, start_pose;
		pose = native_poses_[ poseindex ];
		start_pose = native_poses_[ poseindex ];
		Size const nresidues_pose( pose.size() );

		//Real score2 = (*sfxn2)( pose );
		//Real score1 = (*sfxn)( pose );
		//if ( std::abs( score1 - score2 ) > 1e-8 ) {
		// std::cerr << "Score discrepancy for " << native_pdb_names[ poseindex ] << "score1: " << score1 << " score2: " << score2 << std::endl;
		//}

		/// repack the pose
		pack_mover->apply( pose );

		rotamer_recovery_context_poses_[ poseindex ] = pose; // save this for the next round of optimization

		/// measure rotamer recovery

		// use a 'dummy' PackerTask to determine which residues were repacked in this pose
		// this should reflect the same PackerTask that PackRotamersMover used for repacking
		core::pack::task::PackerTaskOP ptask;
		// If desired make sure the task logic is aligned to the native.
		if ( option[ optE::constant_logic_taskops_file ].user() ) {
			ptask = copy_native_packertask_logic( native_poses_[ poseindex ],
				pose,
				task_factory_for_repacking);
		} else {
			ptask = task_factory_for_repacking->create_task_and_apply_taskoperations( pose ) ;
		}

		for ( Size resi = 1; resi <= nresidues_pose; ++resi ) {
			// do not compile statistics for residues that were not repacked
			if ( ! ptask->being_packed(resi) ) continue;
			if ( ! pose.residue(resi).is_protein() ) continue;
			if ( start_pose.residue(resi).nchi() == 0 ) continue; // don't count gly/ala in stats.
			++nresidues_repacked;

			RotVector original_rotbins, repacked_rotbins;
			rotamer_from_chi( start_pose.residue(resi), original_rotbins );
			rotamer_from_chi( pose.residue(resi), repacked_rotbins );

			bool all_chi_match( true );
			for ( Size chi_index = 1; chi_index <=original_rotbins.size(); ++chi_index ) {
				if ( original_rotbins[ chi_index ] != repacked_rotbins[ chi_index ] ) {
					all_chi_match = false;
					break;
				}
			}
			if ( all_chi_match ) ++nrotamers_recovered;
		}

		/// don't write out new pdb for posterity
		/// pose.dump_pdb( names_for_output_pdbs[ poseindex ] );

	}

	Real recovery( 0.0 );
	if ( nresidues_repacked != 0 ) {
		recovery = ( static_cast< Real > (nrotamers_recovered) ) / nresidues_repacked;
	}
	return recovery;
}

///
Real
IterativeOptEDriver::opte_weight_mixing_factor( Size outer_loop_counter, Size inner_loop_counter )
{
	if ( outer_loop_counter == 1 ) {
		return 1.0;
	} else if ( inner_loop_counter <= 5 ) {
		return ( 1.0 / ( outer_loop_counter + inner_loop_counter) );
	} else {
		return 0.1;
	}
}

///
/// @brief
/// Reads in the files specified by opt_e::free and opt_e::fixed. Figures out what ScoreType the user placed on each line of the file
/// and then sets the free_parameters array with that ScoreType.  If the user does not place a starting weight, a random starting
/// weight is given for that type.  Also sets the fixed terms in fixed_parameters.  Both of these EnergyMap references that are
/// passed in are actually vectors of EnergyMaps?  Either way, the free and fixed params are set in this method.
/// If no fixed or free files are found, then there are some hardcoded defaults that get used.
///
void
IterativeOptEDriver::initialize_free_and_fixed( core::scoring::EnergyMap & free_parameters, core::scoring::EnergyMap & fixed_parameters )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ optE::free ].user() && option[ optE::fixed ].user() ) {
		utility::io::izstream input_free( option[ optE::free ]() );
		if ( !input_free ) utility_exit_with_message("Couldn't find input file for 'free' parameters");
		Size free_line_number = 1;
		while ( input_free ) {
			utility::vector1< std::string > line_tokens = core::pack::task::tokenize_line( input_free );
			if ( line_tokens.size() == 0 ) {
				// noop
			} else if ( line_tokens.size() == 1 && ! option[ optE::design_first ].user() ) {
				// free value randomized
				ScoreType free_score_type = ScoreTypeManager::score_type_from_name( line_tokens[ 1 ] );
				free_parameters[ free_score_type ] = numeric::random::rg().uniform();
			} else if ( line_tokens.size() == 2 ) {
				ScoreType free_score_type = ScoreTypeManager::score_type_from_name( line_tokens[ 1 ] );
				Real free_starting_weight = utility::from_string( line_tokens[ 2 ], Real(0.0) );
				free_parameters[ free_score_type ] = free_starting_weight;
				if ( option[ optE::randomly_perturb_starting_free_weights ].user() && free_parameters[ free_score_type ] != 0.0 ) {
					Real perturb_range = option[ optE::randomly_perturb_starting_free_weights ]();
					free_parameters[ free_score_type ] += 2 * perturb_range * numeric::random::rg().uniform() - perturb_range;
					if ( free_parameters[ free_score_type ] == 0.0 ) {
						free_parameters[ free_score_type ] = 0.0001; // correct if we should accidentally end up here.
					}
				}
			} else {
				if ( option[ optE::design_first ].user()  ) {
					std::cerr << "\n\n";
					std::cerr << "Error reading weight file line: " << free_line_number << " ";
					for ( Size ii = 1; ii <= line_tokens.size(); ++ii ) {
						std::cerr << line_tokens[ ii ] << " ";
					}
					std::cerr << std::endl << "Expected exactly 2 arguments (i.e. you cannot ask for a random starting weight!) since optE:design_first flag found on command line" << std::endl;
				} else {
					std::cerr << "Error reading free weight file line: " << free_line_number << " ";
					for ( Size ii = 1; ii <= line_tokens.size(); ++ii ) {
						std::cerr << ii << ": " << line_tokens[ ii ];
					}
					std::cerr << std::endl << "Expected only 2 tokens" << std::endl;
				}
				utility_exit();
			}
			++free_line_number;
		}

		utility::io::izstream input_fixed( option[ optE::fixed ]() );
		if ( !input_fixed ) utility_exit_with_message("Couldn't find input file for 'fixed' parameters");
		Size fixed_line_number = 1;
		while ( input_fixed ) {
			utility::vector1< std::string > line_tokens = core::pack::task::tokenize_line( input_fixed );
			if ( line_tokens.size() == 0 ) {
				// noop
			} else if ( line_tokens.size() == 2 ) {
				ScoreType fixed_score_type = ScoreTypeManager::score_type_from_name( line_tokens[ 1 ] );
				Real fixed_weight = utility::from_string( line_tokens[ 2 ], Real(0.0) );
				fixed_parameters[ fixed_score_type ] = fixed_weight;
				if ( free_parameters[ fixed_score_type ] != 0 ) {
					std::cerr << "Error reading free weights file.  Term '" << line_tokens[ 1 ] << "' is listed as both free and fixed.";
					utility_exit();
				}
			} else {
				std::cerr << "Error reading free weight file line: " << free_line_number << " ";
				for ( Size ii = 1; ii <= line_tokens.size(); ++ii ) {
					std::cerr << ii << ": " << line_tokens[ ii ];
				}
				std::cerr << std::endl << "Expected exactly 2 tokens" << std::endl;
				utility_exit();
			}
			++fixed_line_number;
		}

		///  HARD CODED DEFAULTS FOR THOSE THAT LIKE RECOMPILING

	} else if ( option[ in::file::centroid_input ] ) {

		free_parameters[ vdw ] = 1; //optE_RG.uniform();
		free_parameters[ pair ] = 1; //optE_RG.uniform();
		free_parameters[ rama ] = 1; //optE_RG.uniform();
		free_parameters[ p_aa_pp ] = 1; //optE_RG.uniform();
		free_parameters[ cenpack ] = 1; //optE_RG.uniform();

		fixed_parameters[ env ] = 0.4;
	} else {
		/*
		free_parameters[ fa_rep ]      = 1.0;
		free_parameters[ fa_sol ]      = 1.0;
		free_parameters[ fa_dun ]      =  1.0;
		free_parameters[ fa_pair ]     =  1.0;
		free_parameters[ p_aa_pp ]     =  1.0;
		free_parameters[ hbond_bb_sc ] =  1.0;
		free_parameters[ hbond_sc ]    = 1.0;
		free_parameters[ envsmooth ]   = 1.0;
		*/
		free_parameters[ envsmooth ]   = 0.001;


		free_parameters[ fa_rep ]      = 0.44;
		//free_parameters[ fa_sol_apo ]      = 0.65;
		//free_parameters[ fa_sol_chr ]      = 0.65;
		//free_parameters[ fa_sol_pol ]      = 0.65;
		free_parameters[ fa_dun ]      =  0.56;
		//free_parameters[ fa_pair ]     =  0.49;
		free_parameters[ p_aa_pp ]     =  0.64;
		//free_parameters[ hbond_chr_chr ]    = 1.1;
		//free_parameters[ hbond_chr_pol ]    = 1.1;
		//free_parameters[ hbond_pol_pol ]    = 1.1;

		fixed_parameters[ fa_atr ] = 0.8;
		fixed_parameters[ hbond_sr_bb ] = 1.17;
		fixed_parameters[ hbond_lr_bb ] = 1.17;
		fixed_parameters[ dslf_ss_dst ] = 1.0;
		fixed_parameters[ dslf_cs_ang ] = 1.0;
		fixed_parameters[ dslf_ss_dih ] = 1.0;
		fixed_parameters[ dslf_ca_dih ] = 1.0;
	}

	if ( MPI_rank_ == 0 ) {
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( free_parameters[ (ScoreType) ii ] != 0 ) {
				TR_VERBOSE << "initialize_free_and_fixed(): initial free_parameters: " << name_from_score_type( (ScoreType) ii )
					<< " " << free_parameters[ (ScoreType) ii ] << std::endl;
			}
		}
		//fixed_parameters[ vdw ] = 1;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( fixed_parameters[ (ScoreType) ii ] != 0 ) {
				TR_VERBOSE << "initialize_free_and_fixed(): initial fixed_parameters: " << name_from_score_type( (ScoreType) ii )
					<< " " << fixed_parameters[ (ScoreType) ii ] << std::endl;
			}
		}
	}

}

///
/// @brief
/// This function is not used.
///
bool
IterativeOptEDriver::converged(
	core::scoring::EnergyMap & free_parameters_prev,
	core::scoring::EnergyMap & free_parameters_curr,
	utility::vector1< Real > const & reference_energies_prev,
	utility::vector1< Real > const & reference_energies_curr
)
{
	using namespace core::scoring;

	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		runtime_assert( reference_energies_prev.size() == reference_energies_curr.size() );
	}

	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( std::abs( free_parameters_prev[ (ScoreType) ii ] - free_parameters_curr[ (ScoreType) ii ]) > 0.001 ) {
			return false;
		}
	}

	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		for ( Size ii = 1; ii <= reference_energies_prev.size(); ++ii ) {
			if ( std::abs( reference_energies_prev[ ii ] - reference_energies_curr[ ii ]) > 0.001 ) {
				return false;
			}
		}
	}

	TR << "Converged: " << std::endl;
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( free_parameters_prev[ (ScoreType) ii ] != 0 || free_parameters_curr[ (ScoreType) ii ] != 0 ) {
			TR << name_from_score_type( ScoreType( ii ) ) << " prev " << free_parameters_prev[ (ScoreType) ii ] << " curr " << free_parameters_curr[ (ScoreType) ii ] << std::endl;
		}
	}

	if ( ! option[ optE::dont_use_reference_energies ].user() ) {
		for ( Size ii = 1; ii <= reference_energies_prev.size(); ++ii ) {
			TR <<  ii << " prev " << reference_energies_prev[ ii ] << " curr " << reference_energies_curr[ ii ] << std::endl;
		}
	}


	return true;
}

void
IterativeOptEDriver::write_parameters_to_std_out( core::scoring::EnergyMap & free_parameters, utility::vector1< Real > const & reference_energies )
{
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( free_parameters[ (ScoreType) ii ] != 0 ) {
			TR << name_from_score_type( ScoreType( ii ) ) << " " << free_parameters[ (ScoreType) ii ] << std::endl;
		}
	}
	for ( Size ii = 1; ii <= reference_energies.size(); ++ii ) {
		TR << "Reference energy for " <<  ii << " " << reference_energies[ ii ] << std::endl;
	}

}

///
/// @detail
/// This method takes the filenames from the current round and creates new filenames based on what iteration of the
/// outer_loop_counter we're on.  So, each iteration we want to use the newest pdbs.  This method ensures the filename
/// get updated.
///
/// @remark
/// Requires workdir_{0..MPI_numprocs_ - 1 } directories already exist.
///
void
IterativeOptEDriver::setup_pdbnames_next_round(
	Size const outer_loop_counter,
	utility::vector1< std::string  > & pdbs_next_round,
	utility::vector1< std::string > const & native_pdb_names
)
{
	// assumption, native_pdb_names end in ".pdb"
	pdbs_next_round.resize( native_pdb_names.size() );
	for ( Size ii = 1; ii <= native_pdb_names.size(); ++ii ) {
		//std::string native_substr = native_pdb_names[ ii ].substr( 0, native_pdb_names[ ii ].size() - 4 );
		utility::file::FileName natfilename( native_pdb_names[ ii ] );
#ifndef USEMPI
		pdbs_next_round[ ii ] = natfilename.base() + "_" + to_string( outer_loop_counter ) + ".pdb";
#else
		// Write to separate directories.
		pdbs_next_round[ ii ] = "workdir_" + to_string( MPI_rank_ ) + "/" + natfilename.base() + "_" + to_string( outer_loop_counter ) + ".pdb";
#endif

	}
}

void
IterativeOptEDriver::repack_and_minimize_pose(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP sfxn
) const
{
	using namespace moves;
	using namespace core::pack::task;
	using core::pack::task::operation::TaskOperationCOP;

	protocols::simple_moves::PackRotamersMover packer( sfxn );
	TaskFactoryOP factory( new TaskFactory );
	factory->push_back( TaskOperationCOP( new operation::RestrictToRepacking ) );
	factory->push_back( TaskOperationCOP( new operation::IncludeCurrent ) );
	factory->push_back( TaskOperationCOP( new operation::InitializeExtraRotsFromCommandline ) );
	packer.task_factory( factory );

	packer.apply( pose );

	protocols::simple_moves::MinMover minmover;
	minmover.min_type( "lbfgs_armijo_nonmonotone_atol" );
	minmover.score_function( sfxn );

	minmover.apply( pose );

	( *sfxn )( pose );
}

///
/// @details input file should be white-space delimited component-name/weight pairs.
///
void
load_component_weights(
	utility::vector1< core::Real > & component_weights
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	component_weights.resize( n_optE_data_types );
	std::fill( component_weights.begin(), component_weights.end(), 1.0 );

	if ( option[ optE::component_weights ].user() ) {
		utility::io::izstream input( option[ optE::component_weights ]() );
		if ( !input ) utility_exit_with_message("Couldn't find input file for 'component_weights_file' parameters");
		Size line_number = 1; Size nread = 0;
		while ( input ) {
			utility::vector1< std::string > line_tokens = core::pack::task::tokenize_line( input );
			if ( line_tokens.size() == 0 ) {
				// noop
			} else if ( line_tokens.size() == 2 ) {
				if ( ! OptEPositionDataFactory::is_optE_type_name( line_tokens[ 1 ] ) ) {
					utility_exit_with_message("Error reading optE component weights file: Token 1 on line " + utility::to_string( line_number ) + " " + line_tokens[ 1 ] + " is not recognized as an optE data type" );
				}
				OptEPositionDataType component_type = OptEPositionDataFactory::optE_type_from_name( line_tokens[ 1 ] );
				Real component_weight = utility::from_string( line_tokens[ 2 ], Real(1.0) );
				component_weights[ component_type ] = component_weight;
				++nread;
			} else {
				std::cerr << "Error reading optE component weights file line: " << line_number << " ";
				for ( Size ii = 1; ii <= line_tokens.size(); ++ii ) {
					std::cerr << ii << ": " << line_tokens[ ii ];
				}
				std::cerr << std::endl << "Expected exactly 2 tokens" << std::endl;
				utility_exit();
			}
			++line_number;
		}
		if ( line_number == 1 ) {
			TR << "WARNING: read no lines from component weight file: " << option[ optE::component_weights ]() << std::endl;
		} else if ( nread == 0 ) {
			TR << "WARNING: only blank lines found in component weight file: " << option[ optE::component_weights ]() << std::endl;
		}
	}
}

///
/// @brief
/// Copies the logic in the native task factory from the native_pose
/// to the context pose  The context pose should be filled with a parsable
/// file that does NOT restrict any residues that might be packable in the native
core::pack::task::PackerTaskOP
IterativeOptEDriver::copy_native_packertask_logic(core::pose::Pose native_pose,
	core::pose::Pose context_pose,
	core::pack::task::TaskFactoryOP native_taskfactory){
	using namespace core::pack::task;
	TaskFactoryOP context_taskfactory( new TaskFactory );
	std::string context_tagfile( option[ optE::constant_logic_taskops_file ]() );

	read_tagfile_to_taskfactory(context_tagfile, context_taskfactory);
	PackerTaskOP context_task = context_taskfactory->create_task_and_apply_taskoperations( context_pose );

	// Lock the task to the starting pose
	operation::TaskOperationOP mimic_nat_task_op( new operation::ReplicateTask(native_pose, native_taskfactory) );
	mimic_nat_task_op->apply( context_pose, *(context_task) );

	return context_task;
}

///
/// @brief for parallel applications.  Wait at a specific point and stay there until
/// you can attach a gdb process (with the --pid <ID> flag in gdb) and internally
/// modify the variable "i" to some non-zero value with a "set var i = 7" command.
void attach_debugger()
{
#ifdef USEMPI
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
        sleep(5);
#endif
}

///
std::string
IterativeOptEDriver::node_name( int rank ) {

	if ( rank == 0 ) {
		return "master node";
	} else {
		std::stringstream r;
		r << "slave node " << rank;
		return r.str();
	}
}

///
void
IterativeOptEDriver::print_energies(
	pose::Pose & pose,
	scoring::ScoreFunctionOP sfxn,
	std::ostream & os /* = std::cout */
)
{
	scoring::EnergyMap const & wts( sfxn->weights() );
	scoring::EnergyMap const & unweighted_scores( pose.energies().total_energies() );

	os << "---------------------------------------------------" << std::endl;

	// for each energy term, print the weighted energy
	float sum_weighted = 0.0;
	for ( int jj = 1; jj <= scoring::n_score_types; ++jj ) {
		Real const weight = wts[ scoring::ScoreType(jj) ];

		switch( scoring::ScoreType( jj ) ) {
		case scoring::fa_atr:
		case scoring::fa_rep:
		case scoring::fa_sol:
		case scoring::fa_intra_rep:
		case scoring::fa_pair:
		case scoring::hbond_sr_bb:
		case scoring::hbond_lr_bb:
		case scoring::hbond_bb_sc:
		case scoring::hbond_sc:
		case scoring::rama:
		case scoring::omega:
		case scoring::fa_dun:
		case scoring::ref:
		case scoring::p_aa_pp:
		case scoring::pro_close:
		case scoring::surface:
		case scoring::unfolded :
			if ( weight != 0.0 ) {
				Real const val = unweighted_scores[ scoring::ScoreType(jj) ];
				//TR << A(18,scoring::ScoreType(jj)) << ": weight:" << F(5,2,weight) << ", rawE:" << F(5,2,val) << ", weightedE: " << F(5,2, weight * val ) << std::endl;
				os << LJ(18,ScoreType(jj)) << F(7,4,weight) << X(4) << F(10,3,val) << X(4) << F(10,3, weight * val ) << std::endl;
				sum_weighted += weight * val;
			}
			break;
		default :
			break;
		}
	}

	os << "---------------------------------------------------\n" << LJ(25, "Total weighted score: ") << X(12) << F(10,3,sum_weighted) << std::endl;

}

}
}
