// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <protocols/rigid/RigidBodyMover.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;


void dump_ddg_pdb(std::string const & tag, Size const & resnum, pose::Pose const & pose, std::string const & prefix ) {
	std::ostringstream outPDB_name;
	outPDB_name << prefix << "_" << tag << "_" << resnum << ".pdb";
	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(outPDB_name.str(), std::ios::out);
	pose.dump_pdb( outPDB_stream );
	outPDB_stream.close();
	outPDB_stream.clear();
	return;
}


/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	devel::init(argc, argv);

	std::cout << "Starting ddG calculations" << std::endl;

	bool const output_all_pdbs = true;
	std::string const output_prefix = "ddg";

	pose::Pose bound_pose, native_pose;

	core::import_pose::pose_from_pdb( bound_pose, "bound.pdb" );
	core::import_pose::pose_from_pdb( native_pose, "native.pdb" );

	if ( bound_pose.total_residue() != native_pose.total_residue() ) {
		std::cout << "ERROR - number of residues in bound PDB and native PDB do not match" << std::endl;
	}

	// Setup the unbound pose
	pose::Pose unbound_pose;
	unbound_pose = bound_pose;
	core::Real const unbound_dist = 40.;
	Size const rb_jump = 1; // use the first jump as the one between partners
	protocols::rigid::RigidBodyTransMover trans_mover( unbound_pose, rb_jump );
	trans_mover.trans_axis( trans_mover.trans_axis() );
	trans_mover.step_size(unbound_dist);
	trans_mover.apply( unbound_pose );

	// Setup for scoring/repacking
	scoring::ScoreFunctionOP scorefxn( get_score_function() );

	(*scorefxn)(bound_pose);
	(*scorefxn)(unbound_pose);
	core::Real const starting_bound_score = bound_pose.energies().total_energy();
	core::Real const starting_unbound_score = unbound_pose.energies().total_energy();
	std::cout << "Starting bound score is: " << starting_bound_score << std::endl;
	std::cout << "Starting unbound score is: " << starting_unbound_score << std::endl;

	pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( bound_pose ));
	base_packer_task->set_bump_check( true );
	base_packer_task->initialize_from_command_line();
	base_packer_task->or_include_current( false );

	// Open output file, generate the header line (save it for printing in the log later), print to file
	std::ostringstream header_string_stream;
	header_string_stream << std::setw(9) << "Mutation";
	header_string_stream << std::setw(15) << "dG_WT_bound";
	header_string_stream << std::setw(15) << "dG_WT_unbound";
	header_string_stream << std::setw(9) << "ddG_WT";
	header_string_stream << std::setw(15) << "dG_Ala_bound";
	header_string_stream << std::setw(15) << "dG_Ala_unbound";
	header_string_stream << std::setw(9) << "ddG_Ala";
	std::string ddg_outfname = output_prefix + ".out";
	utility::io::ozstream ddg_outstream;
	ddg_outstream.open(ddg_outfname, std::ios::out);
	ddg_outstream << header_string_stream.str() << std::endl;

	// Loop over all mutated positions, sequentially mutate each residue in bound_pose
	// to Ala and to the native residue, in both the bound and the unbound
	utility::vector1< Size > mutated_residues;
	for ( int ii = 1, nres = bound_pose.total_residue(); ii <= nres; ++ii ) {
		chemical::AA const bound_aa( bound_pose.residue(ii).aa());
		chemical::AA const native_aa( native_pose.residue(ii).aa());
		if ( bound_aa != native_aa ) {
			mutated_residues.push_back(ii);

			std::ostringstream mutation_name;
			mutation_name << oneletter_code_from_aa(native_aa) << ii << oneletter_code_from_aa(bound_aa);
			std::cout << "Processing mutation: " << mutation_name.str() << std::endl;

			// Figure out which residues should be repacked for this seqpos
			pack::task::PackerTaskOP position_packer_task( base_packer_task->clone() );
			utility::vector1<bool> allow_repacked( nres, false );
			allow_repacked.at(ii) = true;
			// Set the neighbors to repack, everything else fixed
			core::scoring::TenANeighborGraph const & graph = bound_pose.energies().tenA_neighbor_graph();
			for ( core::graph::Graph::EdgeListConstIter
							iter = graph.get_node( ii )->const_edge_list_begin(),
							iter_end = graph.get_node( ii )->const_edge_list_end();
						iter != iter_end; ++iter ) {
				Size const neighbor_id( (*iter)->get_other_ind( ii ) );
				position_packer_task->nonconst_residue_task( neighbor_id ).restrict_to_repacking();
			}
			position_packer_task->restrict_to_residues( allow_repacked );

			// Repack the environment in the bound state, save the energy
			pose::Pose repacked_bound_pose;
			repacked_bound_pose = bound_pose;
			pack::task::PackerTaskOP repacked_packer_task( position_packer_task->clone() );
			utility::vector1< bool > repack_aalist( chemical::num_canonical_aas, false );
			repack_aalist[ bound_aa ] = true;
			repacked_packer_task->nonconst_residue_task(ii).restrict_absent_canonical_aas( repack_aalist );
			pack::pack_rotamers( repacked_bound_pose, *scorefxn, repacked_packer_task );
			(*scorefxn)(repacked_bound_pose);
			core::Real const repacked_bound_score = repacked_bound_pose.energies().total_energy();
			if ( output_all_pdbs ) dump_ddg_pdb("DES_BOUND", ii, repacked_bound_pose, output_prefix);

			// Repack the same environment in the unbound state, save the energy
			pose::Pose repacked_unbound_pose;
			repacked_unbound_pose = unbound_pose;
			pack::pack_rotamers( repacked_unbound_pose, *scorefxn, repacked_packer_task );
			(*scorefxn)(repacked_unbound_pose);
			core::Real const repacked_unbound_score = repacked_unbound_pose.energies().total_energy();
			if ( output_all_pdbs ) dump_ddg_pdb("DES_UNBOUND", ii, repacked_unbound_pose, output_prefix);

			// Mutate to WT in the bound state, save the energy
			pose::Pose reverted_bound_pose;
			reverted_bound_pose = bound_pose;
			pack::task::PackerTaskOP reverted_packer_task( position_packer_task->clone() );
			utility::vector1< bool > native_aalist( chemical::num_canonical_aas, false );
			native_aalist[ native_aa ] = true;
			reverted_packer_task->nonconst_residue_task(ii).restrict_absent_canonical_aas( native_aalist );
			pack::pack_rotamers( reverted_bound_pose, *scorefxn, reverted_packer_task );
			(*scorefxn)(reverted_bound_pose);
			core::Real const reverted_bound_score = reverted_bound_pose.energies().total_energy();
			if ( output_all_pdbs ) dump_ddg_pdb("NAT_BOUND", ii, reverted_bound_pose, output_prefix);

			// Mutate to WT in the unbound state, save the energy
			pose::Pose reverted_unbound_pose;
			reverted_unbound_pose = unbound_pose;
			pack::pack_rotamers( reverted_unbound_pose, *scorefxn, reverted_packer_task );
			(*scorefxn)(reverted_unbound_pose);
			core::Real const reverted_unbound_score = reverted_unbound_pose.energies().total_energy();
			if ( output_all_pdbs ) dump_ddg_pdb("NAT_UNBOUND", ii, reverted_unbound_pose, output_prefix);

			// Mutate to Ala in the bound state, save the energy
			pose::Pose AlaMut_bound_pose;
			AlaMut_bound_pose = bound_pose;
			pack::task::PackerTaskOP AlaMut_packer_task( position_packer_task->clone() );
			utility::vector1< bool > ala_aalist( chemical::num_canonical_aas, false );
			ala_aalist[ chemical::aa_ala ] = true;
			AlaMut_packer_task->nonconst_residue_task(ii).restrict_absent_canonical_aas( ala_aalist );
			pack::pack_rotamers( AlaMut_bound_pose, *scorefxn, AlaMut_packer_task );
			(*scorefxn)(AlaMut_bound_pose);
			core::Real const AlaMut_bound_score = AlaMut_bound_pose.energies().total_energy();
			if ( output_all_pdbs ) dump_ddg_pdb("ALA_BOUND", ii, AlaMut_bound_pose, output_prefix);

			// Mutate to Ala in the unbound state, save the energy
			pose::Pose AlaMut_unbound_pose;
			AlaMut_unbound_pose = unbound_pose;
			pack::pack_rotamers( AlaMut_unbound_pose, *scorefxn, AlaMut_packer_task );
			(*scorefxn)(AlaMut_unbound_pose);
			core::Real const AlaMut_unbound_score = AlaMut_unbound_pose.energies().total_energy();
			if ( output_all_pdbs ) dump_ddg_pdb("ALA_UNBOUND", ii, AlaMut_unbound_pose, output_prefix);

			// Compute differences
			core::Real const reversion_dG_bound = reverted_bound_score - repacked_bound_score;
			core::Real const reversion_dG_unbound = reverted_unbound_score - repacked_unbound_score;
			core::Real const reversion_ddG = reversion_dG_bound - reversion_dG_unbound;
			core::Real const AlaMut_dG_bound = AlaMut_bound_score - repacked_bound_score;
			core::Real const AlaMut_dG_unbound = AlaMut_unbound_score - repacked_unbound_score;
			core::Real const AlaMut_ddG = AlaMut_dG_bound - AlaMut_dG_unbound;

			// Generate string summary for this seqpos
			std::ostringstream data_string_stream;
			data_string_stream << std::setw(9) << mutation_name.str();
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << reversion_dG_bound;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << reversion_dG_unbound;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(9) << reversion_ddG;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << AlaMut_dG_bound;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << AlaMut_dG_unbound;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(9) << AlaMut_ddG;

			// Print to file
			ddg_outstream << data_string_stream.str() << std::endl;

			// Print to log as well, repeat the header line
			std::cout << header_string_stream.str() << std::endl;
			std::cout << data_string_stream.str() << std::endl;

		}
	}

	ddg_outstream.close();
	ddg_outstream.clear();

	std::cout << "Successfully finished computing ddGs" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}

