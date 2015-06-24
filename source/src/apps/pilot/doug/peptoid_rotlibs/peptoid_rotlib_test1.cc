// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/doug/peptoid_rotlibs/peptoid_rotlib_test1.cc
/// @brief Test the peptoid rotlibs
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// core headers
#include <core/init.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResiduePeptoidLibrary.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.tmpl.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

#include <core/conformation/Residue.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/graph/Graph.hh>

// protocols headers
#include <protocols/simple_moves/PackRotamersMover.hh>

// basic headers
#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// utility headers
#include <utility/io/izstream.hh>

// c++ headers
#include <iostream>
#include <string>

// namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int
main( int argc, char * argv [] )
{
	// init options, rng, etc.
	core::init(argc, argv);

// create score function
	core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );
	score_fxn->set_weight( unfolded, 0.0 );

	// get a ResidueType
	std::cout << "RT" << std::endl;
	std::string aa_name( "P28" );
	ResidueTypeSetCAP rsd_type_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
 	ResidueType const & rsd_type( rsd_type_set->name_map( aa_name ) );
	Residue rsd( rsd_type, true );

	std::cout << "MAKING POSES" << std::endl;
	pose::Pose pose;
	pose.append_residue_by_jump( rsd, 1 );
	pose.append_residue_by_bond( rsd, true );
	pose.append_residue_by_bond( rsd, true );

	pose.set_omega( 1, 176 );
	pose.set_phi( 2, 86 );
	pose.set_psi( 2, 176 );

	std::string filename("pose.pdb");
	pose.dump_scored_pdb( filename, *score_fxn );

	// get a rotlib (basically a small slice of RotamerLibrary::get_peptoid_rotamer_library() and )
	std::cout << "DB" << std::endl;
	Size n_rotlib_chi( rsd_type.nchi() - rsd_type.n_proton_chi() );
	std::string dir_name = basic::database::full_name( "/rotamer/peptoid_rotlibs/" );
	std::string file_name = rsd_type.get_peptoid_rotlib_path();
	utility::io::izstream rotlib_in( dir_name + file_name );
	std::cout << "Reading in rot lib " << dir_name + file_name << "..." << std::endl;

	std::cout << "CREATION" << std::endl;
	//pack::rotamers::RotamericSingleResiduePeptoidLibrary< pack::dunbrack::THREE > peptoid_rotlib;
	pack::rotamers::RotamericSingleResiduePeptoidLibrary< pack::dunbrack::TWO > peptoid_rotlib;

	std::cout << "SET_N_CHI_BINS" << std::endl;
	peptoid_rotlib.set_n_chi_bins( rsd_type.get_peptoid_rotlib_n_bin_per_rot() );

	std::cout << "GET_OMGPHIPSI_BINS" << std::endl;
	Size omg_bin_lower, omg_bin_upper, phi_bin_lower, phi_bin_upper, psi_bin_lower, psi_bin_upper;
	Real omg, phi, psi, omg_alpha, phi_alpha, psi_alpha;

	omg_bin_lower = omg_bin_upper = phi_bin_lower = phi_bin_upper = psi_bin_lower = psi_bin_upper = 0;
	omg = phi = psi = omg_alpha = phi_alpha = psi_alpha = 0;

	for ( omg = -360.00; omg <= 360.00; omg += 1 ) {
		peptoid_rotlib.get_omgphipsi_bins( omg, phi, psi,
			omg_bin_lower, phi_bin_lower, psi_bin_lower,
			omg_bin_upper, phi_bin_upper, psi_bin_upper,
			omg_alpha, phi_alpha, psi_alpha );
	}

	std::cout << "READ_FROM_FILE" << std::endl;
	peptoid_rotlib.read_from_file( rotlib_in );

	std::cout << "MEMORY" << std::endl;
	std::cout << peptoid_rotlib.memory_usage_in_bytes() << std::endl;

	std::cout << "ROTAMER PROB LOOKUP" << std::endl;
	std::cout << "180\t90\t180\t" << peptoid_rotlib.get_probability_for_rotamer(-180, 90, -180, 1) << std::endl; // why doesnt 180 work
	std::cout << "-50\t-90\t50\t" << peptoid_rotlib.get_probability_for_rotamer(-50, -90, 50, 1) << std::endl;
	std::cout << "80\t0\t10\t"    << peptoid_rotlib.get_probability_for_rotamer( 80,   0, 10, 1) << std::endl;

	std::cout << "30\t170\t1700\t"    << peptoid_rotlib.get_probability_for_rotamer( 30, 170, 170, 1) << std::endl;
	std::cout << "30\t170\t1700\t"    << peptoid_rotlib.get_probability_for_rotamer( 30, 170, 170, 2) << std::endl;
	std::cout << "30\t170\t1700\t"    << peptoid_rotlib.get_probability_for_rotamer( 30, 170, 170, 3) << std::endl;
	std::cout << "30\t170\t1700\t"    << peptoid_rotlib.get_probability_for_rotamer( 30, 170, 170, 4) << std::endl;
	std::cout << "30\t170\t1700\t"    << peptoid_rotlib.get_probability_for_rotamer( 30, 170, 170, 5) << std::endl;
	std::cout << "30\t170\t1700\t"    << peptoid_rotlib.get_probability_for_rotamer( 30, 170, 170, 6) << std::endl;

	std::cout << "ROTAMER SAMPLE DATA" << std::endl;
	for ( Size j(1); j <= 2; ++j) {
		pack::dunbrack::DunbrackRotamerSampleData drsa( peptoid_rotlib.get_rotamer(-180, 90, -180, j) ); // why doesnt 180 work
		std::cout << drsa.nrchi_sample() << "\t" << drsa.nchi() << "\t" << drsa.probability() << std::endl;
		for ( Size i(1); i <= 4; ++i ) {
			pack::dunbrack::Size4 rw( drsa.rot_well() ); pack::dunbrack::Real4 cm( drsa.chi_mean() ); pack::dunbrack::Real4 cs( drsa.chi_sd() );
			std::cout << "\t" << rw[i] << "\t" << cm[i] << "\t" << cs[i] << std::endl;
		}
	}

	std::cout << "PACKING" << std::endl;
	pack::task::TaskFactoryOP task_factory = new pack::task::TaskFactory;
	task_factory->push_back( new pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		task_factory->push_back( new pack::task::operation::ReadResfile );
	}
	protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
	pack_mover->task_factory( task_factory );
	pack_mover->score_function( score_fxn );
	std::string before_filename( "before.pdb" );
	pose.dump_scored_pdb( before_filename, *score_fxn );
	pack_mover->apply( pose );

	std::string after_filename( "after.pdb" );
	pose.dump_scored_pdb( after_filename, *score_fxn );
	/*
	std::cout << "FILL ROTAMER VECTOR" << std::endl;
	pack::task::TaskFactoryOP task_factory = new pack::task::TaskFactory;
	task_factory->push_back( new pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		task_factory->push_back( new pack::task::operation::ReadResfile );
	}
	pack::task::PackerTaskOP packer_task( task_factory->create_packer_task( pose ) );

	(*score_fxn)( pose );

	pose.update_residue_neighbors();

	score_fxn->setup_for_packing( pose, packer_task->repacking_residues(), packer_task->designing_residues() );

	graph::GraphOP packer_neighbor_graph( pack::create_packer_graph( pose, *score_fxn, packer_task ) );

	pack::rotamer_set::RotamerSetsOP rotsets;
	rotsets->set_task( packer_task );

	utility::vector1< utility::vector1< Real > > extra_chi_steps( rsd_type.nchi() );

	bool buried(true);

	int nneighbs( 10000 );


	// DOUG DOUG DOUG need to iterator over rotamer_set_ objects using rotamer_sets.begin()/end()

	for ( 	utility::vector1< pack::rotamer_set::RotamerSetOP >::const_iterator jj(rotsets->begin() ); jj != rotsets->end(); ++jj ) {
		for ( Size ii = 1; ii <= rsd_type.nchi(); ++ii ) {
			rotsets[jj]->set_extra_samples( packer_task, nneighbs, ii, rsd_type, extra_chi_steps[ ii ] );
		}
	}

	utility::vector1< ResidueOP > suggested_rotamers;

	peptoid_rotlib.fill_rotamer_vector( pose, score_fxn, packer_task, packer_neighbor_graph, rsd_type, pose.residue(2), extra_chi_steps, buried, suggested_rotamers);


template < Size T >
void
RotamericSingleResiduePeptoidLibrary< T >::fill_rotamer_vector(
x	pose::Pose const & pose,
x	scoring::ScoreFunction const & scorefxn,
X	pack::task::PackerTask const & task,
o	graph::GraphCOP packer_neighbor_graph,
x	chemical::ResidueTypeCAP concrete_residue,
x	conformation::Residue const & existing_residue,
o	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
x	bool buried,
o	RotamerVector & rotamers
) const
{
*/

	std::cout << "DONE" << std::endl;

	return 0;
}
