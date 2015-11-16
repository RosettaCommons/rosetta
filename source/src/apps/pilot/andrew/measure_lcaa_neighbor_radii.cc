// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#include <devel/init.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>

#include <core/id/AtomID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pose_io.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/io/pdb/pdb_dynamic_reader.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <utility/LexicographicalIterator.hh>

OPT_1GRP_KEY( Real, measure, step_size )

int main( int argc, char * argv [] )
{
	try {
		using namespace core;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;
		using namespace core::io::pdb;
		using namespace core::pose;
 		using namespace core::scoring;
		using namespace core::graph;
 		using namespace core::pack;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT( measure::step_size, "The size of the step, in degrees, to use for each chi dihedral", 5 );

		devel::init( argc, argv );

		pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AAA", "fa_standard", true );

		ScoreFunctionOP sfxn = get_score_function();

		PackerTaskOP task = TaskFactory::create_packer_task( pose );
		sfxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( pose, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( pose, *sfxn );

		Real step_size = option[ measure::step_size ];
		Size nsteps = core::Size( 360.0 / step_size );

		RotamerSetOP rotset2 = rotsets->rotamer_set_for_residue( 2 );
		for ( Size ii = 1; ii <= rotset2->get_n_residue_types(); ++ii ) {
		  ResidueOP iirot = rotset2->rotamer( rotset2->get_residue_type_begin(ii) )->clone();
			Real maxd2 = 0;
		  Size nchi = iirot->nchi();
			if ( nchi != 0 ) {
				utility::vector1< Size > dims( nchi, nsteps );
				utility::LexicographicalIterator lex( dims );
				utility::vector1< Real > chi( nchi );
				Vector nbat = iirot->xyz( iirot->nbr_atom() );
				while ( ! lex.at_end() ) {
					for ( Size jj = 1; jj <= nchi; ++jj ) {
						chi[ jj ] = step_size * lex[jj];
					}
					iirot->set_all_chi( chi );
					for ( Size jj = iirot->first_sidechain_atom(), jjend = iirot->nheavyatoms(); jj <= jjend; ++jj ) {
						Real d2 = nbat.distance_squared( iirot->xyz( jj ) );
						if ( maxd2 <= d2 ) {
							maxd2 = d2;
						}
					}
					++lex;
				}
			}
			std::cout << "Max D for " << iirot->name() << " " << std::sqrt(maxd2) << std::endl;
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

