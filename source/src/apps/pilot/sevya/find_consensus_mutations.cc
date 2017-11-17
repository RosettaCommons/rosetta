// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/sevya/msd.cc
/// @brief  Experimental multi state design implementation. Written to encourage convergence in protein designs occurring
/// simultaneously, rather than enforce that they have identical sequences.
/// @author Alex Sevy


//core library
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/NegativePackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <devel/init.hh>
#include <protocols/protein_interface_design/design_utils.hh>
// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/io/izstream.hh>
#include <utility/json_spirit/json_spirit_reader.h>
#include <protocols/jobdist/standard_mains.hh>
#include <basic/Tracer.hh>
#include <sstream>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>


namespace basic { namespace options { namespace OptionKeys { namespace msd
{
basic::options::StringVectorOptionKey positive_states( "msd:positive_states" );
basic::options::StringVectorOptionKey negative_states( "msd:negative_states" );
basic::options::StringOptionKey upper_chains( "msd:upper_chains" );
basic::options::StringOptionKey lower_chains( "msd:lower_chains" );
basic::options::StringOptionKey corr( "msd:corr" );
}}}}


///////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR("msd_pilot.main");

std::string
get_out_number(core::Size in) {
	std::stringstream num ("");
	if ( in < 10 ) {
		num << "000" << in;
	} else if ( in < 100 ) {
		num << "00" << in;
	} else if ( in < 1000 ) {
		num << "0" << in;
	} else {
		num << in;
	}
	return num.str();
}

utility::vector1< core::Size >
parse_resfile ( core::pack::task::PackerTaskCOP design_task )
{
	utility::vector1< core::Size > vector;
	utility::vector1<bool> designing = design_task->designing_residues();
	for ( core::Size i = 1; i <= designing.size(); ++i ) {
		if ( designing[ i ] ) {
			vector.push_back( i );
		}
	}
	return vector;
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		//Add local options
		option.add( msd::positive_states, "Positive states for MSD").def("");
		option.add( msd::negative_states, "Negative states for MSD").def("");

		devel::init( argc, argv );

		utility::vector1< std::string > pdb_positive, pdb_negative;
		if ( option[ msd::positive_states ].user() ) {
			pdb_positive = option[ msd::positive_states ];
		}
		if ( option[ msd::negative_states ].user() ) {
			pdb_negative = option[ msd::negative_states ];
		}
		if ( pdb_positive.size() + pdb_negative.size() < 2 ) {
			utility_exit_with_message( "Expected at least two pdbs to be specified from the -msd:positive_states and -msd:negative_states flags" );
		}


		//Poses plus bool representing whether it's a positive state
		utility::vector1<std::pair < core::pose::PoseOP, bool > > poses;

		for ( core::Size i = 1; i <= pdb_positive.size(); ++i ) {
			poses.push_back(std::make_pair( core::import_pose::pose_from_file( pdb_positive[ i ], false), 1 ) , core::import_pose::PDB_file);
		}
		for ( core::Size i = 1; i <= pdb_negative.size(); ++i ) {
			poses.push_back(std::make_pair( core::import_pose::pose_from_file( pdb_negative[ i ], false), 0 ) , core::import_pose::PDB_file);
		}

		core::Size nstruct = 1;
		if ( option[ out::nstruct ].user() ) {
			nstruct = option[ out::nstruct ];
		}
		std::string suffix ("");
		if ( option[ out::suffix ].user() ) {
			suffix = option[ out::suffix ];
		}
		/* Set up task operations and assign to the right task factory */
		core::pack::task::operation::InitializeFromCommandlineOP ifcl = new core::pack::task::operation::InitializeFromCommandline;
		core::pack::task::operation::ReadResfileOP rrf = new core::pack::task::operation::ReadResfile;
		core::pack::task::operation::RestrictToRepackingOP rtr = new core::pack::task::operation::RestrictToRepacking;

		core::pack::task::TaskFactoryOP design_task_factory = new core::pack::task::TaskFactory;
		design_task_factory->push_back( ifcl );
		design_task_factory->push_back( rrf );

		core::pack::task::PackerTaskOP design_task = design_task_factory->create_task_and_apply_taskoperations( *(poses[1].first) );

		core::pack::task::TaskFactoryOP packing_min_task_factory = new core::pack::task::TaskFactory;
		packing_min_task_factory->push_back( ifcl );
		packing_min_task_factory->push_back( rtr );
		packing_min_task_factory->push_back( rrf );

		/* Create score functions */
		core::scoring::ScoreFunctionOP sfxn_clean = core::scoring::getScoreFunction();
		core::scoring::ScoreFunctionOP soft_rep = core::scoring::ScoreFunctionFactory::create_score_function("soft_rep_design");
		core::scoring::ScoreFunctionOP sfxn_design = sfxn_clean->clone();
		soft_rep->set_weight( core::scoring::res_type_constraint, 1.0 );
		sfxn_design->set_weight( core::scoring::res_type_constraint, 1.0 );

		protocols::simple_moves::PackRotamersMoverOP rp = new protocols::simple_moves::PackRotamersMover;
		rp->score_function( sfxn_clean );
		rp->task_factory( packing_min_task_factory );

		utility::vector1< std::pair<std::string, core::Size> > res_links;
		utility::vector1< core::Size > resfile_res_links;
		resfile_res_links = parse_resfile( design_task );


		utility::vector1< std::pair < core::pose::PoseOP, bool > > modifiable_poses;
		for ( core::Size n = 1; n <= nstruct; ++n ) {
			modifiable_poses.clear();

			/* Create copy of each pose so original list stays unchanged */
			for ( core::Size i = 1; i <= poses.size(); ++i ) {
				modifiable_poses.push_back( std::make_pair( poses[i].first->clone(), poses[i].second ) );
			}

			TR << "finding consensus mutations" << std::endl;
			/* Revert mutations to find the best aa at each spot */
			core::Real scaling_factor = ( (core::Real) pdb_positive.size() )/pdb_negative.size();
			TR << "scaling factor: " << scaling_factor << std::endl;
			for ( core::Size i = 1; i <= resfile_res_links.size(); ++i ) {
				core::Size seqpos = resfile_res_links[ i ];
				TR << seqpos << std::endl;
				//check to see if they are all the same
				bool diff = false;
				TR << modifiable_poses[ 1 ].first->residue( seqpos ).name3() << std::endl;
				for ( core::Size j = 2; j <= modifiable_poses.size(); ++j ) {
					TR << modifiable_poses[ j ].first->residue( seqpos ).name3() << std::endl;
					if ( modifiable_poses[ j ].first->residue( seqpos ).aa() != modifiable_poses[ 1 ].first->residue( seqpos ).aa() ) {
						diff = true;
						break;
					}
				}
				TR << diff << std::endl;
				if ( diff ) {
					TR << "position " << seqpos << " is different" << std::endl;
					core::Size min_index = 0;
					core::Real min_score = 0;
					utility::vector1< core::Real > scores ( modifiable_poses.size(), 0 );
					for ( core::Size ref = 1; ref <= modifiable_poses.size(); ++ref ) {
						for ( core::Size comp = 1; comp <= modifiable_poses.size(); ++comp ) {
							TR << "comparing " << modifiable_poses[ ref ].first->residue( seqpos ).name3() << " to " << modifiable_poses[ comp ].first->residue( seqpos ).name3() << std::endl;
							if ( modifiable_poses[ ref ].first->residue( seqpos ).aa() ==
									modifiable_poses[ comp ].first->residue( seqpos ).aa() ) {
								scores[ ref ] += (modifiable_poses[ comp ].second ?
									modifiable_poses[ comp ].first->energies().total_energy() :
									scaling_factor*-modifiable_poses[ comp ].first->energies().total_energy() );
								//         continue;
							} else {
								core::pose::PoseOP mutpose = modifiable_poses[ comp ].first->clone();
								mutpose->replace_residue( seqpos, modifiable_poses[ ref ].first->residue( seqpos ), true );
								rp->apply( *mutpose );
								scores[ ref ] += (modifiable_poses[ comp ].second ?
									mutpose->energies().total_energy() :
									scaling_factor*-mutpose->energies().total_energy() );
							}

						}
						if ( scores[ ref ] < min_score || ref == 1 ) {
							min_index = ref;
							min_score = scores[ ref ];
						}
					}
					for ( core::Size k = 1; k <= scores.size(); ++k ) TR << scores[ k ] << std::endl;
					TR << "best residue at position " << seqpos << ": " << modifiable_poses[ min_index ].first->residue( seqpos ).name3() << std::endl;
					for ( core::Size pose = 1; pose <= modifiable_poses.size(); ++pose ) {
						if ( modifiable_poses[ pose ].first->residue( seqpos ).aa() !=
								modifiable_poses[ min_index ].first->residue( seqpos ).aa() ) {
							modifiable_poses[ pose ].first->replace_residue( seqpos, modifiable_poses[ min_index ].first->residue( seqpos ), true );
							rp->apply( *modifiable_poses[ pose ].first );
						}
					}
				}
			}


			/* Finished - dump out the scored pdbs to the output file */
			core::Size j = 1;
			for ( core::Size i = 1; i <= pdb_positive.size(); ++i ) {
				std::string file_name = pdb_positive[ i ];
				file_name = file_name.substr(0, file_name.length()-4);
				std::string out_name = (suffix=="") ?
					file_name + "_" + get_out_number(n) + ".consensus.pdb" :
					file_name + "_" + suffix + "_" + get_out_number(n) + ".consensus.pdb";
				modifiable_poses[ j ].first->dump_scored_pdb( out_name, *sfxn_clean );
				j++;
			}
			for ( core::Size i = 1; i <= pdb_negative.size(); ++i ) {
				std::string file_name = pdb_negative[i];
				file_name = file_name.substr(0, file_name.length()-4);
				std::string out_name = (suffix=="") ?
					file_name + "_" + get_out_number(n) + ".consensus.pdb" :
					file_name + "_" + suffix + "_" + get_out_number(n) + ".consensus.pdb";
				modifiable_poses[ j ].first->dump_scored_pdb( out_name, *sfxn_clean );
				j++;
			}
		}


	} catch (utility::excn::Exception const & e ) {
		utility_exit_with_message("caught exception " + e.msg());
	}
}
