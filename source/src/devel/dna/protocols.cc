// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <devel/dna/protocols.hh>
#include <devel/dna/util.hh>

#include <protocols/dna/util.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

//#include <protocols/moves/Mover.hh>

// ObjexxFCL Headers

// Utility Headers
#include <utility/vector1.hh>

// c++ headers
#include <iostream>
#include <string>
#include <vector>

#include <utility/vector0.hh>
#include <ObjexxFCL/format.hh>


using utility::vector1;

namespace devel {
namespace dna {

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

// file-scope -- prob bad
static THREAD_LOCAL basic::Tracer tt( "devel.dna.protocols", basic::t_info );


////////////////////////////////////////////////////////////////////////
void
repack_base_pair_neighbors(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn,
	Size const seqpos,
	bool const include_current,
	bool const repack_dna
)
{
	using namespace scoring;
	using namespace scoring::dna;

	Size const nloop( 25 );

	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );
	Size const seqpos_partner( partner[ seqpos ] );

	Size const nres( pose.total_residue() );
	utility::vector1< bool > is_base_pair_neighbor( nres, false );
	for ( Size i=1; i<= nres; ++i ) {
		if ( scorefxn.are_they_neighbors( pose, i, seqpos ) ||
				( seqpos_partner > 0 && scorefxn.are_they_neighbors( pose, i, seqpos_partner ) ) ) {
			is_base_pair_neighbor[i] = true;
		}
	}


	//// setup the packer task
	//
	// pack at protein positions that are neighbors
	// pack at dna positions that are neighbors IF repack_dna is true
	//
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );
	task->initialize_from_command_line();
	if ( include_current ) task->or_include_current( true );
	else                   task->or_include_current( false );

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() ) {
			if ( is_base_pair_neighbor[i] ) {
				task->nonconst_residue_task( i ).restrict_to_repacking();
			} else {
				task->nonconst_residue_task( i ).prevent_repacking();
			}
		} else if ( pose.residue(i).is_DNA() ) {
			if ( repack_dna && is_base_pair_neighbor[i] ) {
				task->nonconst_residue_task( i ).restrict_to_repacking();
			} else {
				task->nonconst_residue_task( i ).prevent_repacking();
			}
		}
	}


	// now call pack_rotamers
	scorefxn(pose);
	vector1< std::pair< Real, std::string > > results;
	pack::pack_rotamers_loop( pose, scorefxn, task, nloop, results );
}

///////////////////////////////////////////////////////////////////////////////////////////////
void
packing_specificity_test_fast(
	pose::Pose const & start_pose,
	scoring::ScoreFunction const & scorefxn,
	Size const motif_begin,
	Size const motif_size,
	Size const nloop,
	std::string const & min_type,
	Real const min_tol,
	bool const postmin,
	std::string const & output_tag,
	bool const add_extra_rotamers, // =false
	bool const dump_pdbs // = true
)
{
	vector1< Size > motif_positions;
	for ( Size i=motif_begin; i< motif_begin+motif_size; ++i ) {
		motif_positions.push_back(i);
	}
	packing_specificity_test_fast( start_pose, scorefxn, motif_positions, nloop, min_type, min_tol, postmin, output_tag, add_extra_rotamers, dump_pdbs );
}
/// @details Try all possible dna basepairs at the motif positions
/// evaluate their energies with a repack
/// requires that the pose already have base partner info set


void
packing_specificity_test_fast(
	pose::Pose const & start_pose,
	scoring::ScoreFunction const & scorefxn,
	vector1< Size > const & motif_positions,
	Size const nloop,
	std::string const & min_type,
	Real const min_tol,
	bool const postmin,
	std::string const & output_tag,
	bool const add_extra_rotamers, // =false
	bool const dump_pdbs // = true
)
{
	using namespace scoring;
	using namespace chemical;
	using namespace pose;
	using namespace conformation;
	using namespace optimization;
	using namespace scoring::dna;

	basic::prof_reset();

	tt << "packing_specificity_test::  nmotif_positions= " << motif_positions.size() <<
		" min_type= " << min_type << " scorefxn: " << scorefxn << std::endl;

	Size const nres( start_pose.total_residue() );
	BasePartner const & partner( retrieve_base_partner_from_pose( start_pose ) );

	// setup the positions
	vector1< Size > pos_list;
	std::string nat_seq;
	{
		for ( Size ii=1; ii<= motif_positions.size(); ++ii ) {
			Size const i( motif_positions[ii] );
			assert( start_pose.residue(i).is_DNA() && partner[i] );
			pos_list.push_back( i );
			pos_list.push_back( partner[ i ] );
			nat_seq += start_pose.residue(i).name1();
		}
	}


	{ // now repack the pose, allowing dna coupled rotamer switches
		Pose pose;
		pose = start_pose;

		// figure out which protein positions are within contact distance of the DNA
		Real const dis2_cutoff( 18 * 18 );
		vector1< bool > allow_repack( nres, false );
		for ( Size i=1; i<= nres; ++i ) {
			Residue const & i_rsd( pose.residue(i) );
			if ( !i_rsd.is_protein() ) continue;
			bool contact( false );
			for ( Size j=1; j<= nres && !contact; ++j ) {
				Residue const & j_rsd( pose.residue(j) );
				if ( !j_rsd.is_DNA() ) continue;

				for ( Size ii=1; ii<= i_rsd.natoms() && !contact; ++ii ) {
					for ( Size jj=1; jj<= j_rsd.natoms() && !contact; ++jj ) {
						if ( i_rsd.xyz(ii).distance_squared( j_rsd.xyz(jj) ) < dis2_cutoff ) contact = true;
					}
				}
			}
			if ( contact ) allow_repack[i] = true;
		}

		//setup task + minimizer movemap
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line();
		tt << "packloop: " << pose.sequence() << " START" << std::endl;

		kinematics::MoveMap mm;
		for ( Size ii = 1; ii <= nres; ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				if ( allow_repack[ii] ) {
					task->nonconst_residue_task( ii ).restrict_to_repacking();
					mm.set_chi( ii, true );
					if ( add_extra_rotamers ) {
						task->nonconst_residue_task(ii).or_ex1( true );
						task->nonconst_residue_task(ii).or_ex2( true );
						task->nonconst_residue_task(ii).or_ex1aro( true );
						task->nonconst_residue_task(ii).or_ex1aro_sample_level( pack::task::EX_THREE_THIRD_STEP_STDDEVS );
						task->nonconst_residue_task(ii).and_extrachi_cutoff( 0 );
					}
				} else {
					task->nonconst_residue_task( ii ).prevent_repacking();
				}
			} else {
				if ( std::find( pos_list.begin(), pos_list.end(), ii ) == pos_list.end() ) {
					task->nonconst_residue_task( ii ).prevent_repacking();
				} else {
					task->nonconst_residue_task( ii ).allow_aa( na_ade );
					task->nonconst_residue_task( ii ).allow_aa( na_thy );
					task->nonconst_residue_task( ii ).allow_aa( na_gua );
					task->nonconst_residue_task( ii ).allow_aa( na_cyt );
					if ( add_extra_rotamers ) {
						// 100 rotamers
						task->nonconst_residue_task( ii ).or_exdna_sample_level( static_cast<pack::task::ExtraRotSample>(4));
					}
				}
			}
		}

		{ // setup residue couplings
			using namespace pack::rotamer_set;
			RotamerCouplingsOP couplings( new RotamerCouplings() );
			couplings->resize( nres );
			for ( Size i=1; i<= nres; ++i ) {
				if ( partner[i] ) {
					(*couplings)[i].first = partner[i];
					(*couplings)[i].second = conformation::ResidueMatcherOP( new conformation::WatsonCrickResidueMatcher() );
				}
			}
			task->rotamer_couplings( couplings );
		}


		// call PACK_ROTAMERS !!!!!!!!!!!!!!!
		scorefxn(pose);
		vector1< std::pair< Real, std::string > > results;
		pack::pack_rotamers_loop( pose, scorefxn, task, nloop, results );
		Energy pack_score = scorefxn( pose );

		//scorefxn.accumulate_residue_total_energies( pose );
		//pose.energies().show( std::cout );

		{ // test for mismatches
			WatsonCrickResidueMatcher m;
			for ( Size i=1; i<= nres; ++i ) {
				if ( partner[i]>i ) {
					assert( m( pose.residue(i), pose.residue(partner[i])));
				}
			}
		}

		// what's the sequence?
		std::string packed_seq;
		for ( Size i=1; i<= motif_positions.size(); ++i ) {
			packed_seq += pose.residue( motif_positions[i] ).name1();
		}

		tt << "packed_score: " << packed_seq << ' ' << nat_seq << ' ' << ObjexxFCL::format::F(9,3,pack_score);

		// calculate frequency correct at each position
		Real total_correct(0.0);
		for ( Size i=1; i<= motif_positions.size(); ++i ) {
			Real correct(0.0);
			for ( Size n=1; n<= results.size(); ++n ) {
				correct += ( results[n].second[ motif_positions[i]-1 ] == nat_seq[i-1] );
			}
			correct /= results.size();
			tt << ObjexxFCL::format::F(6,3,correct);
			total_correct += correct;
		}
		total_correct /= motif_positions.size();

		tt << " total_correct= " << ObjexxFCL::format::F(6,3,total_correct) << ' ' << output_tag << std::endl <<
			"packed_terms: " << output_tag << ' ' << pose.energies().total_energies().show_nonzero() << std::endl;

		if ( dump_pdbs ) io::pdb::dump_pdb( pose, output_tag+"_"+packed_seq+"_"+nat_seq+".pdb" );

		if ( postmin ) {
			//MinimizerOptions min_options( min_type, min_tol, true /*use_nblist*/ );
			AtomTreeMinimizer().run( pose, mm, scorefxn, MinimizerOptions( min_type, min_tol, true ) );//min_options );
			tt << "postmin_score: " << packed_seq << ' ' << scorefxn( pose ) << ' ' << output_tag << std::endl <<
				"postmin_terms: " << output_tag << ' ' << pose.energies().total_energies().show_nonzero() << std::endl;
			if ( dump_pdbs ) io::pdb::dump_pdb( pose, output_tag+"_"+packed_seq+"_"+nat_seq+"_postmin.pdb" );
		}

	} // loop over sequence combinations
	basic::prof_show();
}

//   pack::task::PackerTaskOP
//    pack_task( pack::task::TaskFactory::create_packer_task( pose )),
//    design_task( pack::task::TaskFactory::create_packer_task( pose )),
//    coupled_design_task( pack::task::TaskFactory::create_packer_task( pose ));

//   pack_task->initialize_from_command_line();
//   design_task->initialize_from_command_line();
//   coupled_design_task->initialize_from_command_line();
//   tt << "packloop: " << pose.sequence() << " START" << std::endl;

//   for ( Size ii = 1; ii <= nres; ++ii ) {
//    if ( pose.residue(ii).is_protein() ) {
//     pack_task->nonconst_residue_task( ii ).restrict_to_repacking();
//     design_task->nonconst_residue_task( ii ).restrict_to_repacking();
//     coupled_design_task->nonconst_residue_task( ii ).restrict_to_repacking();
//    } else {
//     pack_task->nonconst_residue_task( ii ).prevent_repacking();
//     design_task->nonconst_residue_task( ii ).allow_aa( na_ade );
//     design_task->nonconst_residue_task( ii ).allow_aa( na_thy );
//     design_task->nonconst_residue_task( ii ).allow_aa( na_gua );
//     design_task->nonconst_residue_task( ii ).allow_aa( na_cyt );
//     coupled_design_task->nonconst_residue_task( ii ).allow_aa( na_ade );
//     coupled_design_task->nonconst_residue_task( ii ).allow_aa( na_thy );
//     coupled_design_task->nonconst_residue_task( ii ).allow_aa( na_gua );
//     coupled_design_task->nonconst_residue_task( ii ).allow_aa( na_cyt );
//    }
//   }

/// @details Try all possible dna basepairs at the motif positions
/// evaluate their energies with a repack
/// requires that the pose already have base partner info set


void
packing_specificity_test(
	pose::Pose const & start_pose,
	scoring::ScoreFunction const & scorefxn,
	Size const motif_begin,
	Size const motif_size,
	std::string const & min_type,
	Real const min_tol,
	bool const postmin,
	std::string const & output_tag,
	bool const repack_DNA// = false
)
{
	using namespace scoring;
	using namespace chemical;
	using namespace pose;
	using namespace conformation;
	using namespace optimization;
	using namespace scoring::dna;

	basic::prof_reset();

	tt << "packing_specificity_test::  motif_begin= " << motif_begin << " motif_size= " << motif_size <<
		" min_type= " << min_type << " min_tol= " << min_tol << " scorefxn: " << scorefxn << std::endl;

	Size const nres( start_pose.total_residue() );
	BasePartner const & partner( retrieve_base_partner_from_pose( start_pose ) );

	// setup the positions
	vector1< Size > pos_list;
	std::string nat_seq;
	{
		for ( Size i=motif_begin; i< motif_begin+motif_size; ++i ) {
			assert( start_pose.residue(i).is_DNA() && partner[i] );
			pos_list.push_back( i );
			pos_list.push_back( partner[ i ] );
			nat_seq += start_pose.residue(i).name1();
		}
	}

	// setup a movemap for protein sidechain minimization
	kinematics::MoveMap mm;

	for ( Size ii = 1; ii <= nres; ++ii ) {
		if ( start_pose.residue(ii).is_protein() ) mm.set_chi( ii, true );
	}


	// setup minimizer options
	MinimizerOptions min_options( min_type, min_tol, true /*use_nblist*/ );


	// score the starting structure
	if ( false ) {
		Pose pose;
		pose = start_pose;

		Real const start_score( scorefxn( pose ) );

		// minimize
		AtomTreeMinimizer().run( pose, mm, scorefxn, min_options );

		Real const start_min_score( scorefxn( pose ) );

		pack::task::PackerTaskOP task
			( pack::task::TaskFactory::create_packer_task( pose ));

		task->initialize_from_command_line();
		for ( Size ii = 1; ii <= nres; ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
				assert( task->pack_residue(ii) );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
		// dont include current
		task->or_include_current( false );

		pack::pack_rotamers( pose, scorefxn, task);
		Energy const start_min_pack_score = scorefxn( pose );

		tt << "start scores: " << nat_seq << ' ' << start_score << ' ' << start_min_score << ' ' <<
			start_min_pack_score << std::endl << pose.energies().total_energies().show_nonzero() << std::endl;
	}


	// setup for the combinatorics
	vector1< vector1< AA > > aa_combinations;
	{

		aa_combinations.push_back( vector1< AA >() );
		while ( aa_combinations[1].size() < pos_list.size() ) {

			vector1< vector1< AA > > new_combinations;
			for ( Size i=1; i<= aa_combinations.size(); ++i ) {
				for ( int j = first_DNA_aa; j<= last_DNA_aa; ++j ) {
					new_combinations.push_back( aa_combinations[i] );
					vector1< AA > & l( new_combinations[ new_combinations.size() ] );
					l.push_back( AA(j) );
					l.push_back( protocols::dna::dna_base_partner( AA(j) ) );
				}
			}
			aa_combinations.swap( new_combinations );
		}
	}

	ResidueTypeSetCOP residue_set( start_pose.residue(1).residue_type_set() );

	tt << "Motif size = " << motif_size << ". Trying " << aa_combinations.size() << " different motif sequences" <<
		std::endl;

	for ( Size ncombo=1; ncombo<= aa_combinations.size(); ++ncombo ) {

		Pose pose;
		pose = start_pose;

		vector1< AA > const & combo( aa_combinations[ncombo] );
		std::string combo_tag, nat_seq;

		// apply the sequence modifications
		for ( Size i=1; i<= pos_list.size(); ++i ) {
			AA const & aa( combo[i] );
			int const seqpos( pos_list[i] );

			// Representative type should have no/minimal variants
			ResidueTypeCOP rsd_type( residue_set->get_representative_type_aa( aa ) );

			Residue const & existing_residue( pose.residue( seqpos ) );
			assert( existing_residue.is_DNA() );

			ResidueOP rsd = ResidueFactory::create_residue( *rsd_type, existing_residue, pose.conformation() );
			rsd->set_chi( 1, existing_residue.chi(1) );

			//tt << "replace rsd: " << aa << ' ' << seqpos << ' ' << rsd->name() << std::endl;
			pose.replace_residue( seqpos, *rsd, false );

			if ( i%2 == 1 ) {
				combo_tag += rsd->name1();
				nat_seq += start_pose.residue( seqpos ).name1();
			}
		}


		tt << "score before packing: " << scorefxn( pose ) << std::endl; // hack -- needed to fill tenA graph

		tt << "packing with motif sequence: " << output_tag << ' ' << combo_tag << std::endl;
		pack::task::PackerTaskOP task
			( pack::task::TaskFactory::create_packer_task( pose ));

		task->set_bump_check( true );

		task->initialize_from_command_line();
		for ( Size ii = 1; ii <= nres; ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
				task->nonconst_residue_task( ii ).restrict_to_repacking();
				//     task->nonconst_residue_task( ii ).or_ex1aro_sample_level ( pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
				//     task->nonconst_residue_task( ii ).or_ex2aro_sample_level ( pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
				assert( task->pack_residue(ii) );
			} else {
				if ( repack_DNA ) {
					task->nonconst_residue_task( ii ).restrict_to_repacking();
				} else {
					task->nonconst_residue_task( ii ).prevent_repacking();
				}
				// assert( !task->pack_residue(ii) );
			}
			assert( !task->design_residue(ii) );
		}
		// dont include current
		task->or_include_current( false );

		Size const nloop( 50 );
		utility::vector1< std::pair< Real, std::string > > results;
		pack::pack_rotamers_loop( pose, scorefxn, task, nloop, results );
		//   pack::pack_rotamers( pose, scorefxn, task);
		Energy pack_score = scorefxn( pose );

		if ( false ) { //debugging
			io::pdb::dump_pdb( pose, combo_tag+"_packed.pdb" );
			pose.energies().show(tt);
			pose.energies().clear();
			Energy pack_rescore( scorefxn( pose ) );
			tt << "pack_rescore: " << pack_score << ' ' << pack_rescore << std::endl;
		}


		{ // test for mismatches -- unnecessary for packing but debugs pose setup
			WatsonCrickResidueMatcher m;
			for ( Size i=1; i<= nres; ++i ) {
				if ( partner[i]>i ) {
					assert( m( pose.residue(i), pose.residue(partner[i])));
				}
			}
		}

		tt << "packed score: " << output_tag << ' ' << combo_tag << ' ' << nat_seq << ' ' << pack_score << std::endl;

		//pose.energies().total_energies().show_nonzero( tt );
		std::ostringstream buf;
		pose.energies().total_energies().show_nonzero( buf );
		tt << buf.str();


		/// PB -- I think SVN MERGE mucked this up a bit. Have to take a closer look.
		io::pdb::dump_pdb( pose, output_tag+"_"+combo_tag+"_"+nat_seq+".pdb" );

		if ( postmin ) {
			AtomTreeMinimizer().run( pose, mm, scorefxn, min_options );
			tt << "postmin score: " << output_tag << ' ' << combo_tag << ' ' << nat_seq << ' ' << scorefxn( pose ) <<
				std::endl;
			//pose.energies().total_energies().show_nonzero( tt );
			pose.energies().total_energies().show_nonzero( buf );
			tt << buf.str();
		}

		// Calculate protein-DNA and DNA RMSDs ...
		Real ca_rmsd, interface_ca_rmsd, interface_allatom_rmsd, dna_bb_rmsd;
		calc_protein_DNA_rmsd( pose, start_pose, ca_rmsd, interface_ca_rmsd, interface_allatom_rmsd );
		calc_DNA_bb_rmsd( pose, start_pose, dna_bb_rmsd );

		std::cout << " RMSDs: " << output_tag << ' ' << combo_tag << ' ' << nat_seq
			<< " interface_allatom=" << interface_allatom_rmsd << ", dna_bb=" << dna_bb_rmsd
			<< std::endl;

		io::pdb::dump_pdb( pose, output_tag+"_"+combo_tag+"_"+nat_seq+".pdb" );

	} // loop over sequence combinations
	basic::prof_show();
}

} // namespace dna
} // namespace devel
