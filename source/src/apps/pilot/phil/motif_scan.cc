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


// libRosetta headers
#include <devel/dna/protocols.hh>
#include <devel/dna/util.hh>
//#include <devel/dna/util.hh>
// AUTO-REMOVED #include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/frags/TorsionFragment.hh>

#include <protocols/viewer/viewers.hh>
#include <utility/excn/Exceptions.hh>
// AUTO-REMOVED #include <protocols/simple_moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/TrialMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/moves/rigid_body_moves.hh>

#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/dna/setup.hh>
// AUTO-REMOVED #include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
// AUTO-REMOVED #include <core/scoring/dna/DNA_BasePotential.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds_geom.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondSet.hh>
// AUTO-REMOVED #include <core/scoring/elec/FA_ElecEnergy.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairAll.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairFunction.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

#include <core/types.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AA.hh>

// AUTO-REMOVED #include <core/conformation/util.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueMatcher.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerCouplings.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/WaterPackingInfo.hh>

// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/kinematics/visualize.hh>
#include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>

#include <basic/options/util.hh>
//#include <basic/options/after_opts.hh>

// AUTO-REMOVED #include <basic/prof.hh> // profiling
// AUTO-REMOVED #include <basic/basic.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sstream>

//silly using/typedef


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/dna.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <protocols/moves/Mover.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end







////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;

static thread_local basic::Tracer tw( "demo.phil.motif_scan", basic::t_warning );
//static basic::Tracer td( "demo.phil.motif_scan", basic::t_debug );
static thread_local basic::Tracer tt( "demo.phil.motif_scan", basic::t_trace );


///////////////////////////////////////////////////////////////////////////////
typedef std::map< char, Real > MotifColumn;

typedef std::map< Size, MotifColumn > Motif;


///////////////////////////////////////////////////////////////////////////////
void
smooth_motif_column( MotifColumn & C )
{
	Real const w( 0.05 );

	for ( MotifColumn::iterator it=C.begin(), ite=C.end(); it != ite; ++it ) {
		it->second = ( it->second + w ) / ( 1 + C.size() * w );
	}

}

///////////////////////////////////////////////////////////////////////////////
Real
motif_column_deviation( MotifColumn P, MotifColumn E )
{
	smooth_motif_column( P );
	smooth_motif_column( E );
	Real dev(0.0);
	for ( MotifColumn::const_iterator it=P.begin(), ite=P.end(); it != ite; ++it ) {
		assert( E.count( it->first ) );
		Real const p( it->second );
		Real const e( E.find( it->first )->second );
		dev += e * std::log( e / p );
	}
	return dev;
}

///////////////////////////////////////////////////////////////////////////////
Real
motif_deviation( Motif const & P, Motif const & E )
{
	Size const N( P.size() );
	assert( N == E.size() );
	Real dev(0.0);

	for ( Motif::const_iterator it = P.begin(), ite= P.end(); it != ite; ++it ) {
		assert( E.count( it->first ) );
		dev += motif_column_deviation( it->second, E.find( it->first )->second );
	}
	return dev / N;
}

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// helper, should go into library

vector1< Size >
parse_pdb_pos( pose::Pose const & pose )
{
	core::pose::PDBPoseMap pose_map(pose.pdb_info()->pdb2pose());

	vector1< string > pdb_pos_list( option[ OK::dna::specificity::pdb_pos ]() );
	vector1< Size >  pose_pos_list;

	for ( Size i=1; i<= pdb_pos_list.size(); ++i ) {
		std::string const & resid( pdb_pos_list[i] );
		int pos;
		char chain;
		std::size_t fpos( resid.find(':') );
		if ( fpos != string::npos ) {
			pos = int_of( resid.substr(0,fpos) );
			if ( fpos == resid.size()-1 ) {
				chain = ' ';
			} else {
				chain = resid[ fpos+1 ];
			}
		} else {
			pos = int_of( resid );
			chain = ' ';
		}
		if ( chain == '_' ) chain = ' ';
		pose_pos_list.push_back( pose_map.find( chain, pos ) );
	}
	return pose_pos_list;
}

///////////////////////////////////////////////////////////////////////////////
// sets a default foldtree
//

void
make_dna_only_pose( pose::Pose const & pose, Pose & dna_pose )//, vector1< Size > & old2new )
{
	dna_pose.clear();
	//old2new.resize( pose.total_residue() );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.is_DNA() ) {
			if ( rsd.is_lower_terminus() || !rsd.is_polymer() ) {
				dna_pose.append_residue_by_jump( rsd, 1 );
			} else {
				dna_pose.append_residue_by_bond( rsd );
			}
// 			old2new[ i ] = dna_pose.total_residue();
// 		} else {
// 			old2new[ i ] = 0;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
class MotifScanMover : public moves::Mover {
public:

// 	MotifScanMover( Size const seqpos, string const & mode, ScoreFunctionCOP scorefxn ):
// 		seqpos_( seqpos ),
// 		mode_( mode ),
// 		scorefxn_( scorefxn )
// 	{
// 		type( "mscan_"+mode );
// 	}

	MotifScanMover( string const & mode, ScoreFunctionCOP scorefxn ):
		seqpos_( 0 ),
		mode_( mode ),
		scorefxn_( scorefxn )
	{
		type( "mscan_"+mode );
	}

	void
	seqpos( Size const setting )
	{
		seqpos_ = setting;
	}

	void
	apply( Pose & pose )
	{
		using devel::dna::repack_base_pair_neighbors;
		if ( mode_ == "only_score" || mode_ == "0" ) {
			// pass
		} else if ( mode_ == "protein_repack" || mode_ == "1" ) {
			// repack protein
			repack_base_pair_neighbors( pose, *scorefxn_, seqpos_, false /*include_current*/, false /* repack_dna */ );
		} else if ( mode_ == "protein_dna_repack" || mode_ == "2" ) {
			// repack protein and dna
			repack_base_pair_neighbors( pose, *scorefxn_, seqpos_, false /*include_current*/, true /* repack_dna */ );
		} else {
			utility_exit_with_message( "unrecognized mode: "+mode_ );
		}
	}

private:
	Size seqpos_;
	string mode_;
	ScoreFunctionCOP scorefxn_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

void
single_position_motif_scan(
													 Pose const & start_pose,
													 ScoreFunction const & scorefxn,
													 Size const seqpos,
													 protocols::moves::Mover & complex_mover,
													 protocols::moves::Mover & dna_mover,
													 std::string const & tag,
													 Real const KT,
													 bool const dump_pdbs,
													 MotifColumn & prob
													 )
{
	prob.clear();

	vector1< Real > base_energies;
	vector1< char > bases;
	Real best_energy(0.0);

	for ( Size ii=1; ii<= 4; ++ii ) {
		AA const na = AA( first_DNA_aa + ii - 1 );
		tt << "mutating position " << seqpos << " to " << na << std::endl;

		// make the mutation
		Pose pose( start_pose );
		scoring::dna::set_base_partner( pose );
		devel::dna::make_base_pair_mutation( pose, seqpos, na );


		tt << "score-after-mutation: " << scorefxn( pose ) << std::endl;

		// optimize complex as requested
		complex_mover.apply( pose );

		Real const complex_score( scorefxn( pose ) );

		// output complex energies
		std::cout << "complex_energies: " << F(9,3,complex_score ) << ' ' <<
			pose.energies().total_energies().weighted_string_of( scorefxn.weights() ) << std::endl;

		if ( dump_pdbs ) {
			pose.dump_scored_pdb( tag+"motif_scan_"+lead_zero_string_of( seqpos, 4 )+pose.residue(seqpos).name1()+".pdb",
														scorefxn );
		}
		// make dna-only pose, optimize as requested
		Pose dna_pose;
		make_dna_only_pose( pose, dna_pose );
		scoring::dna::set_base_partner( dna_pose );

		// optimize dna pose as requested
		dna_mover.apply( dna_pose );

		Real const dna_score( scorefxn( dna_pose ) );
		// output dna-only energies
		std::cout << "dna_energies: " << F(9,3,dna_score ) << ' ' <<
			dna_pose.energies().total_energies().weighted_string_of( scorefxn.weights() ) << std::endl;

		Real const binding_energy( complex_score - dna_score );
		base_energies.push_back( binding_energy );
		bases.push_back( pose.residue( seqpos ).name1() );
		if ( ii == 1 || best_energy > binding_energy ) best_energy = binding_energy;
	}


	// compute boltzmann probs
	Real partition(0.0);
	for ( Size i=1; i<= 4; ++i ) partition += std::exp( ( best_energy - base_energies[i] ) / KT );

	std::cout << "BOLTZ_PROB: " << start_pose.residue(seqpos).name1() << I(4,seqpos) <<
		I( 6, start_pose.pdb_info()->number( seqpos ) ) << ' ' << start_pose.pdb_info()->chain( seqpos );
	for ( Size i=1; i<= 4; ++i ) {
		Real const bp( std::exp( ( best_energy - base_energies[ i ] ) / KT ) / partition );
		std::cout << ' ' << bases[i] << ' ' << F(9,3,base_energies[ i ] - best_energy ) << F(9,6,bp );
		prob[ bases[i] ] = bp;
	}
	std::cout << ' ' << complex_mover.type() << ' ' << dna_mover.type() << ' ' << tag << std::endl;

}






///////////////////////////////////////////////////////////////////////////////
MotifColumn
contact_model_prediction( Pose const & pose, Size const seqpos )
{
	using namespace scoring::dna;

	// from Morozov et al
	Size const Nmax( 20 );
	Real const dis2_cutoff( 4.5 * 4.5 );

	// count contacts between protein heavyatoms and this basepair
	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );
	Size N( 0 );
	for ( int r=1; r<= 2; ++r ) {
		Size const i( r == 1 ? seqpos : partner[ seqpos ] );
		if ( !i ) continue;
		Residue const & rsd1( pose.residue(i) );
		assert( rsd1.is_DNA() );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			Residue const & rsd2( pose.residue(j) );
			if ( rsd2.is_protein( )) {
				for ( Size ii=rsd1.first_sidechain_atom(); ii<= rsd1.nheavyatoms(); ++ii ) {
					for ( Size jj=1; jj<= rsd2.nheavyatoms(); ++jj ) {
						if ( rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) ) < dis2_cutoff ) ++N;
					}
				}
			}
		}
	}

	Real p0,p1;
	if ( N >= Nmax ) {
		p1 = 1.0;
		p0 = 0.0;
	} else {
		p0 = ( 1- Real(N)/Nmax ) / 4;
		p1 = p0 + Real(N)/Nmax;
	}

	assert( std::abs( p1 + 3*p0 - 1.0 ) < 1e-3 );
	MotifColumn C;
	C['a'] = p0;
	C['c'] = p0;
	C['g'] = p0;
	C['t'] = p0;
	C[ pose.residue( seqpos ).name1() ] = p1;
	return C;
}



///////////////////////////////////////////////////////////////////////////////
void
motif_scan()
{
	bool const pre_minimize( option[ OK::dna::specificity::pre_minimize ] );
	bool const pre_pack    ( option[ OK::dna::specificity::pre_pack ] );
	Real const KT( 0.8 );
	bool const dump_pdbs( option[ OK::dna::specificity::dump_pdbs ] );
	std::string const mode( option[ OK::dna::specificity::mode ] );

	// setup the scoring function
	ScoreFunctionOP scorefxn = new ScoreFunction();
	scorefxn->initialize_from_file( option[ OK::dna::specificity::score_function ] );

	// setup the movers
	MotifScanMover complex_mover( mode, scorefxn ), dna_mover( mode, scorefxn );

	// loop over the files
	vector1< string > const files( start_files() );

	for ( Size m=1; m<= files.size(); ++m ) {
		string const & filename( files[m] );

		Pose pose;
		core::import_pose::pose_from_pdb( pose, filename );
		Size const nres( pose.total_residue() );
		core::pose::PDBPoseMap pose_map(pose.pdb_info()->pdb2pose());


		scoring::dna::set_base_partner( pose );



		//// optionally pre-optimize the complex
		if ( pre_minimize ) {
			// minimize the protein sidechains
			kinematics::MoveMap mm;
			mm.set_chi( false );
			for ( Size i=1; i<= nres; ++i ) {
				if ( pose.residue(i).is_protein() ) mm.set_chi(i, true );
			}
			Real const score1( (*scorefxn)(pose) );
			optimization::AtomTreeMinimizer().run( pose, mm, *scorefxn, optimization::MinimizerOptions("dfpmin",0.001,true ) );
			Real const score2( (*scorefxn)(pose) );
			tt << "pre-min scores: " << F(9,3,score1) << F(9,3,score2) << std::endl;
		}

		if ( pre_pack ) {
			// repack the protein sidechains, with include_current = true
			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
			task->initialize_from_command_line();
			task->or_include_current( true );

			for ( Size ii = 1; ii <= nres; ++ii ) {
				if ( pose.residue(ii).is_protein() ) {
					task->nonconst_residue_task( ii ).restrict_to_repacking();
				} else {
					task->nonconst_residue_task( ii ).prevent_repacking();
				}
			}

			Real const score1( (*scorefxn)(pose) );
			pack::rotamer_trials( pose, *scorefxn, task );
			Real const score2( (*scorefxn)(pose) );
			pack::pack_rotamers_loop( pose, *scorefxn, task, 15 );
			Real const score3( (*scorefxn)(pose) );
			tt << "pre-pack scores: " << F(9,3,score1) << F(9,3,score2) << F(9,3,score3) << std::endl;
		}


		// sequence mapping from pose to dna-only pose
		vector1< Size > dna_seqpos_mapping( pose.total_residue(), 0 );
		{
			for ( Size i=1, n=0; i<= pose.total_residue(); ++i ) {
				if ( pose.residue(i).is_DNA() ) dna_seqpos_mapping[ i ] = ++n;
			}
		}


		// read the motif positions
		vector1< Motif > matches;
		if ( option[ OK::dna::specificity::pdb_pos ].user() ) {
			vector1< Size > pos_list( parse_pdb_pos( pose ) );
			Motif M;
			MotifColumn flat;
			flat['a'] = flat['c'] = flat['g'] = flat['t'] = 0.25;
			for ( Size i=1; i<= pos_list.size(); ++i ) {
				M[ pos_list[i] ] = flat;
			}
			matches.push_back( M );

		} else {
			// read from motif pdb file
			std::ifstream data( filename.c_str() );
			Motif motif; // motif currently  being parsed
			std::string line, tag1, tag2;
			while ( getline( data, line ) ) {
				std::istringstream l( line );
				l >> tag1 >> tag2;
				if ( !l.fail() && tag1 == "REMARK" ) {
					if ( tag2 == "MOTIF_CHAIN" ) {
						// start a new motif match:
						if ( !motif.empty() ) matches.push_back( motif );
						motif.clear();
					} else if ( tag2 == "MOTIF" ) {
						Size pos;
						string chain;
						l >> tag1 >> tag1 >> tag1 >> tag1 >> pos >> chain >>
							tag2 >> tag2 >> tag2 >> tag2 >> tag2 >> tag2 >> tag2 >> tag2 >> tag2;
						if ( !l.fail() && tag1 == "pdb_pos=" && tag2 == "motif=" && chain.size() == 1 ) {
							MotifColumn C;
							for ( Size i=1; i<= 4; ++i ) {
								Real freq;
								l >> tag1 >> freq;
								if ( l.fail() || tag1.size()!= 2 || tag1[1] != '=' ) break;
								C[ tag1[0] ] = freq;
							}
							if ( l.fail() || C.size() != 4 ) {
								tw << "bad line! " << line << std::endl;
								utility_exit();
							} else {
								// success
								motif[ pose_map.find( chain[0], pos ) ] = C;
							}
						} else {
							tw << "bad line! " << line << std::endl;
							utility_exit();
						}
					}
				}
			}
			data.close();
			if ( !motif.empty() ) matches.push_back( motif );
		}


		for ( Size m=1; m<= matches.size(); ++m ) {
			Motif const & motif( matches[m] );
			Motif prediction, contact_model_motif;
			for ( Motif::const_iterator it= motif.begin(), ite= motif.end(); it != ite; ++it ) {
				Size const seqpos( it->first );
				Size const dna_seqpos( dna_seqpos_mapping[ seqpos ] );
				complex_mover.seqpos( seqpos );
				dna_mover.seqpos( dna_seqpos );
				MotifColumn prob;
				single_position_motif_scan( pose, *scorefxn, seqpos, complex_mover, dna_mover, filename, KT, dump_pdbs, prob );
				prediction[ seqpos ] = prob;
				MotifColumn const & experimental( it->second );
				MotifColumn const contact_model( contact_model_prediction( pose, seqpos ) );
				contact_model_motif[ seqpos ] = contact_model;
				std::cout << "COLUMN_DEV " << pose.residue(seqpos).name1() << " P:" <<
					F(9,3,motif_column_deviation( prob, experimental )) << " C:" <<
					F(9,3,motif_column_deviation( contact_model, experimental )) << " motifs(P,C,E): ";
				for ( MotifColumn::const_iterator it=prob.begin(), ite=prob.end(); it != ite; ++it ) {
					char const na( it->first );
					std::cout << ' ' << na <<
						F(5,2, it->second ) <<
						F(5,2, contact_model.find( na )->second ) <<
						F(5,2,  experimental.find( na )->second );
				}
				std::cout << ' ' << m << ' ' << seqpos << ' ' << filename << std::endl;
			}

			// now compare predicted and experimental:
			std::cout << "MOTIF_DEV " << F(9,3,motif_deviation( prediction, motif )) <<
				F(9,3,motif_deviation( contact_model_motif, motif )) <<
				' ' << m << ' ' << filename << std::endl;


		} // motif matches for this pdb

	} // pdbfiles




}

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	motif_scan();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
