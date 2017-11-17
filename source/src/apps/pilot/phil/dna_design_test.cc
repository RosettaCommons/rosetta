// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers

#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>


#include <protocols/viewer/viewers.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/ScoreFunction.hh>


#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <utility/excn/Exceptions.hh>


#include <core/pose/Pose.hh>

#include <basic/options/util.hh>

#include <basic/prof.hh> // profiling
#include <basic/Tracer.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>


#include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/format.hh>


// // C++ headers
#include <fstream>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


// //silly using/typedef


using basic::Error;
using basic::Warning;


using namespace core;
//using namespace protocols;

using utility::vector1;
using std::string;
using std::cout;
using std::endl;

using namespace ObjexxFCL::format;

static basic::Tracer tt( "demo.phil.dna_design_test", basic::t_info );

///////////////////////////////////////////////////////////////////////////////////////////
// rescale phosphate charges by desired amount
void
rescale_phosphate_charges(
	Real const scale_factor
)
{
	using namespace chemical;

	vector1< string > atoms;
	atoms.push_back( "P" );
	atoms.push_back( "OP2" );
	atoms.push_back( "OP1" );
	atoms.push_back( "O5'" );
	atoms.push_back( "O3'" );

	ResidueTypeSet & rsd_set( ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ) );

	for ( int i= first_DNA_aa; i<= last_DNA_aa; ++i ) {
		ResidueTypeCOPs const rsd_types( rsd_set.aa_map( AA(i) ) );
		for ( Size k=1; k<= rsd_types.size(); ++k ) {
			ResidueType & rsd_type( rsd_set.nonconst_name_map( rsd_types[k]->name() ) );
			assert( rsd_type.is_DNA() );
			for ( Size j=1; j<= atoms.size(); ++j ) {
				Real const old_charge( rsd_type.atomic_charge( rsd_type.atom_index( atoms[j] ) ) );
				Real const new_charge( old_charge * scale_factor );
				rsd_type.set_atomic_charge( atoms[j], new_charge );
				tt << "rescaling phosphate charge on atom " << rsd_type.name() << ' ' << atoms[j] << " from " <<
					old_charge << " to " << new_charge << '\n';
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////
vector1< Real >
get_boltzmann_probabilities( vector1< Real > const & E, Real const KT )
{
	vector1< Real > P;
	if ( !E.empty() ) {
		Real const mn( utility::min( E ) );
		Real total( 0.0 );
		Size const n( E.size() );
		for ( Size i=1; i<= n; ++i ) {
			Real const prob( std::exp( ( mn - E[i] ) / KT ) );
			P.push_back( prob );
			total += prob;
		}
		for ( Size i=1; i<= n; ++i ) {
			P[i] /= total;
		}
	}

	return P;
}

///////////////////////////////////////////////////////////////////////////////////////////

void
design_test()
{
	using namespace basic::options;
	namespace OK = OptionKeys;
	// using namespace basic::options::OptionKeys;
	// using namespace basic::options::OptionKeys::dna::specificity;
	using namespace scoring;
	using namespace conformation;
	using namespace chemical;
	using namespace kinematics;
	using namespace optimization;
	//using namespace scoring::dna;
	using namespace pose;

	bool const dry_run( option[ OK::run::dry_run ] );
	bool const only_repack( option[ OK::dna::specificity::only_repack ] );
	bool const  dna_design( option[ OK::dna::specificity::design_DNA  ] );
	bool const verbose_chi_stats( false );

	Size const nloop( 25 );
	Real const boltzmann_KT( 0.4 );
	Size const mincore( 3 ); // at least 3 "core" basepairs per structure
	Size const maxcore( 6 ); // at  most 6 "core" basepairs per structure
	Size const ncontacts_core_min1( 4 ); // stop if we fall below this value
	Size const ncontacts_core_min2( 8 ); // stop if we fall below this value and n > mincore

	// setup the scoring function
	assert( option[ OK::dna::specificity::score_function ].user() );
	std::string const scorefxn_file( option[ OK::dna::specificity::score_function ] );
	ScoreFunctionOP scorefxn( new ScoreFunction() );
	scorefxn->initialize_from_file( scorefxn_file );

	// HACK -- fiddle with phosphate charge parameters
	{
		std::ifstream data( scorefxn_file.c_str() );
		std::string line,tag;
		bool found( false );
		Real rescale( 1.0 );
		while ( getline( data,line ) ) {
			std::istringstream l( line );
			l >> tag >> rescale;
			if ( !l.fail() && tag == "#PCR" ) { // Phosphate Charge Rescale
				found = true;
				break;
			}
		}
		data.close();
		if ( found ) {
			// rescale phosphate charges by desired amount
			tt << "rescaling phosphate charges by " << rescale << std::endl;
			rescale_phosphate_charges( rescale );
		}
	}

	vector1< string > files( start_files() );

	// permute files so that all 8 processes won't be on the large PDBs at the same time
	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size n=1; n<= files.size(); ++n ) {
		std::string const & filename( files[n] );

		tt << "read file: " << filename << std::endl;

		// read structure
		Pose pose;
		core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
		Size const nres( pose.size() );

		tt << "read file: " << filename << ' ' << nres << std::endl;

		Pose const start_pose( pose );

		// for profiling
		basic::prof_reset();

		{ // score the starting structure
			(*scorefxn)(pose);
			//std::cout << "start score: " << (*scorefxn)( pose ) << std::endl;
			scorefxn->show( tt, pose );
			basic::prof_show();
		}

		protocols::dna::RestrictDesignToProteinDNAInterfaceOP dna_inter_restrictor =
			new protocols::dna::RestrictDesignToProteinDNAInterface;

		pack::task::TaskFactory taskfactory;
		taskfactory.push_back( dna_inter_restrictor );
		pack::task::PackerTaskOP task = taskfactory.create_task_and_apply_taskoperations( pose );

		if ( only_repack || dna_design ) {
			// convert designed protein positions to repacked
			for ( Size i=1; i<= nres; ++i ) {
				if ( task->design_residue( i ) ) {
					task->nonconst_residue_task(i).restrict_to_repacking();
				}
			}
		}

		// a mask indicating the "core" dna positions
		// only filled if dna_design is true
		vector1< bool > core_dna( nres, false );
		vector1< int > contacts( nres, 0 );
		scoring::dna::set_base_partner( pose );

		if ( dna_design ) {

			scoring::dna::BasePartner const & partner( scoring::dna::retrieve_base_partner_from_pose( pose ) );
			// "design" dna positions that contact the moving aa's
			Real const dis2_threshold( 4.5 * 4.5 );
			vector1< bool > interface_dna( nres, false );
			for ( Size i=1; i<= nres; ++i ) {
				if ( task->pack_residue( i ) ) {
					Residue const & rsd1( pose.residue(i) );
					assert( rsd1.is_protein() );
					// look for positions that this guy contacts, if any
					for ( Size j=1; j<= nres; ++j ) {
						Residue const & rsd2( pose.residue(j) );
						if ( !rsd2.is_DNA() || partner[j] <1 ) continue; // skip unpaired DNA or non-DNA for rsd2

						Real const nbr_dis( rsd1.nbr_atom_xyz().distance( rsd2.nbr_atom_xyz() ) );
						if ( nbr_dis < rsd1.nbr_radius() + rsd2.nbr_radius() + 5.5 ) {
							++contacts[i]; // count DNA nbrs for protein residues
						}
						assert ( partner[j] );
						for ( Size ii=1; ii<= rsd1.nheavyatoms(); ++ii ) {
							for ( Size jj=rsd2.first_sidechain_atom(); jj<= rsd2.nheavyatoms(); ++jj ) {
								if ( rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) ) < dis2_threshold ) {
									interface_dna[j] = true;
									interface_dna[ partner[j] ] = true;
									++contacts[j]; // count heavyatom contacts for DNA residues
								}
							}
						}
					} // j : DNA residues potentially contacted by rsd i
				} // rsd i is being packed?
			} // i : protein residues

			// now allow these positions to be "designed"
			for ( Size i=1; i<= nres; ++i ) {
				if ( interface_dna[i] ) {
					assert( !task->pack_residue( i ) );
					task->temporarily_set_pack_residue( i, true );
					assert(  task->pack_residue( i ) );
					assert( !task->design_residue( i ) );
					pack::task::ResidueLevelTask & rtask( task->nonconst_residue_task(i) );
					rtask.allow_aa( na_ade );
					rtask.allow_aa( na_cyt );
					rtask.allow_aa( na_gua );
					rtask.allow_aa( na_thy );
					assert( task->design_residue( i ) );
				}
			}


			{ // setup residue couplings
				using namespace pack::rotamer_set;
				RotamerCouplingsOP couplings( new RotamerCouplings() );
				couplings->resize( nres );
				for ( Size i=1; i<= nres; ++i ) {
					if ( partner[i] ) {
						(*couplings)[i].first = partner[i];
						(*couplings)[i].second = new conformation::WatsonCrickResidueMatcher();
					}
				}
				task->rotamer_couplings( couplings );
			}

			// define a set of core DNA positions using the contact counts
			{
				vector1< std::pair< int, int > > l;
				for ( Size i=1; i<= nres; ++i ) {
					if ( partner[i] > i ) {
						assert( i == partner[ partner[i] ] );
						l.push_back( std::make_pair( contacts[i] + contacts[ partner[i] ], i ) );
					}
				}
				std::sort( l.begin(), l.end() );
				std::reverse( l.begin(), l.end() );
				Size prev_ncontacts( 0 );
				for ( Size i=1; i<= l.size(); ++i ) {
					Size const ncontacts( l[i].first );
					if ( ncontacts < ncontacts_core_min1 ) break;
					if ( i > mincore && ncontacts < ncontacts_core_min2 ) break;
					if ( i > maxcore && ncontacts < prev_ncontacts ) break;

					Size const pos( l[i].second );
					core_dna[          pos   ] = true;
					core_dna[ partner[ pos ] ] = true;
					prev_ncontacts = ncontacts;
				}
			}
		} // if ( dna_design )

		// stats
		{
			int repack(0), design(0);
			for ( Size i=1; i<= nres; ++i ) {
				if ( task->pack_residue(i) ) ++repack;
				if ( task->design_residue(i) ) ++design;
			}
			std::cout << "START: " << filename << " nres: " << nres << " n_repack " << repack << " n_design " << design <<
				'\n';
		}

		if ( dry_run ) continue;

		// call pack rotamers /////////////////////////////////////////////////////////////////////////////////////////
		utility::vector1< std::pair< Real, std::string > > results;
		pack::pack_rotamers_loop( pose, (*scorefxn), task, nloop, results );
		scorefxn->show( tt, pose );

		// write the final statistics for this PDB ////////////////////////////////////////////////////////////////////

		{ // FIRST THE POSITIONS BEING REPACKED:
			for ( Size i=1; i<= nres; ++i ) {
				if ( task->pack_residue( i ) && !task->design_residue( i ) ) {
					using namespace pack::dunbrack;
					using namespace ObjexxFCL::format;
					Residue const &       rsd(       pose.residue(i) );
					Residue const & start_rsd( start_pose.residue(i) );
					assert( rsd.is_protein() );
					Size const nchi( rsd.nchi() );

					// correct by chi dev
					ChiVector final_chi(nchi), start_chi(nchi), chidev(nchi);
					for ( Size j=1; j<= nchi; ++j ) {
						start_chi[j] = start_rsd.chi(j);
						final_chi[j] =       rsd.chi(j);
					}

					for ( Size j=1; j<= nchi; ++j ) {
						chidev[j] = subtract_chi_angles( final_chi[j], start_chi[j], rsd.aa(), j );
					}

					// correct by rot number
					RotVector final_rot(nchi), start_rot(nchi);
					rotamer_from_chi( rsd.type(), final_chi, final_rot );
					rotamer_from_chi( rsd.type(), start_chi, start_rot );

					std::cout << "REPACKED " << I(4,i) << I(4,contacts[i]) << ' ' << start_pose.residue(i).name() << I(4,nchi);
					// rot correct
					std::cout << " rc:";
					bool all_correct( true );
					for ( Size j=1; j<= nchi; ++j ) {
						if ( final_rot[j] != start_rot[j] ) all_correct = false;
						std::cout << ' ' << ( all_correct && ( final_rot[j] == start_rot[j] ) );
					}
					// chi dev
					std::cout << " cd:";
					for ( Size j=1; j<= nchi; ++j ) {
						std::cout << F(6,1,std::abs( chidev[j]) );
					}
					if ( !verbose_chi_stats ) {
						std::cout << '\n';
						continue;
					}

					// start rot
					std::cout << " start_rot:";
					for ( Size j=1; j<= nchi; ++j ) {
						std::cout << ' ' << start_rot[j];
					}
					// start chi
					std::cout << " start_chi:";
					for ( Size j=1; j<= nchi; ++j ) {
						std::cout << ' ' << start_chi[j];
					}
					// final rot
					std::cout << " final_rot:";
					for ( Size j=1; j<= nchi; ++j ) {
						std::cout << ' ' << final_rot[j];
					}
					// final chi
					std::cout << " final_chi:";
					for ( Size j=1; j<= nchi; ++j ) {
						std::cout << ' ' << final_chi[j];
					}
					std::cout << std::endl;
				}
			}
		} // scope

		{ // NOW POSITIONS REDESIGNED
			scoring::dna::BasePartner const & partner( scoring::dna::retrieve_base_partner_from_pose( pose ) );
			assert( results.size() == nloop && results[1].second.size() == nres );

			vector1< Real > energies;
			for ( Size i=1; i<= nloop; ++i ) {
				energies.push_back( results[i].first );
			}
			vector1< Real > const boltz_prob( get_boltzmann_probabilities( energies, boltzmann_KT ) );

			for ( Size i=1; i<= nres; ++i ) {
				if ( task->design_residue( i ) ) {
					Size const basepair_contacts( ( partner[i] > 0 ) ? contacts[i] + contacts[ partner[i] ] : contacts[i] );
					std::cout << "DESIGNED " << I(4,i) << I(4,partner[i]) << I(4,basepair_contacts) << " from " <<
						start_pose.residue(i).aa() << " to " << pose.residue(i).aa() << '\n';

					// get boltz-weighted frequencies for each aa that occurred at this position
					std::map< char, Real > frequency;
					for ( Size j=1; j<= nloop; ++j ) {
						char const aa( results[j].second[ i-1 ] );
						frequency[ aa ] += boltz_prob[j];
					}

					// show results
					char const start_aa( start_pose.sequence()[i-1] );
					std::cout << "BOLTZ_PROB " << I(4,i) << I(4,partner[i]) << I(4,basepair_contacts) << " from " <<
						start_aa << " to";
					for ( std::map< char, Real >::const_iterator it=frequency.begin(), ite=frequency.end(); it!=ite; ++it ) {
						std::cout << ' ' << it->first << F(5,2,it->second);
					}
					std::cout << '\n';
					if ( core_dna[i] ) {
						std::cout << "BP_CORE " << I(4,i) << I(4,partner[i]) << I(4,basepair_contacts) << " from " <<
							start_aa << " to";
						for ( std::map< char, Real >::const_iterator it=frequency.begin(), ite=frequency.end(); it!=ite; ++it ) {
							std::cout << ' ' << it->first << F(5,2,it->second);
						}
						std::cout << '\n';
					}
				}
			}
		}

		basic::prof_show();

		if ( option[ OK::out::file::o ].user() ) {
			pose.dump_pdb( option[ OK::out::file::o ] );
		}
	}
}


///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	design_test();

	exit(0);

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
		// initialize option and random number system
		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
