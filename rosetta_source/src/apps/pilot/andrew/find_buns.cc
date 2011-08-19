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

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

// C++ headers
#include <fstream>

int main( int argc, char * argv [] )
{
	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::io::pdb;
	using namespace core::graph;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	using namespace devel::vardist_solaccess;

	using namespace basic::options;

	devel::init( argc, argv );

	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

	if ( input_jobs.size() != 1 ) {
		utility_exit_with_message( "Expected exactly one pdb to be specified from the -s or -l flags" );
	}

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, input_jobs[ 1 ]->input_tag() );
	ScoreFunctionOP sfxn = getScoreFunction();
	(*sfxn)(pose);

	HBondSet hbset;
	fill_hbond_set( pose, false, hbset );

	utility::vector1< VarSolDRotamerDotsOP > rotamer_dots( pose.total_residue() );
	Vector com( 0.0 );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		rotamer_dots[ ii ] = new VarSolDRotamerDots( new core::conformation::Residue( pose.residue( ii )), true );
		rotamer_dots[ ii ]->increment_self_overlap();
		com += pose.residue( ii ).xyz( pose.residue( ii ).nbr_atom() );
	}
	com /= pose.total_residue();
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ii)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(ii)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			rotamer_dots[ ii ]->intersect_residues( *rotamer_dots[ (*iru)->get_second_node_ind() ] );
		}
	}

	std::list< AtomID > buns;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Size jj = 1; jj <= pose.residue(ii).nheavyatoms(); ++jj) {
			AtomID atid( jj, ii );

			if ( pose.residue(ii).atom_type( jj ).is_acceptor() || pose.residue(ii).atom_type( jj ).is_donor() ) {

				if ( pose.residue(ii).atom_type( jj ).is_donor() &&
						pose.residue(ii).type().attached_H_begin( jj ) > pose.residue(ii).type().attached_H_end( jj )) {
						// i.e. proline backbone N; not really a donor
						continue;
				}

				if ( rotamer_dots[ ii ]->any_exposed_dots( jj ) ) continue;
				bool exposed_hydrogens = false;
				for( Size kk = pose.residue(ii).type().attached_H_begin( jj ); kk<= pose.residue(ii).type().attached_H_end( jj ); kk++){
					if ( rotamer_dots[ ii ]->any_exposed_dots( kk ) ) { exposed_hydrogens = true; break; }
				}
				if ( exposed_hydrogens ) continue;

				utility::vector1< HBondCOP > const & jjhbonds( hbset.atom_hbonds( atid ));
				bool attached_h_form_hbonds = false;
				if ( jjhbonds.empty() ) {
					for( Size kk = pose.residue(ii).type().attached_H_begin( jj ); kk<= pose.residue(ii).type().attached_H_end( jj ); ++kk ){
						utility::vector1< HBondCOP > const & kkhbonds( hbset.atom_hbonds( AtomID( kk, ii ) ));
						if ( ! kkhbonds.empty() ) {
							attached_h_form_hbonds = true;
							break;
						}
					}
					if ( attached_h_form_hbonds ) continue;
					buns.push_back( atid );
				}
			}
		}
	}

	ConformationKinWriter writer;
	std::ofstream fout( "test.kin" );
	write_kinemage_header( fout, 1, "test", com );
	writer.write_coords( fout, pose.conformation() );
	fout << "@group { Dots } dominant off\n";
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) rotamer_dots[ ii ]->write_dot_kinemage( fout );
	fout << "@group { BUns } dominant\n";
	fout << "@balllist radius= 0.6 color= white\n";
	for ( std::list< AtomID >::const_iterator iter = buns.begin(); iter != buns.end(); ++iter ) {
		fout << "{" << iter->rsd() << " " << pose.residue( iter->rsd() ).atom_name( iter->atomno() ) << "} P " << pose.xyz( *iter ).x() << " " << pose.xyz( *iter ).y() << " " << pose.xyz( *iter ).z() << "\n";
	}
}

