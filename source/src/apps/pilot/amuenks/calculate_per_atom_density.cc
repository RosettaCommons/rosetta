#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/Atom.hh>
#include <core/chemical/AtomType.hh> 
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/DensitySymmInfo.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>
#include <cmath>
#include <limits>
#include <numeric>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

//for debugging coordinates
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

static basic::Tracer TR("amuenks.per_atom_density");

using namespace basic::options;
using namespace basic::options::OptionKeys;

int
main( int argc, char * argv[] ){
    try {
        
        devel::init(argc, argv);

        core::pose::Pose p;
        utility::vector1< utility::file::FileName > fname;
        fname = basic::options::option[in::file::s]();
        core::import_pose::pose_from_file(p, fname[1], core::import_pose::PDB_file);
        core::Size nres = p.size();

        core::scoring::electron_density::ElectronDensity density;
        density.readMRCandResize( basic::options::option[edensity::mapfile](),
                                  basic::options::option[edensity::mapreso](),
                                  basic::options::option[edensity::grid_spacing]() );
        
	core::conformation::Residue const &lig( p.residue( nres ) );
				
	utility::vector1<core::Real> dens_values;
	for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
		core::Real atom_dens = core::scoring::electron_density::getDensityMap().matchAtomFast( nres, iatm, lig, p, nullptr );
		dens_values.push_back(atom_dens);
		TR << "Atom density at " << lig.atom_name(iatm) << " " << iatm << ": " << atom_dens  << std::endl; 		
	}

	//std::sort(dens_values.begin(), dens_values.end());
	
	core::Real total_sum = std::accumulate( dens_values.begin(), dens_values.end(), 0.0 );
	core::Real mean = total_sum/dens_values.size();

	core::Real distance_from_mean = 0.0;
	for ( core::Size i = 1; i <= dens_values.size(); ++i ) {
		distance_from_mean += std::pow(dens_values[i] - mean, 2);
	}

	core::Real std_dev = std::sqrt( distance_from_mean / dens_values.size() );
	TR << "Mean: " << mean << std::endl;
	TR << "Standard deviation: " << std_dev << std::endl;

	utility::vector1< core::Size > centrality_scores( lig.nheavyatoms(), 0 );

	core::Size loop_counter = 0;
	while ( std::find( centrality_scores.begin(), centrality_scores.end(), 0 ) != centrality_scores.end() ) {
		if ( lig.nheavyatoms() <= 1 ) { break; }
		
		for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
			if ( centrality_scores[ iatm ] == 0 ) {  
				if ( lig.get_adjacent_heavy_atoms( iatm ).size() == 1 ) {
					centrality_scores[iatm] = 1;
					//TR << "Centrality score at atom " << lig.atom_name(iatm) << " " << iatm << ": " << centrality_scores[iatm] << std::endl;
				}

				//handles rings or another scenario where an atom has 2 unassigned neighbors
				//hacky, will potentially leave 1 ring atom and its non-ring neighbors with too high pf a value
				else if ( loop_counter > lig.nheavyatoms() ) {
					core::Size smallest_assigned_score = lig.nheavyatoms();
					for ( core::Size inbr = 1; inbr <= lig.get_adjacent_heavy_atoms( iatm ).size(); ++inbr ) {
                                                core::Size nbr_atm_index = lig.get_adjacent_heavy_atoms( iatm )[inbr];

						if ( centrality_scores[ nbr_atm_index ] < smallest_assigned_score && centrality_scores[ nbr_atm_index ] != 0) {
                                                        smallest_assigned_score = centrality_scores[ nbr_atm_index ];
						}
					}

					centrality_scores[ iatm ] = smallest_assigned_score + 1;
				}

				else {
					core::Size unassigned_centrality = 0;
					core::Size largest_assigned_score = 0;
					for ( core::Size inbr = 1; inbr <= lig.get_adjacent_heavy_atoms( iatm ).size(); ++inbr ) {
						core::Size nbr_atm_index = lig.get_adjacent_heavy_atoms( iatm )[inbr];
						
						if ( std::find( lig.cut_bond_neighbor( iatm ).begin(), lig.cut_bond_neighbor( iatm ).end(), nbr_atm_index ) != lig.cut_bond_neighbor( iatm ).end() ) {
							TR << "Atom " << iatm << " has a cut bond with " << nbr_atm_index << std::endl;
							continue;
						}

						if ( iatm == 1 ) { TR << nbr_atm_index << std::endl; }
						
						if ( centrality_scores[ nbr_atm_index ] == 0 ) {
							unassigned_centrality++;
						}
						else if ( centrality_scores[ nbr_atm_index ] > largest_assigned_score ) {
							largest_assigned_score = centrality_scores[ nbr_atm_index ];
						}
					}

					if ( unassigned_centrality == 1 ) {
						centrality_scores[ iatm ] = largest_assigned_score + 1;
						//TR << "Centrality score at atom " << lig.atom_name(iatm) << " " << iatm << ": " << centrality_scores[iatm] << std::endl;
					}

					if (unassigned_centrality == 0 ) { //last atom to be assigned
						utility::vector1< core::Size > neighbor_centralities;
						for ( core::Size inbr = 1; inbr <= lig.get_adjacent_heavy_atoms( iatm ).size(); ++inbr ) {
                                                	core::Size nbr_atm_index = lig.get_adjacent_heavy_atoms( iatm )[inbr];
                                        		neighbor_centralities.push_back( centrality_scores[ nbr_atm_index ] );
						}
	       					
						std::sort(neighbor_centralities.rbegin(), neighbor_centralities.rend());
						
						centrality_scores[ iatm ] = neighbor_centralities[2] + 1;
					}
				}
			}
			//TR << "Centrality score at atom " << lig.atom_name(iatm) << " " << iatm << ": " << centrality_scores[iatm] << std::endl;
		}
		
		loop_counter++;
	}

	for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
		TR << "Centrality score at atom " << lig.atom_name(iatm) << " " << iatm << ": " << centrality_scores[iatm] << std::endl;
	}

	core::Real penalty = 0.0;	
	for ( core::Size i = 1; i <= dens_values.size(); ++i ) {
                if ( dens_values[i] < mean - std_dev ) {
			TR << "Atom " << i << " " << dens_values[i] << " will be penalized" << std::endl;
			penalty += ( mean - std_dev - dens_values[i] ) * std::log( /*2 */ centrality_scores[i] + 1 );
		}
        }
		
	TR << "Penalty: " << penalty << std::endl;//std_dev * 10.0 << std::endl;
	//core::Real bottom_sum = 0.0;
	//for ( core::Size i = 1; i <= dens_values.size()*0.75; ++i ) {
	//	TR << i << std::endl;
	//	bottom_sum += dens_values[i];
	//}
	//TR << "bottom sum: " << bottom_sum << std::endl;

	//TR << "All atom per atom density: " << total_sum/lig.nheavyatoms() << std::endl;
	//TR << "Lowest 80 percent of atoms per atom density: " << bottom_sum/(std::floor(dens_values.size()*0.75)) << std::endl;

    } catch ( utility::excn::Exception const & e ){
        e.display();
        return -1;
    }
    return 0;
}
