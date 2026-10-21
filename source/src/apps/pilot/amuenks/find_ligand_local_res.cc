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

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

//for debugging coordinates
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

static basic::Tracer TR("amuenks.ligand_local_res");

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
        
	ObjexxFCL::FArray3D< float > const & densdata = core::scoring::electron_density::getDensityMap().get_data();

	core::conformation::Residue const &lig( p.residue( nres ) );
	core::Real distance_to_atom = 0.0;
	core::Real total_local_res = 0.0;
	core::Size nearby_voxel_count = 0;
	core::Size zero_voxel_count = 0;
	for ( int z=1; z<=(int)densdata.u3(); z+=1 ) {
                for ( int y=1; y<=(int)densdata.u2(); y+=1 ) {
                        for ( int x=1; x<=(int)densdata.u1(); x+=1 ) {
                                numeric::xyzVector< core::Real > x_idx(x,y,z);
                                numeric::xyzVector< core::Real > x_cart( x_idx[0], x_idx[1], x_idx[2] );
                                core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
				//core::Real dens_value = -densdata(x,y,z);
				

                                for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
					distance_to_atom = sqrt( pow( x_cart[0] - lig.xyz(iatm)[0], 2 ) + pow( x_cart[1] - lig.xyz(iatm)[1], 2 ) + pow( x_cart[2] - lig.xyz(iatm)[2], 2 ) );

                                	if ( distance_to_atom <= 5.0 ) {//10.0 ) { //default 10.0
						core::Real density = core::scoring::electron_density::getDensityMap().get(x,y,z);
						TR << "HETATM" << std::setw(5) << "1" << " " <<  std::setw(4) << std::left << "O1" << " " <<  std::setw(3) << std::right << "ALA" << " A" << std::setw(4) << "1" << "    " << std::setw(8) << x_cart[0] << std::setw(8) << x_cart[1] << std::setw(8) << x_cart[2] << std::setw(6) << "1.00" << std::setw(6) << std::fixed << std::setprecision(2) << density << std::endl;
						core::Real local_res = density;
						if ( local_res > 0.0 ) {
							total_local_res += local_res;
							nearby_voxel_count++;
						}
						else {
							zero_voxel_count++;
						}
						break;
                                	}
				}
                        }
                }
        }

	TR << "Average local resolution around ligand: " << total_local_res/nearby_voxel_count << std::endl;
	TR << "number of voxels: " << nearby_voxel_count << std::endl;
	TR << "Percentage of voxels with zero value: " << float(zero_voxel_count/(nearby_voxel_count+zero_voxel_count)) << std::endl;
	/*core::scoring::electron_density::ElectronDensity localres;
        localres.readMRCandResize( basic::options::option[map_localres](),
                                   basic::options::option[edensity::mapreso](),
                                   basic::options::option[edensity::grid_spacing]() );
        */
        // smoothing window, adjustable
        /*density.setScoreWindowContext( true );
        density.setWindow( basic::options::option[edens_window]() );

        TR << "PERRESCC, rosettaNum, resName, perResCC, perResLocReso" << std::endl;
        for (core::Size resi=1; resi <= nres; resi++ ){
            std::string resn = p.residue(resi).name();
            // can just call the function for zscores that ray used - calculate_nbr_dens - probably ok
            // need to score the pose before doing this
            core::Real per_res_dens = density.matchRes( resi, p.residue(resi), p, nullptr, false );
            core::Real per_res_locres = get_localreso( p.residue(resi), localres );

            TR << "PERRESCC," << resi << "," << resn << "," << per_res_dens << "," << per_res_locres << std::endl;
        }*/

    } catch ( utility::excn::Exception const & e ){
        e.display();
        return -1;
    }
    return 0;
}
