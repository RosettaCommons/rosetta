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

static basic::Tracer TR("reggiano.val_metrics.density_zscores");

using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY(File, map_localres)
OPT_KEY(Integer, edens_window)
OPT_KEY(Real, local_radius)
OPT_KEY(Boolean, use_bfactors)
OPT_KEY(Boolean, bfac_as_local)
OPT_KEY(Real, max_dist_bfac)

// for debugging
void
dump_coords(utility::vector1 < numeric::xyzVector< core::Real > > all_xyz,
	std::string pdb_name){

	core::pose::Pose temppose;

	for ( core::Size i=1; i<=all_xyz.size(); i++ ) {

		core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::conformation::ResidueOP mgres = core::conformation::ResidueFactory::create_residue( *(rsd_set->get_representative_type_name3(" MG")));
		mgres->set_xyz("MG",all_xyz[i]);
		temppose.append_residue_by_jump( *mgres, 1 );

	}

	temppose.dump_pdb(pdb_name+".pdb");

}


// this has to exist somewhere else
bool
nearlyEqual(float a, float b, float epsilon){
	float absA = std::abs(a);
	float absB = std::abs(b);
	float diff = std::abs(a - b);

	if ( a == b ) {
		return true;
	} else if ( a == 0 || b == 0 || (absA+absB) < std::numeric_limits<float>::lowest() ) {
		return diff < (epsilon * std::numeric_limits<float>::lowest());
	} else {
		return diff / std::min((absA + absB), std::numeric_limits<float>::max()) < epsilon;
	}

}


core::Real
get_localreso( core::conformation::Residue const & rsd,
	core::scoring::electron_density::ElectronDensity & locres){

	// may be inefficient to call this below for each res...
	ObjexxFCL::FArray3D<float> used = locres.get_data();
	float used_val = 9999.0;

	utility::vector1< std::string > atoms( 4 );
	atoms[1] = "C"; atoms[2] = "CA"; atoms[3] = "N"; atoms[4] = "O";

	std::vector < float > resos;

	//get indices for map
	numeric::xyzVector< int > ixyz1(0,0,0);
	numeric::xyzVector< int > ixyz2(0,0,0);

	// get the ix, iy, iz for all the points we're interested in
	numeric::xyzVector< core::Real > tmp_vec(0,0,0);

	for ( core::Size a = 1; a <= atoms.size(); a++ ) {
		core::Size atmn = rsd.atom_index(atoms[a]);
		numeric::xyzVector< core::Real > nbr_xyz = rsd.xyz(atmn);
		core::Real nbr_radius = rsd.atom_type(atmn).lj_radius();

		for ( int i = 0; i <= 2; i++ ) {

			tmp_vec[i] = nbr_radius;
			numeric::xyzVector< core::Real > cart1 = nbr_xyz - tmp_vec;
			numeric::xyzVector< core::Real > cart2 = nbr_xyz + tmp_vec;

			numeric::xyzVector< core::Real > frac1 = locres.get_c2f()*cart1;
			numeric::xyzVector< core::Real > frac2 = locres.get_c2f()*cart2;

			ixyz1[i] = frac1[i]*locres.getGrid()[i]-locres.getOrigin()[i]+1;
			ixyz2[i] = frac2[i]*locres.getGrid()[i]-locres.getOrigin()[i]+1;

			tmp_vec[0] = 0.0;
			tmp_vec[1] = 0.0;
			tmp_vec[2] = 0.0;
		}

		for ( int ix = std::min(ixyz1[0], ixyz2[0]); ix <= std::max(ixyz1[0],ixyz2[0]); ix++ ) {
			for ( int iy = std::min(ixyz1[1], ixyz2[1]); iy <= std::max(ixyz1[1], ixyz2[1]); iy++ ) {
				for ( int iz = std::min(ixyz1[2], ixyz2[2]); iz <= std::max(ixyz1[2], ixyz2[2]); iz++ ) {
					core::Real reso = locres.get(ix, iy, iz);
					// skip if we already used this pixel for this residue
					if ( nearlyEqual(used(ix, iy, iz), used_val, 0.0001) ) { continue; }
					resos.push_back(reso);
					used( ix, iy, iz ) = used_val;
				} //iz
			} //iy
		} //ix

	}

	std::sort( resos.begin(), resos.end() );
	if ( resos.size() % 2 == 0 ) {
		return (resos[resos.size()/2]+resos[resos.size()/2+1])/2;
	} else {
		return resos[resos.size()/2];
	}

}

core::Real
get_localbfacs( core::pose::Pose & pose, core::Vector x, core::Real const MAX_DIST){
	core::Real avg_bfac = 0.0;
	core::Size count = 0;
	for ( core::Size ri=1; ri <= pose.total_residue(); ri ++ ) {
		if ( pose.residue(ri).is_virtual_residue() ) { continue; }
		core::Vector nbr_x = pose.residue(ri).nbr_atom_xyz();
		if ( (x - nbr_x).length() < MAX_DIST ) {
			for ( core::Size ai=1; ai<=pose.residue(ri).nheavyatoms(); ai++ ) {
				if ( pose.residue(ri).atom_type(ai).is_virtual() ) { continue; }
				avg_bfac += pose.pdb_info()->temperature( ri, ai );
				count++;
			}
		}
	}
	avg_bfac /= count;
	return avg_bfac;
}


int
main( int argc, char * argv[] ){
	try {

		NEW_OPT(max_dist_bfac, "distance for calculating nearby bfactors", 8.0);

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

		density.setScoreWindowContext( true );
		density.setWindow( 3 );

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scorefxn->set_weight( core::scoring::elec_dens_window, 1.0 );
		(*scorefxn)(p);

		core::Real per_res_locres = 0.0;
		TR << "PERRESCC, rosettaNum, resName, pdbNum, chID, win3_ResCC, win1_ResCC, perResLocReso" << std::endl;
		for ( core::Size resi=1; resi <= nres; resi++ ) {
			if ( !p.residue(resi).is_protein() ) continue;
			std::string resn = p.residue(resi).name3();
			density.setWindow(1);
			core::Real res1_cc = density.matchRes( resi, p.residue(resi), p, nullptr,
				false, false );
			density.setWindow(3);
			core::Real res3_cc = density.matchRes( resi, p.residue(resi), p, nullptr,
				false, false );
			core::Vector x = p.residue(resi).nbr_atom_xyz();
			per_res_locres = get_localbfacs(p, x, option[ max_dist_bfac ]() );

			core::pose::PDBInfoOP pdbinfo = p.pdb_info();
			TR << "PERRESCC," << resi << "," << resn << "," << pdbinfo->number(resi) << pdbinfo->icode(resi);
			TR << "," << pdbinfo->chain(resi) << "," << res3_cc << "," << res1_cc << "," << per_res_locres << std::endl;

		}

	} catch ( utility::excn::Exception const & e ){
		e.display();
		return -1;
	}
	return 0;
}
