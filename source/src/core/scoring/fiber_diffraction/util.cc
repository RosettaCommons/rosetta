// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffractiobn/FiberDiffractionEnergyDens.cc
/// @brief  FiberDiffraction utility functions
/// @author Wojciech Potrzebowski and Ingemar Andre


// Project headers
#include <core/scoring/fiber_diffraction/util.hh>
#include <core/scoring/fiber_diffraction/CentroidScatter.hh>
#include <core/scoring/fiber_diffraction/FAScatter.hh>
#include <core/pose/util.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/io/pdb/build_pose_as_is.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


// ObjexxFCL headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <numeric/statistics/functions.hh>
#include <numeric/util.hh>

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

namespace core {
namespace scoring {
namespace fiber_diffraction {

static basic::Tracer TR("core.scoring.fiber_diffraction.util");

void setup_cylindrical_coords(
	pose::Pose const & pose,
	core::Size & natoms,
	utility::vector1< Size > & atom_type_number,
	std::map<  core::id::AtomID, core::Size > & AtomID_to_atomnbr,
	utility::vector1< Real > & phi,
	utility::vector1< Real > & z,
	utility::vector1< Real > & r,
	utility::vector1< Real > & bfactors
) {

	// Are we symmetric?
	const core::conformation::symmetry::SymmetryInfo *symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >(
			pose.conformation()).Symmetry_Info().get();
	}

	if ( !symminfo ) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}

	natoms = 0;
	Size total_atoms(0);
	for ( Size res=1; res <= pose.total_residue(); ++res ) {
		conformation::Residue const &rsd_i (pose.residue(res));
		for ( Size j=1 ; j<=rsd_i.natoms(); ++j ) {
			total_atoms++;
		}
	}

	for ( Size res=1; res <= pose.total_residue(); ++res ) {
		if ( ! symminfo->bb_is_independent( res ) ) continue;
		conformation::Residue const &rsd_i (pose.residue(res));
		for ( Size j=1 ; j<=rsd_i.natoms(); ++j ) {
			id::AtomID id( j, res );
			chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			std::string elt_i = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			if ( elt_i == "H" || elt_i == "X" ) continue;
			natoms++;
			AtomID_to_atomnbr.insert( std::make_pair( id, natoms ) );
			Real B = pose.pdb_info()->temperature( res, j );
			bfactors.push_back(B);
			if ( elt_i != "C" &&
					elt_i != "O" &&
					elt_i != "N" &&
					elt_i != "S" &&
					elt_i != "P" ) {
				TR.Error<<" Unsupported atom type for fiber diffraction " << elt_i<<std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}

			if ( elt_i == "C" ) atom_type_number.push_back(1);
			if ( elt_i == "O" ) atom_type_number.push_back(2);
			if ( elt_i == "N" ) atom_type_number.push_back(3);
			if ( elt_i == "S" ) atom_type_number.push_back(4);
			if ( elt_i == "P" ) atom_type_number.push_back(5);
			// calculate phi and z
			numeric::xyzVector< Real > coord ( pose.xyz(id) );
			core::Real z_j(coord(3));
			numeric::xyzVector< Real > z_axis (0,0,z_j);
			numeric::xyzVector< Real > x_axis (1,0,0);
			numeric::xyzVector< Real > transformed_coord( coord - z_axis );
			core::Real r_j( transformed_coord.length());

			r.push_back( r_j );
			z.push_back( z_j );
			core::Real angle ( atan2( transformed_coord(2), transformed_coord(1) ));
			if ( angle < 0 ) angle += 2*M_PI;
			if ( fabs(r_j)<1e-2 ) angle=0;
			phi.push_back( angle );
		}
	}
	if ( natoms < 1 ) {
		utility_exit_with_message("No atoms loaded for fiber diffraction calculation!");
	}
}

void find_pitch(
	pose::Pose const & pose,
	core::Real & pitch
) {
	const core::conformation::symmetry::SymmetryInfo *symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >(
			pose.conformation()).Symmetry_Info().get();
	}

	if ( !symminfo ) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}
	Size master_res(1);
	for ( Size res=1; res <= pose.total_residue(); ++res ) {
		if ( symminfo->bb_is_independent( res ) ) {
			master_res = res;
			break;
		}
	}
	pitch = 1e10;
	numeric::xyzVector< Real > start ( pose.xyz( id::AtomID( 1, master_res ) ) );
	for ( utility::vector0 < Size>::const_iterator
			clone     = symminfo->bb_clones( master_res ).begin(),
			clone_end = symminfo->bb_clones( master_res ).end();
			clone != clone_end; ++clone ) {
		numeric::xyzVector< Real > stop ( pose.xyz( id::AtomID( 1, *clone ) ) );
		Real z_diff( fabs( start(3) -stop(3) ) );
		if ( z_diff < pitch ) pitch = z_diff;
	}
	TR << " Pitch is calculated to " << pitch << std::endl;
	if ( pitch == 1e10 ) {
		utility_exit_with_message("Pitch could not be calculated. Aborting...");
	}
}

void
find_min_xyz(
	pose::Pose const & pose,
	core::Real &minX,
	core::Real &minY,
	core::Real &minZ,
	core::Real &maxX,
	core::Real &maxY,
	core::Real &maxZ
) {

	// Are we symmetric?
	const core::conformation::symmetry::SymmetryInfo *symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >(
			pose.conformation()).Symmetry_Info().get();
	}

	if ( !symminfo ) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}

	maxX = -1e5;
	maxY = -1e5;
	maxZ = -1e5;
	minX = 1e5;
	minY = 1e5;
	minZ = 1e5;

	int nres = pose.total_residue();

	for ( int i=1 ; i<=nres; ++i ) {

		if ( ! symminfo->bb_is_independent( i ) ) continue;

		conformation::Residue const &rsd_i (pose.residue(i));

		// skip vrts & masked reses
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

		int nheavyatoms = rsd_i.nheavyatoms();

		for ( int j=1 ; j<=nheavyatoms; ++j ) {
			conformation::Atom const &atm_i( rsd_i.atom(j) ); //(pose.residue(i).atom("CA"));

			Real x = atm_i.xyz()[0];
			Real y = atm_i.xyz()[1];
			Real z = atm_i.xyz()[2];


			if ( x > maxX ) maxX = x;
			if ( y > maxY ) maxY = y;
			if ( z > maxZ ) maxZ = z;

			if ( x < minX ) minX = x;
			if ( y < minY ) minY = y;
			if ( z < minZ ) minZ = z;
		}
	}
}

void
find_max_r(
	pose::Pose const & pose,
	core::Real & maxR
) {
	// Are we symmetric?
	const core::conformation::symmetry::SymmetryInfo *symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >(
			pose.conformation()).Symmetry_Info().get();
	}

	if ( !symminfo ) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}

	maxR = 0;
	core::Real maxR2 = 1e-5;
	int nres = pose.total_residue();
	for ( int i=1 ; i<=nres; ++i ) {
		if ( ! symminfo->bb_is_independent( i ) ) continue;
		conformation::Residue const &rsd_i (pose.residue(i));
		// skip vrts & masked reses
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		int nheavyatoms = rsd_i.nheavyatoms();
		for ( int j=1 ; j<=nheavyatoms; ++j ) {
			conformation::Atom const &atm_i( rsd_i.atom(j) ); //(pose.residue(i).atom("CA"));
			Real x = atm_i.xyz()[0];
			Real y = atm_i.xyz()[1];
			core::Real r2( x*x+y*y );
			if ( r2 > maxR2 ) maxR2 = r2;
		}
	}
	maxR = sqrt(maxR2);
}

void
find_num_scattering_atoms(
	pose::Pose & pose,
	core::Size & nscatterers
)
{
	nscatterers = 0;
	// Are we symmetric?
	const core::conformation::symmetry::SymmetryInfo *symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >(
			pose.conformation()).Symmetry_Info().get();
	}

	int nres = pose.total_residue();
	for ( int i=1 ; i<=nres; ++i ) {
		if ( ! symminfo->bb_is_independent( i ) ) continue;
		conformation::Residue const &rsd_i (pose.residue(i));
		// skip vrts & masked reses
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		nscatterers +=  rsd_i.nheavyatoms();
	}
}

void
centroid_scatter(
	const std::string & res_name,
	OneGaussianScattering & sig_centroid
) {
	using namespace core::conformation;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	ResidueOP fa_res( ResidueFactory::create_residue( rsd_set->name_map(res_name ) ) );

	Atoms::const_iterator sc_begin( fa_res->sidechainAtoms_begin() );
	Atoms::const_iterator sc_end( fa_res->atom_end() );
	core::Real sigma_tot(0), weight_tot(0);
	core::Size num_atoms(0);
	for ( Atoms::const_iterator atom=sc_begin; atom!=sc_end; ++atom ) {
		chemical::AtomTypeSet const & atom_type_set( fa_res->atom_type_set() );
		std::string elt_i = atom_type_set[ atom->type() ].element();
		if ( elt_i == "H" ) continue;
		OneGaussianScattering sig_atom = get_A( elt_i );
		sigma_tot += sig_atom.s();
		weight_tot += sig_atom.a();
		num_atoms++;
	}
	if ( num_atoms == 0 ) ++num_atoms;
	//core::Real sigma_av ( sigma_tot/num_atoms );
	//core::Real weight_av ( weight_tot/num_atoms );
	sig_centroid = OneGaussianScattering( weight_tot,sigma_tot );
}

utility::vector1< OneGaussianScattering > setup_centroid_scatter(
	pose::Pose & pose
) {
	utility::vector1< OneGaussianScattering > sig_centroid_;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_CEN_SCATTERING ) ) {
		CentroidScatter & centroid_sct_factor( retrieve_centroid_scatter_from_pose( pose ) );
		sig_centroid_ = centroid_sct_factor.getValues();
	} else {
		sig_centroid_.resize( chemical::num_canonical_aas );
		for ( int aa = 1; aa <= chemical::num_canonical_aas; ++aa ) {
			std::string name_3 = core::chemical::name_from_aa( core::chemical::AA(aa) );
			OneGaussianScattering scatter;
			centroid_scatter( name_3, scatter );
			sig_centroid_[aa] = scatter;
		}
		pose.data().set(core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_CEN_SCATTERING,
			basic::datacache::DataCache_CacheableData::DataOP( new CentroidScatter( sig_centroid_ ) ) );
	}
	return sig_centroid_;
}

bool isPowerOfTwo(int n)
{
	if ( n == 0 ) return false;
	return !((n-1) & n);
}

utility::vector0< utility::vector1< utility::vector1< core::Real > > > setup_form_factors (
	pose::Pose & pose,
	core::Size const & lmax,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	core::Real const & c,
	core::Real const & B_factor,
	core::Real const & B_factor_solv,
	core::Real const & Ksolv
) {
	utility::vector0< utility::vector1< utility::vector1< core::Real > > > form_factors_;
	TR << "Calculating form factors..." << std::endl;
	std::string atoms[] = {"C", "O", "N", "S", "P" };

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_FA_SCATTERING ) ) {
		FAScatter & fa_sct_factor( retrieve_fa_scatter_from_pose( pose ) );
		form_factors_  = fa_sct_factor.getValues();
	} else {
		form_factors_.resize(lmax+1);
		for ( Size i=0; i <= lmax; ++i ) {
			form_factors_[i].resize(5);
		}

		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R = 1; R <= layer_lines_R[l].size(); ++R ) {
				Real Rinv=layer_lines_R[l][R];
				Real D_2=Rinv*Rinv + l*l/(c*c);
				Real solv_factor ( 1-Ksolv*exp( -B_factor_solv*D_2/4 ) );
				Real B_factor_exp( exp( -B_factor*D_2/4 ) );
				for ( Size atom =1; atom <= 5; ++atom ) {
					core::scoring::fiber_diffraction::KromerMann scatter(
						core::scoring::fiber_diffraction::get_km( atoms[atom-1] ) );
					form_factors_[l][atom].push_back( B_factor_exp*solv_factor*scatter.f0( D_2 ) );
				}
			}
		}
		pose.data().set(core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_FA_SCATTERING,
			basic::datacache::DataCache_CacheableData::DataOP( new FAScatter( form_factors_ ) ) );
	}
	return form_factors_;
}

//generating Shannon sampling points
void generate_shannon_points(
	core::Size const & lmax,
	core::Real const & dmax,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0 < core::Size > ::iterator & sampling_points_l,
	utility::vector0 < utility::vector1< core::Real > >::iterator & shannon_points_lS,
	utility::vector0 < core::Size > ::iterator & lowest_bessel_orders_l,
	utility::vector0 < utility::vector0 < int > >::iterator const & nvals
) {
	utility::vector0 < core::Size > sampling_points_l_;
	utility::vector0 < utility::vector1< core::Real > > shannon_points_lS_;
	utility::vector0 < core::Size > lowest_bessel_orders_l_;
	core::Real RinvMax(0.0);
	core::Real RinvMin(0.0);
	if ( sampling_points_l_.size() == 0 ) {
		sampling_points_l_.resize(lmax+1);
		shannon_points_lS_.resize(lmax+1);
		lowest_bessel_orders_l_.resize(lmax+1);
		for ( Size l=0; l <= lmax; ++l ) {
			Size max_R_values( layer_lines_R[l].size() );
			if ( max_R_values==0 ) continue;
			lowest_bessel_orders_l_[l] = 1000;
			Size max_b_order( nvals[l].size() );
			for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {
				int n( nvals[l][b_order-1] );
				core::Size abs_n( abs(n) );
				if ( abs_n<lowest_bessel_orders_l_[l] ) {
					lowest_bessel_orders_l_[l] = abs_n;
				}
			}
			RinvMax = layer_lines_R[l][max_R_values];
			RinvMin = layer_lines_R[l][1];
			sampling_points_l_[l]  = (RinvMax - RinvMin)*4*dmax;
			shannon_points_lS_[l].resize(sampling_points_l_[l]);
			for ( Size sn = 1; sn<=sampling_points_l_[l]; ++sn ) {
				shannon_points_lS_[l][sn] = RinvMin + sn/(4*dmax);
			}
		}
	}
	sampling_points_l = sampling_points_l_.begin();
	shannon_points_lS = shannon_points_lS_.begin();
	lowest_bessel_orders_l = lowest_bessel_orders_l_.begin();
	TR<<"Shannon done"<<std::endl;
}

void bessel_roots(
	core::Size const & lmax,
	core::Real const & c,
	core::Real const & res_cutoff_high,
	core::Real const & res_cutoff_low,
	core::Real & structure_cutoff,
	utility::vector0 < utility::vector1< core::Real > > & bessel_roots_lE_,
	utility::vector0 < core::Size >  & sampling_points_l_,
	utility::vector0 < core::Size >  & lowest_bessel_orders_l_,
	utility::vector0 < core::Size >  & highest_resolution_l_,
	utility::vector0 < core::Size >  & lowest_resolution_l_,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0 < utility::vector0 < int > >::iterator const & nvals
) {
	core::Size max_bessel_roots = 400;
	core::Real Rinv, Rinv_lc;
	core::Real res_cutoff_high2(res_cutoff_high*res_cutoff_high);
	core::Real res_cutoff_low2(res_cutoff_low*res_cutoff_low);
	core::Real structure_cutoff_(0.0);

	sampling_points_l_.resize(lmax+1);
	lowest_bessel_orders_l_.resize(lmax+1);
	highest_resolution_l_.resize(lmax+1);
	lowest_resolution_l_.resize(lmax+1);
	bessel_roots_lE_.resize(lmax+1);
	for ( Size l=0; l <= lmax; ++l ) {
		Size max_R_values( layer_lines_R[l].size() );
		if ( max_R_values==0 ) continue;
		lowest_bessel_orders_l_[l] = 1000;
		Size max_b_order( nvals[l].size() );
		for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {
			int n( nvals[l][b_order-1] );
			core::Size abs_n( abs(n) );
			if ( abs_n<lowest_bessel_orders_l_[l] ) {
				lowest_bessel_orders_l_[l] = abs_n;
			}
		}
		core::Real high_res_value_ = 1000.0;
		core::Real low_res_value_ = -1000.0;
		for ( Size R=1; R<= max_R_values; ++R ) {
			Rinv=layer_lines_R[l][R];
			Rinv_lc = Rinv*Rinv+(l/c)*(l/c);
			//Resolution cutoff is set to be on the safe side
			//When using trim data it should be ok, without
			if ( 1/Rinv<high_res_value_ && Rinv_lc <= res_cutoff_high2 ) {
				high_res_value_=1/Rinv;
				highest_resolution_l_[l]=R;
			}
			if ( 1/Rinv>low_res_value_ && Rinv_lc >= res_cutoff_low2 ) {
				low_res_value_=1/Rinv;
				lowest_resolution_l_[l]=R;
			}
		}
		utility::vector1< core::Real > zeroj;
		zeroj.resize(max_bessel_roots);
		core::Size npoints;
		structure_cutoff_ = structure_cutoff/high_res_value_;
		rootj(2*lowest_bessel_orders_l_[l],structure_cutoff_, zeroj, npoints);
		bessel_roots_lE_[l].resize(npoints);
		core::Size sm_points(1);
		for ( Size e=1; e<=npoints; ++e ) {
			//Rounding to five digist that it then is in the same format as Rinv
			core::Real zeroj_tmp = zeroj[e]/structure_cutoff;
			zeroj_tmp = ((int)(zeroj_tmp * 10000 + .5) / 10000.0);
			if ( zeroj_tmp*zeroj_tmp+(l/c)*(l/c) < res_cutoff_low2 || zeroj_tmp*zeroj_tmp+(l/c)*(l/c) > res_cutoff_high2 ) continue;
			bessel_roots_lE_[l][sm_points] = zeroj_tmp;
			sm_points++;
		}
		sampling_points_l_[l] = sm_points-1;
	}
}

//Used by FreeSet app
void interpolate_sampled_to_grid(
	core::Size const & lmax,
	utility::vector0 < utility::vector1< core::Real > > const & bessel_roots_lE,
	utility::vector0 < core::Size > const & sampling_points_l,
	utility::vector0 < core::Size > const & highest_resolution_l,
	utility::vector0 < core::Size > const & lowest_resolution_l,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0 < utility::vector1< core::Size > > & selected_R_l_,
	utility::vector0 < utility::vector1< core::Real > > & selected_Rinv_l_
) {
	core::Size selected_R(1);
	selected_Rinv_l_.resize(lmax+1);
	selected_R_l_.resize(lmax+1);
	for ( Size l=0; l <= lmax; ++l ) {
		selected_R_l_[l].resize(sampling_points_l[l]);
		selected_Rinv_l_[l].resize(sampling_points_l[l]);
		for ( core::Size b=1; b<=sampling_points_l[l]; ++b ) {
			core::Real min_dist(1000);
			selected_R = 1;
			for ( core::Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				if ( R<lowest_resolution_l[l] || R>highest_resolution_l[l] ) continue;
				if ( fabs(bessel_roots_lE[l][b]-layer_lines_R[l][R]) < min_dist ) {
					min_dist = fabs(bessel_roots_lE[l][b]-layer_lines_R[l][R]);
					selected_R = R;
				}
			}
			selected_R_l_[l][b] = selected_R;
			selected_Rinv_l_[l][b] = layer_lines_R[l][selected_R];
		}
	}
}

void calculate_I_of_E(
	core::Size const & lmax,
	core::Size const & k_iteration,
	utility::vector0 < utility::vector1< core::Real > > const & sampling_points_lE,
	core::Size const & natoms,
	core::Size const & c_,
	utility::vector0< utility::vector0 < int > >::iterator const & nvals,
	utility::vector1< Size > const & atom_type_number,
	utility::vector1< Real > const & phi,
	utility::vector1< Real > const & z,
	utility::vector1< Real > const & r,
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator const & form_factors,
	utility::vector0< utility::vector1 < Real > > & I_E
) {
	core::Real Rinv(0.0);
	for ( Size l=0; l <= lmax; ++l ) {
		if ( l==k_iteration ) continue;
		Size max_b_order( nvals[l].size() );
		for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {
			int n( nvals[l][b_order-1] );
			int abs_n( abs(n) );
			// Precalculate phases. Initialize vector
			utility::vector1< utility::vector1 < Real > > phases;
			phases.resize(natoms);
			for ( Size atom=1; atom <= natoms; ++atom ) phases[atom].resize(natoms);
			// store phases for each value of n
			for ( Size atom1=1; atom1 <= natoms; ++atom1 ) {
				for ( Size atom2=atom1+1; atom2 <= natoms; ++atom2 ) {
					Real phase ( cos( n*(phi[atom2]-phi[atom1] )+2*M_PI*l/c_*( z[atom1]-z[atom2] ) ) );
					phases[atom1][atom2] = phase;
				}
			}
			// loop over all sampling points
			Size max_R_values( sampling_points_lE[l].size() );
			I_E[l].resize(max_R_values);
			for ( Size R=1; R<= max_R_values; ++R ) {
				Rinv=sampling_points_lE[l][R];
				if ( Rinv==0 ) continue;
				Real x_factor( 2*M_PI*Rinv );
				// precalculate jn for each n,R
				utility::vector1< Real > jn_vec;
				jn_vec.resize(natoms);
				// Calculate scattering...
				for ( Size atom1=1; atom1 <= natoms; ++atom1 ) {
					// Calculate only the jn values once.
					if ( atom1 == 1 ) {
						for ( Size at=1; at <= natoms; ++at ) {
							jn_vec[at] = jn(n,x_factor*r[at]);
						}
					}
					Real X1 (x_factor*r[atom1]);
					// calculate R_min
					if ( abs_n > X1 +2 ) continue;
					Real jn1 (jn_vec[atom1]);
					Size atom_type_index1 ( atom_type_number[atom1] );
					Real dummy( form_factors[l][atom_type_index1][R]*jn1 );
					Real dummy_sum( 0 );
					I_E[l][R] += dummy*form_factors[l][atom_type_index1][R]*jn1;
					for ( Size atom2=atom1+1; atom2 <= natoms; ++atom2 ) {
						Real X2 (x_factor*r[atom2]);
						if ( abs_n > X2 +2 ) continue;
						Real jn2 (jn_vec[atom2]);
						Size atom_type_index2 ( atom_type_number[ atom2 ] );
						dummy_sum += form_factors[l][atom_type_index2][R]*jn2*phases[atom1][atom2];
					} // atoms2
					I_E[l][R] +=2*dummy*dummy_sum;
				} // atom1
			} // R
		} //b_order
	} // l
	//I_E = I_E_.begin();
}

core::Real calculate_chi_k(
	core::Size const & lmax,
	core::Size const & k_iteration,
	utility::vector0 < utility::vector1< core::Size > > const & selected_R_l,
	utility::vector0 < utility::vector1< core::Real > > const & selected_Rinv_l,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_I,
	utility::vector0 < utility::vector1< core::Real > > const & I_interpolated
) {
	core::Real prod(0);
	core::Real square_obs_ (0);
	core::Real sum_obs_ (0);
	core::Size R;
	for ( Size l=0; l <= lmax; ++l ) {
		if ( l==k_iteration ) continue;
		for ( Size bin=1; bin<=selected_Rinv_l[l].size(); ++bin ) {
			R = selected_R_l[l][bin];
			if ( layer_lines_I[l][R] == 0.0 ) continue;
			core::Real I_obs_ ( layer_lines_I[l][R]*layer_lines_I[l][R] );
			prod += I_interpolated[l][bin]*I_obs_;
			square_obs_ += I_obs_*I_obs_;
			sum_obs_ +=  I_obs_;
		}
	}

	core::Real scale_factor_ = square_obs_/prod;
	core::Real chi2_(0);
	for ( Size l=0; l <= lmax; ++l ) {
		if ( l==k_iteration ) continue;
		for ( Size bin=1; bin<=selected_Rinv_l[l].size(); ++bin ) {
			R = selected_R_l[l][bin];
			chi2_ +=  (scale_factor_*I_interpolated[l][bin]-layer_lines_I[l][R]*layer_lines_I[l][R])*\
				(scale_factor_*I_interpolated[l][bin]-layer_lines_I[l][R]*layer_lines_I[l][R]);
		}
	}
	chi2_ /=square_obs_;
	return chi2_;
}

core::Real calculate_chi2_free(
	pose::Pose & pose,
	core::Size const & chi_free_iterations_,
	core::Size const & lmax,
	//utility::vector0 < utility::vector1< core::Real > >::iterator const & sampling_points_lE,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_I,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	core::Size const & natoms,
	core::Size const & c_,
	utility::vector0 < utility::vector0 < int > >::iterator const & nvals,
	utility::vector1< Size > const & atom_type_number,
	utility::vector1< Real > const & phi,
	utility::vector1< Real > const & z,
	utility::vector1< Real > const & r,
	core::Real const b_factor_,
	core::Real const b_factor_solv,
	core::Real const b_factor_solv_K
) {
	core::Real chi2free(0);
	core::Real chi2_k(0);
	utility::vector1<core::Real> chi2_k_vec;
	chi2_k_vec.resize(chi_free_iterations_+1);
	//utility::vector0 < utility::vector1< utility::vector1< core::Real > > >::iterator form_factors;
	//setup_form_factors( form_factors, lmax, sampling_points_lE, c_, b_factor_, b_factor_solv, b_factor_solv_K );
	//TODO: Check this function
	utility::vector0< utility::vector1< utility::vector1< core::Real > > > form_factors_(
		setup_form_factors( pose, lmax, layer_lines_R, c_, b_factor_, b_factor_solv, b_factor_solv_K ));
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator form_factors(form_factors_.begin());


	for ( Size k=0; k < chi_free_iterations_; ++k ) {
		utility::vector0 < utility::vector1< core::Real > > I_interpolated;
		utility::vector0 < utility::vector1< core::Real > > selected_Rinv_l;
		utility::vector0 < utility::vector1< core::Size > > selected_R_l;
		selected_Rinv_l.resize(lmax+1);
		selected_R_l.resize(lmax+1);
		I_interpolated.resize(lmax+1);

		sample_layer_lines_from_bins(lmax, k, layer_lines_R, selected_R_l, selected_Rinv_l);
		for ( Size l=0; l <= lmax; ++l ) {
			core::Size overzero(0);
			for ( Size R=1; R<=selected_Rinv_l[l].size(); ++R ) {
				if ( selected_Rinv_l[l][R]>0.0 ) overzero++;
			}
		}

		calculate_I_of_E(lmax, k, selected_Rinv_l, natoms,\
			c_, nvals, atom_type_number, phi, z, r, form_factors, I_interpolated);

		chi2_k = calculate_chi_k(lmax, k, selected_R_l,\
			selected_Rinv_l, layer_lines_I, I_interpolated);
		TR<<"k chi2 "<<k<<" : "<<chi2_k<<std::endl;
		chi2_k_vec[k+1]=chi2_k;
		chi2free += chi2_k;
	}
	core::Real m1 = numeric::statistics::mean( chi2_k_vec.begin(),chi2_k_vec.end(),0.0 );
	core::Real sd1 = numeric::statistics::std_dev_with_provided_mean( chi2_k_vec.begin(), chi2_k_vec.end(), m1 );
	core::Real md1 = numeric::median(chi2_k_vec);
	//Calculating average
	TR << "Median, mean and standard dev "<<md1<<" "<<m1<<" "<<sd1<<std::endl;
	chi2free /= chi_free_iterations_;
	return chi2free;
}

//Omitting one layer line at the time
void sample_layer_lines_from_bins(
	core::Size const & lmax,
	core::Size const & k_iteration,
	utility::vector0< utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0< utility::vector1< core::Size > > & selected_R_l_,
	utility::vector0< utility::vector1< core::Real > > & selected_Rinv_l_
) {
	core::Size max_R_values;
	for ( Size l=0; l <= lmax; ++l ) {
		max_R_values = layer_lines_R[l].size();
		selected_Rinv_l_[l].resize(max_R_values);
		selected_R_l_[l].resize(max_R_values);
		for ( Size R=1; R<=max_R_values; ++R ) {
			if ( l==k_iteration ) {
				selected_R_l_[l][R] = 0;
				selected_Rinv_l_[l][R] = 0.0;
			} else {
				selected_R_l_[l][R] = R;
				selected_Rinv_l_[l][R] = layer_lines_R[l][R];
			}
		}
	}
}

void rootj(
	int const & N,
	core::Real const & CUTOFF,
	utility::vector1< core::Real > & zeros,
	core::Size & npoints
)  {
	core::Real zeroJ,b0,b1,b2,b3,b5,b7;
	core::Real t0,t1,t3,t5,t7,fN,fK;
	core::Real c1,c2,c3,c4,c5;
	core::Real f1,f2,f3,tol,errJ,pi;
	int ierror,nitmx;
	int k=1;
	tol=1E-8; nitmx=10; pi=4.0*atan(1.0);
	c1=1.8557571; c2=1.033150; c3=0.00397; c4=0.0908; c5=0.043;
	fN = 1.0*N;

	//first zero
	if ( N==0 ) {
		zeroJ = c1+c2-c3-c4+c5;
		secant(N,nitmx,tol,&zeroJ,&ierror);
		zeros[k] = zeroJ;
	} else {
		f1 = pow(fN,(1.0/3.0));
		f2 = f1*f1*fN;
		f3 = f1*fN*fN;
		zeroJ = fN+c1*f1+(c2/f1)-(c3/fN)-(c4/f2)+(c5/f3);
		secant(N,nitmx,tol,&zeroJ,&ierror);
		zeros[k]=zeroJ;
	}

	t0 = 4.0*fN*fN;
	t1 = t0-1.0;
	t3 = 4.0*t1*(7.0*t0-31.0);
	t5 = 32.0*t1*((83.0*t0-982.0)*t0+3779.0);
	t7 = 64.0*t1*(((6949.0*t0-153855.0)*t0+1585743.0)*t0-6277237.0);

	//next zeros
	k=2;
	while ( zeroJ<=CUTOFF ) {
		fK = 1.0*k;
		//Mac Mahon series
		b0 = (fK+0.5*fN-0.25)*pi;
		b1 = 8.0*b0;
		b2 = b1*b1;
		b3 = 3.0*b1*b2;
		b5 = 5.0*b3*b2;
		b7 = 7.0*b5*b2;

		zeroJ = b0-(t1/b1)-(t3/b3)-(t5/b5)-(t7/b7);

		errJ=fabs(jn(N,zeroJ));

		if ( errJ > tol ) secant(N,nitmx,tol,&zeroJ,&ierror);
		zeros[k]=zeroJ;
		k++;
	}
	npoints = k-1;
}

void secant(
	int const & N,
	int const & nitmx,
	core::Real tol,
	core::Real *zeroJ,
	int *ier)  {

	core::Real p0,p1,q0,q1;
	core::Real dP(0.0),p(0.0);
	core::Real c[3];
	core::Size nev,ntry;
	int it;
	c[1]=0.95; c[2]=0.999;
	*ier=0;

	for ( ntry=1; ntry<=2; ntry++ ) {
		p0 = c[ntry]*(*zeroJ);
		p1 = *zeroJ;
		nev = 2;
		q0 = jn(N,p0);
		q1 = jn(N,p1);
		for ( it = 1; it <= nitmx; it++ ) {
			if ( q1==q0 ) break;
			p = p1-q1*(p1-p0)/(q1-q0);
			dP = p-p1;
			if ( it!=1 ) {
				if ( fabs(dP) < tol ) break;
			}
			nev = nev+1;
			p0 = p1;
			q0 = q1;
			p1 = p;
			q1 = jn(N,p1);
		}
		if ( fabs(dP) < tol ) break;
		*ier = ntry;
	}
	*zeroJ = p;
}


} // fiber_diffraction
} // scoring
} // core
