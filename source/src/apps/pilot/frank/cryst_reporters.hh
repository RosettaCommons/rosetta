/// @file
/// @brief


// libRosetta headers
#include <protocols/cryst/refinable_lattice.hh>

#include <protocols/moves/PyMOLMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <utility/graph/Graph.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/vdwaals/VDW_Energy.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <devel/init.hh>

// neil stuff
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <basic/MetricValue.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight.hh>

#include <utility/tools/make_map.hh>

#include <numeric/xyz.io.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/options/keys/optimization.OptionKeys.gen.hh> // to confirm -new_sym_min

using namespace basic;
using namespace core;
using namespace core::pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TS("cryst.reporters");

typedef utility::vector1< core::Size > Sizes;
typedef numeric::xyzMatrix< core::Real > Matrix;

template <class T>
bool has_element (utility::vector1<T> const &list, T elt) {
	return (std::find( list.begin(), list.end(), elt ) != list.end() );
}

template <class T>
bool has_element (T elt, utility::vector1<T> const &list) {
	return (has_element(list,elt));
}



Vector
random_unit_vector() {
	// phi and theta both in radians
	Real const phi( 2 * numeric::constants::d::pi * numeric::random::uniform() );
	Real const theta
		( std::acos(numeric::sin_cos_range( 1.0 - 2.0 * numeric::random::uniform() ) ) );
	assert( 0 <=   phi &&   phi <= 2 * numeric::constants::d::pi );
	assert( 0 <= theta && theta <=     numeric::constants::d::pi );
	return core::kinematics::Stub().spherical( phi, theta, 1.0 );
}
//
void
randomly_shift_and_tilt_pose( Pose & pose ) { // randomly re-orient
	runtime_assert( !pose::symmetry::is_symmetric( pose ) ); // could cause trouble
	Vector const x( random_unit_vector() );
	Vector ytmp( random_unit_vector() );
	while ( x.cross( ytmp ).length_squared() < 0.1 ) ytmp = random_unit_vector();
	Vector const z( x.cross( ytmp ).normalized() ), y( z.cross(x) );
	Matrix const R( Matrix::cols( x, y, z ) );
	Vector const v(numeric::random::uniform(), numeric::random::uniform(), numeric::random::uniform());
	pose.apply_transform_Rx_plus_v( R,v);
}

void
randomly_shift_and_tilt_poses( Pose & pose1, Pose & pose2 ) { // randomly re-orient, w SAME TRANSFORM
	runtime_assert( !pose::symmetry::is_symmetric( pose1 ) ); // could cause trouble
	runtime_assert( !pose::symmetry::is_symmetric( pose2 ) ); // could cause trouble
	Vector const x( random_unit_vector() );
	Vector ytmp( random_unit_vector() );
	while ( x.cross( ytmp ).length_squared() < 0.1 ) ytmp = random_unit_vector();
	Vector const z( x.cross( ytmp ).normalized() ), y( z.cross(x) );
	Matrix const R( Matrix::cols( x, y, z ) );
	Vector const v(numeric::random::uniform(), numeric::random::uniform(), numeric::random::uniform());
	pose1.apply_transform_Rx_plus_v( R,v);
	pose2.apply_transform_Rx_plus_v( R,v);
}

struct SingleInterface {
	SingleInterface( numeric::xyzMatrix<core::Real> R_in, numeric::xyzVector<core::Real> T_in, Real cb_overlap_in ) {
		cb_overlap_ = cb_overlap_in;
		R_ = R_in;
		T_ = T_in;
	}

	Real cb_overlap_;
	numeric::xyzMatrix<core::Real> R_;
	numeric::xyzVector<core::Real> T_;
};

inline bool transforms_equiv(
	numeric::xyzMatrix<Real> const &S1, numeric::xyzVector<Real> const &T1,
	numeric::xyzMatrix<Real> const &S2, numeric::xyzVector<Real> const &T2
) {
	Real const threshold( 1e-4 );
	bool equiv( true );
	for ( Size i=1; i<= 3 && equiv ; ++i ) {
		equiv = equiv && ( fabs( T1(i)-T2(i) )< threshold );
		for ( Size j=1; j<= 3 && equiv; ++j ) {
			equiv = equiv && ( fabs( S1(i,j)-S2(i,j) )< threshold );
		}
	}
	return equiv;
}


core::Real
get_max_interface(
	Pose /*pose*/,
	utility::vector1<SingleInterface> &allInterfaces,
	utility::vector1<core::kinematics::RT> const & /*rts_all*/ ) {

	core::Real retval=0.0;
	for ( int k=1; k<=(int)allInterfaces.size(); ++k ) {
		retval = std::max( retval, allInterfaces[k].cb_overlap_ );
	}
	return retval;
}


void
reconstruct_lattice_pose_info(
	core::Size const nres_monomer,
	core::pose::Pose const & pose, // already a symmetric pose
	core::Size & base_monomer,
	core::Size & num_monomers,
	utility::vector1<core::Size> & monomer_jumps,
	utility::vector1<core::Size> & monomer_anchor_virtuals) {
	monomer_jumps.clear();
	monomer_anchor_virtuals.clear();

	conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
	core::Size nres_protein( pose.size() );
	while ( !pose.residue( nres_protein ).is_protein() ) --nres_protein;
	runtime_assert( nres_protein%nres_monomer == 0  );
	num_monomers = nres_protein/nres_monomer;
	base_monomer = 0;

	{
		core::kinematics::FoldTree const & f( pose.fold_tree() );
		for ( Size imon=1; imon<= num_monomers; ++imon ) {
			core::Size const offset( (imon-1)*nres_monomer );
			if ( symminfo->bb_is_independent( offset+1 ) ) {
				runtime_assert( !base_monomer );
				base_monomer = imon;
			}
			core::Size jumpno(0);
			for ( Size j=offset+1; j<= offset+nres_monomer; ++j ) {
				kinematics::Edge const e( f.get_residue_edge(j) );
				if ( e.is_jump() ) {
					runtime_assert( !jumpno );
					jumpno = e.label();
				}
			}
			if ( jumpno==0 ) {
				TS.Fatal << imon << " " << offset << " " << nres_monomer << std::endl;
				utility_exit_with_message("Bad jump number.");
			}
			monomer_jumps.push_back( jumpno );
			monomer_anchor_virtuals.push_back( f.upstream_jump_residue( jumpno ) );
			runtime_assert( pose.residue( monomer_anchor_virtuals.back() ).name() == "VRT" );
		}
	}
	runtime_assert( base_monomer );
}


///
///
core::Real
get_interface_score(
	Pose /*pose*/,
	utility::vector1<SingleInterface> &allInterfaces,
	utility::vector1<core::kinematics::RT> const &rts_all,
	bool verbose = false ) {
	// [STAGE 2] see if we form a connected lattice
	//    * idea: expand all contacting symmops
	//    *       ensure we generate all symmops
	//    *       ensure we generate (+/-1,0,0), (0,+/-1,0), and (0,0,+/-1)
	int EXPAND_ROUNDS=20;  // i think this is enough

	int nxformsOrig = (int)allInterfaces.size();

	// for (int j=1; j<=(int)allInterfaces.size(); ++j) {  // 1 is identity
	// numeric::xyzMatrix<Real> const &Sj=allInterfaces[j].R_;
	// numeric::xyzVector<Real> const &Tj=allInterfaces[j].T_;
	// std::cerr << "At " << j << ": " << std::endl;
	// std::cerr << "   " << Sj.xx() << " " << Sj.xy() << " " << Sj.xz() << " " << Tj[0] << std::endl;
	// std::cerr << "   " << Sj.yx() << " " << Sj.yy() << " " << Sj.yz() << " " << Tj[1] << std::endl;
	// std::cerr << "   " << Sj.zx() << " " << Sj.zy() << " " << Sj.zz() << " " << Tj[2] << std::endl;
	// }
	//
	// for (int i=2; i<=(int)rts_all.size(); ++i) {  // 1 is identity
	// numeric::xyzMatrix<Real> const &Si= rts_all[i].get_rotation();
	// numeric::xyzVector<Real> const &Ti= rts_all[i].get_translation();
	// std::cerr << "*** At " << i << ": " << std::endl;
	// std::cerr << "   " << Si.xx() << " " << Si.xy() << " " << Si.xz() << " " << Ti[0] << std::endl;
	// std::cerr << "   " << Si.yx() << " " << Si.yy() << " " << Si.yz() << " " << Ti[1] << std::endl;
	// std::cerr << "   " << Si.zx() << " " << Si.zy() << " " << Si.zz() << " " << Ti[2] << std::endl;
	// }

	bool done = false;
	for ( int rd=1; rd<=EXPAND_ROUNDS && !done; ++rd ) {
		if ( verbose ) {
			std::cerr << "Round " << rd << ":" << std::endl;
			for ( int j=1; j<=(int)allInterfaces.size(); ++j ) {  // 1 is identity
				numeric::xyzMatrix<Real> const &Sj=allInterfaces[j].R_;
				numeric::xyzVector<Real> const &Tj=allInterfaces[j].T_;
				std::cerr << "At " << j << ", score = " << allInterfaces[j].cb_overlap_ << std::endl;
				std::cerr << "  R=[ " << Sj.xx() << " " << Sj.xy() << " " << Sj.xz() << " ; "
					<< Sj.yx() << " " << Sj.yy() << " " << Sj.yz() << " ; "
					<< Sj.zx() << " " << Sj.zy() << " " << Sj.zz() << "]; T=[ "
					<< Tj[0] << " " <<  Tj[1] << " " <<Tj[2] << "]" << std::endl;
			}
		}

		done = true;
		int nxforms = (int)allInterfaces.size();
		for ( int i=1; i<=nxforms; ++i ) {
			for ( int j=1; j<=nxformsOrig; ++j ) {
				// keep in inner loop
				numeric::xyzMatrix<Real> const &Si=allInterfaces[i].R_;
				numeric::xyzVector<Real> const &Ti=allInterfaces[i].T_;
				numeric::xyzMatrix<Real> const &Sj=allInterfaces[j].R_;
				numeric::xyzVector<Real> const &Tj=allInterfaces[j].T_;

				numeric::xyzMatrix<Real> Sij = Si*Sj;
				numeric::xyzVector<Real> Tij = Ti + Si*Tj;
				Real overlap_ij = std::min( allInterfaces[i].cb_overlap_, allInterfaces[j].cb_overlap_ );

				if ( (Tij[0])>2 || (Tij[1])>2 || (Tij[2])>2 ) continue;
				if ( (Tij[0])<-2 || (Tij[1])<-2 || (Tij[2])<-2 ) continue;

				// check if it is in the set
				// we could hash these for a small speed increase...
				bool unique_xform = true;
				for ( int k=1; k<=(int)allInterfaces.size(); ++k ) {
					if ( transforms_equiv(Sij, Tij, allInterfaces[k].R_, allInterfaces[k].T_) ) {
						// check if min interface is better
						if ( allInterfaces[k].cb_overlap_ < overlap_ij ) {
							allInterfaces[k].cb_overlap_ = overlap_ij;
							done = false;
						}
						unique_xform = false;
					}
				}
				if ( unique_xform ) {
					allInterfaces.push_back( SingleInterface( Sij, Tij, overlap_ij ) );
					done = false;
				}
			}
		}
	}


	// check if we are fully connected
	//   1--- all symmops
	Real score = 1e30;
	bool connected = true;
	for ( int i=2; i<=(int)rts_all.size() && connected; ++i ) {  // 1 is identity
		bool contains_j = false;
		numeric::xyzMatrix<Real> const &Si= rts_all[i].get_rotation();
		numeric::xyzVector<Real> const &Ti= rts_all[i].get_translation();
		for ( int j=1; j<=(int)allInterfaces.size() && !contains_j; ++j ) {
			numeric::xyzMatrix<Real> const &Sj=allInterfaces[j].R_;
			numeric::xyzVector<Real> const &Tj=allInterfaces[j].T_;
			if ( transforms_equiv( Si, Ti, Sj, Tj) ) {
				contains_j = true;
				score = std::min( score, allInterfaces[j].cb_overlap_ );
			}
		}
		connected &= contains_j;
	}
	if ( !connected ) return 0.0;


	//   2--- +/- 1 in each direction
	Real score_00p=0, score_00m=0, score_0m0=0, score_0p0=0, score_p00=0, score_m00=0;
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	for ( int i=1; i<=(int)allInterfaces.size(); ++i ) {
		numeric::xyzMatrix<Real> const &Si=allInterfaces[i].R_;
		numeric::xyzVector<Real> const &Ti=allInterfaces[i].T_;
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,0, 1)) ) {
			score_00p = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,0,-1)) ) {
			score_00m = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0, 1,0)) ) {
			score_0p0 = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,-1,0)) ) {
			score_0m0 = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>( 1,0,0)) ) {
			score_p00 = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(-1,0,0)) ) {
			score_m00 = allInterfaces[i].cb_overlap_;
		}
	}
	score = std::min( score, std::min(score_00p, score_00m));
	score = std::min( score, std::min(score_0m0, score_0p0));
	score = std::min( score, std::min(score_p00, score_m00));

	return score;
}





///////////////////////////////////////////////////////////////////////////////

Real
calc_solvent_content( Real const vm ) {
	return 1.0 - ( 1.23 / vm );
}


Real
calc_vm( Real const volume, std::string const & sequence, Size const ncopies )
{
	std::map< char, int > daltons ( utility::tools::make_map(
		'A', 71,  'R', 156, 'N', 114, 'D', 115, 'C', 103,
		'Q', 128, 'E', 129, 'G', 57,  'H', 137, 'I', 113,
		'L', 113, 'K', 128, 'M', 131, 'F', 147, 'P', 97,
		'S', 87,  'T', 101, 'W', 186, 'Y', 163, 'V', 99 ) );

	Real weight(0.0);
	for ( char const c : sequence ) weight += daltons[ c ];
	return ( volume / Real( ncopies * weight ) );
}

///////////////////////////////////////////////////////////////////////////////

void
get_sc( protocols::cryst::MakeLatticeMover &setup, core::pose::Pose & pose, core::Real & total_sc, core::Real & weakest_sc, core::Real &total_sasa, core::Real &weakest_sasa) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;

	core::scoring::sc::ShapeComplementarityCalculator scc;
	scc.Init();

	SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info = ( SymmConf.Symmetry_Info() );
	core::Size nsubunits = symm_info->subunits();
	core::Size nres_asu = symm_info->num_independent_residues();

	// Figure out which chains touch chain A, and add the residues from those chains
	// into the sc surface objects
	utility::vector1<SingleInterface> sc_interfaces, sasa_interfaces;

	for ( Size i=1; i<=nres_asu; ++i ) scc.AddResidue(0, pose.residue(i));
	for ( Size i=2; i<=nsubunits; ++i ) {
		core::scoring::sc::ShapeComplementarityCalculator scc_i;
		scc_i.Init();
		for ( Size j=1; j<=nres_asu; ++j ) scc_i.AddResidue(0, pose.residue(j));

		bool contact = false;
		Size start = (i-1)*nres_asu;
		for ( Size ir=1; ir<=nres_asu; ir++ ) {
			if ( pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0 ) {
				contact = true;
				break;
			}
		}
		if ( contact ) {
			for ( Size ir=1; ir<=nres_asu; ir++ ) {
				scc.AddResidue(1, pose.residue(ir+start));
				scc_i.AddResidue(1, pose.residue(ir+start));
			}

			bool result = false;
			try {
				result = scc_i.Calc();
			} catch ( utility::excn::EXCN_Base const & e ) {
				result = false;
			}

			if ( result ) {
				numeric::xyzMatrix<core::Real> Ri;
				numeric::xyzVector<core::Real> Ti;
				setup.getRT( i, Ri, Ti );
				sc_interfaces.push_back( SingleInterface( Ri, Ti, scc_i.GetResults().sc ) );
				sasa_interfaces.push_back( SingleInterface( Ri, Ti, scc_i.GetResults().surface[2].trimmedArea ) );
			}
		}
	}
	weakest_sc = get_interface_score (pose, sc_interfaces, setup.symmops() );
	weakest_sasa = get_interface_score (pose, sasa_interfaces, setup.symmops() );

	if ( scc.Calc() ) {
		total_sc = scc.GetResults().sc;
		total_sasa = scc.GetResults().surface[2].trimmedArea;
	}
}

void
single_sasa_calc( core::pose::Pose & pose, core::Real probe_radius, core::Real & total_sasa, core::Real &total_nonpolar_sasa ) {
	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;

	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius );

	total_nonpolar_sasa = total_sasa = 0.0;

	for ( Size i=1; i<= pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( !rsd.is_protein() ) continue;
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			core::chemical::AtomType const & at( rsd.atom_type(j) );
			Real const sasa_this_atom( atom_sasa[ core::id::AtomID( j, i ) ] );
			if ( !at.is_donor() && !at.is_acceptor() && !at.is_polar_hydrogen() ) {
				total_nonpolar_sasa += sasa_this_atom;
			}
			total_sasa += sasa_this_atom;
		}
	}
}


//
void
get_sasa( protocols::cryst::MakeLatticeMover &setup, core::pose::Pose & pose, core::Real & total_sasa, core::Real & weakest_sasa, core::Real &total_nonpolar_sasa, core::Real &weakest_nonpolar_sasa) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;

	SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info = ( SymmConf.Symmetry_Info() );
	core::Size nsubunits = symm_info->subunits();
	core::Size nres_asu = symm_info->num_independent_residues();

	utility::vector1<SingleInterface> sasa_interfaces, nonpolar_sasa_interfaces;

	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	// per subunit
	for ( Size i=2; i<=nsubunits; ++i ) {
		core::pose::Pose pose_i = pose_asu;

		bool contact = false;
		Size start = (i-1)*nres_asu;
		for ( Size ir=1; ir<=nres_asu; ir++ ) {
			if ( pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0 ) {
				contact = true;
				break;
			}
		}

		if ( contact ) {
			Size const anchor( nres_asu );
			pose_i.append_residue_by_jump( pose.residue(start+1), anchor );
			for ( Size ir=2; ir<=nres_asu; ir++ ) {
				pose_i.append_residue_by_bond( pose.residue(start+ir) );
			}

			core::pose::Pose pose_i_apart = pose_i;
			pose_i_apart.apply_transform_Rx_plus_v( Matrix::identity(), Vector(-1000,-1000,-1000) );
			for ( Size ir=1; ir<= nres_asu; ++ir ) {
				pose_i_apart.replace_residue( ir, pose_i.residue(ir), false );
			}

			core::Real sasa_i=0, nonpolar_sasa_i=0, sasa_i_apart=0, nonpolar_sasa_i_apart=0;
			single_sasa_calc( pose_i, 1.4, sasa_i, nonpolar_sasa_i );
			single_sasa_calc( pose_i_apart, 1.4, sasa_i_apart, nonpolar_sasa_i_apart );

			numeric::xyzMatrix<core::Real> Ri;
			numeric::xyzVector<core::Real> Ti;
			setup.getRT( i, Ri, Ti );

			sasa_interfaces.push_back( SingleInterface( Ri, Ti, sasa_i_apart-sasa_i ) );
			nonpolar_sasa_interfaces.push_back( SingleInterface( Ri, Ti, nonpolar_sasa_i_apart-nonpolar_sasa_i ) );
		}
	}
	weakest_sasa = get_interface_score (pose, sasa_interfaces, setup.symmops() );
	weakest_nonpolar_sasa = get_interface_score (pose, nonpolar_sasa_interfaces, setup.symmops() );

	// whole complex
	core::pose::Pose pose_apart( pose );
	pose::symmetry::make_asymmetric_pose( pose_apart );
	core::pose::Pose pose_together( pose_apart );
	pose_apart.apply_transform_Rx_plus_v( Matrix::identity(), Vector(-1000,-1000,-1000) );
	for ( Size i=1; i<= nres_asu; ++i ) {
		pose_apart.replace_residue( i, pose_together.residue(i), false );
	}

	core::Real sasa_i=0, nonpolar_sasa_i=0, sasa_i_apart=0, nonpolar_sasa_i_apart=0;
	single_sasa_calc( pose_together, 1.4, sasa_i, nonpolar_sasa_i );
	single_sasa_calc( pose_apart, 1.4, sasa_i_apart, nonpolar_sasa_i_apart );
	total_sasa = sasa_i_apart-sasa_i;
	total_nonpolar_sasa = nonpolar_sasa_i_apart-nonpolar_sasa_i;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// skips non-protein
void
find_buried_unsatisfied_polars(
	Pose const & pose_in,
	core::Size nres_asu,
	utility::vector1< id::AtomID > & buried_unsatisfied_donors,
	utility::vector1< id::AtomID > & buried_unsatisfied_acceptors,
	Real const probe_radius = 1.0 )
{
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;

	bool const new_unsats( true );
	Pose pose( pose_in );

	pose::symmetry::make_asymmetric_pose( pose );
	Real const hbond_energy_threshold( -0.01 );
	Real const burial_threshold( 0.01 );

	id::AtomID_Map< Real > atom_sasa, total_hbond_energy;
	utility::vector1< Real > rsd_sasa;
	bool const use_big_polar_H( new_unsats ? false : true );
	scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, use_big_polar_H );

	if ( new_unsats ) { // accumulate sasas of dependent polar hydrogens into their parents' sasa
		for ( Size i=1; i<= pose.size(); ++i ) {
			conformation::Residue const & rsd( pose.residue(i) );
			for ( chemical::AtomIndices::const_iterator
					hnum  = rsd.Hpos_polar().begin(),
					hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
				id::AtomID const hatm( *hnum, i ), datm( rsd.atom_base( *hnum ), i );
				atom_sasa[ datm ] += atom_sasa[ hatm ];
			}
		}
	}

	// HACK!!!!!!!!!!!!!! TO TRIGGER 10A nbr graph update
	scoring::hbonds::HBondSet hbond_set;
	{
		using namespace scoring;
		ScoreFunction sf;
		sf.set_weight( hbond_sc, 1.0 );
		sf(pose);
		pose.update_residue_neighbors();
		scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
	}

	core::pose::initialize_atomid_map( total_hbond_energy, pose, 0.0 );

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		if ( true || hbond_set.allow_hbond(i) ) { // EVEN UNALLOWED HBONDS...
			scoring::hbonds::HBond const & hb( hbond_set.hbond(i) );
			id::AtomID const hatm( hb.don_hatm(), hb.don_res() );
			id::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
			total_hbond_energy[ hatm ] += hb.energy(); // unweighted
			total_hbond_energy[ aatm ] += hb.energy(); // unweighted
		}
	}

	// now find unsat+buried
	buried_unsatisfied_donors.clear();
	buried_unsatisfied_acceptors.clear();

	for ( Size i=1; i<= nres_asu; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( !rsd.is_protein() ) continue; /// SKIPPING NON-PROTEIN RIGHT NOW !!!!!!!!!!!!!!!!!!!!!!!!

		// donors
		for ( chemical::AtomIndices::const_iterator
				hnum  = rsd.Hpos_polar().begin(),
				hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			id::AtomID const hatm( *hnum, i ), datm( rsd.atom_base( *hnum ), i );
			Real const hatm_sasa( new_unsats ? atom_sasa[ datm ] : atom_sasa[ hatm ] );
			if ( total_hbond_energy[ hatm ] >= hbond_energy_threshold && hatm_sasa <= burial_threshold ) {
				buried_unsatisfied_donors.push_back( hatm );
			}
		}

		// acceptors
		for ( chemical::AtomIndices::const_iterator
				anum  = rsd.accpt_pos().begin(),
				anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
			id::AtomID const aatm( *anum, i );
			if ( total_hbond_energy[ aatm ] >= hbond_energy_threshold &&
					atom_sasa[ aatm ] <= burial_threshold ) {
				buried_unsatisfied_acceptors.push_back( aatm );
			}
		}
	} // i
}



core::Real
unsatisfied_buried_polars( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sf, core::Real probe, utility::vector1< bool > &scs_buried ) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;
	using id::AtomID;

	Size nrepeat = 5;

	Real unsats(0);
	//Real unsatdonbb(0), unsatdonsc(0), unsataccbb(0), unsataccsc(0);

	SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info = ( SymmConf.Symmetry_Info() );
	core::Size nres_asu = symm_info->num_independent_residues();

	core::pose::Pose pose_apart( pose );
	pose::symmetry::make_asymmetric_pose( pose_apart );
	core::pose::Pose pose_together( pose_apart );
	pose_apart.apply_transform_Rx_plus_v( Matrix::identity(), Vector(-1000,-1000,-1000) );
	for ( Size i=1; i<= nres_asu; ++i ) {
		pose_apart.replace_residue( i, pose_together.residue(i), false );
	}

	ScoreFunctionOP asym_sf( scoring::symmetry::asymmetrize_scorefunction( *sf ) );

	id::AtomID_Map<Real> unsat_frequency;
	core::pose::initialize_atomid_map( unsat_frequency, pose_together, 0.0 );

	Real increment(1.0/nrepeat);
	for ( Size R=1; R<= nrepeat; ++R ) {
		core::pose::Pose pose_together_i( pose );
		core::pose::Pose pose_apart_i( pose_apart );
		randomly_shift_and_tilt_poses( pose_together, pose_apart );

		utility::vector1< id::AtomID > unsatdons1, unsataccs1, unsatdons2, unsataccs2;

		find_buried_unsatisfied_polars( pose_apart_i   , nres_asu, unsatdons1, unsataccs1, probe );
		find_buried_unsatisfied_polars( pose_together_i, nres_asu, unsatdons2, unsataccs2, probe );

		for ( id::AtomID const & id : unsatdons2 ) {
			if ( !has_element( id, unsatdons1 ) ) {
				unsat_frequency[id] += increment;
				unsats += increment;
			}
		}
		//for ( id::AtomID const & id : unsatdons1 ) {
		//if ( !has_element( id, unsatdons2 ) ) TS << "Unsatisfied in unbound only: " << id << std::endl;
		//}
		for ( id::AtomID const & id : unsataccs2 ) {
			if ( !has_element( id, unsataccs1 ) ) {
				unsat_frequency[id] += increment;
				unsats += increment;
			}
		}
		//for ( id::AtomID const & id : unsataccs1 ) {
		//if ( !has_element( id, unsataccs2 ) ) TS << "Unsatisfied in unbound only: " << id << std::endl;
		//}
	}


	scs_buried.clear();
	scs_buried.resize(nres_asu, false);
	for ( Size i=1; i<= nres_asu; ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			if ( rsd.atom_is_backbone( j ) ) continue;
			if ( unsat_frequency[ AtomID( j,i ) ] >= 0.5 ) {
				scs_buried[i]=true;
			}
		}
	}

	return unsats;
}


core::Real
unsatisfied_buried_polars( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sf, core::Real probe ) {
	utility::vector1< bool > dummy;
	return unsatisfied_buried_polars( pose, sf, probe, dummy );
}

void
get_perresE( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sf,utility::vector1< core::Real > &perresE, bool int_only ) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;

	//
	core::Size nres_asu = pose.size();
	while ( pose.residue(nres_asu).aa() == core::chemical::aa_vrt ) nres_asu--;

	SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		symm_info = ( SymmConf.Symmetry_Info() );
		nres_asu = symm_info->num_independent_residues();
	}

	// per interface energies
	//   with bbhbond
	core::scoring::ScoreFunctionOP sfc = sf->clone();
	//sfc->set_weight( core::scoring::cart_bonded, 0.0 );
	//sfc->set_weight( core::scoring::res_type_constraint, 0.0 );

	core::Real Etot=0.0, E1b=0.0;
	core::scoring::methods::EnergyMethodOptions opts = sfc->energy_method_options();
	core::scoring::hbonds::HBondOptions hbo = opts.hbond_options();
	hbo.decompose_bb_hb_into_pair_energies( true );
	opts.hbond_options( hbo );
	sfc->set_energy_method_options( opts );

	perresE.clear();
	perresE.resize(nres_asu,0.0);

	(*sfc)(pose);

	Energies & energies( pose.energies() );
	EnergyGraph & energy_graph( energies.energy_graph() );
	for ( Size i=1, i_end = nres_asu; i<= i_end; ++i ) {
		perresE[i] = energies.onebody_energies(i).dot(sfc->weights());
		E1b += perresE[i];
	}

	for ( Size i=1, i_end = nres_asu; i<= i_end; ++i ) {
		conformation::Residue const & rsd1( pose.residue( i ) );
		for ( utility::graph::Graph::EdgeListIter
				iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge & edge( static_cast< EnergyEdge & > (**iru) );

			Size const j( edge.get_second_node_ind() );
			conformation::Residue const & rsd2( pose.residue( j ) );

			if ( i==j ) continue;
			if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;
			if ( int_only && (!symm_info || symm_info->bb_is_independent(rsd2.seqpos())) ) continue;

			Real thisEdgeE = edge.dot( sfc->weights() );
			if ( j<=nres_asu ) {
				perresE[i] += 0.5*thisEdgeE;
				perresE[j] += 0.5*thisEdgeE;
			} else {
				perresE[i] += thisEdgeE;
			}
			Etot += thisEdgeE;
		}
	}


	(*sf)(pose);
}

void do_sc_relax( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sf, utility::vector1<bool> residues, bool rb, int rpts=1 ) {
	// movemap (for minimization)
	core::kinematics::MoveMapOP movemap (new core::kinematics::MoveMap);
	movemap->set_bb(false); movemap->set_jump(rb);

	for ( core::Size i=1; i<=residues.size(); ++i ) {
		movemap->set_chi(i, residues[i]);
	}

	// "sc relax" w/ angles + RB
	protocols::relax::FastRelax relax_prot( sf, rpts );
	relax_prot.min_type("lbfgs_armijo_nonmonotone");
	relax_prot.set_movemap( movemap );
	relax_prot.apply(pose);
}

void do_sc_relax( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sf, bool rb, int rpts=1 ) {
	// movemap (for minimization)
	core::kinematics::MoveMapOP movemap (new core::kinematics::MoveMap);
	movemap->set_bb(false); movemap->set_chi(true); movemap->set_jump(rb);

	// "sc relax" w/ angles + RB
	protocols::relax::FastRelax relax_prot( sf, rpts );
	relax_prot.min_type("lbfgs_armijo_nonmonotone");
	relax_prot.set_movemap( movemap );
	relax_prot.apply(pose);
}

void do_sc_min( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sf, utility::vector1<bool> residues, bool rb, int rpts=1 ) {
	// movemap (for minimization)
	core::kinematics::MoveMapOP movemap (new core::kinematics::MoveMap);
	movemap->set_bb(false); movemap->set_jump(rb);

	for ( core::Size i=1; i<=residues.size(); ++i ) {
		movemap->set_chi(i, residues[i]);
	}

	// "sc relax" w/ angles + RB
	protocols::relax::FastRelax relax_prot( sf, rpts );
	relax_prot.min_type("lbfgs_armijo_nonmonotone");
	relax_prot.set_movemap( movemap );
	relax_prot.apply(pose);
}

// fast ddg
core::Real
ddg (
	core::pose::Pose &pose , core::scoring::ScoreFunctionOP sf,
	utility::vector1< core::Real > &perresddg,
	core::Real &esymm,
	core::Real &easymm,
	int ncyc=1,
	int nrepeats=3) {
	core::Real ddg=0;
	core::Real esymmPR=0.0,easymmPR=0.0;
	esymm = easymm = 0.0;
	perresddg.clear();

	core::scoring::ScoreFunctionOP sf_asu = core::scoring::symmetry::asymmetrize_scorefunction(*sf), sfc = sf->clone();
	sfc->set_weight( core::scoring::res_type_constraint , 0.0 );

	(*sfc)(pose);
	//sfc->show( pose );

	for ( int rpt = 1; rpt<=nrepeats; ++rpt ) {
		core::pose::Pose pose_asu;
		core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

		// asymm relax
		do_sc_relax( pose_asu, sf_asu, false, ncyc );

		utility::vector1< core::Real > perresAsymm, perresSymm;
		get_perresE (pose_asu, sf_asu, perresAsymm, false );
		get_perresE (pose, sf, perresSymm, false );
		if ( rpt==1 ) perresddg.resize( perresSymm.size(), 0.0 );

		ddg += (0.5*(*sfc)(pose) - (*sf_asu)(pose_asu)) / nrepeats;
		esymm += 0.5*(*sfc)(pose) / nrepeats;
		easymm += (*sf_asu)(pose_asu)/ nrepeats;
		for ( int i=1; i<=(int)perresSymm.size(); ++i ) {
			perresddg[i] += (0.5*perresSymm[i] - perresAsymm[i]) / nrepeats;
			esymmPR += 0.5*perresSymm[i] / nrepeats;
			easymmPR += perresAsymm[i] / nrepeats;
		}
	}

	return ddg;
}

void
get_interface_energies( protocols::cryst::MakeLatticeMover &setup, core::pose::Pose & pose , core::scoring::ScoreFunctionOP sf, core::Real &totalE, core::Real &intE, core::Real &weakestE ) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;

	//
	SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	core::Size nsubunits = symm_info->subunits();
	core::Size nres_asu = symm_info->num_independent_residues();

	// per interface energies
	//   with bbhbond
	core::scoring::ScoreFunctionOP sfc = sf->clone();
	sfc->set_weight( core::scoring::res_type_constraint, 0.0 );
	totalE = (*sfc)(pose);  // score before setting decompose_bb_hb_into_pair_energies to get overall energy right

	core::scoring::methods::EnergyMethodOptions opts = sfc->energy_method_options();
	core::scoring::hbonds::HBondOptions hbo = opts.hbond_options();
	hbo.decompose_bb_hb_into_pair_energies( true );
	opts.hbond_options( hbo );
	sfc->set_energy_method_options( opts );
	(*sfc)(pose);          // need to rescore with this to get hbond energies assigned to residues correctly

	utility::vector1<core::Real> score_perinterface(nsubunits,0.0);
	intE=(0.0);
	{
		(*sf)(pose);
		Energies & energies( pose.energies() );
		EnergyGraph & energy_graph( energies.energy_graph() );
		for ( Size i=1, i_end = nres_asu; i<= i_end; ++i ) {
			conformation::Residue const & rsd1( pose.residue( i ) );

			for ( utility::graph::Graph::EdgeListIter
					iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->upper_edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge & edge( static_cast< EnergyEdge & > (**iru) );

				Size const j( edge.get_second_node_ind() );
				conformation::Residue const & rsd2( pose.residue( j ) );

				if ( i==j ) continue;
				if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;
				if ( symm_info->bb_is_independent(rsd2.seqpos()) ) continue;

				Real thisEdgeE = edge.dot( sf->weights() );
				score_perinterface[symm_info->subunit_index(j)] += thisEdgeE;
				intE += thisEdgeE;
			}
		}
	}

	utility::vector1<SingleInterface> motif_interfaces;
	for ( int i=2; i<=(int)score_perinterface.size(); ++i ) {
		numeric::xyzMatrix<core::Real> Ri;
		numeric::xyzVector<core::Real> Ti;
		setup.getRT( i, Ri, Ti );
		motif_interfaces.push_back( SingleInterface( Ri, Ti, -score_perinterface[i] ) );  // negative since get_int_score looks for max
	}
	weakestE = -get_interface_score (pose, motif_interfaces, setup.symmops() ); // negative since get_int_score looks for max
}

