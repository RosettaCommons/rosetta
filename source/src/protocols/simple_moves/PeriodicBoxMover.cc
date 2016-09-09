// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/PeriodicBoxMover.cc
/// @brief Mover for running liquid simulation (or related other)
/// @author Frank DiMaio & Hahnbeom Park

// Unit headers
#include <protocols/simple_moves/PeriodicBoxMover.hh>
#include <protocols/simple_moves/PeriodicBoxMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator_Private.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/CrystInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/graph/Graph.fwd.hh>
#include <utility/graph/Graph.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/kinematics/DomainMap.fwd.hh>
//#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/EnergyEdge.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/packer_neighbors.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/types.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizerMap.hh>

#include <numeric/random/random.hh>
#include <cmath>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#include <basic/Tracer.hh>

using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.moves.PeriodicBoxMover");

#include <utility/tag/Tag.hh>

#include <core/pose/util.hh>

#include <fstream>
#include <cmath>

namespace protocols {
namespace simple_moves {


// rotate about a random axis by the specified angle
// this may duplicate other function in Rosetta but I could not locate it
numeric::xyzMatrix<core::Real> random_rotation(core::Real magnitude) {
	core::Real u1=numeric::random::rg().uniform();
	core::Real u2=numeric::random::rg().uniform();
	core::Real theta = std::acos(numeric::sin_cos_range( 1.0 - 2.0*u1 ));
	core::Real phi = 2*numeric::constants::d::pi * u2;
	core::Vector axis ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) );
	return (numeric::rotation_matrix_degrees( axis, magnitude ));
}

std::string
PeriodicBoxMoverCreator::keyname() const {
	return PeriodicBoxMoverCreator::mover_name();
}

protocols::moves::MoverOP
PeriodicBoxMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PeriodicBoxMover );
}

std::string
PeriodicBoxMoverCreator::mover_name() {
	return "PeriodicBoxMover";
}

PeriodicBoxMover::PeriodicBoxMover()
: Mover("PeriodicBoxMover") {

	// 125 molecule sim, 1g/ml
	nmol_side_ = 5;
	initial_density_ = 1;

	// total guess
	temp_ = 160;

	// from paper
	vol_step_ = 40.0;
	trans_step_ = 0.15;
	rot_step_ = 15.0;
	tor_step_ = 15.0;
	probability_rigid_ = 0.9;
	resize_vol_every_ = 600;
	nsteps_equilibrate_ = 500000;
	nsteps_sim_ = 1000000;
	istart_ = 0;

	dump_every_ = 0;
	report_every_ = 0;
	report_silent_ = "";
	report_scorefile_ = "";

	origin_ = 0;

	P0_ = 1.0; // unit in atmosphere
	rg_atom1_ = "";
	rg_atom2_ = "";

	central_molecule_pdb_ = "";
	central_resno_ = 0;

	thermodynamic_data_ = ThermodynamicData();
	correct_LJtruncation_ = false;
	ljcorrection_factor_ = 0.0;

	sf_ = core::scoring::get_score_function();
}

std::string
PeriodicBoxMover::get_name() const {
	return PeriodicBoxMoverCreator::mover_name();
}

protocols::moves::MoverOP
PeriodicBoxMover::clone() const {
	return protocols::moves::MoverOP( new PeriodicBoxMover( *this ) );
}

protocols::moves::MoverOP
PeriodicBoxMover::fresh_instance() const {
	return protocols::moves::MoverOP( new PeriodicBoxMover );
}

void
PeriodicBoxMover::dump_ASU( Pose & pose, core::Size &lattice_jump, std::string filename ) {
	Pose pose_asu;
	// lets skip this for now
	//core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);
	pose_asu = pose;

	core::Real lattice = std::abs( pose.jump(lattice_jump).get_translation()[0] );
	core::io::CrystInfo ci;
	ci.A(lattice); ci.B(lattice); ci.C(lattice);
	ci.alpha(90.0); ci.beta(90.0); ci.gamma(90.0);
	ci.spacegroup("P 1");

	//pose_asu.conformation().delete_residue_slow( pose_asu.size() );
	pose_asu.pdb_info()->set_crystinfo(ci);

	pose_asu.dump_pdb( filename );
}

void
PeriodicBoxMover::setup_pose( Pose & pose, core::Real &mweight, core::Size &lattice_jump ) {
	using namespace core;
	using namespace core::conformation;
	using core::conformation::symmetry::SymDof;

	Pose pose_asu, pose_periodic;

	// get com of input pose
	Size nres = pose.size();
	mweight=0;
	Vector com(0,0,0);
	for ( Size i=1; i<= nres; ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		for ( Size j=1; j<= pose.residue(i).natoms(); ++j ) {
			std::string elt = pose.residue(i).atom_type(j).element();
			core::Real wt = 0.0;
			if ( elt=="H" ) wt = 1.008;
			else if ( elt=="C" ) wt = 12.011;
			else if ( elt=="N" ) wt = 14.007;
			else if ( elt=="O" ) wt = 15.999;
			else if ( elt=="P" ) wt = 30.974;
			else if ( elt=="S" ) wt = 32.060;
			else {
				TR << "ERROR: " << elt << std::endl;
			}

			com += wt*pose.residue(i).xyz(j);
			mweight+=wt;
		}
	}
	com /= mweight;

	// build box at specified density
	core::Size nmolecules = nmol_side_*nmol_side_*nmol_side_;
	core::Real sidelen = std::pow ( (nmolecules*mweight)/(0.60220*initial_density_) , 1.0/3.0 );

	// for now...
	if ( pose.size() == 1 ) {
		pose.apply_transform_Rx_plus_v( numeric::xyzMatrix<core::Real>::identity(),-com );

		TR << "Creating a box " << sidelen << "A on each side with " << nmolecules << " molecules." << std::endl;
		for ( int i=1; i<=nmol_side_; ++i ) {
			for ( int j=1; j<=nmol_side_; ++j ) {
				for ( int k=1; k<=nmol_side_; ++k ) {
					core::Vector center (
						sidelen * (2.0*i-1.0)/(2.0*nmol_side_),
						sidelen * (2.0*j-1.0)/(2.0*nmol_side_),
						sidelen * (2.0*k-1.0)/(2.0*nmol_side_) );


					// replace a position if running liqsim for two systems
					if ( i == (nmol_side_+1)/2 && j == (nmol_side_+1)/2 && k == (nmol_side_+1)/2
							&& central_molecule_pdb_ != "" ) {
						core::chemical::ResidueTypeSetCOP rsd_set =
							core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
						core::pose::Pose pose_central;
						core::import_pose::pose_from_file( pose_central, *rsd_set, central_molecule_pdb_ , core::import_pose::PDB_file);

						runtime_assert( pose_central.size() == 1 );
						core::conformation::Residue rsd_new = pose_central.residue(1);
						for ( Size z=1; z<= rsd_new.natoms(); ++z ) {
							rsd_new.set_xyz( z, rsd_new.xyz(z) + center );
						}
						pose_asu.append_residue_by_jump( rsd_new, 1 );

						central_resno_ = pose_asu.size();
						TR << "central resno: " << central_resno_ << std::endl;

					} else {
						core::conformation::Residue rsd_new = pose.residue(1);
						for ( Size z=1; z<= rsd_new.natoms(); ++z ) {
							rsd_new.set_xyz( z, rsd_new.xyz(z) + center );
						}
						pose_asu.append_residue_by_jump( rsd_new, 1 );
					}
				}
			}
		}
	} else {
		pose_asu = pose;
		io::CrystInfo ci = pose.pdb_info()->crystinfo();
		runtime_assert(ci.A()*ci.B()*ci.C() != 0 && ci.A()==ci.B() && ci.A()==ci.C());

		sidelen = ci.A();

		pose.apply_transform_Rx_plus_v( numeric::xyzMatrix<core::Real>::identity(),-com+0.5*sidelen );


		nmolecules = pose.size();
		nmol_side_ = std::pow(nmolecules, 1.0/3.0);
		std::cerr << nmolecules << " " << nmol_side_ << std::endl;
		runtime_assert( nmolecules == std::floor( nmol_side_*nmol_side_*nmol_side_ +0.5 ) );

		core::Real density = (nmolecules*mweight)/(0.60220*sidelen*sidelen*sidelen);

		TR << "Loading configuration " << sidelen << "A on each side (density = " << density << ") with " << nmolecules << " molecules." << std::endl;
	}

	pose_asu.dump_pdb( "asu.pdb" );
	//dump_ASU( pose_asu, "asu.pdb" );

	// make periodic box
	Vector Ax(1,0,0), Bx(0,1,0), Cx(0,0,1), Ay(0,1,0), By(0,0,1), Cy(1,0,0);

	ObjexxFCL::FArray3D<int> vrtX, vrtY, vrtZ;
	vrtX.dimension(3,1,1); vrtX=0;
	vrtY.dimension(3,3,1); vrtY=0;
	vrtZ.dimension(3,3,3); vrtZ=0;

	utility::vector1<Size> Ajumps, Bjumps, Cjumps, SUBjumps;

	origin_ = -0.5*sidelen;

	// 1 expand A (j==k==1)
	for ( int i=1; i<=3; ++i ) {
		core::Vector O ((i-1.5)*sidelen, -0.5*sidelen, -0.5*sidelen);
		ResidueOP vrt_x = make_vrt(O,Ax,Ay);
		ResidueOP vrt_y = make_vrt(O,Bx,By);
		ResidueOP vrt_z = make_vrt(O,Cx,Cy);

		if ( i==1 ) {
			pose_periodic.append_residue_by_bond( *vrt_x );    vrtX(1,1,1) = 1;
			pose_periodic.append_residue_by_jump( *vrt_y, 1);  vrtY(1,1,1) = 2;
			pose_periodic.append_residue_by_jump( *vrt_z, 2);  vrtZ(1,1,1) = 3;
		} else {
			pose_periodic.append_residue_by_jump( *vrt_x, vrtX(i-1,1,1));  vrtX(i,1,1) = pose_periodic.size();
			Ajumps.push_back(pose_periodic.fold_tree().num_jump());
			pose_periodic.append_residue_by_jump( *vrt_y, vrtX(i  ,1,1));  vrtY(i,1,1) = pose_periodic.size();
			pose_periodic.append_residue_by_jump( *vrt_z, vrtY(i  ,1,1));  vrtZ(i,1,1) = pose_periodic.size();
		}
	}

	// expand B (k==1)
	for ( int i=1; i<=3; ++i ) {
		for ( int j=2; j<=3; ++j ) {
			core::Vector O ((i-1.5)*sidelen, (j-1.5)*sidelen, -0.5*sidelen);
			ResidueOP vrt_y = make_vrt(O,Bx,By);
			ResidueOP vrt_z = make_vrt(O,Cx,Cy);

			pose_periodic.append_residue_by_jump( *vrt_y, vrtY(i,j-1,1)); vrtY(i,j,1) = pose_periodic.size();
			Bjumps.push_back(pose_periodic.fold_tree().num_jump());
			pose_periodic.append_residue_by_jump( *vrt_z, vrtY(i,j,1));   vrtZ(i,j,1) = pose_periodic.size();
		}
	}

	// expand C
	for ( int i=1; i<=3; ++i ) {
		for ( int j=1; j<=3; ++j ) {
			for ( int k=2; k<=3; ++k ) {
				core::Vector O ((i-1.5)*sidelen, (j-1.5)*sidelen, (k-1.5)*sidelen);
				ResidueOP vrt_z = make_vrt(O,Cx,Cy);

				pose_periodic.append_residue_by_jump( *vrt_z, vrtZ(i,j,k-1));  vrtZ(i,j,k) = pose_periodic.size();
				Cjumps.push_back(pose_periodic.fold_tree().num_jump());
			}
		}
	}

	Size nvirtuals = pose_periodic.size();

	// add subunits
	core::Size nres_monomer = pose_asu.size();
	core::Size nres_protein=0;
	for ( int n=1; n<=(int)nres_monomer; ++n ) {
		pose_periodic.insert_residue_by_jump( pose_asu.residue(n), nres_protein+1, vrtZ(2,2,2)+nres_protein ); ++nres_protein;
		//pose_periodic.append_residue_by_jump( pose_asu.residue(n), vrtZ(2,2,2) );
		SUBjumps.push_back(pose_periodic.fold_tree().num_jump());
	}

	for ( int i=1; i<=3; ++i ) {
		for ( int j=1; j<=3; ++j ) {
			for ( int k=1; k<=3; ++k ) {
				if ( i==2 && j==2 && k==2 ) continue;
				for ( int n=1; n<=(int)nres_monomer; ++n ) {
					pose_periodic.insert_residue_by_jump( pose_asu.residue(n), nres_protein+1, vrtZ(i,j,k)+nres_protein ); ++nres_protein;
					//pose_periodic.append_residue_by_jump( pose_asu.residue(n), vrtZ(i,j,k) );
				}
			}
		}
	}

	// make proper symmetric pose
	conformation::symmetry::SymmetryInfo symminfo;
	core::Size njumps_monomer = SUBjumps.size();
	for ( int i=2; i<=3*3*3; ++i ) {
		for ( Size j=1; j<= nres_monomer; ++j ) {
			symminfo.add_bb_clone ( j, (i-1)*nres_monomer+j );
			symminfo.add_chi_clone( j, (i-1)*nres_monomer+j );
		}
		for ( Size j=1; j<=njumps_monomer; ++j ) {
			symminfo.add_jump_clone( SUBjumps[j], (i-1)*njumps_monomer+SUBjumps[j], 0.0 );
		}
	}

	Size cellMaster=Ajumps[1];
	for ( Size i=2; i<=Ajumps.size(); ++i ) symminfo.add_jump_clone( cellMaster, Ajumps[i], 0.0 );
	for ( Size i=1; i<=Bjumps.size(); ++i ) symminfo.add_jump_clone( cellMaster, Bjumps[i], 0.0 );
	for ( Size i=1; i<=Cjumps.size(); ++i ) symminfo.add_jump_clone( cellMaster, Cjumps[i], 0.0 );

	std::map< Size, SymDof > symdofs;
	SymDof symdof_a;
	symdof_a.read( "x" );
	symdofs[ Ajumps[1] ] = symdof_a;
	symminfo.set_dofs( symdofs );

	// jump names
	TR << "Initializing " << pose_periodic.num_jump() << " jumps." << std::endl;
	for ( Size v=1; v<=pose.num_jump(); ++v ) symminfo.set_jump_name(v, "v_"+utility::to_string(v));
	symminfo.set_jump_name(Ajumps[1], "A");

	symminfo.num_virtuals( nvirtuals );
	symminfo.set_use_symmetry( true );
	symminfo.set_nres_subunit( nres_monomer );

	symminfo.set_flat_score_multiply( pose_periodic.size(), 0 );
	for ( Size i=1; i<= nres_protein; ++i ) {
		if ( symminfo.bb_is_independent( i ) ) symminfo.set_score_multiply( i, 2 );
		else symminfo.set_score_multiply( i, 1 );
	}
	symminfo.update_score_multiply_factor();

	pose::symmetry::make_symmetric_pose( pose_periodic, symminfo );

	for ( int i=1; i<=(int)SUBjumps.size(); ++i ) {
		core::kinematics::Jump j = pose_periodic.jump( SUBjumps[i] );
		pose_periodic.set_jump( SUBjumps[i], j );
	}

	// and we're done
	lattice_jump = Ajumps[1];
	pose = pose_periodic;
}

// function to check if pressure is being preserved
void
PeriodicBoxMover::check_virial_pressure( Pose & pose, core::Real L ) const
{
	using namespace core::optimization;
	using namespace ObjexxFCL::format;

	//TR << "check virial" << std::endl;

	core::kinematics::MoveMapOP movemap = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	movemap->set_chi( true ); movemap->set_bb( true );  movemap->set_jump( true );
	movemap->set( core::id::THETA, true ); movemap->set( core::id::D, true );

	// First get virial

	// just a way to get force
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	}

	CartesianMinimizerMap min_map;
	min_map.setup( pose, *movemap );

	sf_->score( pose );

	Multivec force( min_map.ndofs(), 0.0 ), xyz( min_map.ndofs(), 0.0 );
	//Multivec force( min_map.ndofs(), 0.0 ), xyz( min_map.ndofs(), 0.0 );
	min_map.copy_dofs_from_pose( pose, xyz ); // get xyz

	CartesianMultifunc f_ros( pose, min_map, *sf_, false, false );

	//TR << "to get force2: " << std::endl;
	f_ros.dfunc( xyz, force ); // get force

	/*
	core::Vector F1, F2;
	core::kinematics::DomainMap dommap = pose.energies().domain_map();
	for( core::Size imol = 1; imol <= pose.size(); ++imol ){
	for( core::Size iatm = 1; iatm <= pose.residue(imol).natoms(); ++iatm ){
	sf_->eval_npd_atom_derivative( core::id::AtomID( iatm, imol ),
	pose, dommap, F1, F2 );
	TR << F1[0] << " " << F1[1] << " " << F1[2] << std::endl;
	}
	}
	*/

	// just for check
	core::Vector forcesum( 0.0 );
	for ( core::Size i = 1; i <= force.size(); ++i ) forcesum[i%3] += force[i];

	TR  << forcesum[0] << " " << forcesum[1] << " " << forcesum[2] << std::endl;

	core::Real vir( 0.0 );
	for ( core::Size i = 1; i <= force.size(); ++i ) vir += force[i]*xyz[i];
	vir /= 3.0;

	// Then
	core::Real const V = L*L*L;
	//core::Size N( nmol_side_*nmol_side_*nmol_side_ );
	core::Size const N( force.size()/3 );
	core::Real const rho = (core::Real)(N)/V;

	core::Real const Boltzmann( 0.0019872065 ); // Boltmann constant in kcal
	// unit of P is kcal/molAng^3, rho: 1/Ang^3,
	core::Real Pvir = rho*temp_*Boltzmann + vir/V;
	Pvir /= 0.0001426; // convert unit to atm

	core::Real const dP(Pvir - P0_ );
	TR << "Pressure assigned / Pressure from virial / delta: ";
	TR << F(8,5,P0_) << " " << F(8,5,Pvir) << " " << F(8,5,dP) << std::endl;
}

void
PeriodicBoxMover::report_thermodynamics( Pose & pose, core::Size lattice_jump )
{
	using namespace ObjexxFCL::format;
	using namespace core::scoring;

	core::Size const nres( pose.size() );

	core::Real lattice = std::abs( pose.jump(lattice_jump).get_translation()[0] );
	core::Real volume = lattice*lattice*lattice;

	// 1. Density
	core::Real density( 0.0 );
	for ( core::Size ires = 1; ires <= nmol_side_*nmol_side_*nmol_side_; ++ires ) { // mol1
		if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;

		for ( Size j=1; j<= pose.residue(ires).natoms(); ++j ) {
			std::string elt = pose.residue(ires).atom_type(j).element();
			core::Real wt = 0.0;
			if ( elt=="H" ) wt = 1.008;
			else if ( elt=="C" ) wt = 12.011;
			else if ( elt=="N" ) wt = 14.007;
			else if ( elt=="O" ) wt = 15.999;
			else if ( elt=="P" ) wt = 30.974;
			else if ( elt=="S" ) wt = 32.060;
			else wt = 0.0;
			// TR << "ERROR: " << elt << std::endl;
			density += wt;
		}
	}
	density /= volume; // (g/mole)/Ang^3
	density /= 0.6023; // -> g/cm^3

	// 2. Intra/Inter energy decomposition
	// Assumption - just use two-body shortranged only
	EnergyMap weights( sf_->weights() );
	Energies const & pose_energies( pose.energies() );
	EnergyGraph const & egraph( pose_energies.energy_graph() );

	core::Real intraE( 0.0 ), interE( 0.0 ), centralE( 0.0 );
	for ( core::Size ires = 1; ires <= nres; ++ires ) { // mol1

		EnergyEdge const * edge(egraph.find_energy_edge(ires, ires));
		EnergyMap pair_energies;

		if ( edge != nullptr ) {
			EnergyMap unwt_pair_energies( edge->fill_energy_map() );
			pair_energies += unwt_pair_energies * weights;
		}
		EnergyMap onebody_energies = pose_energies.onebody_energies( ires );

		intraE += pair_energies.sum() + (onebody_energies * weights).sum();
		interE += pose_energies.residue_total_energy( ires );

		if ( ires == central_resno_ ) {
			centralE = pose_energies.residue_total_energy( ires );
		}
	}

	intraE /= (core::Real)(nres);
	interE /= (core::Real)(nres);

	// report
	TR << "THERMO: (Volume/Density/<IntraEnergy>/<InterEnergy>) " << F(8,3,volume);
	TR << " " << F(8,5,density);
	TR << " " << F(8,3,intraE) << " " << F(8,3,interE) << std::endl;

	// cache current data
	thermodynamic_data_.volume = volume;
	thermodynamic_data_.density = density;
	thermodynamic_data_.intraE = intraE;
	thermodynamic_data_.interE = interE;
	thermodynamic_data_.centralE = centralE;

	// 3. Radial distribution
	// skip if rg_atoms undefined
	if ( rg_atom1_ == "" || rg_atom2_ == "" ) return;

	core::Real const rgmin = 1.0, rgmax = 9.0, rgbin=0.2;
	core::Size const nbins = (core::Size)(( rgmax - rgmin )/rgbin) + 1;
	utility::vector1< core::Size > rgcounts( 0, nbins );

	// count RG
	for ( core::Size ires = 1; ires <= pose.size(); ++ires ) { // mol1
		if ( !pose.residue(ires).has( rg_atom1_ ) ) continue;
		core::Vector const &crd1 = pose.residue(ires).xyz( rg_atom1_ );

		for ( core::Size jres = 1; jres <= pose.size(); ++jres ) { // mol1
			if ( !pose.residue(ires).has( rg_atom1_ ) ) continue;
			core::Vector const &crd2 = pose.residue(jres).xyz( rg_atom2_ );
			core::Real const dist = crd1.distance(crd2);

			if ( dist > rgmin || dist < rgmax ) continue;

			core::Size ibin = (dist - rgmin)/rgbin + 1;
			rgcounts[ibin]++;
		}
	}

	TR << "RG: ";
	for ( core::Size ibin = 1; ibin <= nbins; ++ibin ) {
		TR << " " << rgcounts[ibin];
	}
	TR << std::endl;

}

void
PeriodicBoxMover::change_volume_move( Pose & pose, core::Size lattice_jump, bool &accept ) {
	using namespace ObjexxFCL::format;

	core::Real lattice = std::abs( pose.jump(lattice_jump).get_translation()[0] );
	core::Real volume = lattice*lattice*lattice;

	Pose pose_trial( pose );

	core::Real del_volume = vol_step_ * (1.0 - 2.0*numeric::random::rg().uniform() );
	core::Real new_lattice = std::pow( volume+del_volume, 1.0/3.0 );

	// 0.5: correction for symmetry setup
	core::Real const Eprv = 0.5*(sf_->score( pose_trial ) + ljcorrection_factor_/volume);

	core::Size const N( nmol_side_*nmol_side_*nmol_side_ );

	// change the lattice size
	core::kinematics::Jump j = pose_trial.jump( lattice_jump );
	j.set_translation( core::Vector(-new_lattice,0,0) );
	pose_trial.set_jump( lattice_jump, j );

	core::Real const volume2 = volume + del_volume;
	recenter_pose( pose_trial, lattice, new_lattice );

	// 0.5: correction for symmetry setup
	core::Real const Eaft = 0.5*(sf_->score( pose_trial ) + ljcorrection_factor_/volume2);

	// Boltzmann probability to change volume
	core::Real const beta = 1.0/((temp_+0.00100)*0.0019872065); // unit: (kcal/mole)-1, to prevent from being infinity at 0K
	core::Real const convert = 0.000014586; // 1 atm*Ang^3 = 0.00014586 kcal/mole

	// convert units for all the vars inside parenthesis to kcal/mole; multiplying beta makes unitless
	core::Real const arg1 = beta*((Eaft - Eprv) + P0_*del_volume*convert);
	core::Real const arg2 = -((core::Real)(N))*log( volume2/volume );

	core::Real probability = exp( -(arg1+arg2) ) ;
	if ( probability > 1.0 ) probability = 1.0;

	// Metropolis criteria
	if ( probability >= numeric::random::rg().uniform() ) {
		accept = true;
		pose = pose_trial;
	}

	TR << "Accept/Prob/dE/N/V/dV/PdV(kcal/mole)/arg1/arg2: " << accept << " " <<  F(8,5,probability);
	TR << " " << F(8,5,Eaft-Eprv) << " " << N << " " << F(10,5,volume2) << " " << F(10,5,del_volume);
	TR << " " << F(8,5,P0_*del_volume*convert) << " " << F(8,5,ljcorrection_factor_*(1.0/volume2 - 1.0/volume) )
		<< " " << F(10,5,arg1) << " " << F(10,5,arg2) << std::endl;

	// Check consistency by comparing to virial pressure
	//check_virial_pressure( pose, new_lattice );
}

void
PeriodicBoxMover::recenter_pose( Pose & pose, core::Real const lattice, core::Real const new_lattice) const {

	// now recenter
	for ( int mol_num=1; mol_num<=nmol_side_*nmol_side_*nmol_side_; ++mol_num ) {
		core::Vector com(0,0,0);
		for ( int i=1; i<=(int)pose.residue(mol_num).natoms(); ++i ) {
			com += pose.residue(mol_num).xyz(i);
		}
		com /= pose.residue(mol_num).natoms();

		//origin_ = new_lattice*0.5;
		com = com-origin_-0.5*new_lattice;
		core::Vector offset = com*new_lattice/lattice - com;

		for ( int i=1; i<=(int)pose.residue(mol_num).natoms(); ++i ) {
			pose.set_xyz( core::id::AtomID( i, mol_num ),  pose.residue(mol_num).xyz(i) + offset );
		}
	}
}

void
PeriodicBoxMover::perturb_molecule_move( Pose & pose, core::Size lattice_jump, bool &accept ) {
	core::Real lattice = std::abs( pose.jump(lattice_jump).get_translation()[0] );

	core::Real const Eprv = 0.5*sf_->score( pose ); // 0.5: correction for symmetry setup

	//1 pick a random molecule
	core::Size mol_num = numeric::random::rg().random_range(1,nmol_side_*nmol_side_*nmol_side_);
	core::conformation::Residue res_old (pose.residue(mol_num));

	if ( pose.residue(mol_num).nchi() == 0 || numeric::random::rg().uniform() <= probability_rigid_ ) {
		//2 random rotation & translation
		numeric::xyzMatrix<core::Real> R = random_rotation(  rot_step_ * numeric::random::rg().gaussian() );
		core::Vector T (
			trans_step_ * numeric::random::gaussian(),
			trans_step_ * numeric::random::gaussian(),
			trans_step_ * numeric::random::gaussian()
		);

		//3 apply
		core::Vector com(0,0,0);
		for ( int i=1; i<=(int)pose.residue(mol_num).natoms(); ++i ) {
			com += pose.residue(mol_num).xyz(i);
		}
		com /= pose.residue(mol_num).natoms();

		// boundry check
		for ( int i=0; i<3; ++i ) {
			if ( com[i]+T[i]-origin_ > 1.5*lattice ) T[i]-=lattice;
			if ( com[i]+T[i]-origin_ < 0.5*lattice ) T[i]+=lattice;
		}

		for ( core::Size i=1; i<=pose.residue(mol_num).natoms(); ++i ) {
			core::Vector Rb = pose.residue(mol_num).xyz(i);
			core::Vector Rc = R*( Rb - com) +T+com;
			pose.set_xyz( core::id::AtomID( i, mol_num ), R*( pose.residue(mol_num).xyz(i) - com) +T+com );
		}

	} else {
		// 2-2. torsional change
		core::Real dtor = tor_step_*(1.0 - 2.0*numeric::random::rg().uniform());
		core::Size ichi = numeric::random::rg().random_range( 1, pose.residue(mol_num).nchi() );
		core::Real newchi = pose.chi( ichi, mol_num ) + dtor;
		pose.set_chi( ichi, mol_num, newchi );
	}

	core::Real const Eaft = 0.5*sf_->score( pose ); // 0.5: correction for symmetry setup

	// Boltzmann probability
	core::Real kT( 0.008315*0.239*temp_ );
	core::Real const arg = (Eaft - Eprv)/kT;
	core::Real probability = exp( -(arg) ) ;
	if ( probability > 1.0 ) probability = 1.0;

	// Metropolis criteria
	if ( probability >= numeric::random::rg().uniform() ) {
		accept = true;
	} else {
		pose.replace_residue( mol_num, res_old, false );
	}
}

/*
void
PeriodicBoxMoverWaterReporter::report( PeriodicBoxMover & mover, int step, core::pose::Pose & pose, core::Real mweight, core::Size lattice_jump ) {
// calculate and write density
core::Real sidelen = std::abs( pose.jump(lattice_jump).get_translation()[0] );
core::Real box_vol = std::pow( sidelen, 3.0 );
core::Real nres = mover.nmol_side_*mover.nmol_side_*mover.nmol_side_;
core::Real box_den = (mweight*nres)/(box_vol*0.6022); // box volume in gm/cm3
volout_eq_ << step << " " << sidelen << " " << box_vol << " " << box_den << std::endl;


// calculate and write RDF data
core::Real cut = 10.0; // cutoff for RDF
core::Real bw = 0.1; // bin width for RDF
if (print_RDF_header_) {
RDF_eq_ << "#RDF cutoff distance = "<< cut <<"; bin width = "<< bw <<"\n";
RDF_eq_ << "#step: RDF bins\n";
print_RDF_header_ = false;
}

int nbins = (core::Size)std::ceil(cut/bw);
utility::vector1<int> RDF_data(nbins);

// create normalization vector for given bin width and cutoff and number density
core::Real nden = nres/std::abs((std::pow(pose.jump(lattice_jump).get_translation()[0],3.0)));
utility::vector1<core::Real> RDF_norm(std::ceil(cut/bw));
for (std::string::size_type ni = 1; ni<=RDF_norm.size(); ni++)
RDF_norm[ni] = nden*(4/3*M_PI*std::pow((bw*(ni-1))+bw, 3.0) - 4/3*M_PI*std::pow(bw*(ni-1), 3.0));

// find all interacting pairs and sort them into RDF bins
for (int ires=1; ires<=nres; ++ires) {
core::scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );
for (utility::graph::Graph::EdgeListConstIter iru = energy_graph.get_node(ires)->const_edge_list_begin(),
irue = energy_graph.get_node(ires)->const_edge_list_end(); iru != irue; ++iru) {
core::scoring::EnergyEdge const * edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
core::Size const e1( edge->get_first_node_ind() );
core::Size const e2( edge->get_second_node_ind() );
numeric::xyzVector< core::Real > p1(pose.residue(e1).xyz("O"));
numeric::xyzVector< core::Real > p2(pose.residue(e2).xyz("O"));
int bin = (int)std::floor( p1.distance(p2)/bw ) +1;
if (bin>=1 && bin<=nbins)
++RDF_data[bin]; // increment RDF bins
}
}
// normalize simulation data with ideal gas case
utility::vector1<core::Real> RDF_out(std::ceil(cut/bw));
for (std::string::size_type k=1; k <= RDF_data.size(); k++) {
RDF_out[k] = RDF_data[k]/nres/RDF_norm[k];
}

RDF_eq_ << step <<": ";
for (std::string::size_type rdfi=1; rdfi <= RDF_out.size(); rdfi++) {
RDF_eq_ << RDF_out[rdfi] <<" ";
}
RDF_eq_ << std::endl;

// calculate hydrogen bonds
core::scoring::hbonds::HBondSet hbond_set;
core::scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
Size nhbond_tot( 0 );
for ( Size ihb = 1; ihb <= Size(hbond_set.nhbonds()); ++ihb ) {
core::scoring::hbonds::HBond const & hb( hbond_set.hbond(ihb) );
if ((hb.don_res() <= nres) && (hb.acc_res() <= nres)) {
nhbond_tot += 2;
}
else {
nhbond_tot++;
}
}
hbond_eq_ << step << " "<<  hbond_set.nhbonds() <<" "<< nhbond_tot <<" " <<ObjexxFCL::format::F(6,3,(float)nhbond_tot/(float)nres)<< std::endl;
}
*/

void
PeriodicBoxMover::setup_LJcorrection( core::pose::Pose const &pose )
{
	using namespace ObjexxFCL::format;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// first count numbers
	std::map< core::Size, core::Size > atommap;
	utility::vector1< core::Size > Natoms( 0 );
	utility::vector1< core::Size > atmtypes( 0 );
	ljcorrection_factor_ = 0.0;

	// For now, assume that we're running only for single type molecule simulation
	for ( core::Size iatm = 1; iatm <= pose.residue( 1 ).natoms(); ++iatm ) {
		core::Size atmtype( pose.residue(1).atom_type_index(iatm) );
		if ( atommap.count( atmtype ) == 0 ) {
			Natoms.push_back( 0 );
			atmtypes.push_back( atmtype );
			atommap[ atmtype ] = Natoms.size();
		}
		Natoms[ atommap.at( atmtype ) ] ++;
	}

	//core::chemical::AtomTypeSet const atom_set = pose.residue( 1 ).atom_type_set();
	core::chemical::AtomTypeSetCOP atom_set =
		core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");

	TR << "LJ tail correction (energy times 1000 Ang^3)" << std::endl;
	TR << "Atom1 Atom2   Ni   Nj    Ecorr       arg1       arg2      eps      sig" <<  std::endl;

	core::Size Ntot( nmol_side_*nmol_side_*nmol_side_ );
	bool const score_Hatr( option[ score::fa_Hatr]() );

	for ( core::Size i = 1; i <= atmtypes.size(); ++i ) {
		core::Size const Ni = Natoms[i]*Ntot;
		core::Real const eps1 = (*atom_set)[atmtypes[i]].lj_wdepth();
		core::Real const sig1 = (*atom_set)[atmtypes[i]].lj_radius();

		bool const is_iH = (*atom_set)[atmtypes[i]].is_hydrogen();
		if ( !score_Hatr && is_iH ) continue;

		for ( core::Size j = 1; j <= atmtypes.size(); ++j ) {
			core::Size const Nj = Natoms[j]*Ntot;
			core::Real const eps2 = (*atom_set)[atmtypes[j]].lj_wdepth();
			core::Real const sig2 = (*atom_set)[atmtypes[j]].lj_radius();

			core::Real sig = sig1+sig2;
			core::Real eps = std::sqrt(eps1*eps2);

			bool const is_jH = (*atom_set)[atmtypes[j]].is_hydrogen();

			if ( !score_Hatr && is_jH ) continue;

			core::Real factor1(0.0), factor2(0.0);
			core::Real rc = option[ score::fa_max_dis ]();
			core::Real rc3 = rc*rc*rc;
			core::Real sig6 = std::pow( sig, 6 );

			// correction for max_dis to infinity
			factor1 = (sig6/(3.0*rc3*rc3) - 2.0) * (sig6/rc3);
			factor1 *= 4.0*Ni*Nj*3.141592*eps/3.0;

			// correction for max_dis-1.5 to max_dis (truncated part)
			// using apprx function of (r-max_dis+1.5)/1.5, integrating over max_dis-1.5 to max_dis
			core::Real rt = rc - 1.5;
			/*
			core::Real rt3 = rt*rt*rt;
			core::Real factor2 = (-rc/18.0+rc/3.0)*rc3 - (-rt/18.0+rc/3.0)*rt3;
			*/
			core::Real ra = 0.5*(rc + rt);
			factor2 = 4.0*3.141592*0.5*1.5*ra*ra; // ~20 at rc = 6.0
			factor2 *= Ni*Nj*(-2.0)*(sig6/(rc3*rc3))*eps;

			// correction factor for i,j pair
			// correction factor * 1/V gives correction energy
			ljcorrection_factor_ += factor1 + factor2;
			TR << std::setw(5) << (*atom_set)[atmtypes[i]].name() << " " << std::setw(5) << (*atom_set)[atmtypes[j]].name()
				<< " " << I(4,Ni) << " " << I(4,Nj) << " " << F(8,3,(factor1 + factor2)/1000.0)
				<< " " << F(10,3,factor1/1000.0) << " " << F(10,3,factor2/1000.0)
				<< " " << F(8,5,eps) << " " << F(8,3,sig) << std::endl;
		}
	}

	// this number has same aspect as total score; should be divided by 2 to consider symmetry
	TR << "Total: " << ljcorrection_factor_/1000.0 << " times 1000 Ang^3 " << std::endl;
}

void
PeriodicBoxMover::add_thermodynamic_data_to_silent( core::io::silent::SilentStructOP ss ) const
{
	ss->add_energy( "volume", thermodynamic_data_.volume );
	ss->add_energy( "density", thermodynamic_data_.density );
	ss->add_energy( "intraE", thermodynamic_data_.intraE );
	ss->add_energy( "interE", thermodynamic_data_.interE );
	ss->add_energy( "LJcorr", ljcorrection_factor_/thermodynamic_data_.volume );
	if ( central_molecule_pdb_ != "" ) {
		ss->add_energy( "centralE", thermodynamic_data_.centralE );
	}
}

void
PeriodicBoxMover::apply( Pose & pose ) {

	using namespace ObjexxFCL::format;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// build the water grid
	core::Real mweight;
	core::Size lattice_jump;

	setup_pose( pose, mweight, lattice_jump );
	//pose.dump_pdb( "setup.pdb" );
	//dump_ASU( pose, lattice_jump, "setup2.pdb" );

	sf_ = core::scoring::symmetry::symmetrize_scorefunction( *sf_ );
	core::io::silent::SilentFileData sfd, sfd2;

	if ( correct_LJtruncation_ ) setup_LJcorrection( pose );

	// burn in
	core::Size nvolume_change( 0 ), nvolume_change_accept( 0 );
	core::Size nperturb( 0 ), nperturb_accept( 0 );
	core::Real starttime, currtime;
	starttime = time(nullptr);

	// get initial thermodynamics properties

	for ( int i=istart_+1; i<=(int)nsteps_equilibrate_; ++i ) {
		if ( i%resize_vol_every_ == 0 ) {
			nvolume_change++;
			bool accept( false );
			change_volume_move(pose, lattice_jump, accept );
			if ( accept ) nvolume_change_accept++;
		} else {
			nperturb++;
			bool accept( false );
			perturb_molecule_move(pose, lattice_jump, accept);
			if ( accept ) nperturb_accept++;
		}

		// this is for scorefile
		currtime = time(nullptr);
		if ( report_every_>0 && i%report_every_ == 0 ) {
			core::Real dt_in_min = ( currtime - starttime ) / 60.0;
			TR << "[equilibrate cycle " << i << ", time(min) " << dt_in_min << "]" << std::endl;
			core::Real ratio_volume = (core::Real)( nvolume_change_accept )/(core::Real)( nvolume_change );
			core::Real nperturb_volume = (core::Real)( nperturb_accept )/(core::Real)( nperturb );
			TR << "change_volume    trials= " << I(6,nvolume_change) << ";  accepts= " << F(8,4,ratio_volume) << std::endl;
			TR << "perturb_molecule trials= " << I(6,nperturb) << ";  accepts= " << F(8,4,nperturb_volume) << std::endl;

			if ( report_scorefile_ != "" ) {
				core::io::silent::SilentStructOP ss =
					core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

				ss->fill_struct( pose );
				ss->add_energy( "sim_hours", ( currtime - starttime )/3600.0 );
				ss->set_decoy_tag( "equilibrate_"+utility::to_string(i) );
				add_thermodynamic_data_to_silent( ss );
				sfd2.write_silent_struct( *ss, report_scorefile_, true );
			}
		}

		if ( report_thermodynamics_>0 && ( i == 1 || i%report_thermodynamics_ == 0) ) {
			report_thermodynamics( pose, lattice_jump );
		}

		//if (water_reporter_ && report_water_>0 && i%report_water_==0)
		// water_reporter_->report( *this, i, pose, mweight, lattice_jump );

		// this is for structure
		if ( dump_every_>0 && i%dump_every_==0 ) {

			if ( report_silent_ == "pdb" ) {
				//pose.dump_pdb( "equilibrate_"+utility::to_string(i)+".pdb" );
				dump_ASU( pose, lattice_jump, "equilibrate_"+utility::to_string(i)+".pdb" );

			} else if ( report_silent_ != "" ) {
				core::io::silent::SilentStructOP ss =
					core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

				core::pose::Pose pose_asu;
				core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);
				core::Real lattice = std::abs( pose.jump(lattice_jump).get_translation()[0] );
				core::io::CrystInfo ci;
				ci.A(lattice); ci.B(lattice); ci.C(lattice);
				ci.alpha(90.0); ci.beta(90.0); ci.gamma(90.0);
				ci.spacegroup("P 1");

				pose_asu.conformation().delete_residue_slow( pose_asu.size() );
				pose_asu.pdb_info()->set_crystinfo(ci);

				ss->fill_struct( pose_asu );
				ss->add_energy( "sim_hours", ( currtime - starttime )/3600.0 );
				ss->set_decoy_tag( "equilibrate_"+utility::to_string(i) );
				add_thermodynamic_data_to_silent( ss );
				sfd.write_silent_struct( *ss, report_silent_, false );
			}
			// dump_ASU has a bug; lets replace it to dump_pose instead for now
			dump_ASU( pose, lattice_jump, option[ out::prefix ]()+"restart"+option[ out::suffix ]()+".pdb" );
			//pose.dump_pdb( option[ out::prefix ]()+"restart"+option[ out::suffix ]()+".pdb" );
		}

	}

	//if (water_reporter_) // reset water reporter for sim run after eq completes
	// water_reporter_ = PeriodicBoxMoverWaterReporterOP( new PeriodicBoxMoverWaterReporter("sim") );

	nvolume_change_accept = 0; nvolume_change = 0;
	nperturb = 0; nperturb_accept = 0;
	for ( int i=1; i<=(int)nsteps_sim_; ++i ) {
		if ( i%resize_vol_every_ == 0 ) {
			nvolume_change++;
			bool accept( false );
			change_volume_move(pose, lattice_jump, accept );
			if ( accept ) nvolume_change_accept++;
		} else {
			nperturb++;
			bool accept( false );
			perturb_molecule_move(pose, lattice_jump, accept);
			if ( accept ) nperturb_accept++;
		}

		currtime = time(nullptr);
		if ( report_every_>0 && i%report_every_ == 0 ) {
			core::Real dt_in_min = ( currtime - starttime ) / 60.0;
			TR << "[simulate cycle " << i << ", time(min) " << dt_in_min << "]" << std::endl;
			core::Real ratio_volume = (core::Real)( nvolume_change_accept )/(core::Real)( nvolume_change );
			core::Real nperturb_volume = (core::Real)( nperturb_accept )/(core::Real)( nperturb );
			TR << "change_volume    trials= " << I(6,nvolume_change) << ";  accepts= " << F(8,4,ratio_volume) << std::endl;
			TR << "perturb_molecule trials= " << I(6,nperturb) << ";  accepts= " << F(8,4,nperturb_volume) << std::endl;

			if ( report_scorefile_ != "" ) {
				core::io::silent::SilentStructOP ss =
					core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

				ss->fill_struct( pose );
				ss->add_energy( "sim_hours", ( currtime - starttime )/3600.0 );
				ss->set_decoy_tag( "simulate_"+utility::to_string(i) );
				add_thermodynamic_data_to_silent( ss );
				sfd2.write_silent_struct( *ss, report_scorefile_, true );
			}
		}

		if ( report_thermodynamics_>0 && i%report_thermodynamics_ == 0 ) {
			report_thermodynamics( pose, lattice_jump );
		}

		//if (water_reporter_ && report_water_>0 && i%report_water_==0)
		// water_reporter_->report( *this, i, pose, mweight, lattice_jump );

		if ( dump_every_>0 && i%dump_every_==0 ) {
			if ( report_silent_ != "" ) {
				core::io::silent::SilentStructOP ss =
					core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

				ss->fill_struct( pose );
				ss->add_energy( "sim_hours", ( currtime - starttime )/3600.0 );
				ss->set_decoy_tag( "simulate_"+utility::to_string(i) );
				add_thermodynamic_data_to_silent( ss );
				sfd.write_silent_struct( *ss, report_silent_, false );
			}
			dump_ASU( pose, lattice_jump, option[ out::prefix ]()+"restart"+option[ out::suffix ]()+".pdb" );
		}
	}
}

void
PeriodicBoxMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &  ) {
	using namespace core::conformation;
	using namespace core::pack::task;

	sf_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );


	if ( tag->hasOption("nmol_side") ) {
		nmol_side_ = tag->getOption<core::Real>("nmol_side");
	}
	if ( tag->hasOption("initial_density") ) {
		initial_density_ = tag->getOption<core::Real>("initial_density");
	}
	if ( tag->hasOption("temp") ) {
		temp_ = tag->getOption<core::Real>("temp");
	}
	if ( tag->hasOption("vol_step") ) {
		vol_step_ = tag->getOption<core::Real>("vol_step");
	}
	if ( tag->hasOption("trans_step") ) {
		trans_step_ = tag->getOption<core::Real>("trans_step");
	}
	if ( tag->hasOption("rot_step") ) {
		rot_step_ = tag->getOption<core::Real>("rot_step");
	}
	if ( tag->hasOption("tor_step") ) {
		tor_step_ = tag->getOption<core::Real>("tor_step");
	}
	if ( tag->hasOption("probability_rigid") ) {
		probability_rigid_ = tag->getOption<core::Real>("probability_rigid");
	}
	if ( tag->hasOption("resize_vol_every") ) {
		resize_vol_every_ = tag->getOption<core::Real>("resize_vol_every");
	}
	if ( tag->hasOption("report_every") ) {
		report_every_ = tag->getOption<core::Real>("report_every");
	}
	if ( tag->hasOption("report_silent") ) {
		report_silent_ = tag->getOption<std::string>("report_silent");
	}
	if ( tag->hasOption("report_scorefile") ) {
		report_scorefile_ = tag->getOption<std::string>("report_scorefile");
	}
	if ( tag->hasOption("dump_every") ) {
		dump_every_ = tag->getOption<core::Real>("dump_every");
	}
	if ( tag->hasOption("nsteps_equilibrate") ) {
		nsteps_equilibrate_ = tag->getOption<core::Real>("nsteps_equilibrate");
	}
	if ( tag->hasOption("nsteps_sim") ) {
		nsteps_sim_ = tag->getOption<core::Real>("nsteps_sim");
	}
	if ( tag->hasOption("istart") ) {
		istart_ = tag->getOption<core::Real>("istart");
	}
	if ( tag->hasOption("central_molecule_pdb") ) {
		central_molecule_pdb_ = tag->getOption<std::string>("central_molecule_pdb");
	}

	// correction factor for LJ truncation
	if ( tag->hasOption("correct_LJtruncation") ) {
		correct_LJtruncation_ = tag->getOption< bool>("correct_LJtruncation");
	}

	// added for thermodynamics stuffs
	if ( tag->hasOption("report_thermodynamics") ) {
		report_thermodynamics_ = tag->getOption< core::Real>("report_thermodynamics");
	}
	if ( tag->hasOption("pressure") ) {
		P0_ = tag->getOption<core::Real>("pressure");
	}
	if ( tag->hasOption("rg_atom1") ) {
		rg_atom1_ = tag->getOption<std::string>("rg_atom1");
	}
	if ( tag->hasOption("rg_atom2") ) {
		rg_atom2_ = tag->getOption<std::string>("rg_atom1");
	}

	// reporters for water
	/*
	if (tag->hasOption("report_water")) {
	report_water_ = tag->getOption<core::Real>("report_water");
	water_reporter_ = PeriodicBoxMoverWaterReporterOP( new PeriodicBoxMoverWaterReporter("eq") );
	}
	*/
}

} // moves
} // protocols
