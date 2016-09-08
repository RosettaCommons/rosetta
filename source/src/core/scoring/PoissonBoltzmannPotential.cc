// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// Compute electrostatic.
///
/// Options:
///
/// vector1<Integer> pb_potential::charged_chain
///    The chain numbers (>=1) of which charge is non-zero.
///    The electrostatic will be computed for the atoms in these chains.
///
/// @file   core/scoring/PoissonBoltzmannPotential.cc
/// @brief  PoissonBoltzmann potential class implementation
/// @author Yifan Song (yfsong@uw.edu)

// basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pb_potential.OptionKeys.gen.hh>

// Unit Headers
#include <core/scoring/PoissonBoltzmannPotential.hh>
#include <core/scoring/APBSWrapper.hh>

// Package Headers
#include <core/scoring/electron_density/util.hh>
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Numeric Headers

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// utility
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <vector>
#include <ctime>

// option key includes


#include <fstream>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.PoissonBoltzmannPotential" );

static std::string const chr_chains("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcba");

namespace core {
namespace scoring {

typedef PoissonBoltzmannPotential PB;

const std::string PB::APBS_CONFIG_EXT = ".in";
const std::string PB::APBS_PQR_EXT = ".pqr";
const std::string PB::APBS_DX_EXT = ".dx";
const std::string PB::DEFAULT_APBS_PATH = "apbs";

// @brief Auto-generated virtual destructor
PB::~PoissonBoltzmannPotential() = default;

PB::PoissonBoltzmannPotential()
:config_filename_("Unknown.in"),
pqr_filename_("Unknown.pqr"),
dx_filename_("Unknown.dx"),
apbs_path_(DEFAULT_APBS_PATH),
calcenergy_(false)
{
}

core::Real
PB::get_potential(ObjexxFCL::FArray3D< core::Real > const & potential,
	numeric::xyzVector<core::Real> const & cartX) const {
	if ( out_of_bounds(cartX) ) return 0.0;

	numeric::xyzVector<core::Real> idxX;
	cart2idx(cartX, idxX);
	return core::scoring::electron_density::interp_linear(potential, idxX);
}
void
PB::eval_PB_energy_residue(
	core::conformation::Residue const & rsd,
	Real & PB_energy_residue,
	Real & PB_energy_backbone,
	Real & PB_energy_sidechain,
	Real const & PB_burial_weight
) const {
	PB_energy_residue = 0.0;
	PB_energy_backbone = 0.0;
	PB_energy_sidechain = 0.0;

	// Compute potential. Linear sum.
	for ( Size iatom=1; iatom<=rsd.natoms(); ++iatom ) {
		if ( rsd.atomic_charge(iatom) > -1e-6 && rsd.atomic_charge(iatom) < 1e-6 ) continue;
		core::Real iatom_potential = get_potential(potential_, rsd.xyz(iatom));

		core::Real atom_energy = rsd.atomic_charge(iatom)*iatom_potential * PB_burial_weight;

		if ( rsd.atom_is_backbone(iatom) ) {
			PB_energy_backbone += atom_energy;
		} else {
			PB_energy_sidechain += atom_energy;
		}
		PB_energy_residue += atom_energy;
	}
}

//============================================================================
#ifdef LINK_APBS_LIB   // APBS libraries are linked

void
PB::solve_pb(
	core::pose::Pose const & pose,
	std::string const & tag,
	std::map<std::string, bool> const & charged_residues
) {
	using namespace std;
	time_t begin;
	time(&begin);

	config_filename_ = tag + APBS_CONFIG_EXT;
	pqr_filename_ = tag + APBS_PQR_EXT;
	dx_filename_ = tag + APBS_DX_EXT;

	int apbs_dbg = basic::options::option[ basic::options::OptionKeys::pb_potential::apbs_debug ];
	calcenergy_ = basic::options::option[ basic::options::OptionKeys::pb_potential::calcenergy ];
	APBSWrapper apbs(pose, charged_residues, apbs_dbg, calcenergy_);

	APBSResultCOP result = 0;

	result = apbs.exec();

	if ( result == 0 ) {
		TR.Error << "APBS failed!  Terminating the program..." << std::endl;
		TR.flush();
		runtime_assert(false);
	}

	TR.Debug << "Solved PB. Loading potential..." << std::endl;
	const double * meta = result->grid_meta;
	const double * data = &(result->grid_data[0][0]);
	load_potential(meta, data);

	time_t end;
	time(&end);
	TR << "PB took " << end-begin << " seconds" << std::endl;
}

void
PB::load_potential(
	const double grid_meta[],
	const double pot[]
) {
	int nx = grid_meta[1];
	int ny = grid_meta[2];
	int nz = grid_meta[3];
	n_grid_[0] = nx;
	n_grid_[1] = ny;
	n_grid_[2] = nz;
	grid_spacing_[0] = grid_meta[4];
	grid_spacing_[1] = grid_meta[5];
	grid_spacing_[2] = grid_meta[6];
	lower_bound_[0] = grid_meta[10];
	lower_bound_[1] = grid_meta[11];
	lower_bound_[2] = grid_meta[12];

	i2c_ = numeric::xyzMatrix <core::Real>::rows(
		grid_spacing_[0], 0., 0.,
		0., grid_spacing_[1], 0.,
		0., 0., grid_spacing_[2] );
	c2i_ = numeric::xyzMatrix <core::Real>::rows(
		1./grid_spacing_[0], 0., 0.,
		0., 1./grid_spacing_[1], 0.,
		0., 0., 1./grid_spacing_[2] );

	potential_.dimension(nx, ny, nz);

	double cap =  basic::options::option[ basic::options::OptionKeys::pb_potential::potential_cap ];

	int icol=0;
	int u;
	int lines = 0;
	//using namespace ObjexxFCL::format;
	//std::ofstream ofs(dx_filename_.c_str());
	//ofs << "object 1 class gridpositions counts " << nx << " " << ny << " " << nz << std::endl;
	//ofs << "origin " << lower_bound_[0] << " " << lower_bound_[1] << " " << lower_bound_[2] << std::endl;
	//ofs << "object 3 array type double rank 0 items " << nx*ny*nz << " data follows" << std::endl;
	for (int i=1; i<=nx; i++) {
		for (int j=1; j<=ny; j++) {
			for (int k=1; k<=nz; k++) {
				u = (k-1)*(nx)*(ny)+(j-1)*(nx)+(i-1);
				//ofs << E(12,6,pot[u]) << " ";
				icol++;
				if (icol == 3) {
				    icol = 0;
					lines++;
					//ofs << std::endl;
				}
				if ( pot[u] > cap ) {
					potential_(i,j,k) = cap;
				}
				else if ( pot[u] < -cap ) {
					potential_(i,j,k) = -cap;
				}
				else{
					potential_(i,j,k) = pot[u];
				}
			}
		}
	}
	//ofs.close();
	TR.Debug << "PB potential is successfully loaded." << std::endl;

	idx2cart(n_grid_, upper_bound_);
	TR << "Convertion of PB potential to Certesian coordinates is completed" << std::endl;
}

#else  // APBS libraries are not linked.  Use system call.

void
PB::solve_pb( core::pose::Pose const & pose,
	std::string const & tag,
	std::map<std::string, bool> const &charged_residues )
{
#ifndef __native_client__
	using namespace std;
	time_t begin;
	time(&begin);

	if ( basic::options::option[basic::options::OptionKeys::pb_potential::apbs_path].user() ) {
		apbs_path_ = basic::options::option[basic::options::OptionKeys::pb_potential::apbs_path];
	}

	calcenergy_ = basic::options::option[ basic::options::OptionKeys::pb_potential::calcenergy ];

	// Generate filenames based on the given tag.
	config_filename_ = tag + APBS_CONFIG_EXT;
	pqr_filename_ = tag + APBS_PQR_EXT;
	dx_filename_ = tag + APBS_DX_EXT;

	write_pqr(pose, charged_residues );
	write_config(pose);
	std::string command_line(apbs_path_ + " " + config_filename_);
	if ( system(command_line.c_str()) == -1 ) {
		TR.Error << "Shell command failed to run!" << std::endl;
	}

	// Check if APBS succeeded.  If not, get out.
	std::ifstream dxstream(dx_filename_.c_str());
	if ( ! dxstream.good() ) {
		TR << "APBS failed to generate the result file.  Terminating the program." << std::endl;
		TR.flush();
		runtime_assert(false);
	}
	dxstream.close();

	// load the result
	load_APBS_potential();

	time_t end;
	time(&end);
	TR << "PB took " << end-begin << " seconds" << std::endl;
#endif

}
void
PB::load_APBS_potential()
{
	utility::io::izstream p_stream( dx_filename_ );
	runtime_assert(p_stream);

	std::string line;
	while ( p_stream ) {
		// # Comments
		char first_char = p_stream.peek();
		if ( first_char == '#' ) {
			p_stream.getline( line );
			continue;
		}

		//object 1 class gridpositions counts nx ny nz
		std::string buff;
		for ( Size i=0; i<5; ++i ) p_stream >> buff;
		p_stream >> n_grid_[0] >> n_grid_[1] >> n_grid_[2];

		//origin xmin ymin zmin
		p_stream >> buff;
		p_stream >> lower_bound_[0];
		p_stream >> lower_bound_[1];
		p_stream >> lower_bound_[2];

		//delta hx 0.0 0.0found
		//delta 0.0 hy 0.0
		//delta 0.0 0.0 hz

		p_stream >> buff >> grid_spacing_[0] >> buff >> buff;
		p_stream >> buff >> buff >> grid_spacing_[1] >> buff;
		p_stream >> buff >> buff >> buff >> grid_spacing_[2];
		i2c_ = numeric::xyzMatrix <core::Real>::rows(
			grid_spacing_[0],0.,0.,
			0.,grid_spacing_[1],0.,
			0.,0.,grid_spacing_[2]);
		c2i_ = numeric::xyzMatrix <core::Real>::rows(
			1./grid_spacing_[0],0.,0.,
			0.,1./grid_spacing_[1],0.,
			0.,0.,1./grid_spacing_[2]);

		//object 2 class gridconnections counts nx ny nz
		//object 3 class array type double rank 0 items n data follows
		p_stream.getline( line );
		p_stream.getline( line );
		p_stream.getline( line );

		//u(0,0,0) u(0,0,1) u(0,0,2)
		potential_.dimension(n_grid_[0],n_grid_[1],n_grid_[2]);

		for ( int i=1; i<=potential_.u1(); i++ ) {
			for ( int j=1; j<=potential_.u2(); j++ ) {
				for ( int k=1; k<=potential_.u3(); k++ ) {
					p_stream >> buff;
					potential_(i,j,k) = atof(buff.c_str());
					if ( potential_(i,j,k) > basic::options::option[ basic::options::OptionKeys::pb_potential::potential_cap ] ) {
						potential_(i,j,k) = basic::options::option[ basic::options::OptionKeys::pb_potential::potential_cap ];
					}
					if ( potential_(i,j,k) < -basic::options::option[ basic::options::OptionKeys::pb_potential::potential_cap ] ) {
						potential_(i,j,k) = -basic::options::option[ basic::options::OptionKeys::pb_potential::potential_cap ];
					}
				}
			}
		}
		break;
	}

	TR << "PB potential is successfully loaded from: " << dx_filename_ << std::endl;

	idx2cart(n_grid_, upper_bound_);
	TR << "Convertion of PB potential to Certesian coordinates is completed" << std::endl;

	//attribute "dep" string "positions"
	//object "regular positions regular connections" class field
	//component "positions" value 1
	//component "connections" value 2
	//component "data" value 3
}

#endif // LINK_APBS_LIB
//==============================================================================
void
PB::write_pqr(
	core::pose::Pose const & pose,
	std::map<std::string, bool> const & is_residue_charged_by_name_
) const {
	// Generate .pqr
	std::ofstream pqr_ostr(pqr_filename_.c_str());

	using namespace basic::options::OptionKeys;
	utility::vector1<Size> charged_chains;
	charged_chains.push_back(1);

	Size const nres( pose.size() );
	Size number( 0 );

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		// TODO: Sachko 11/19/2012
		// This search result should be cached.
		bool residue_charged = const_cast<std::map<std::string,bool>&>(is_residue_charged_by_name_)[rsd.type().name()];
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

			//skip outputing virtual atom unless specified.
			//fixed so that the last atom in atom type set can be something other than a virtual atom --steven combs
			if ( rsd.atom_type(j).is_virtual() ) continue;

			++number;

			char const chain( chr_chains[ (rsd.chain()-1)%chr_chains.size() ] );
			if ( residue_charged ) {
				using namespace ObjexxFCL::format;
				pqr_ostr << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
					rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
					F(8,3,atom.xyz()(1)) <<
					F(8,3,atom.xyz()(2)) <<
					F(8,3,atom.xyz()(3)) <<
					F(8,3,pose.residue_type(i).atom(j).charge()) <<
					F(8,3,rsd.atom_type(j).lj_radius()) << '\n';
			} else {
				using namespace ObjexxFCL::format;
				pqr_ostr << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
					rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
					F(8,3,atom.xyz()(1)) <<
					F(8,3,atom.xyz()(2)) <<
					F(8,3,atom.xyz()(3)) <<
					F(8,3,0.0) <<
					F(8,3,rsd.atom_type(j).lj_radius()) << '\n';
			}
		}
	}
	pqr_ostr.close();

	TR << pqr_filename_ << " is successfully written." << std::endl;
}

///
/// Write out the configurati
///
void
PB::write_config(
	core::pose::Pose const & pose
) const {

	// Generate .in
	std::ofstream config_ostr(config_filename_.c_str());

	numeric::xyzVector <core::Real> min_r(9999,9999,9999);
	numeric::xyzVector <core::Real> max_r(-9999,-9999,-9999);
	// Find the min & max coords within the moleculer system to define the grid.
	for ( core::Size ires=1; ires<=pose.size(); ++ires ) {
		for ( core::Size iatom=1; iatom<=pose.residue(ires).natoms(); ++iatom ) {
			for ( core::Size i=0; i<3; ++i ) {
				if ( pose.residue(ires).xyz(iatom)[i] < min_r[i] ) {
					min_r[i] = pose.residue(ires).xyz(iatom)[i];
				}
				if ( pose.residue(ires).xyz(iatom)[i] > max_r[i] ) {
					max_r[i] = pose.residue(ires).xyz(iatom)[i];
				}
			}
		}
	}

	numeric::xyzVector <core::Real> length = max_r - min_r;       // grid widths
	numeric::xyzVector <core::Real> center = (min_r + max_r)/2.;  // grid center coords
	//APBS psize.py paramters:
	core::Real cfac(1.7);
	core::Real fadd(20);
	core::Real space(0.5);
	numeric::xyzVector <core::Real> dimension_fine = length + numeric::xyzVector <core::Real> (fadd,fadd,fadd);
	numeric::xyzVector <core::Real> dimension_coarse = length * cfac;

	numeric::xyzVector <core::Real> tmpv = dimension_fine / space + numeric::xyzVector <core::Real> (1.,1.,1.); // bazzoli: to avoid warning by my compiler
	numeric::xyzVector <core::Size> n_grid;
	n_grid.x() = (core::Size)tmpv.x();
	n_grid.y() = (core::Size)tmpv.y();
	n_grid.z() = (core::Size)tmpv.z();

	using namespace ObjexxFCL::format;
	config_ostr << "#" << std::endl;
	config_ostr << "# Note that most of the comments here were taken from sample" << std::endl;
	config_ostr << "# input files that came with APBS.  You can find APBS at" << std::endl;
	config_ostr << "# http://agave.wustl.edu/apbs/" << std::endl;
	config_ostr << "# Note that APBS is BSD & MIT'd code." << std::endl;
	config_ostr << "#" << std::endl;
	config_ostr << "read" << std::endl;
	config_ostr << "mol pqr "<< pqr_filename_ <<"       # read molecule 1" << std::endl;
	config_ostr << "end" << std::endl;
	config_ostr << "elec" << std::endl;
	config_ostr << "mg-auto" << std::endl;
	config_ostr << "dime   " << I(5,n_grid.x()) << " "<< I(5,n_grid.y()) << " "<< I(5,n_grid.z()) << "  # number of find grid points" << std::endl;
	config_ostr << "# calculated by psize.py" << std::endl;
	config_ostr << "cglen   " << F(11,6,dimension_coarse.x()) << " " << F(11,6,dimension_coarse.y()) << " " << F(11,6,dimension_coarse.z()) << " # coarse mesh lengths (A)" << std::endl;
	config_ostr << "fglen   " << F(11,6,dimension_fine.x()) << " " << F(11,6,dimension_fine.y()) << " " << F(11,6,dimension_fine.z()) << " # fine mesh lengths (A)" << std::endl;
	config_ostr << "# calculated by psize.py" << std::endl;
	config_ostr << "cgcent " << F(11,6,center.x()) << " " << F(11,6,center.y()) << " " << F(11,6,center.z()) << "  # (could also give (x,y,z) form psize.py) #known center" << std::endl;
	config_ostr << "fgcent " << F(11,6,center.x()) << " " << F(11,6,center.y()) << " " << F(11,6,center.z()) << "  # (could also give (x,y,z) form psize.py) #known center" << std::endl;
	config_ostr << "npbe               # solve the full nonlinear PBE with npbe" << std::endl;
	config_ostr << "#lpbe            # solve the linear PBE with lpbe" << std::endl;
	config_ostr << "bcfl sdh          # Boundary condition flag" << std::endl;
	config_ostr << "#  0 => Zero" << std::endl;
	config_ostr << "#  1 => Single DH sphere" << std::endl;
	config_ostr << "#  2 => Multiple DH spheres" << std::endl;
	config_ostr << "#  4 => Focusing" << std::endl;
	config_ostr << "#" << std::endl;
	config_ostr << "#ion 1 0.000 2.0 # Counterion declaration:" << std::endl;
	config_ostr << "ion  1 0.150000 2.000000     # Counterion declaration:" << std::endl;
	config_ostr << "ion -1 0.150000 2.000000     # ion <charge> <conc (M)> <radius>" << std::endl;
	config_ostr << "ion  2 0.000000 2.000000     # ion <charge> <conc (M)> <radius>" << std::endl;
	config_ostr << "ion -2 0.000000 2.000000     # ion <charge> <conc (M)> <radius>" << std::endl;
	config_ostr << "pdie 4.000000          # Solute dielectric" << std::endl;
	config_ostr << "sdie 80.000000          # Solvent dielectric" << std::endl;
	config_ostr << "chgm spl2          # Charge disc method" << std::endl;
	config_ostr << "# 0 is linear splines" << std::endl;
	config_ostr << "# 1 is cubic b-splines" << std::endl;
	config_ostr << "mol 1            # which molecule to use" << std::endl;
	config_ostr << "srfm smol        # Surface calculation method" << std::endl;
	config_ostr << "#  0 => Mol surface for epsilon;" << std::endl;
	config_ostr << "#       inflated VdW for kappa; no" << std::endl;
	config_ostr << "#       smoothing" << std::endl;
	config_ostr << "#  1 => As 0 with harmoic average" << std::endl;
	config_ostr << "#       smoothing" << std::endl;
	config_ostr << "#  2 => Cubic spline " << std::endl;
	config_ostr << "srad 1.400000          # Solvent radius (1.4 for water)" << std::endl;
	config_ostr << "swin 0.3         # Surface cubic spline window .. default 0.3" << std::endl;
	config_ostr << "temp 310.000000          # System temperature (298.15 default)" << std::endl;
	config_ostr << "sdens 10.000000         # Specify the number of grid points per square-angstrom to use in Vacc object. Ignored when srad is 0.0 (see srad) or srfm is spl2 (see srfm). There is a direct correlation between the value used for the Vacc sphere density, the accuracy of the Vacc object, and the APBS calculation time. APBS default value is 10.0." << std::endl;
	config_ostr << "gamma 0.105      # Surface tension parameter for apolar forces (in kJ/mol/A^2)" << std::endl;
	config_ostr << "# only used for force calculations, so we don't care, but" << std::endl;
	config_ostr << "# it's always required, and 0.105 is the default" << std::endl;
	config_ostr << "calcenergy " << (calcenergy_? "yes" : "no") << "    # Energy I/O to stdout" << std::endl;
	config_ostr << "#  0 => don't write out energy" << std::endl;
	config_ostr << "#  1 => write out total energy" << std::endl;
	config_ostr << "#  2 => write out total energy and all" << std::endl;
	config_ostr << "#       components" << std::endl;
	config_ostr << "calcforce no     # Atomic forces I/O (to stdout)" << std::endl;
	config_ostr << "#  0 => don't write out forces" << std::endl;
	config_ostr << "#  1 => write out net forces on molecule" << std::endl;
	config_ostr << "#  2 => write out atom-level forces" << std::endl;
	config_ostr << "write pot dx " << dx_filename_.substr(0,dx_filename_.size()-APBS_DX_EXT.size()) << "  # What to write .. this says write the potential in dx" << std::endl;
	config_ostr << "# format to a file." << std::endl;
	config_ostr << "end" << std::endl;
	config_ostr << "quit" << std::endl;

	config_ostr.close();

	TR << config_filename_ << " is successfully written." << std::endl;
}

}
}
