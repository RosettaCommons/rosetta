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
/// @detailed
/// @author Yifan Song

// basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// core
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/types.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/PoissonBoltzmannPotential.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

// protocols
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

// utility
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <ostream>

//Option( 'PB_charged_chains', 'String', desc="Poisson Boltzmann potential is generated with the charges of given chains" , default="A"),
//Option( 'PB_boundary_chains', 'String', desc="Poisson Boltzmann potential is generated with the boundary of given chains" , default="AB"),

namespace PB_potential {
	basic::options::FileOptionKey potential_file("PB_potential:potential_file");
	basic::options::StringOptionKey apbs_exe("PB_potential:apbs_exe");
	basic::options::IntegerVectorOptionKey no_charge_on_chain("PB_potential:no_charge_on_chain");
}

static basic::Tracer TR("pilot.yfsong.PB_potential");


class PBPotentialMover : public protocols::moves::Mover {
public:

void write_APBS_config(core::pose::Pose & pose, std::ostream & config_file, std::string pqr_fn) {
	std::string tag = core::pose::extract_tag_from_pose( pose );
	numeric::xyzVector <core::Real> min_r(9999,9999,9999);
	numeric::xyzVector <core::Real> max_r(-9999,-9999,-9999);
	for (core::Size ires=1; ires<=pose.total_residue(); ++ires) {
		for (core::Size iatom=1; iatom<=pose.residue(ires).natoms(); ++iatom) {
			for (core::Size i=0; i<3; ++i) {
				if (pose.residue(ires).xyz(iatom)[i] < min_r[i]) {
					min_r[i] = pose.residue(ires).xyz(iatom)[i];
				}
				if (pose.residue(ires).xyz(iatom)[i] > max_r[i]) {
					max_r[i] = pose.residue(ires).xyz(iatom)[i];
				}
			}
		}
	}

	//min_r -= numeric::xyzVector <core::Real> (20,20,20);
	//max_r += numeric::xyzVector <core::Real> (20,20,20);
	numeric::xyzVector <core::Real> length = max_r - min_r;
	numeric::xyzVector <core::Real> center = (min_r + max_r)/2.;
	//APBS psize.py paramters:
	core::Real cfac(1.7);
	core::Real fadd(20);
	core::Real space(0.5);
	numeric::xyzVector <core::Real> dimension_fine = length + numeric::xyzVector <core::Real> (fadd,fadd,fadd);
	numeric::xyzVector <core::Real> dimension_coarse = length * cfac;
	numeric::xyzVector <core::Size> n_grid = static_cast< numeric::xyzVector <core::Size> > (dimension_fine / space + numeric::xyzVector <core::Real> (1.,1.,1.));

	// Shrink fine dimensions if they exceed coarse dimensions
	for (core::Size i = 0; i < 3 ; i++)
	{
		if(dimension_fine[i] > dimension_coarse[i])
		{
			dimension_fine[i] = dimension_coarse[i];
		}
	}

	using namespace ObjexxFCL::format;
	config_file << "#" << std::endl;
	config_file << "# Note that most of the comments here were taken from sample" << std::endl;
	config_file << "# input files that came with APBS.  You can find APBS at" << std::endl;
	config_file << "# http://agave.wustl.edu/apbs/" << std::endl;
	config_file << "# Note that APBS is GPL'd code." << std::endl;
	config_file << "#" << std::endl;
	config_file << "read" << std::endl;
	config_file << "mol pqr "<< pqr_fn <<"       # read molecule 1" << std::endl;
	config_file << "end" << std::endl;
	config_file << "elec" << std::endl;
	config_file << "mg-auto" << std::endl;
	config_file << "dime   " << I(5,n_grid.x()) << " "<< I(5,n_grid.y()) << " "<< I(5,n_grid.z()) << "  # number of find grid points" << std::endl;
	config_file << "# calculated by psize.py" << std::endl;
	config_file << "cglen   " << F(11,6,dimension_coarse.x()) << F(11,6,dimension_coarse.y()) << F(11,6,dimension_coarse.z()) << " # coarse mesh lengths (A)" << std::endl;
	config_file << "fglen   " << F(11,6,dimension_fine.x()) << F(11,6,dimension_fine.y()) << F(11,6,dimension_fine.z()) << " # fine mesh lengths (A)" << std::endl;
	config_file << "# calculated by psize.py" << std::endl;
	config_file << "cgcent " << F(11,6,center.x()) << F(11,6,center.y()) << F(11,6,center.z()) << "  # (could also give (x,y,z) form psize.py) #known center" << std::endl;
	config_file << "fgcent " << F(11,6,center.x()) << F(11,6,center.y()) << F(11,6,center.z()) << "  # (could also give (x,y,z) form psize.py) #known center" << std::endl;
	config_file << "npbe               # solve the full nonlinear PBE with npbe" << std::endl;
	config_file << "#lpbe            # solve the linear PBE with lpbe" << std::endl;
	config_file << "bcfl sdh          # Boundary condition flag" << std::endl;
	config_file << "#  0 => Zero" << std::endl;
	config_file << "#  1 => Single DH sphere" << std::endl;
	config_file << "#  2 => Multiple DH spheres" << std::endl;
	config_file << "#  4 => Focusing" << std::endl;
	config_file << "#" << std::endl;
	config_file << "#ion 1 0.000 2.0 # Counterion declaration:" << std::endl;
	config_file << "ion  1 0.150000 2.000000     # Counterion declaration:" << std::endl;
	config_file << "ion -1 0.150000 2.000000     # ion <charge> <conc (M)> <radius>" << std::endl;
	config_file << "ion  2 0.000000 2.000000     # ion <charge> <conc (M)> <radius>" << std::endl;
	config_file << "ion -2 0.000000 2.000000     # ion <charge> <conc (M)> <radius>" << std::endl;
	config_file << "pdie 4.000000          # Solute dielectric" << std::endl;
	config_file << "sdie 80.000000          # Solvent dielectric" << std::endl;
	config_file << "chgm spl2          # Charge disc method" << std::endl;
	config_file << "# 0 is linear splines" << std::endl;
	config_file << "# 1 is cubic b-splines" << std::endl;
	config_file << "mol 1            # which molecule to use" << std::endl;
	config_file << "srfm smol        # Surface calculation method" << std::endl;
	config_file << "#  0 => Mol surface for epsilon;" << std::endl;
	config_file << "#       inflated VdW for kappa; no" << std::endl;
	config_file << "#       smoothing" << std::endl;
	config_file << "#  1 => As 0 with harmoic average" << std::endl;
	config_file << "#       smoothing" << std::endl;
	config_file << "#  2 => Cubic spline " << std::endl;
	config_file << "srad 1.400000          # Solvent radius (1.4 for water)" << std::endl;
	config_file << "swin 0.3         # Surface cubic spline window .. default 0.3" << std::endl;
	config_file << "temp 310.000000          # System temperature (298.15 default)" << std::endl;
	config_file << "sdens 10.000000         # Specify the number of grid points per square-angstrom to use in Vacc object. Ignored when srad is 0.0 (see srad) or srfm is spl2 (see srfm). There is a direct correlation between the value used for the Vacc sphere density, the accuracy of the Vacc object, and the APBS calculation time. APBS default value is 10.0." << std::endl;
	config_file << "gamma 0.105      # Surface tension parameter for apolar forces (in kJ/mol/A^2)" << std::endl;
	config_file << "# only used for force calculations, so we don't care, but" << std::endl;
	config_file << "# it's always required, and 0.105 is the default" << std::endl;
	config_file << "calcenergy no    # Energy I/O to stdout" << std::endl;
	config_file << "#  0 => don't write out energy" << std::endl;
	config_file << "#  1 => write out total energy" << std::endl;
	config_file << "#  2 => write out total energy and all" << std::endl;
	config_file << "#       components" << std::endl;
	config_file << "calcforce no     # Atomic forces I/O (to stdout)" << std::endl;
	config_file << "#  0 => don't write out forces" << std::endl;
	config_file << "#  1 => write out net forces on molecule" << std::endl;
	config_file << "#  2 => write out atom-level forces" << std::endl;
	config_file << "write pot dx " << tag << "  # What to write .. this says write the potential in dx" << std::endl;
	config_file << "# format to a file." << std::endl;
	config_file << "end" << std::endl;
	config_file << "quit" << std::endl;
}


void
apply ( core::pose::Pose & pose )
{
	std::string tag = core::pose::extract_tag_from_pose( pose );
	std::string pqr_fn = tag + ".pqr";
	std::ofstream out_pqr(pqr_fn.c_str());
	utility::vector1<Size> zero_charge_chains;
	if (basic::options::option[PB_potential::no_charge_on_chain].user()) {
		zero_charge_chains = basic::options::option[PB_potential::no_charge_on_chain]();
	}
	core::scoring::dump_pqr(pose, out_pqr, tag, zero_charge_chains);

	std::string config_fn = tag + ".in";
	std::ofstream out_config(config_fn.c_str());
	write_APBS_config(pose, out_config, pqr_fn);
	out_config.close();

	if (!basic::options::option[PB_potential::apbs_exe].user()) {
		std::cerr << "Need to know where the APBS executable is to proceed." << std::endl;
		std::cerr << "Or, run apbs from a command line." << std::endl;
		return;
	}
	else {
		std::string command_line(basic::options::option[PB_potential::apbs_exe]() + " " + config_fn);
		system(command_line.c_str());
	}
}

void initialize() {
}

std::string
get_name() const {
	return "PBPotentialMover";
}
private:
};

int
main( int argc, char * argv [] )
{
	try {
	basic::options::option.add( PB_potential::apbs_exe, "APBS executable position" );
	basic::options::option.add( PB_potential::no_charge_on_chain, "chain with zero charge" );
	//basic::options::option.add( PB_potential::chain, "Only print given chains" );
	//basic::options::option.add( PB_potential::chain, "Only print given chains" );
	devel::init( argc, argv );
	PBPotentialMover PB_potential_mover;
	PB_potential_mover.initialize();

	protocols::jobdist::universal_main( PB_potential_mover );
	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
