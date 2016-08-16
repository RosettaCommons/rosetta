// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PeriodicBoxMover.hh
/// @brief Mover for running liquid simulation (and related others)
/// @author Frank DiMaio & Hahnbeom Park (hahnbeom@gmail.com)
/// @details implementation of MC liquid simulation in Rosetta.
/// Reference: William Jorgensen et al. "Development and testing of the OPLS all-atom
/// force field on conformational energetics and properties", JACS 118 (1996), 11225-11235.
/// The purpose of the implementation is to run "liquid simulation" to get optimized parameters
/// for LJ, hbond, electrostatic terms in Rosetta. Please refer to the reference to see how this works.
/// Converting output scorefile into thermodynamic data requires additional script.
/// Please e-mail hahnbeom@gmail.com to request for analysis script.

#ifndef INCLUDED_protocols_moves_PeriodicBoxMover_hh
#define INCLUDED_protocols_moves_PeriodicBoxMover_hh

#include <protocols/simple_moves/PeriodicBoxMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// C++ Headers
#include <string>
#include <fstream>

namespace protocols {
namespace simple_moves {

/// @brief structure that stores data during simulation
/// @note intraE and interE are incorrect for now;
struct ThermodynamicData {
public:
	core::Real volume;
	core::Real density;
	core::Real intraE;
	core::Real interE;
	core::Real centralE;
};

// WaterReporter commented out for now because of PyRosetta compilation issue...
/*
class PeriodicBoxMoverWaterReporter : public utility::pointer::ReferenceCount {
public:
PeriodicBoxMoverWaterReporter(std::string suffix) {
volout_eq_.open((std::string("volume_")+suffix+std::string(".txt")).c_str());
volout_eq_ << "#step sidelength volume density\n";
RDF_eq_.open((std::string("RDF_")+suffix+std::string(".txt")).c_str());
hbond_eq_.open((std::string("hbond_")+suffix+std::string(".txt")).c_str());
hbond_eq_ << "#step unique_hbond total_hbond hbond_per_residue\n";
print_RDF_header_ = true;
}
void
report( PeriodicBoxMover & mover,
int step,
core::pose::Pose & pose,
core::Real mweight,
core::Size lattice_jump );
private:
std::ofstream volout_eq_, RDF_eq_, hbond_eq_;
bool print_RDF_header_;
};
*/

//typedef utility::pointer::shared_ptr< PeriodicBoxMoverWaterReporter > PeriodicBoxMoverWaterReporterOP;

class PeriodicBoxMover : public protocols::moves::Mover {
	//friend class PeriodicBoxMoverWaterReporter;
public:
	PeriodicBoxMover();
	//PeriodicBoxMover( std::string const & );

	virtual void apply( Pose & pose );

	void setup_pose( Pose & pose, core::Real &, core::Size & );
	void change_volume_move( Pose & pose, core::Size, bool& );
	void perturb_molecule_move( Pose & pose, core::Size, bool& );

	//  void show_residue_hbonds( Pose const & pose, core::Size );
	//utility::vector1<core::Real> RDF( Pose & pose, core::Real, core::Real, core::Size );
	//core::Real dist( Pose & pose, core::Size, std::string, core::Size, std::string);

	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & pose );

private:

	// helper
	core::conformation::ResidueOP make_vrt( core::Vector O, core::Vector X, core::Vector Y ) {
		core::conformation::ResidueOP vrtrsd
			( core::conformation::ResidueFactory::create_residue(
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( "VRT" ) ) );
		vrtrsd->set_xyz("ORIG",O);
		vrtrsd->set_xyz("X",O+X);
		vrtrsd->set_xyz("Y",O+Y);
		return vrtrsd;
	}

	void
	report_thermodynamics( Pose & pose,
		core::Size lattice_jump );
	void
	recenter_pose( Pose & pose,
		core::Real const lattice,
		core::Real const new_lattice
	) const;

	void
	dump_ASU( Pose & pose, core::Size &lattice_jump, std::string filename );

	void
	add_thermodynamic_data_to_silent( core::io::silent::SilentStructOP ss ) const;

	void
	check_virial_pressure( Pose & pose,
		core::Real L
	) const;


	/// @brief correction to energy that accounts for missing vdw interaction from cut-off to infinity
	void
	setup_LJcorrection( core::pose::Pose const & pose );

	// scorefunction
	core::scoring::ScoreFunctionOP sf_;

	// setup
	core::Real nmol_side_, initial_density_, origin_;

	// parameters
	core::Real temp_, vol_step_, trans_step_, rot_step_, tor_step_;
	core::Size resize_vol_every_;
	core::Real probability_rigid_;

	// simulation length
	core::Size nsteps_equilibrate_, nsteps_sim_, istart_;

	// outputs
	core::Size dump_every_, report_every_;
	core::Size report_thermodynamics_;
	core::Size report_water_;
	std::string rg_atom1_, rg_atom2_;
	std::string report_silent_, report_scorefile_;
	ThermodynamicData thermodynamic_data_;

	// data structure for LJ correction
	core::Real ljcorrection_factor_;
	bool correct_LJtruncation_;

	// constant pressure setups
	core::Real P0_; // unit in atmosphere

	//PeriodicBoxMoverWaterReporterOP water_reporter_;

	// for two molecule simulations
	core::Size central_resno_;
	std::string central_molecule_pdb_;

};

} // moves
} // protocols

#endif
