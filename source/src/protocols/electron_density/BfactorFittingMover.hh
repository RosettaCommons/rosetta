// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file

#ifndef INCLUDED_protocols_electron_density_BfactorFittingMover_hh
#define INCLUDED_protocols_electron_density_BfactorFittingMover_hh

#include <protocols/moves/Mover.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>


#include <core/optimization/Multifunc.hh>
#include <ObjexxFCL/FArray3D.hh>

namespace protocols {
namespace electron_density {


/// Bfactor multifunc
class BfactorMultifunc : public core::optimization::Multifunc {
public:
	BfactorMultifunc(
		core::pose::Pose & pose_in,
		core::Real wt_adp,
		core::Real wt_dens,
		core::Real rmax,
		core::Real radius_exp,
		core::Real scorescale,
		bool exact, bool verbose, bool deriv_check );

	~BfactorMultifunc() override = default;

	core::Real
	operator ()( core::optimization::Multivec const & vars ) const override;

	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const override;

	void
	dump( core::optimization::Multivec const & x1, core::optimization::Multivec const & x2 ) const override;

	// helper functions
	void
	poseBfacts2multivec( core::pose::Pose & pose, core::optimization::Multivec &y ) const;

	void
	multivec2poseBfacts( core::optimization::Multivec const &vars, core::pose::Pose & pose ) const;

private:
	core::pose::Pose & pose_;

	core::Real scorescale_;  // overall weighing on score
	core::Real wt_adp_;      // weight on constraints
	core::Real wt_dens_;     // weight on density
	core::Real rmax_;        // max radius (def 5)
	core::Real radius_exp_;  // weigh neighbors by 1/r^thisval (def 1)
	bool exact_;
	bool verbose_, deriv_check_;

	// handles AtomID <-> vector index mapping
	core::id::AtomID_Map< core::Size > atom_indices_;
	utility::vector1< core::id::AtomID > moving_atoms_;

	mutable ObjexxFCL::FArray3D< core::Real > rhoC_, rhoMask_;  // temp storage avoids reallocation


	core::Real B_EPS;
};


/// mover to fit B factors
class BfactorFittingMover : public moves::Mover {
public:
	BfactorFittingMover();
	~BfactorFittingMover() override = default;

	void init();

	void apply( core::pose::Pose & ) override;


	moves::MoverOP clone() const override { return utility::pointer::make_shared< BfactorFittingMover >( *this ); }
	moves::MoverOP fresh_instance() const override { return utility::pointer::make_shared< BfactorFittingMover >(); }

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void set_scorescale( core::Real scorescale ){ scorescale_ = scorescale; };
	void set_wt_adp( core::Real wt_adp ){ wt_adp_ = wt_adp; };
	void set_wt_dens( core::Real wt_dens ){ wt_dens_ = wt_dens; };
	void set_rmax( core::Real rmax ){ rmax_ = rmax; };
	void set_radius_exp( core::Size radius_exp ){ radius_exp_ = radius_exp; };
	void set_max_iter( core::Size max_iter ){ max_iter_ = max_iter; };
	void set_minimizer( std::string minimizer ){ minimizer_ = minimizer; };
	void set_init( bool init ){ init_ = init; };
	void set_exact( bool exact ){ exact_ = exact; };

	core::Real get_scorescale(){ return scorescale_; };
	core::Real get_wt_adp(){ return wt_adp_; };
	core::Real get_wt_dens(){ return wt_dens_; };
	core::Real get_rmax(){ return rmax_; };
	core::Size get_radius_exp(){ return radius_exp_; };
	core::Size get_max_iter(){ return max_iter_; };
	std::string get_minimizer(){ return minimizer_; };
	bool get_init(){ return init_; };
	bool get_exact(){ return exact_; };


private:
	core::Real scorescale_;  // overall weighing on score
	core::Real wt_adp_;      // weight on constraints
	core::Real wt_dens_;     // weight on density
	core::Real rmax_;        // max radius (def 5)
	core::Size radius_exp_;  // weigh neighbors by 1/r^thisval (def 1)
	core::Size max_iter_;    // minimizer iterations
	std::string minimizer_;
	bool init_;
	bool exact_;
	bool verbose_, deriv_check_;

	core::Real res_low_,res_high_;
	core::Size nresbins_;

	utility::vector1< core::Real > weights_to_scan_;
};

} // moves
} // protocols

#endif
