// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Relax Baseclass
///
///
/// @author Mike Tyka


#ifndef INCLUDED_protocols_relax_RelaxProtocolBase_hh
#define INCLUDED_protocols_relax_RelaxProtocolBase_hh

// Unit headers
#include <protocols/relax/RelaxProtocolBase.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>


namespace protocols {
namespace relax {

class RelaxProtocolBase : public moves::Mover {
public:
	typedef moves::Mover parent;

	RelaxProtocolBase( core::scoring::ScoreFunctionOP );
	RelaxProtocolBase( std::string const & movername = "RelaxProtocol" );
	RelaxProtocolBase( std::string const & movername, core::scoring::ScoreFunctionOP );
	~RelaxProtocolBase() override;

	protocols::moves::MoverOP fresh_instance() const override { return clone(); };

	static void register_options();

	void apply_disulfides( core::pose::Pose & pose );

	// Default options -------------------------------------
	void set_defaults();
	void set_default_minimization_settings();
	void set_default_coordinate_settings();
	void set_default_movemap();

	// Public accessors
	core::kinematics::MoveMapCOP get_movemap() const;
	core::kinematics::MoveMapOP get_movemap();
	const core::scoring::ScoreFunctionCOP get_scorefxn() const;
	core::pack::task::TaskFactoryOP const & get_task_factory() const;

	bool cartesian() const { return cartesian_; }
	std::string min_type() const { return min_type_; }
	core::Size max_iter() const { return max_iter_; }
	bool dry_run() const { return dry_run_; }

	bool constrain_relax_to_native_coords() const { return constrain_relax_to_native_coords_; }
	bool constrain_relax_to_start_coords() const { return constrain_relax_to_start_coords_; }
	bool constrain_coords() const { return constrain_coords_; }
	bool coord_constrain_sidechains() const { return coord_constrain_sidechains_; }
	bool ramp_down_constraints() const { return ramp_down_constraints_; }
	bool constrain_relax_segments() const { return constrain_relax_segments_; }

	bool minimize_bond_lengths() const { return minimize_bond_lengths_; }
	bool minimize_bond_angles() const { return minimize_bond_angles_; }

	// Public mutators
	void set_movemap( core::kinematics::MoveMapOP movemap );
	void set_movemap_factory( core::select::movemap::MoveMapFactoryOP mm_factory );
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void set_task_factory( core::pack::task::TaskFactoryOP task_factory );

	/// @brief Use cartesian (minimization step)
	/// @details
	/// Sets to use the lbfgs_armijo_nonmonotone if true or FR default if false
	/// Recommended to set max_iter to 200.
	/// Requires scorefunction setup for non-ideal minimization.
	void cartesian( bool newval );
	void min_type( std::string min_type );
	void max_iter( core::Size max_iter );
	void dry_run( bool setting );

	void constrain_relax_to_native_coords( bool constrain_relax_to_native_coords );
	void constrain_relax_to_start_coords(  bool constrain_relax_to_start_coords );
	void constrain_coords( bool constrain_coords );
	void coord_constrain_sidechains( bool coord_constrain_sidechains );
	void constrain_relax_segments( bool constrain_relax_segments );
	void ramp_down_constraints( bool ramp_down_constraints );

	void minimize_bond_lengths( bool minimize_bond_lengths );
	void minimize_bond_angles( bool minimize_bond_angles );

public: // CitationManager fxns:

	/// @brief Does this mover provide information about how to cite it?
	/// @details Defaults to false.  Derived relax protocols may override this to provide citation info.  If set to
	/// true, the provide_citation_info() override should also be provided.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	bool mover_provides_citation_info() const override;

	/// @brief Provide the citation.
	/// @returns A vector of citation collections.  This allows the mover to provide citations for
	/// itself and for any modules that it invokes.
	/// @details The default implementation of this function provides citations only for the task operations
	/// in the task factory.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const override;

	/// @brief Does this mover indicate that it is unpublished (and, by extension, that the author should be
	/// included in publications resulting from it)?
	/// @details Defaults to false.  Derived relax protocols may override this to provide authorship info.  If set to
	/// true, the provide_authorship_info_for_unpublished() override should also be provided.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	bool mover_is_unpublished() const override;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	/// @returns A list of pairs of (author, e-mail address).  The default version only provides authorship information
	/// for gthe task operations in the task factory.  Empty list if not unpublished.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const override;


protected:

	core::scoring::ScoreFunctionOP get_scorefxn();

	// Accessors -------------------------------------
	bool fix_omega() const { return fix_omega_; }
	int minimize_bondlength_subset() const { return minimize_bondlength_subset_; }
	int minimize_bondangle_subset() const { return minimize_bondangle_subset_; }
	bool limit_aroma_chi2() const { return limit_aroma_chi2_; }
	std::string cst_files( core::Size const i ) const { return cst_files_[i]; }

	// Mutators -------------------------------------
	void fix_omega( bool fix_omega ) { fix_omega_ = fix_omega; }

	void minimize_bondangle_subset( int minimize_bondangle_subset ) { minimize_bondangle_subset_ = minimize_bondangle_subset; }
	void minimize_bondlength_subset( int minimize_bondlength_subset ) { minimize_bondlength_subset_ = minimize_bondlength_subset; }


	void initialize_movemap( core::pose::Pose const & pose, core::kinematics::MoveMap & movemap );
	void set_up_constraints( core::pose::Pose &pose,  core::kinematics::MoveMap & local_movemap );
	void output_debug_structure( core::pose::Pose & pose, std::string prefix );
	void add_cst_files( std::string const & cstfile ) { cst_files_.push_back( cstfile ); }


private:
	// Essentially MoveMap settings
	bool fix_omega_;
	bool minimize_bond_lengths_;
	bool minimize_bond_angles_;
	int minimize_bondangle_subset_;
	int minimize_bondlength_subset_;

	// Constraint settings
	bool constrain_relax_to_native_coords_;
	bool constrain_relax_to_start_coords_;
	bool constrain_coords_;
	bool coord_constrain_sidechains_;
	bool ramp_down_constraints_;
	bool constrain_relax_segments_;
	utility::vector1< std::string > cst_files_;

	// The minimizer algorithm
	std::string min_type_;

	/// Do cartesian-space minimization?
	bool cartesian_;

	/// maximum minimizer iterations
	core::Size max_iter_;

	bool limit_aroma_chi2_;

	// The MoveMap and the MoveMapFactory
	// Derived relax protocol
	core::kinematics::MoveMapOP movemap_;
	core::select::movemap::MoveMapFactoryOP movemap_factory_;

	// Fullatom scoring function used
	core::scoring::ScoreFunctionOP scorefxn_;

	/// task factory for packer ( used only FastRelax as of 11/05/2010 )
	core::pack::task::TaskFactoryOP task_factory_;

	/// @brief Is this a dry run (i.e. do no cycles)
	bool dry_run_;
};

}
} // protocols

#endif
