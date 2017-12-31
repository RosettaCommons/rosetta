// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/import_pose/options/RNA_MinimizerOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_import_pose_options_RNA_MinimizerOptions_HH
#define INCLUDED_core_import_pose_options_RNA_MinimizerOptions_HH

#include <core/import_pose/options/RNA_BasicOptions.hh>
#include <core/import_pose/options/RNA_MinimizerOptions.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

namespace core {
namespace import_pose {
namespace options {

class RNA_MinimizerOptions: public virtual RNA_BasicOptions {

public:

	//constructor
	RNA_MinimizerOptions();

	RNA_MinimizerOptions( RNA_MinimizerOptions const & src );

	//destructor
	~RNA_MinimizerOptions();

public:

	RNA_MinimizerOptionsOP clone() const;

	void
	initialize_from_command_line();
	void
	initialize_from_options( utility::options::OptionCollection const & opts );
	static void
	list_options_read( utility::options::OptionKeyList & opts );

	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag( utility::tag::TagCOP ){}

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const{ return "RNA_MinimizerOptions";}

	void set_vary_bond_geometry( bool const & setting ){ vary_bond_geometry_ = setting; }
	bool vary_bond_geometry() const { return vary_bond_geometry_; }

	void set_extra_minimize_res( utility::vector1< Size > const & extra_minimize_res ) { extra_minimize_res_ = extra_minimize_res; }
	utility::vector1< Size > const & extra_minimize_res() const { return extra_minimize_res_; };

	void set_extra_minimize_chi_res( utility::vector1< Size > const & extra_minimize_chi_res ) { extra_minimize_chi_res_ = extra_minimize_chi_res; }
	utility::vector1< Size > const & extra_minimize_chi_res() const { return extra_minimize_chi_res_; };

	void set_minimizer_use_coordinate_constraints( bool const & setting ){ minimizer_use_coordinate_constraints_ = setting; }
	bool minimizer_use_coordinate_constraints() const { return minimizer_use_coordinate_constraints_; }

	void set_min_type( std::string const & setting ){ min_type_ = setting; }
	std::string min_type() const { return min_type_; }

	void set_minimize_bps( bool const & setting ){ minimize_bps_ = setting; }
	bool minimize_bps() const { return minimize_bps_; }

	void set_minimize_all_protein( bool const & setting ){ minimize_all_protein_ = setting; }
	bool minimize_all_protein() const { return minimize_all_protein_; }

	void set_minimize_protein_sc( bool const & setting ){ minimize_protein_sc_ = setting; }
	bool minimize_protein_sc() const { return minimize_protein_sc_; }

	void set_protein_packing( bool const & setting ){ protein_packing_ = setting; }
	bool protein_packing() const { return protein_packing_; }

	void set_protein_pack_all( bool const & setting ){ protein_pack_all_ = setting; }
	bool protein_pack_all() const { return protein_pack_all_; }

	void set_protein_packing_distance( core::Real const & setting ){ protein_packing_distance_ = setting; }
	core::Real protein_packing_distance() const { return protein_packing_distance_; }

	void set_deriv_check( bool const & setting ){ deriv_check_ = setting; }
	bool deriv_check() const { return deriv_check_; }

	void set_skip_o2prime_trials( bool const & setting ){ skip_o2prime_trials_ = setting; }
	bool skip_o2prime_trials() const { return skip_o2prime_trials_; }

	void set_minimize_rounds( core::Size const & setting ){ minimize_rounds_ = setting; }
	core::Size minimize_rounds() const { return minimize_rounds_; }

private:

	core::Size minimize_rounds_;
	bool deriv_check_;
	bool skip_o2prime_trials_;
	bool vary_bond_geometry_;

	utility::vector1< core::Size > extra_minimize_res_;
	utility::vector1< core::Size > extra_minimize_chi_res_;

	bool minimizer_use_coordinate_constraints_;
	std::string min_type_;
	bool minimize_bps_;
	bool minimize_all_protein_;
	bool minimize_protein_sc_;
	bool protein_packing_;
	bool protein_pack_all_;
	core::Real protein_packing_distance_;

};

} //options
} //denovo
} //rna

#endif
