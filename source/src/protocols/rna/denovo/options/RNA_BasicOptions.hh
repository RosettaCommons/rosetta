// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/options/RNA_BasicOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_BasicOptions_HH
#define INCLUDED_protocols_farna_RNA_BasicOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/rna/denovo/options/RNA_BasicOptions.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

namespace protocols {
namespace rna {
namespace denovo {
namespace options {


class RNA_BasicOptions: public virtual basic::resource_manager::ResourceOptions {

public:

	//constructor
	RNA_BasicOptions();

	RNA_BasicOptions( RNA_BasicOptions const & src );

	//destructor
	~RNA_BasicOptions();

public:

	RNA_BasicOptionsOP clone() const;

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
	type() const{ return "RNA_BasicOptions";}

public:

	void set_dump_pdb( bool const setting ){ dump_pdb_ = setting; };
	bool dump_pdb() const { return dump_pdb_; }

	void set_move_first_rigid_body( bool const & setting ){ move_first_rigid_body_ = setting; }
	bool move_first_rigid_body() const { return move_first_rigid_body_; }

	void set_dock_into_density( bool const & setting ){ dock_into_density_ = setting; }
	bool dock_into_density() const { return dock_into_density_; }

	void set_model_with_density( bool const & setting ){ model_with_density_ = setting; }
	bool model_with_density() const { return model_with_density_; }

	void set_verbose( bool const & setting ){ verbose_ = setting; }
	bool verbose() const { return verbose_; }

private:

	bool dump_pdb_;
	bool move_first_rigid_body_;
	bool dock_into_density_;
	bool model_with_density_;
	bool verbose_;

};

} //options
} //denovo
} //rna
} //protocols

#endif
