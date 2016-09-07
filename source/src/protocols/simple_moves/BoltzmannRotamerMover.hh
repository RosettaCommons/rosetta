// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/BoltzmannRotamerMover.hh
/// @brief definition of BoltzmannRotamerMover class and functions
/// @author Noah Ollikainen (nollikai@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_BoltzmannRotamerMover_hh
#define INCLUDED_protocols_simple_moves_BoltzmannRotamerMover_hh

// Unit headers
#include <protocols/simple_moves/BoltzmannRotamerMover.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class BoltzmannRotamerMover : public protocols::moves::Mover {
public:

	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;
	typedef protocols::moves::MoverOP MoverOP;

public:

	// default constructor
	BoltzmannRotamerMover();

	/// @brief constructor with PackerTask
	BoltzmannRotamerMover(
		ScoreFunctionCOP scorefxn_in,
		PackerTaskOP & task_in
	);

	/// @brief constructor with TaskFactory
	BoltzmannRotamerMover(
		ScoreFunctionCOP scorefxn_in,
		TaskFactoryCOP factory_in
	);

	/// @brief copy constructor
	BoltzmannRotamerMover( BoltzmannRotamerMover const & rval );

	/// @brief destructor
	~BoltzmannRotamerMover() override;

	/// @brief clone this object
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of object
	protocols::moves::MoverOP fresh_instance() const override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	void show(std::ostream & output=std::cout) const override;

	// helpers
	core::Size select_rotamer(
		utility::vector1<std::pair<core::Size, core::Real> > const & boltzmann_factors,
		core::Real const & partition_function);

	// setters
	void set_score_function( core::scoring::ScoreFunctionCOP sf );
	void set_task_factory( core::pack::task::TaskFactoryCOP tf );
	void set_resnum( core::Size resnum );
	void set_ligand_resnum( core::Size ligand_resnum );
	void set_ligand_weight( core::Real ligand_weight );
	void set_temperature( core::Real temperature );
	void set_bias_sampling( bool bias_sampling );
	void set_randomize_resnum( bool randomize_resnum );
	void set_bump_check( bool bump_check );

	// getters
	core::Size get_resnum() const;
	core::Size get_ligand_resnum() const;
	core::Real get_ligand_weight() const;
	core::Real get_temperature() const;
	bool get_bias_sampling() const;
	bool get_randomize_resnum() const;
	bool get_bump_check() const;

public:

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) override;

protected:

	/// @brief read access for derived classes
	ScoreFunctionCOP
	scorefxn() const;

	/// @brief read access for derived classes, pose needed to run TaskFactory
	PackerTaskCOP
	task( core::pose::Pose const & pose ) const;

private:

	ScoreFunctionCOP scorefxn_;

	PackerTaskOP task_;

	TaskFactoryCOP factory_;

	bool show_packer_task_;

	/// @brief residue number specifying the next residue to move
	core::Size resnum_;

	/// @brief residue number specifying the ligand residue (0 if no ligand)
	core::Size ligand_resnum_;

	/// @brief weight for interaction between resnum_ and ligand_resnum_
	core::Real ligand_weight_;

	/// @brief kT value used for Boltzmann probability calculation
	core::Real temperature_;

	/// @brief if true, bias rotamer selection based on energy
	bool bias_sampling_;

	/// @brief if true, choose a random residue for the next move
	bool randomize_resnum_;

	/// @brief if true, use bump check when generating rotamers
	bool bump_check_;

};  // class BoltzmannRotamerMover

std::ostream &operator<< (std::ostream &os, BoltzmannRotamerMover const &mover);

} // moves
} // protocols

#endif
