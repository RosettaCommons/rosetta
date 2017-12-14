// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockingEnsemble.cc
/// @brief container class for ensemble docking information to be used with ConformerSwitchMover
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_docking_DockingEnsemble_HH
#define INCLUDED_protocols_docking_DockingEnsemble_HH

// Unit headers
#include <protocols/docking/DockingEnsemble.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

//Option key includes

// ObjexxFCL Headers

// C++ Headers
#include <map>
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace docking {

/// @brief this mover does the conformer swap in RosettaDock's ensemble docking.  It takes
/// in a multi-model PDB file as an ensemble, and does swaps conformers by superpositioning
/// over interface residues, and selects a conformer based on a partition function using
/// a ScoreFunction.
class DockingEnsemble : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~DockingEnsemble() override;

	//default constructor
	DockingEnsemble() :
		start_res_(1),
		end_res_(1),
		jump_id_(1),
		ensemble_size_(0),
		ensemble_file_path_(""),
		partner_("")
	{}

	//constructor with arguments
	DockingEnsemble(
		Size start_res,
		Size end_res,
		Size jump_id,
		std::string const & ensemble_file_path,
		std::string const & partner,
		core::scoring::ScoreFunctionCOP scorefxn_low,
		core::scoring::ScoreFunctionCOP scorefxn_high
	);

	// @brief recover sidechains for a conformer
	void recover_conformer_sidechains( core::pose::Pose & pose );

	// @brief calculate reference energies for prepacking
	void calculate_lowres_ref_energy( core::pose::Pose & pose );
	void calculate_highres_ref_energy( core::Size conf_num );

	// @brief update pdblist file with reference energies
	void update_pdblist_file();

	void set_packer( protocols::moves::SequenceMoverOP packer );

	void set_current_confnum( core::Size conf_num ) { conf_num_ = conf_num; }
	core::Size get_current_confnum() { return conf_num_; }

	// @brief return a pose from the ensemble (Either fullatom or centroid)
	core::pose::Pose & get_conformer( core::Size conf_num ) { return ensemble_list_[conf_num]; }
	core::pose::Pose & get_conformer_cen( core::Size conf_num ) { return ensemble_list_cen_[conf_num]; } // Add by DK

	// @brief simple getters
	core::Size size() { return ensemble_size_; }
	core::Size jump_id() { return jump_id_; }
	core::Size start_res() { return start_res_; }
	core::Size end_res() { return end_res_; }
	core::Size conf_size() { return conf_size_; }
	std::string partner() { return partner_; }

	// return reference energies.
	// Can either be done with the current conformer (for example at the end of docking)
	// or by passing the conf_num of the conformer you want an energy for (for example when calculating Probability Tables)
	core::Real lowres_reference_energy( core::Size conf_num ) { return lowres_reference_energies_[ conf_num ]; }
	core::Real lowres_reference_energy() { return lowres_reference_energies_[ conf_num_ ]; }
	core::Real highres_reference_energy( core::Size conf_num ) { return highres_reference_energies_[ conf_num ]; }
	core::Real highres_reference_energy() { return highres_reference_energies_[ conf_num_ ]; }

	core::scoring::ScoreFunctionCOP scorefxn_low() { return scorefxn_low_; }

private:
	void load_ensemble();

	Size start_res_, end_res_, conf_size_, jump_id_, ensemble_size_, conf_num_;
	std::string ensemble_file_path_;
	std::string partner_; // which partner is this a conformer for?
	core::scoring::ScoreFunctionCOP scorefxn_low_, scorefxn_high_;
	utility::vector1< std::string > pdb_filenames_;
	utility::vector1< core::pose::Pose > ensemble_list_;
	utility::vector1< core::pose::Pose > ensemble_list_cen_;
	utility::vector1< core::Real > lowres_reference_energies_, highres_reference_energies_;

	protocols::moves::SequenceMoverOP pack_operations_;
}; //mover

} // docking
} // rosetta


#endif
