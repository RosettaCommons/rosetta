// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/minimization_packing/PertMinMover.hh
/// @brief Roberto Chica inspired random perturbation followed by minimization
/// @details For original, see Davey & Chica, Proteins 82:771-784 doi:10.1002/prot.24457
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_protocols_minimization_packing_PertMinMover_hh
#define INCLUDED_protocols_minimization_packing_PertMinMover_hh

// Unit headers
#include <protocols/minimization_packing/PertMinMover.fwd.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace minimization_packing {

///@brief Random perturbation, followed by minimization
class PertMinMover : public protocols::moves::Mover {

public:

	PertMinMover();

	// copy constructor
	PertMinMover( PertMinMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~PertMinMover();

	core::Real pert_size() const { return pert_size_; }
	bool uniform() const { return uniform_; }
	bool sc_only() const { return sc_only_; }
	core::select::residue_selector::ResidueSelectorCOP residues() const { return residues_; }
	core::scoring::ScoreFunctionCOP scorefxn() const { return scorefxn_; }
	core::select::movemap::MoveMapFactoryCOP movemap() const { return movemap_factory_; }

	void pert_size(core::Real setting) { pert_size_ = setting; }
	void uniform(bool setting) { uniform_ = setting; }
	void sc_only(bool setting) { sc_only_ = setting; }
	void residues(core::select::residue_selector::ResidueSelectorCOP setting) { residues_ = setting; }
	void scorefxn(core::scoring::ScoreFunctionCOP setting) { scorefxn_ = setting; }
	void movemap_factory(core::select::movemap::MoveMapFactoryCOP setting) { movemap_factory_ = setting; }

	virtual void
	apply( core::pose::Pose & pose );

	void
	pert( core::pose::Pose & pose, utility::vector1< bool > const & resi ) const;

	void
	min( core::pose::Pose & pose ) const;

public:
	virtual void
	show( std::ostream & output=std::cout ) const;

	std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//PertMinMover & operator=( PertMinMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;

	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	/// @brief The size of the random perturbation
	core::Real pert_size_;
	/// @brief Do a uniform perturbation (default is Gaussian)
	bool uniform_;
	/// @brief Don't perturb the backbone coordinates
	bool sc_only_;
	/// @brief Only perturb the selected residues
	core::select::residue_selector::ResidueSelectorCOP residues_;

	/// @brief the score function to use in minimization
	core::scoring::ScoreFunctionCOP scorefxn_;
	/// @brief The movemap factory to be used in minimization
	core::select::movemap::MoveMapFactoryCOP movemap_factory_;
};

std::ostream &operator<< (std::ostream &os, PertMinMover const &mover);


} //protocols
} //minimization_packing


#endif //protocols/minimization_packing_PertMinMover_hh







