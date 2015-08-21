// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/ResidueVicinityRCG.hh
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009

#ifndef INCLUDED_protocols_forge_remodel_ResidueVicinityRCG_hh
#define INCLUDED_protocols_forge_remodel_ResidueVicinityRCG_hh

// protocol headers
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>

//utility headers
#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace remodel {

class ResidueVicinityRCG;

typedef utility::pointer::shared_ptr< ResidueVicinityRCG > ResidueVicinityRCGOP;
typedef utility::pointer::weak_ptr< ResidueVicinityRCG const > ResidueVicinityRCGCAP;

class ResidueVicinityInfo;
typedef utility::pointer::shared_ptr< ResidueVicinityInfo > ResidueVicinityInfoOP;
typedef utility::pointer::shared_ptr< ResidueVicinityInfo const > ResidueVicinityInfoCOP;

/// @brief small helper class for the ResidueVicinityRCG
class
	ResidueVicinityInfo : public utility::pointer::ReferenceCount
{

public:
	ResidueVicinityInfo(
		core::Size old_seqpos,
		utility::vector1< core::Size > const & residue_atoms,
		utility::vector1< core::Size > const & loopres_atoms,
		core::Size desired_remodelres_in_vicinity
	);

	virtual ~ResidueVicinityInfo();

	core::Size
	old_seqpos() const {
		return old_seqpos_; }

	utility::vector1< core::Size > const &
	residue_atoms() const {
		return residue_atoms_; }

	utility::vector1< core::Size > const &
	residue_base_atoms() const {
		return residue_base_atoms_; }

	utility::vector1< core::Size > const &
	residue_base2_atoms() const {
		return residue_base2_atoms_; }

	utility::vector1< core::Size > const &
	loopres_atoms() const {
		return loopres_atoms_; }

	utility::vector1< core::Size > const &
	loopres_base_atoms() const {
		return loopres_base_atoms_; }

	utility::vector1< core::Size > const &
	loopres_base2_atoms() const {
		return loopres_base2_atoms_; }

	void
	set_residue_base_atoms(
		utility::vector1< core::Size > const & residue_base_atoms ) {
		residue_base_atoms_ = residue_base_atoms; }

	void
	set_residue_base2_atoms(
		utility::vector1< core::Size > const & residue_base2_atoms ) {
		residue_base2_atoms_ = residue_base2_atoms; }

	void
	set_loopres_base_atoms(
		utility::vector1< core::Size > const & loopres_base_atoms ) {
		loopres_base_atoms_ = loopres_base_atoms; }

	void
	set_loopres_base2_atoms(
		utility::vector1< core::Size > const & loopres_base2_atoms ) {
		loopres_base2_atoms_ = loopres_base2_atoms; }


	core::Size
	desired_remodelres_in_vicinity() const {
		return desired_remodelres_in_vicinity_; }

	core::scoring::func::FuncOP
	dis() const;

	core::scoring::func::FuncOP
	loop_ang() const;

	core::scoring::func::FuncOP
	targ_ang() const;

	core::scoring::func::FuncOP
	loop_dih() const;

	core::scoring::func::FuncOP
	targ_dih() const;

	core::scoring::func::FuncOP
	lt_dih() const;

	void
	set_dis( core::scoring::func::FuncOP dis );

	void
	set_loop_ang( core::scoring::func::FuncOP loop_ang );

	void
	set_targ_ang( core::scoring::func::FuncOP targ_ang );

	void
	set_loop_dih( core::scoring::func::FuncOP loop_dih );

	void
	set_targ_dih( core::scoring::func::FuncOP targ_dih );

	void
	set_lt_dih( core::scoring::func::FuncOP lt_dih );

private:

	core::Size old_seqpos_;
	utility::vector1< core::Size > residue_atoms_, residue_base_atoms_, residue_base2_atoms_;
	utility::vector1< core::Size > loopres_atoms_, loopres_base_atoms_, loopres_base2_atoms_;

	core::scoring::func::FuncOP dis_, loop_ang_, targ_ang_, loop_dih_, targ_dih_, lt_dih_;

	core::Size desired_remodelres_in_vicinity_;

}; //ResidueVicinityInfo


/// @brief a RemodelConstraintGenerator that creates AmbiguousMultiConstraints for all positions
/// @brief in the remodeled loop to the desired positions, such that during remodeling, the remodeled
/// @brief region will be driven towards the vicinity of these residues
class ResidueVicinityRCG : public RemodelConstraintGenerator
{

public:
	// Constructors and virtual functions
	ResidueVicinityRCG();

	/// @brief copy construtor
	ResidueVicinityRCG( ResidueVicinityRCG const & rval );

	ResidueVicinityRCG(
		core::Size lstart,
		core::Size lstop,
		utility::vector1< ResidueVicinityInfoOP > const & rv_infos
	);

	virtual ~ResidueVicinityRCG();

	virtual
	void
	generate_remodel_constraints(
		core::pose::Pose const & pose );

	virtual void
	parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	virtual std::string
	get_name() const;

	virtual protocols::moves::MoverOP
	fresh_instance() const;

	virtual protocols::moves::MoverOP
	clone() const;

public:
	// Public member functions
	void
	clear_rv_infos(){
		rv_infos_.clear(); }


	void
	add_rv_info(
		ResidueVicinityInfoOP rv_info ){
		rv_infos_.push_back( rv_info ); }

	void
	lstart( core::Size const lstart );

	void
	lstop( core::Size const lstop );

	void
	set_rv_infos( utility::vector1< ResidueVicinityInfoOP > const & rv_infos );

protected:

	void
	generate_remodel_constraints_for_respair(
		core::pose::Pose const & pose,
		core::Size const loopres,
		ResidueVicinityInfo const & rv_info,
		utility::vector1< core::scoring::constraints::ConstraintCOP > & csts
	);


private:

	core::Size lstart_, lstop_;

	utility::vector1< ResidueVicinityInfoOP > rv_infos_;


}; //class ResidueVicinityRCG


} //namespace remodel
} //namespace forge
} //namespace protocols


#endif // INCLUDED_protocols_forge_remodel_ResidueVicinityRCG_HH
