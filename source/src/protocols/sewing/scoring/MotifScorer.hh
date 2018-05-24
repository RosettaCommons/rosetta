// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MotifScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_MotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_MotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/MotifScorer.fwd.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.fwd.hh>

//Package headers
#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.fwd.hh>
//Core headers
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/types.hh>

//Utility headers
#include <utility/vector1.hh>
#include <numeric/xyzTransform.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class MotifScorer : public AssemblyScorer {

public:

	///@brief default construct
	MotifScorer();

	virtual ~MotifScorer()=default;
	MotifScorer( MotifScorer const & );

	static std::string type_name();

	//Virtual because derived classes override this as well
	virtual
	core::Real
	score(
		data_storage::SmartAssemblyCOP assembly
	) override;

	core::Real
	norm_motif_score(
		data_storage::SmartAssemblyCOP assembly
	);

	core::Real
	full_motif_score(
		data_storage::SmartAssemblyCOP assembly
	);

	numeric::xyzTransform<core::Real>
	get_stub(
		data_storage::SmartSewingResidueCOP residue
	) const;

	core::Real
	get_score(
		numeric::xyzTransform<core::Real> stub1,
		char ss1,
		char aa1,
		numeric::xyzTransform<core::Real> stub2,
		char ss2,
		char aa2
	) const;

	virtual std::string
	get_name() const override;

	virtual core::Real
	get_weight() const override;

	virtual void
	set_weight( core::Real ) override;


	virtual void
	set_options_from_tag(
		utility::tag::TagCOP scorer_tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose) override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	/*
	void
	dump_motif(
	AssemblyCOP assembly
	) const;
	*/
	core::Real
	get_last_score() const override;
	void
	set_last_score(core::Real score) override;
	core::Real
	get_old_last_score() const override;
	void
	set_old_last_score( core::Real score ) override;

protected:

	core::scoring::motif::MotifHashManager & mman_;
	core::chemical::ResidueTypeSetCOP res_type_set_;
	//private:
	core::Real last_score_=1000;
	core::Real old_last_score_=1000;

private:
	core::Real weight_=1;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
