// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file EnvClaim
/// @brief Virtual base class for all claims used with environment claiming system.
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_claims_EnvClaim_hh
#define INCLUDED_protocols_environment_claims_EnvClaim_hh

// Unit Headers
#include <protocols/environment/claims/EnvClaim.fwd.hh>

// Package Headers
#include <core/environment/LocalPosition.fwd.hh>
#include <core/environment/FoldTreeSketch.hh>
#include <core/environment/SequenceAnnotation.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.hh>

#include <protocols/environment/claims/BrokerElements.hh>
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/ProtectedConformation.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

#include <basic/datacache/DataMap.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <utility/vector1.hh>

// option key includes


namespace protocols {
namespace environment {
namespace claims {

class EnvClaim : public utility::pointer::ReferenceCount {
	typedef core::select::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;
	typedef core::environment::FoldTreeSketch FoldTreeSketch;
	typedef std::map< std::string, ResidueSelectorCOP > AnnotatingSelectors;

public:
	/// @brief factory method for claims.
	/// @note I chose not to make a full-on factory implementation with registrators and creators and a factory
	///       because there aren't that many different kinds of claims, I think, and that level of complexity is
	///       likely superfluous.
	static
	EnvClaimOP make_claim( std::string const& name,
		ClientMoverOP owner,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap& datamap );

	static
	bool is_claim( std::string const& name );


	static
	void
	define_envclaim_schema_group( utility::tag::XMLSchemaDefinition & xsd );

	static
	std::string
	envclaim_ct_namer( std::string );

	static
	std::string
	envclaim_group_name();


	/// @brief Virtual destructor
	virtual ~EnvClaim();

	EnvClaim( ClientMoverOP );

	virtual EnvClaimOP clone() const = 0;

	/// @brief A clone used by the EnvClaimFactory to instantiate new EnvClaims using an XML tag.
	ClientMoverOP owner() const;

	void set_owner( ClientMoverOP owner ) { claim_source_ = owner; }

	/// @brief allow the claim to use any internally queued ResidueSelectors to create sequence annotations.
	void annotate( core::pose::Pose const&, core::environment::SequenceAnnotationOP ) const;

	/// @brief build ResidueElements that indicate the introduction of a new peptide edge into the fold tree.
	virtual void yield_elements( FoldTreeSketch const&, ResidueElements& ) const {};

	/// @brief build the JumpElements that represent the inclusion of a jump in the nascent FoldTree
	virtual void yield_elements( FoldTreeSketch const&, JumpElements& ) const {};

	/// @brief build and export the CutElements that represent the inclusion of a cut in the tree.
	virtual void yield_elements( FoldTreeSketch const&, CutElements& ) const {};

	/// @brief build and export the CutElements that represent the inclusion of a cut in the tree.
	virtual void yield_elements( FoldTreeSketch const&, CutBiasElements& ) const {};

	/// @brief build and export DOFElements, which represent control over non-jump dofs (torsions, bond lengths, angles) final conformation.
	virtual void yield_elements( core::pose::Pose const&, DOFElements& ) const {};

	virtual void show( std::ostream& os ) const;

	virtual std::string type() const = 0;

protected:
	virtual
	DOFElement wrap_dof_id( core::id::DOF_ID const& id ) const;

	ControlStrength parse_ctrl_str( std::string const& str ) const;

	void queue_for_annotation( std::string const& label, ResidueSelectorCOP selector );

private:

	AnnotatingSelectors selector_list_;
	ClientMoverOP claim_source_;

}; //class EnvClaim

extern std::ostream& operator<<( std::ostream& os, EnvClaim const& );
extern std::ostream& operator<<( std::ostream& os, EnvClaims const& );

}
}
}

#endif
