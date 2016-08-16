// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/EnvClaimBroker.hh
/// @brief A class responsible for managing claims from ClientMovers, and retaining important information about the result of the brokering process.
///
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_environment_EnvClaimBroker_hh
#define INCLUDED_protocols_environment_EnvClaimBroker_hh

// Unit Headers
#include <protocols/environment/EnvClaimBroker.fwd.hh>

// Package headers
#include <core/environment/DofPassport.fwd.hh>
#include <core/environment/SequenceAnnotation.fwd.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <protocols/environment/Environment.fwd.hh>
#include <protocols/environment/ClientMover.fwd.hh>
#include <protocols/environment/ProtectedConformation.fwd.hh>

#include <protocols/environment/claims/BrokerElements.hh>
#include <protocols/environment/claims/EnvClaim.fwd.hh>

// Project headers
#include <basic/datacache/WriteableCacheableMap.fwd.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <map>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class EnvClaimBroker : public utility::pointer::ReferenceCount {
	typedef std::map< ClientMoverOP, core::environment::DofPassportOP > MoverPassMap;
	typedef core::environment::SequenceAnnotationCOP SequenceAnnotationCOP;
	typedef core::environment::SequenceAnnotationOP SequenceAnnotationOP;
	typedef core::environment::FoldTreeSketch FoldTreeSketch;

	typedef core::conformation::Conformation Conformation;
	typedef core::conformation::ConformationOP ConformationOP;
	typedef core::conformation::ConformationCOP ConformationCOP;

	typedef utility::vector1< std::pair< claims::ResidueElement, ClientMoverOP > > ResElemVect;
	typedef utility::vector1< std::pair< claims::JumpElement, ClientMoverOP > > JumpElemVect;
	typedef utility::vector1< std::pair< claims::CutElement, ClientMoverOP > > CutElemVect;
	typedef utility::vector1< std::pair< claims::CutBiasElement, ClientMoverOP > > CutBiasElemVect;
	typedef utility::vector1< std::pair< claims::DOFElement, ClientMoverOP > > DOFElemVect;

	typedef std::map< core::Size, std::string > SizeToStringMap;
	typedef std::map< std::string, std::pair< Size, Size > > StringToSizePairMap;
	typedef std::map< std::string, std::pair< std::string, std::string > > StringToStringPairMap;
	typedef utility::vector1< core::Real > BiasVector;

public:
	EnvClaimBroker( EnvironmentCAP env,
		MoverPassMap const& movers_and_passes,
		core::pose::Pose const& in_pose,
		SequenceAnnotationOP ann );

	virtual ~EnvClaimBroker();

	struct BrokerResult {
		core::pose::Pose pose;
		SizeToStringMap new_vrts;
		core::kinematics::FoldTreeOP closer_ft;
		std::set< core::Size > auto_cuts;
		basic::datacache::WriteableCacheableMapCOP cached_data;
		SequenceAnnotationCOP ann;
	};

	BrokerResult const& result() const;

private:

	/// @brief an inner class for tracking the properties of jumps that've been brokered.
	class BrokeredJumpData : public utility::pointer::ReferenceCount {
	public:
		BrokeredJumpData( std::pair< core::Size, core::Size > const& positions,
			std::pair< std::string, std::string > const& atoms,
			bool put_jump_stub_intra_residue  );
		bool operator==( BrokeredJumpData const& ) const;
		bool operator!=( BrokeredJumpData const& other ) const { return !this->operator==( other ); }
		std::ostream& operator<<( std::ostream& ) const;

		std::pair< core::Size, core::Size > const pos;
		std::pair< std::string, std::string > const atoms;
		bool const put_jump_stub_intra_residue;
	};

	typedef utility::pointer::shared_ptr< BrokeredJumpData > BrokeredJumpDataOP;
	typedef utility::pointer::shared_ptr< BrokeredJumpData const > BrokeredJumpDataCOP;
	typedef std::map< std::string, BrokeredJumpDataCOP > JumpDataMap;

	void broker_fold_tree( Conformation&, basic::datacache::BasicDataCache& );

	void broker_dofs( core::pose::Pose& );

	core::kinematics::FoldTreeOP render_fold_tree( FoldTreeSketch& fts,
		utility::vector1< core::Real > const& bias,
		basic::datacache::BasicDataCache&,
		core::conformation::Conformation const& input_conf );

	void annotate_fold_tree( core::kinematics::FoldTreeOP,
		JumpDataMap const& new_jumps,
		SequenceAnnotationOP = NULL );

	void add_virtual_residues( Conformation&,
		SizeToStringMap const& new_vrts,
		SequenceAnnotationOP );

	void process_elements( ResElemVect const& elems,
		FoldTreeSketch& fts,
		SizeToStringMap& new_vrts );

	void process_elements( JumpElemVect const& elems,
		FoldTreeSketch& fts,
		JumpDataMap& new_jumps );

	void process_elements( CutElemVect const& elems,
		FoldTreeSketch& fts );

	void process_elements( CutBiasElemVect const& elems, BiasVector& bias );

	void grant_access( claims::DOFElement const& e, ClientMoverOP owner ) const;

	void setup_passports( DOFElemVect& elements,
		claims::ControlStrength const& (*str_access)( std::pair< claims::DOFElement, ClientMoverOP > const& ) );

	template < typename T, typename I >
	utility::vector1< std::pair< T, ClientMoverOP > > collect_elements( I const& info ) const;

	claims::EnvClaims collect_claims( MoverPassMap const& movers_and_passes,
		core::pose::Pose& pose );

	/// @brief broker new residues
	void build_new_residues( claims::EnvClaims const& claims, FoldTreeSketch& fts, SequenceAnnotationOP ann );

	/// @brief use accepted claims to build DofPassport objects for movers.
	void assign_passports( claims::EnvClaims const&, ProtectedConformation const& );

	void add_chainbreak_variants( core::Size cut, core::conformation::Conformation& conf ) const;

	MoverPassMap const& movers_and_passes_;
	SequenceAnnotationOP ann_;
	claims::EnvClaims claims_;
	BrokerResult result_;
	EnvironmentCAP env_;

}; // EnvClaimBroker

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_EnvClaimBroker_hh
