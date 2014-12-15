// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange
#ifndef INCLUDED_protocols_abinitio_Templates_hh
#define INCLUDED_protocols_abinitio_Templates_hh

// Unit Headers
#include <protocols/abinitio/Templates.fwd.hh>

// Package Headers
#include <protocols/abinitio/PairingStatistics.fwd.hh>
#include <protocols/abinitio/TemplateJumpSetup.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/dssp/PairingsList.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
// AUTO-REMOVED #include <core/scoring/dssp/StrandPairing.hh>
#ifdef __clang__
#include <core/pose/Pose.hh>
#include <core/fragment/SecondaryStructure.hh>
#endif
// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

//// C++ headers
#include <string>
#include <map>

#include <protocols/abinitio/Template.fwd.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <utility/vector1.hh>

#ifdef WIN32
	#include <protocols/abinitio/Template.hh>
#endif


namespace protocols {
namespace abinitio {

class Templates : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Templates >
{
public:
  typedef std::map< std::string, TemplateOP > TemplateMap;

public:
	typedef TemplateMap::const_iterator const_iterator;
	typedef utility::vector1< TemplateCOP > TemplateList;

	static void register_options();

	Templates( std::string const& config_file, core::pose::PoseCOP native = NULL );
	virtual ~Templates();
	
	core::fragment::FragSetOP pick_frags( core::fragment::FragSetOP, core::fragment::FragDataCOP frag_type, Size min_nr_frags, Size ncopies = 1 ) const;
	Size pick_frags( core::fragment::FragSet&, core::fragment::FragDataCOP frag_type, Size ncopies = 1 ) const;

	void read_pairings( std::string const& filename, core::scoring::dssp::PairingsList& ) const;

	TemplateJumpSetupOP create_jump_def( core::fragment::SecondaryStructureCOP = NULL ) const;

	bool has_pairings() const {
		return pairings_.size();
	}

	core::scoring::dssp::PairingsList const& pairings() const {
		return pairings_;
	}

	std::string const& target_sequence() const {
		return target_sequence_;
	}

	std::string& target_sequence() {
		return target_sequence_;
	}

	core::Size target_total_residue() const {
		return target_sequence().size();
	}

	const_iterator begin() const {
		return templates_.begin();
	}

	const_iterator end() const {
		return templates_.end();
	}

	TemplateList const& helixjump_picks() const {
		return helixjump_pick_list_;
	}

	void
	add_target_constraints( core::scoring::constraints::ConstraintSetOP, core::pose::Pose const& ) const;

	bool has_template( std::string const& ModelID ) const {
		return ( templates_.find( ModelID ) != templates_.end() );
	}

	Template const& get_template( std::string const& ModelID ) const {
		TemplateMap::const_iterator iter ( templates_.find( ModelID ) );
		if ( iter == templates_.end() ) {
			utility_exit_with_message("Unknown model name: " + ModelID );
		}
		return *(iter->second);
	}

	bool
	is_good() const { return good_; };


	Size pick_large_frags(
			core::fragment::FragSet& frag_set,
			core::fragment::SingleResidueFragDataOP frag_type,
			core::Size ncopies = 1
	) const;

	void set_native( core::pose::PoseCOP native );

	PairingStatistics const& strand_pairing_stats() {
		return *strand_stats_;
	}

	inline TemplatesCOP get_self_ptr() const { return shared_from_this(); }
	inline TemplatesOP  get_self_ptr() { return shared_from_this(); }
	inline TemplatesCAP get_self_weak_ptr() const { return TemplatesCAP( shared_from_this() ); }
	inline TemplatesAP  get_self_weak_ptr() { return TemplatesAP( shared_from_this() ); }

private:
	void get_cst_list( TemplateList& cst_list, TemplateList& cull_list ) const;
	void scored_fragpick_list( TemplateList& frag_pick_list ) const;
	void _get_scored_list( TemplateList& cst_list, core::Size topN, core::Real wTopol, core::Real wExtern) const;

	bool good_;

	TemplateMap templates_;
	core::scoring::dssp::PairingsList pairings_;
	std::string target_sequence_;

	PairingStatisticsCOP strand_stats_;
	core::pose::PoseCOP native_;

	TemplateList cull_list_;
	TemplateList fragpick_list_;
	TemplateList helixjump_pick_list_;
};

} //abinitio
} //protocols

#endif
