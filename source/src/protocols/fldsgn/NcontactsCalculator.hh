// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/NcontactsCalculator.hh
/// @brief header file for NcontactsCalculator class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_NcontactsCalculator_hh
#define INCLUDED_protocols_fldsgn_NcontactsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

#include <utility/vector1.hh>


// Utility headers

//// C++ headers

namespace protocols {
namespace fldsgn {

class NcontactsCalculator : public core::pose::metrics::StructureDependentCalculator {
public:


	typedef core::pose::metrics::StructureDependentCalculator Super;
	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Atom Atom;
	typedef core::conformation::Residue Residue;
	typedef core::pose::metrics::PoseMetricCalculatorOP PoseMetricCalculatorOP;
  typedef basic::MetricValueBase MetricValueBase;
	//typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	//typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;


public:// constructor/destructor


	/// @brief default constructor
	NcontactsCalculator();

	/// @brief default constructor
	NcontactsCalculator( Real const condist, Size const isep_sep );

	/// @brief copy constructor
	NcontactsCalculator( NcontactsCalculator const & rval );

	/// @brief destructor
	virtual ~NcontactsCalculator();


public:// virtual constructor


	/// @brief make clone
  PoseMetricCalculatorOP clone() const { return PoseMetricCalculatorOP( new NcontactsCalculator( *this ) ); }


public:// mutator

	/// @brief set contact distance
	void contact_distance( Real const r ) { condist_ = r; }

	/// @brief ignore loops for calculation
	void ignore_loops( bool const b ) { ignore_loops_ = b; };

	/// @brief ignore residue pairs of which residue belong to same ss element
	void ignore_same_sselement( bool const b ) { ignore_same_sselement_ = b; }

	/// @brief ignore residue pairs of which residue belong to same beta sheet ( default true )
	void ignore_same_sheet( bool const b ) { ignore_same_sheet_ = b; }

	/// @brief use only calpha for calculation
	void use_only_calpha( bool const b ) { use_only_calpha_ = b; }

protected:


  virtual void lookup( String const & key, MetricValueBase * valptr ) const;
  virtual std::string print( String const & key ) const;
  virtual void recompute( Pose const & this_pose );


private:

	/// @brief distact used for juding contact pair ( default 6.0 )
	Real condist_;

	/// @brief residue pairs of i < i+isep_residue_ are used for counting #countacts ( default 4 )
	Size isep_residue_;

	/// @brief ignore loops for calculation  ( default false )
	bool ignore_loops_;

	/// @brief ignore residue pairs of which residue belong to same ss element ( default false )
	bool ignore_same_sselement_;

	/// @brief ignore residue pairs of which residue belong to same beta sheet ( default false )
	bool ignore_same_sheet_;

	/// @brief use only calpha atoms for calculation ( default false )
	bool use_only_calpha_;

	/// @brief  #atom-contacts among all heavy atoms
	Real nc_allatm_;

  /// @brief  #atom-contacts among hydrophobic heavy atoms
	Real nc_hpatm_;

	/// @brief #atom-contacts among heavy atoms of sidechains of hydrophobic residues
	Real nc_hpres_;

	/// @brief
	Real ss_entrpy_;

}; //NontactCalculator


} // ns toolbox
} // ns protocols

#endif
