// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/LimitAromaChi2Operation.hh
/// @brief  eliminate aromatic rotamers, of which chi2 are around 0, 180 degree.
/// @detail Chi2=0, 180 rotamers of aromatic residues ( PHE, TYR, HIS ) are not observed in nature very much,
/// however Rosetta really like them. This is really pathology. For design purpose, we don't need them actually.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_toolbox_task_operations_LimitAromaChi2Operation_hh
#define INCLUDED_protocols_toolbox_task_operations_LimitAromaChi2Operation_hh


// Unit Headers
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class LimitAromaChi2_RotamerSetOperation : public core::pack::rotamer_set::RotamerSetOperation
{
public:


	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::rotamer_set::RotamerSet RotamerSet;
	typedef core::pack::rotamer_set::Rotamers Rotamers;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::graph::GraphCOP GraphCOP;
	typedef core::pack::rotamer_set::RotamerSetOperationOP RotamerSetOperationOP;


public:


	LimitAromaChi2_RotamerSetOperation();

	LimitAromaChi2_RotamerSetOperation( Real const chi2max, Real const chi2min );

	virtual ~LimitAromaChi2_RotamerSetOperation();

	virtual
	RotamerSetOperationOP
	clone() const;

	virtual
	void
	alter_rotamer_set(
		Pose const &,
		ScoreFunction const &,
		PackerTask const & ptask,
		GraphCOP,
		RotamerSet & rotamer_set
	);


public:

	/// @brief max chi2 for picking rotamers of YFH
	void chi2max( Real const r )
	{
		chi2max_ = r;
	}

	/// @brief min chi2 for picking rotamers of YFH
	void chi2min( Real const r )
	{
		chi2min_ = r;
	}

	/// @brief include TRP ?
	void include_trp( bool const b )
	{
		include_trp_ = b;
	}

private: // data

	/// @brief max chi2 for picking rotamers of YFH
	Real chi2max_;

	/// @brief min chi2 for picking rotamers of YFH
	Real chi2min_;

	/// @brief include TRP ? ( default false )
	bool include_trp_;

};


//////////////////////////////////////////////////////////////////////////////////////////////
class LimitAromaChi2Operation : public core::pack::task::operation::TaskOperation {
public:


	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;


public:


	/// @brief default constructor
	LimitAromaChi2Operation();

	/// @brief destructor
	virtual ~LimitAromaChi2Operation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;

public:


	/// @brief max chi2 for picking rotamers of YFH
	Real chi2max() const
	{
		return chi2max_;
	}

	/// @brief min chi2 for picking rotamers of YFH
	Real chi2min() const
	{
		return chi2min_;
	}

	/// @brief include TRP ?
	bool include_trp() const
	{
		return include_trp_;
	}


public:


	/// @brief max chi2 for picking rotamers of YFH
	void chi2max( Real const r )
	{
		chi2max_ = r;
	}

	/// @brief min chi2 for picking rotamers of YFH
	void chi2min( Real const r )
	{
		chi2min_ = r;
	}

	/// @brief include TRP ?
	void include_trp( bool const b )
	{
		include_trp_ = b;
	}


public:


	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;


public:


	void parse_tag( TagCOP tag , DataMap & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "LimitAromaChi2"; }


private: // data


	/// @brief max chi2 for picking rotamers of YFH
	Real chi2max_;

	/// @brief min chi2 for picking rotamers of YFH
	Real chi2min_;

	/// @brief include TRP ? ( default false )
	bool include_trp_;

};


} // TaskOperations
} // toolbox
} // protocols


#endif
