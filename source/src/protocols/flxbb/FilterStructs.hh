// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/flxbb/FilterStructs.hh
/// @brief filter structures
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_flxbb_FilterStructs_hh
#define INCLUDED_protocols_flxbb_FilterStructs_hh

// Unit Header
#include <protocols/flxbb/FilterStructs.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace flxbb {


////////////////////////////////////////////////////////////////////////////////////////
class FilterStructs : public utility::pointer::ReferenceCount {
public:

	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;


public:// constructor/destructor


	/// @brief default constructor
	FilterStructs();

	/// @brief value constructor
	FilterStructs( String name );

	/// @brief value constructor
	FilterStructs( String name, Size const ntrial );

	/// @brief value constructor
	FilterStructs( String name, Pose const & pose, Size const ntrial );

	/// @brief copy constructor
	FilterStructs( FilterStructs const & rval );

	/// @brief destructor
	~FilterStructs() override;


public:// virtural constructor


	/// @brief
	virtual FilterStructsOP clone() const;

	/// @brief
	virtual FilterStructsOP fresh_instance() const;


public:// virtual operations


	/// @brief
	virtual void apply( Pose const & ) = 0;

	/// @brief
	virtual void reset( Pose const & ) = 0;


public:// accessors


	/// @brief
	inline String name() const { return name_; }

	/// @brief
	inline bool filter_on() const { return filter_on_; }

	/// @brief
	inline Size current_trial() const { return current_trial_; }

	/// @brief
	PoseOP get_bestpose() const;


public:// mutators


	/// @brief
	void name( String const & name );

	/// @brief
	void set_ntrial( Size const ntrial );


protected://


	/// @brief
	void initialize( Pose const & pose );

	/// @brief
	void set_filter_off() { filter_on_ = false; }

	/// @brief
	void set_filter_on() { filter_on_ = true; }

	/// @brief
	void count_ntrial();

	/// @brief
	bool filter_is_over();

	/// @brief
	void set_bestpose( Pose const & pose );


private:

	String name_;
	bool filter_on_;
	Size ntrial_;
	Size current_trial_;
	PoseOP best_pose_;

};

////////////////////////////////////////////////////////////////////////////////////////
class FilterStructs_Packstat : public FilterStructs {
public:

	typedef FilterStructs Super;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;


public:

	/// @brief default constructor
	FilterStructs_Packstat( Size const ntrial=10 );

	/// @brief constructor
	FilterStructs_Packstat( Pose const & pose, Size const ntrial=10 );

	/// @brief copy constructor
	FilterStructs_Packstat( FilterStructs_Packstat const & rval );

	/// @brief destructor
	~FilterStructs_Packstat() override;


public:// virtural constructor


	/// @brief
	FilterStructsOP clone() const override;

	/// @brief
	FilterStructsOP fresh_instance() const override;


public:// virtual main operation


	/// @brief
	void apply( Pose const & ) override;

	/// @brief
	void reset( Pose const & ) override;


private:

	Real best_packscore_;

};

////////////////////////////////////////////////////////////////////////////////////////
class FilterStructs_TotalCharge : public FilterStructs {
public:

	typedef FilterStructs Super;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;


public:


	/// @brief constructor
	FilterStructs_TotalCharge( Size const ntrial=20 );

	/// @brief constructor
	FilterStructs_TotalCharge( Pose const & pose, Size const ntrial=20 );

	/// @brief copy constructor
	FilterStructs_TotalCharge( FilterStructs_TotalCharge const & rval );

	/// @brief destructor
	~FilterStructs_TotalCharge() override;


public:// virtural constructor


	/// @brief
	FilterStructsOP clone() const override;

	/// @brief
	FilterStructsOP fresh_instance() const override;


public:// virtual main operation


	/// @brief
	void apply( Pose const & ) override;

	/// @brief
	void reset( Pose const & ) override;


private:

	Size disallowed_value_;

};


} // flxbb
} // protocols


#endif
