// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/FragData.hh
/// @brief  A fragment as list of SingleResidue Data
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_FragData_HH
#define INCLUDED_core_fragment_FragData_HH

// Unit Headers
#include <core/fragment/FragData.fwd.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.hh>

// C/C++ headers
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

typedef utility::vector1 < Size > PositionList;

class FragData : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< FragData >
{
	typedef utility::vector1 < SingleResidueFragDataOP > SRFD_List;

public:
	FragData () : valid_( false ), score_( 0.0 ) {};

	//@brief convience constructor to create FragData object that contains n SRFD of same type
	// argument SRFD will only be used for cloning
	FragData( SingleResidueFragDataOP, core::Size n);

	//@brief
	~FragData() override = default;

	virtual
	FragDataOP clone() const;

	/// self pointers
	inline FragDataCOP get_self_ptr() const { return shared_from_this(); }
	inline FragDataOP get_self_ptr() { return shared_from_this(); }
	//inline FragDataCAP get_self_weak_ptr() const { return FragDataCAP( shared_from_this() ); }
	//inline FragDataAP get_self_weak_ptr() { return FragDataAP( shared_from_this() ); }

	//might be overwritten to shortcut for continous changes
	virtual Size apply( kinematics::MoveMap const&, pose::Pose&, Size start, Size end ) const; // continous application to length residues
	virtual Size apply( kinematics::MoveMap const&, pose::Pose&, Frame const& ) const; //application to any set of residues specified by the Fr

	//without movemap --- just apply the full fragment
	virtual Size apply( pose::Pose&, Frame const& ) const; //application to any set of residues specified by the Fr
	virtual Size apply( pose::Pose&, Size start, Size end ) const; // continous application to length residues

	// @brief apply framents sec-struct info to ss-string
	virtual Size apply_ss( kinematics::MoveMap const&, std::string& ss, Frame const& ) const;

	// @brief check weather the dofs changed by this fragment are allowed according to movemap
	// return the number of allowed residues before a disabled one is met
	virtual Size is_applicable( kinematics::MoveMap const&, Size start, Size end ) const;
	virtual Size is_applicable( kinematics::MoveMap const&, Frame const& ) const;

	// set FragData from the given pose ( inverse of apply )
	virtual bool steal( pose::Pose const&, Size start, Size end ); // continous application to length residues
	virtual bool steal( pose::Pose const&, Frame const& ); //application to any set of residues specified by the Frame

	// copy FragData from the given FragData
	virtual void copy(FragData const& frag_data);

	// check if both frag_data contain the same number and types of SRFD
	bool is_compatible( FragData const& frag_data ) const;

	bool is_valid() const {
		return valid_;
	}

	Size size() const {
		return data_.size();
	}

	char secstruct( Size const pos) const {
		return data_[pos]->secstruct();
	}

	char sequence( Size const pos ) const {
		return data_[ pos ] -> sequence();
	}

	std::string secstruct() const {
		std::string str;
		for ( Size pos=1; pos<=size(); pos++ ) {
			str.push_back( secstruct( pos) );
		};
		return str;
	};

	std::string sequence() const {
		std::string str;
		for ( Size pos=1; pos<=size(); pos++ ) {
			str.push_back( sequence( pos) );
		};
		return str;
	};

	void set_sequence( Size const pos, char const setting ) {
		data_[ pos ] -> set_sequence( setting );
	}

	void set_residue( Size pos, SingleResidueFragDataOP res ) {
		if ( pos < data_.size() ) {
			data_.resize( pos );
		};
		data_[ pos ] = res;
	};

	void add_residue( SingleResidueFragDataOP res ) {
		data_.push_back ( res );
	};

	SingleResidueFragDataCOP get_residue( core::Size pos ) const {
		return data_[ pos ];
	}

	//@brief Generates a fragment referencing the residues start..stop, does not copy SRFD.
	virtual FragDataOP generate_sub_fragment( Size start, Size stop ) const;

	virtual void show( std::ostream &out ) const;
	virtual void show( std::ostream &out, Frame const& ) const;
	virtual void show_classic( std::ostream &out ) const;

	//@brief notify FragData that it now contains some proper numbers ( would like to have this protected )
	// but IO needs access to this.
	void set_valid( bool setting = true ) {
		valid_ = setting;
	}

	/// @brief Returns the PDB identifier if it has been specified, "no_pdb" otherwise.
	virtual std::string pdbid() const {
		return "no_pdb";
	}

	/// @brief Returns the chain if it has been specified, '_' otherwise.
	virtual char chain() const {
		return '_';
	}

	virtual
	Size pdbpos() const {
		return 0;
	}

	/// @brief Set a score value for this fragment
	void set_score( Real score ) {
		score_ = score;
	}

	/// @brief Returns the score for this fragment
	Real score() const {
		return score_;
	}

protected: // make private
	FragData( Size nr_res ):
		data_( nr_res ),
		valid_(false), // Does this make sense?
		score_(0) // Does this make sense?
	{}

	SRFD_List data_;

private:

	// this will be true if dofs are set ( via steal, or IO )
	bool valid_;

	Real score_;

};


/// @brief FragData that contains additional information
class AnnotatedFragData : public FragData {
	typedef FragData Parent;

public:
	AnnotatedFragData(const std::string& pdb_id, Size start_pos, char chain='_') {
		// initialize w/ known chain
		initialize(pdb_id, chain, start_pos);
	}

	AnnotatedFragData(const std::string& pdb_id, Size start_pos, const FragData& frag, char chain='_')
	: FragData(frag) {
		initialize(pdb_id, chain, start_pos);
	}

	FragDataOP clone() const override;

	void copy(FragData const& frag_data) override;

	//@brief Generates a fragment referencing a subset of the given fragment, does not copy SRFD.
	FragDataOP generate_sub_fragment( Size start, Size stop ) const override;

	std::string pdbid() const override  {
		return pdbid_;
	}

	Size pdbpos() const override {
		return startpos_;
	}

	char chain() const override {
		return chain_;
	}


private:
	/// @brief common initialization routine
	void initialize(const std::string& pdb_id, char chain, Size start_pos) {
		pdbid_ = pdb_id;
		chain_ = chain;
		startpos_ = start_pos;
	}

	std::string pdbid_;

	/// @brief 1-letter chain identifier or '' if it was not specified.
	char chain_;

	Size startpos_; //or list of indices
};


inline std::ostream& operator<<( std::ostream& out, FragData const& fr ) {
	fr.show( out );
	return out;
}

} // fragment
} // core

#endif
