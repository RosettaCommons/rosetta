// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author gsmurphy
#ifndef INCLUDED_devel_denovo_protein_design_SSClass_hh
#define INCLUDED_devel_denovo_protein_design_SSClass_hh
#include <devel/denovo_protein_design/SSClass.fwd.hh>
#include <core/types.hh>
// C++ Headers
//#include <iosfwd>
// Utility Headers

#include <utility/vector1_bool.hh>
#include <ostream>
#include <functional>

///////////////////////////////////////////////////////////////////////////////
namespace devel {
namespace denovo_protein_design {

/// single secondary structure element definition
class SS {
public:
	/// default constructor
	SS():
		start_(0),
		stop_(0),
		sstype_('L') // by default everything will be loop
	{}
	/// input constructor
	SS(
		 core::Size const start_in, core::Size const stop_in, char const sstype_in
			):
		start_( start_in ),
		stop_( stop_in ),
		sstype_( sstype_in )
	{}


	inline core::Size start() const { return start_; }
	inline core::Size stop() const { return stop_; }
	inline char sstype() const { return sstype_;}
  inline core::Size size() const { return stop_ - start_ + 1; }

	inline void set_start( core::Size input ) { start_ = input; }
	inline void set_stop( core::Size input  ) { stop_ = input; }
	inline void set_sstype( char input ) { sstype_ = input;}

	// not sure what operators to make

	friend std::ostream & operator<<( std::ostream & os, const SS & ss_ );

private:
	core::Size start_;
	core::Size stop_;
	char sstype_;
};


//////////////////////////////////////////////////////////////////////
inline std::ostream & operator<<( std::ostream & os, const SS & ss_ ) {
	// it's super-annoying that this returns a std::endl
	os << "SecondaryStructureElement " << ss_.start_ << " " << ss_.stop_ << " "
		 << ss_.sstype_;
	return os;
}


/// @brief used to sort SS elements by start-res
class SS_lt : public std::binary_function<double, double, bool> {
public:
	bool operator()(SS x, SS y) {
		return (x.start() < y.stop());
	}
};


///////////////////////////////////////////////////////////////////////////
// a list of SS elements
class SSs {

public:
	typedef utility::vector1< SS > SSList;
	typedef SSList::iterator iterator;
	typedef SSList::const_iterator const_iterator;


public:
	inline core::Size num_ss() const { return sss_.size(); }
	inline const_iterator begin() const { return sss_.begin(); }
	inline const_iterator end() const { return sss_.end(); }
	inline iterator v_begin() { return sss_.begin(); }
	inline iterator v_end() { return sss_.end(); }

	//constructor
	SSs(){};
	//copy constructor
	SSs( const SSs & src ):sss_(src.sss_) {}
	//operator
	SSs & operator =( SSs const & src ) {
		sss_ = src.sss_;
		return *this;
	}

	friend std::ostream & operator<<( std::ostream & os, const SSs & sss_ );
	// friend basic::VTracer & operator<<( basic::VTracer & os, const SSs & sss );


	void read_ss_file(
		std::string filename
	);

	void
	write_ss_to_file(
		std::string const & filename
	);

	void
	add_ss( SS ss_ );

	void
	add_ss(
		core::Size const start,
		core::Size const stop,
		char const sstype
	);

	void
	add_ss(
		const SSs::const_iterator & it
	);

	void
	add_ss(
		const SSs::iterator & it
	);

	void push_back( SS ss ) {
		add_ss( ss );
	}

	void
	delete_ss(
		core::Size const start,
		core::Size const stop
	);

	const_iterator one_random_ss_element() const;

	core::Size
	ss_size(
		core::Size const ss_num
	) const;

	core::Size
	ss_size() const;


	core::Size size() const {
	 		return sss_.size();
	}


	// write size() function but it will be the actual number of residues in the secondary structure element
	// write number_elements - this will be the number of secondary structure elements in the list

	void clear();

 	SSList const& sss() const { return sss_; }

	inline
	const SS&
	operator[] ( core::Size const i ) const
	{
		return sss_[i];
	}

	inline
	SS&
	operator[] ( core::Size const i )
	{
		return sss_[i];
	}

private:
	SSList sss_;
	SS ss_;
};


} //namespace denovo_protein_design
} //namespace devel

#endif //INCLUDED_protocols_loops_SSClass_HH
