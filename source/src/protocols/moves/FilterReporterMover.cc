// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  A FilterMover that also calls report() on apply()
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Jul. 30, 2014

// Unit headers
#include <protocols/moves/FilterReporterMover.hh>

namespace protocols {
namespace moves {

FilterReporterMover::FilterReporterMover() :
	Mover(),
	mover_( nullptr ),
	filter_( nullptr ),
	out_( /*streambuf=*/ nullptr )
{
	filter_mover_ = protocols::moves::FilterMoverOP(
		new protocols::moves::FilterMover() );
}

FilterReporterMover::FilterReporterMover( protocols::moves::MoverOP const & mover,
	protocols::filters::FilterOP const & filter, core::Size const max_tries,
	std::ostream & out, protocols::moves::MoverStatus const mover_status ) :
	Mover( mover->type() ),
	mover_( mover ),
	filter_( filter ),
	out_( /*streambuf=*/ nullptr ) // default ctor is protected
{
	filter_mover_ = protocols::moves::FilterMoverOP(
		new protocols::moves::FilterMover( mover, filter, max_tries, mover_status ) );
	set_out( out );
}

FilterReporterMover::FilterReporterMover( FilterReporterMover const & rhs ) :
	Mover( rhs ),
	out_( /*streambuf=*/ nullptr ) {
	FilterReporterMover::assign( *this, rhs );
}

FilterReporterMover::~FilterReporterMover() {
	filter_mover_ = nullptr; // it is an OP, so should be destroyed here once assigning
}

void
FilterReporterMover::apply( core::pose::Pose & pose ) {
	filter_mover_->apply( pose );

	if ( out_.rdbuf() != nullptr ) {
		filter_->report( out_, pose );
	}
}

FilterReporterMover & FilterReporterMover::operator=( FilterReporterMover const & rhs ) {
	if ( this == &rhs ) return *this;

	Mover::operator=(rhs);
	FilterReporterMover::assign( *this, rhs );
	return *this;
}

void
FilterReporterMover::assign(FilterReporterMover & lhs, FilterReporterMover const & rhs) {
	lhs.mover_ = rhs.mover_;
	lhs.filter_ = rhs.filter_;
	lhs.filter_mover_ = rhs.filter_mover_;
	lhs.set_out( rhs.out_ );
}


protocols::moves::MoverOP
FilterReporterMover::clone() const {
	return protocols::moves::MoverOP( new FilterReporterMover(*this) );
}

std::string
FilterReporterMover::get_name() const {
	return "FilterReporterMover";
}

protocols::moves::MoverOP
FilterReporterMover::fresh_instance() const {
	return protocols::moves::MoverOP( new FilterReporterMover );
}

protocols::moves::MoverOP
FilterReporterMover::get_mover() const {
	return mover_;
}

void
FilterReporterMover::set_mover( protocols::moves::MoverOP const & mover ) {
	mover_ = mover; filter_mover_->set_mover(mover);
}

protocols::filters::FilterOP
FilterReporterMover::get_filter() const {
	return filter_;
}

void
FilterReporterMover::set_filter( protocols::filters::FilterOP const & filter ) {
	filter_ = filter;
	filter_mover_->add_filter(filter);
}

std::ostream & FilterReporterMover::get_out/*ta here*/() {
	return out_;
}

void
FilterReporterMover::set_out( std::ostream const & out ) {
	// this seems to be the way to make two ostreams point to the
	// same stream.
	// e.g. http://infocenter.arm.com/help/index.jsp?topic=/com.arm.doc.dui0729a/stdug/loc-io/13-2.htm
	// see also, https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.3/a00869.html
	out_.copyfmt(out); // copy everything except rdstate and rdbuf
	out_.clear(out.rdstate()); // copy rdstate
	out_.rdbuf(out.rdbuf()); // share the buffer
}

} // namespace peptide_deriver
} // namespace protocols
