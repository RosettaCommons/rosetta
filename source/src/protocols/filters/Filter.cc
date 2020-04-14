// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/Filter.cc
/// @brief
/// @details
///   Contains currently:
///
///
/// @author Florian Richter, Sarel Fleishman (sarelf@uw.edu)

// Unit Headers
#include <protocols/filters/Filter.hh>


#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>

// Project Headers
#include <basic/Tracer.hh>

#include <utility>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>


static basic::Tracer TR( "protocols.filters.Filter" );

namespace protocols {
namespace filters {

using namespace core;
typedef std::pair< std::string const, FilterCOP > StringFilter_pair;
using TagCOP = utility::tag::TagCOP;

FilterCollection::~FilterCollection() = default;

bool
FilterCollection::apply( core::pose::Pose const & pose ) const
{
	for ( protocols::filters::FilterCOP filter : filters_ ) {
		if ( ! filter->apply( pose ) ) {
			return false;
		}
	}

	return true;
}

void
FilterCollection::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	for ( protocols::filters::FilterCOP filter : filters_ ) {
		filter->report( out, pose );
	}
}

Filter::Filter()
: utility::VirtualBase(),
	type_( "UNDEFINED TYPE" ),
	scorename_("defaultscorename")
{}

Filter::Filter( std::string const & type )
: utility::VirtualBase(),
	type_( type ),
	scorename_("defaultscorename")
{}

Filter::Filter( Filter const & init )
: utility::VirtualBase(),
	type_( init.type_ ),
	user_defined_name_( init.user_defined_name_ ),
	scorename_("defaultscorename")

{}

Filter::~Filter() = default;

void
Filter::parse_my_tag(
	TagCOP,
	basic::datacache::DataMap &
) {
	TR.Warning << "The parse_my_tag method has been invoked for filter " << name() << " but it hasn't been defined. Are you sure this is appropriate?" << std::endl;
}

core::Real Filter::score( core::pose::Pose & pose ) {
	core::Real score = report_sm( pose );
	core::pose::setPoseExtraScore( pose, scorename_, score );
	return score;
}

/// @brief Does this filter provide information about how to cite it?
/// @details Defaults to false.  Derived classes may override this to provide citation info.  If set to
/// true, the provide_citation_info() override should also be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
Filter::filter_provides_citation_info() const {
	return false;
}

/// @brief Provide the citation.
/// @returns A vector of citation collections.  This allows the filter to provide citations for
/// itself and for any modules that it invokes.
/// @details The default implementation of this function provides an empty vector.  It may be
/// overriden by filters wishing to provide citation information.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
Filter::provide_citation_info() const {
	return utility::vector1< basic::citation_manager::CitationCollectionCOP >();
}

/// @brief Does this filter indicate that it is unpublished (and, by extension, that the author should be
/// included in publications resulting from it)?
/// @details Defaults to false.  Derived classes may override this to provide authorship info.  If set to
/// true, the provide_authorship_info_for_unpublished() override should also be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
Filter::filter_is_unpublished() const {
	return false;
}

/// @brief Provide a list of authors and their e-mail addresses, as strings.
/// @details The citation should be in MLA format.  Only the citation needs to be in the output
/// string. No newlines should be included.
/// @returns A list of pairs of (author, e-mail address).  Empty list if not unpublished.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
Filter::provide_authorship_info_for_unpublished() const {
	return utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >();
}

} // filters
} // protocols
