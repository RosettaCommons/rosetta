// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/dna/WatsonCrickRotamerCouplings.hh>
#include <protocols/dna/WatsonCrickRotamerCouplingsCreator.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/util.hh> // find_basepairs

#include <core/conformation/ResidueMatcher.hh> // WatsonCrickResidueMatcher
#include <core/chemical/ResidueType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/Tracer.hh>

// option key includes

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif

#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


namespace protocols {
namespace dna {

using namespace core;
using namespace chemical;
using namespace conformation;
using namespace basic::options;
using namespace pack;
using namespace rotamer_set;
using namespace task;
using namespace operation;

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static THREAD_LOCAL basic::Tracer TR( "protocols.dna.WatsonCrickRotamerCouplings", t_info );

TaskOperationOP WatsonCrickRotamerCouplingsCreator::create_task_operation() const
{
	return TaskOperationOP( new WatsonCrickRotamerCouplings );
}

void WatsonCrickRotamerCouplingsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	WatsonCrickRotamerCouplings::provide_xml_schema( xsd );
}

std::string WatsonCrickRotamerCouplingsCreator::keyname() const
{
	return WatsonCrickRotamerCouplings::keyname();
}

WatsonCrickRotamerCouplings::~WatsonCrickRotamerCouplings() {}

TaskOperationOP WatsonCrickRotamerCouplings::clone() const
{
	return TaskOperationOP( new WatsonCrickRotamerCouplings( *this ) );
}

void
WatsonCrickRotamerCouplings::parse_tag( TagCOP /*tag*/ , DataMap & )
{}

void WatsonCrickRotamerCouplings::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

void
WatsonCrickRotamerCouplings::apply(
	Pose const & pose,
	PackerTask & ptask
) const
{
	// find DNA chains
	DnaChains dna_chains;
	find_basepairs( pose, dna_chains );
	Size const nres( pose.total_residue() );
	// setup residue couplings
	RotamerCouplingsOP couplings( new RotamerCouplings );
	couplings->resize( nres );
	for ( DnaPositions::const_iterator it( dna_chains.begin() ); it != dna_chains.end(); ++it ) {
		//mjo commenting out 'resid' because it is unused and causes a warning
		//Size const resid( it->first );
		DnaPosition const & dnapos( it->second );
		if ( ! dnapos.paired() ) continue;
		Size top( dnapos.top() );
		Size bot( dnapos.bottom() );
		TR << "base pair " << pose.pdb_info()->chain(top) << "." << pose.pdb_info()->number(top) << "."
			<< dna_full_name3( pose.residue_type(top).name3() ) << " - " << pose.pdb_info()->chain(bot)
			<< "." << pose.pdb_info()->number(bot) << "." << dna_full_name3( pose.residue_type(bot).name3() )
			<< std::endl;
		(*couplings)[ top ].first = bot;
		(*couplings)[ top ].second = ResidueMatcherCOP( ResidueMatcherOP( new conformation::WatsonCrickResidueMatcher() ) );
		(*couplings)[ bot ].first = top;
		(*couplings)[ bot ].second = ResidueMatcherCOP( ResidueMatcherOP( new conformation::WatsonCrickResidueMatcher() ) );
	}
	ptask.rotamer_couplings( couplings );
}

} // namespace dna
} // namespace protocols

