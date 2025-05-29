// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/SAScoreFilter.cc
/// @brief A filter to calculate the Novartis Sythetic Accesibility score using RDKit
/// @author Rocco Moretti (rmorettiase@gmail.com)

/// @ details This filter is a translation of Contrib/SA_Score/sascorer.py from RDKit
/// http://www.rdkit.org/
///
/// Copyright information for sascorer.py:
///
////////////////////////////////////////////////////////////////////////////////////////
///  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
///  All rights reserved.
///
/// Redistribution and use in source and binary forms, with or without
/// modification, are permitted provided that the following conditions are
/// met:
///
///     * Redistributions of source code must retain the above copyright
///       notice, this list of conditions and the following disclaimer.
///     * Redistributions in binary form must reproduce the above
///       copyright notice, this list of conditions and the following
///       disclaimer in the documentation and/or other materials provided
///       with the distribution.
///     * Neither the name of Novartis Institutes for BioMedical Research Inc.
///       nor the names of its contributors may be used to endorse or promote
///       products derived from this software without specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/////////////////////////////////////////////////////////////////////////////////////////
///
/// calculation of synthetic accessibility score as described in:
///
/// Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
/// Peter Ertl and Ansgar Schuffenhauer
/// Journal of Cheminformatics 1:8 (2009)
/// http://www.jcheminf.com/content/1/1/8
///
/// several small modifications to the original paper are included
/// particularly slightly different formula for marocyclic penalty
/// and taking into account also molecule symmetry (fingerprint density)
///
/// for a set of 10k diverse molecules the agreement between the original method
/// as implemented in PipelinePilot and this implementation is r2 = 0.97
///
/// peter ertl & greg landrum, september 2013

#include <protocols/drug_design/SAScoreFilter.hh>
#include <protocols/drug_design/SAScoreFilterCreator.hh>

#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>

#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/conformation/Residue.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <map>
#include <cmath>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RingInfo.h>
#include <rdkit/GraphMol/Fingerprints/MorganFingerprints.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/DataStructs/SparseIntVect.h>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.SAScoreFilter");

protocols::filters::FilterOP
SAScoreFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SAScoreFilter ); }

std::string
SAScoreFilterCreator::keyname() const { return SAScoreFilter::class_name(); }

void SAScoreFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SAScoreFilter::provide_xml_schema( xsd );
}

std::string SAScoreFilter::class_name() {
	return "SAScore";
}

void SAScoreFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "residue", xsct_refpose_enabled_residue_number, "Residue to calculate SAScore on" )
		+ XMLSchemaAttribute( "threshold", xsct_real, "Fail if the SAScore is greater than or equal to this value." );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Calculate the synthetic accessiblity score (of Ertl and Schuffenhauer doi:10.1186/1758-2946-1-8) for a residue.",
		attlist );
}

void
SAScoreFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	residue_ = tag->getOption<std::string>( "residue" );
	threshold_ = tag->getOption<core::Real>( "threshold" );

	TR << "SAScore filter on residue "<<residue_<< " with threshold " << threshold_ << std::endl;
}

bool
SAScoreFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const value( compute( pose ) );
	if ( value < threshold_ ) {
		return true;
	} else {
		return false;
	}
}

void
SAScoreFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size resnum( core::pose::parse_resnum( residue_, pose ) );
	if ( resnum == 0 ) {
		TR.Error << "Cannot computer SAScore for non-existent residue " << residue_ << std::endl;
		utility_exit_with_message("Cannot apply SAScore filter on non-existant residue!");
	}
	out << "SAScore for residue " << pose.residue(resnum).name() << " is " << compute(pose) << std::endl;
}

core::Real
SAScoreFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
SAScoreFilter::compute( core::pose::Pose const & pose ) const {
	core::Size resnum( core::pose::parse_resnum( residue_, pose ) );
	if ( resnum == 0 || resnum > pose.total_residue() ) {
		TR.Error << "Attempted to access residue " << residue_ << " in pose with " <<  pose.total_residue() << " residues. Failing filter. " << std::endl;
		utility_exit_with_message("Cannot apply SAScore filter on non-existant residue!");
	}

	core::chemical::MutableResidueTypeOP restype( utility::pointer::make_shared< core::chemical::MutableResidueType >(pose.residue(resnum).type()) );
	// Want to have neutralized, no-hydro molecules: this should be the default.
	core::chemical::rdkit::RestypeToRDMol to_converter(*restype);
	::RDKit::RWMOL_SPTR rdmol( to_converter.Mol() ); // Convert

	return calculate_rdkit( *rdmol );

}

core::Size
SAScoreFilter::num_spiro(RDKit::ROMol const & mol) const {
	using namespace RDKit;
	RingInfo *ri = mol.getRingInfo(); // Do NOT delete pointer
	std::vector< std::set< int > > arings;
	for ( VECT_INT_VECT::const_iterator itr( ri->atomRings().begin() ), itr_end( ri->atomRings().end() );
			itr != itr_end; ++itr ) {
		std::set< int > aset( itr->begin(), itr->end() );
		arings.push_back( aset );
	}
	std::set< int > spiros;
	for ( std::vector< std::set< int > >::const_iterator itr( arings.begin() ), itr_end( arings.end() );
			itr != itr_end; ++itr ) {
		for ( std::vector< std::set< int > >::const_iterator itr2( itr + 1 ); itr2 != itr_end; ++itr2 ) {
			std::vector<int> out( std::min( itr->size(), itr2->size() ) );
			std::vector<int>::iterator out_end;
			out_end = std::set_intersection( itr->begin(), itr->end(), itr2->begin(), itr2->end(), out.begin() );
			if ( (out_end - out.begin()) == 1 ) {
				spiros.insert( *out.begin() );
			}
		}
	}
	return spiros.size();
}

core::Size
SAScoreFilter::num_bridgeheads(RDKit::ROMol const & mol) const {
	using namespace RDKit;
	RingInfo *ri = mol.getRingInfo(); // Do NOT delete pointer

	// find bonds that are shared between rings that share at least 2 bonds:
	std::vector< std::set< int > > brings;
	for ( VECT_INT_VECT::const_iterator itr( ri->bondRings().begin() ), itr_end( ri->bondRings().end() );
			itr != itr_end; ++itr ) {
		std::set< int > bset( itr->begin(), itr->end() );
		brings.push_back( bset );
	}

	std::set< int > bridges;

	for ( std::vector< std::set< int > >::const_iterator itr( brings.begin() ), itr_end( brings.end() );
			itr != itr_end; ++itr ) {
		for ( std::vector< std::set< int > >::const_iterator itr2( itr + 1 ); itr2 != itr_end; ++itr2 ) {
			std::vector<int> out( std::min( itr->size(), itr2->size() ) );
			std::vector<int>::iterator out_end;
			out_end = std::set_intersection( itr->begin(), itr->end(), itr2->begin(), itr2->end(), out.begin() );
			if ( (out_end - out.begin()) > 1 ) {
				std::map< int, core::Size > atomCounts;

				for ( std::vector<int>::const_iterator out_itr( out.begin() ); out_itr != out_end; ++out_itr ) {
					Bond const * bond( mol.getBondWithIdx(*out_itr) );
					// std::map<int>::operator[] will auto create an entry of 0 if one does not already exist
					atomCounts[bond->getBeginAtomIdx()] +=1;
					atomCounts[bond->getEndAtomIdx()] +=1;
				}

				for ( std::map< int, core::Size >::const_iterator ac_itr( atomCounts.begin() ), ac_end( atomCounts.end() );
						ac_itr != ac_end; ++ac_itr ) {
					if ( ac_itr->second == 1 ) {
						bridges.insert(ac_itr->first);
					}
				}
			}
		}
	}
	return bridges.size();
}

core::Size
SAScoreFilter::num_chiral(RDKit::ROMol & mol, bool includeUnassigned ) const {
	using namespace RDKit;
	MolOps::assignStereochemistry(mol,false,true,false);
	core::Size nchiral(0);
	for ( ROMol::AtomIterator itr( mol.beginAtoms() ), itr_end( mol.endAtoms() ); itr != itr_end; ++itr ) {
		if ( (*itr)->hasProp("_CIPCode") ) {
			++nchiral;
		} else if ( includeUnassigned && (*itr)->hasProp("_ChiralityPossible") ) {
			++nchiral;
		}
	}
	return nchiral;
}

/// @details The caller takes ownership of the returned item
RDKit::SparseIntVect<boost::uint32_t> *
SAScoreFilter::get_morgan_fingerprint(const RDKit::ROMol &mol, int radius) const {

	//int nBits = -1;
	bool useChirality = false;
	bool useBondTypes = true;
	//bool useFeatures = false;
	bool useCounts = true;

	std::vector<boost::uint32_t> *invars = nullptr;
	// Ignore invariants & useFeatures
	std::vector<boost::uint32_t> *froms = nullptr;
	// Ignore fromAtoms
	RDKit::MorganFingerprints::BitInfoMap *bitInfoMap = nullptr;
	// Ignore bitInfo
	RDKit::SparseIntVect<boost::uint32_t> *res;
	// Assume nBits < 0
	res = RDKit::MorganFingerprints::getFingerprint(mol,
		static_cast<unsigned int>(radius),
		invars,froms,useChirality,
		useBondTypes,useCounts,false,bitInfoMap);
	if ( bitInfoMap ) {
		delete bitInfoMap;
	}
	if ( invars ) delete invars;
	if ( froms ) delete froms;
	return res;
}

core::Real
SAScoreFilter::calculate_rdkit(RDKit::ROMol & mol) const {
	using namespace RDKit;
	SAScoreData const & sas_data( *SAScoreData::get_instance() );

	// fragment score - 2 is the *radius* of the circular fingerprint
	RDKit::SparseIntVect<boost::uint32_t> *fp( get_morgan_fingerprint(mol, 2) ); // Delete me!
	typedef std::map<boost::uint32_t,int> StorageType;
	StorageType const & fps( fp->getNonzeroElements() );
	core::Real score1 = 0.0;
	core::Size nf = 0;
	for ( StorageType::const_iterator itr(fps.begin()), itr_end(fps.end()); itr != itr_end; ++itr ) {
		nf += itr->second;
		score1 += sas_data[ itr->first ] * itr->second;
	}
	score1 /= nf;

	// features score
	core::Real nAtoms = mol.getNumAtoms();
	core::Real nChiralCenters = num_chiral(mol,true);
	core::Real nBridgeheads = num_bridgeheads(mol);
	core::Real nSpiro = num_spiro(mol);
	core::Real nMacrocycles=0;
	VECT_INT_VECT const & atomRings( mol.getRingInfo()->atomRings() );
	for ( VECT_INT_VECT::const_iterator itr( atomRings.begin() ), itr_end( atomRings.end() );
			itr != itr_end; ++itr ) {
		if ( itr->size() > 8 ) {
			++nMacrocycles;
		}
	}

	core::Real sizePenalty = std::pow(nAtoms,1.005) - nAtoms;
	core::Real stereoPenalty = std::log10(nChiralCenters+1);
	core::Real spiroPenalty = std::log10(nSpiro+1);
	core::Real bridgePenalty = std::log10(nBridgeheads+1);
	// ---------------------------------------
	// This differs from the paper, which defines:
	// macrocyclePenalty = math.log10(nMacrocycles+1)
	// This form generates better results when 2 or more macrocycles are present
	core::Real macrocyclePenalty = 0.;
	if ( nMacrocycles > 0 ) {
		macrocyclePenalty = std::log10(2);
	}

	core::Real score2 = 0.0 -sizePenalty -stereoPenalty -spiroPenalty -bridgePenalty -macrocyclePenalty;

	// correction for the fingerprint density
	// not in the original publication, added in version 1.1
	// to make highly symmetrical molecules easier to synthetise
	core::Real score3 = 0.0;
	if ( nAtoms > fps.size() ) {
		score3 = std::log( nAtoms / fps.size() ) * 0.5;
	}

	core::Real sascore = score1 + score2 + score3;

	// need to transform "raw" value into scale between 1 and 10
	core::Real const min = -4.0;
	core::Real const max = 2.5;
	sascore = 11. - (sascore - min + 1) / (max - min) * 9.;
	// smooth the 10-end
	if ( sascore > 8.0 ) {
		sascore = 8. + std::log(sascore+1.-9.);
	}
	if ( sascore > 10.0 ) {
		sascore = 10.0;
	}
	if ( sascore < 1. ) {
		sascore = 1.0;
	}

	delete fp;
	return sascore;
}

///////////////// SAScoreData
SAScoreData::SAScoreData():
	default_(-4.0)
{
	load_data_from_file("fpscores.gz");
}

void
SAScoreData::load_data_from_file(std::string const & filename ) {
	//Right now, just load from database - in future may want some sort of local loading.
	fscores_.clear();

	utility::io::izstream input;
	if ( ! basic::database::open( input, "chemical/rdkit/" + filename ) ) {
		utility_exit_with_message("Cannot open SAScore database file.");
	}
	boost::uint32_t item;
	float value;
	std::string line;
	while ( getline( input, line ) ) {
		std::istringstream l( line );
		l >> value;
		l >> item;
		while ( l ) {
			fscores_[item] = value;
			l >> item;
		}
	}
	TR << "Loaded " << fscores_.size() << " fragment scores for SAScoring. (Should be 705292)" << std::endl;
}


float
SAScoreData::operator[] ( boost::uint32_t index ) const {
	fscores_t::const_iterator entry( fscores_.find(index) );
	if ( entry == fscores_.end() ) {
		return default_;
	} else {
		return entry->second;
	}
}

} // namespace drug_design
} // namespace protocols
