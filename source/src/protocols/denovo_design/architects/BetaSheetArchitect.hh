// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/BetaSheetArchitect.hh
/// @brief Architect that creates a beta sheet
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_architects_BetaSheetArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_BetaSheetArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/BetaSheetArchitect.fwd.hh>
#include <protocols/denovo_design/architects/StructureArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/SheetDB.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueVector.fwd.hh>

// Utility headers
#include <utility/excn/EXCN_Base.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace denovo_design {
namespace architects {

///@brief Architect that creates a beta sheet
class BetaSheetArchitect : public protocols::denovo_design::architects::DeNovoArchitect {
public:
	typedef protocols::denovo_design::architects::DeNovoArchitect DeNovoArchitect;
	typedef protocols::denovo_design::architects::DeNovoArchitectOP DeNovoArchitectOP;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataOP StructureDataOP;
	typedef utility::vector1< components::StructureDataCOP > StructureDataCOPs;
	typedef utility::vector1< StrandArchitectOP > StrandArchitectOPs;
	typedef std::map< std::string, core::Size > StrandExtensionsMap;

	typedef components::StrandOrientation StrandOrientation;
	typedef components::StrandOrientations StrandOrientations;
	typedef components::RegisterShift RegisterShift;
	typedef components::RegisterShifts RegisterShifts;
	typedef core::select::residue_selector::ResidueVector ResidueVector;

	typedef std::vector< components::StrandPairingCOP > PairingsInfo;
	typedef std::vector< PairingsInfo >  PairingsInfoVector;

public:
	BetaSheetArchitect( std::string const & id_value );

	virtual ~BetaSheetArchitect();

	static std::string
	class_name() { return "BetaSheetArchitect"; }

	virtual std::string
	type() const;

	DeNovoArchitectOP
	clone() const;

	virtual StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	Lengths
	retrieve_lengths( components::StructureData const & perm ) const;

	StrandOrientations
	retrieve_orientations( components::StructureData const & perm ) const;

	RegisterShifts
	retrieve_register_shifts( components::StructureData const & perm ) const;

	void
	set_strand_extensions( std::string const & extensions_str );

	void
	add_strand_extension( std::string const & strand_name, core::Size const length );

private:
	void
	setup_strand_pairings();

	void
	add_strand( StrandArchitect const & strand );

	/// @brief generates and stores a vector of permutations based on strands
	void
	enumerate_permutations();

	/// @brief merges a list of permutations
	StructureDataOP
	combine_permutations( components::StructureDataCOPs const & chain ) const;

	/// @brief combines the given set of permutations with the current set
	void
	combine_permutations_rec(
		components::StructureDataCOPs const & chain,
		utility::vector1< components::StructureDataCOPs > const & plist );

	/// @brief modifies/stores data into a permutation and adds it
	void
	modify_and_add_permutation( StructureData const & perm );

	/// @brief checks permutations
	/// @throws EXCN_PreFilterFailed if something goes wrong
	void
	check_permutation( StructureData const & perm ) const;

	core::Size
	extension_length( std::string const & strand ) const;

private:
	void
	needs_update();

	void
	store_sheet_idx( StructureData & sd, core::Size const sheet_idx ) const;

	/// @brief set allowed register shifts from a string
	void
	add_register_shifts( std::string const & val );

	void
	add_orientations( std::string const & orientations_str );

	components::StructureDataCOPs
	filter_permutations( components::StructureDataCOPs const & perms ) const;

	components::StructureDataCOPs
	add_pairings( components::StructureDataCOPs const & perms ) const;

private:
	StructureDataCOPs permutations_;
	StrandArchitectOPs strands_;
	utility::vector1< StrandOrientations > orientations_;
	utility::vector1< RegisterShifts > shifts_;
	StrandExtensionsMap extensions_;
	components::SheetDBOP sheetdb_;
	bool use_sheetdb_;
	bool updated_;
};

class EXCN_PreFilterFailed : public utility::excn::EXCN_Base {
public:
	EXCN_PreFilterFailed( std::string const & msg ) :
		utility::excn::EXCN_Base(), msg_( msg )
	{}

	virtual void
	show( std::ostream & os ) const { os << msg_ << std::endl; }

private:
	std::string const msg_;
};

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_BetaSheetArchitect_hh

