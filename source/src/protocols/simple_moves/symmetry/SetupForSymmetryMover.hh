// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_simple_moves_symmetry_SetupForSymmetryMover_hh
#define INCLUDED_protocols_simple_moves_symmetry_SetupForSymmetryMover_hh

// Unit headers
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <utility/vector1.hh>


// Utility Headers

namespace protocols {
namespace simple_moves {
namespace symmetry {

///////////////////////////////////////////////////////////////////////////////

class SetupForSymmetryMover : public protocols::moves::Mover
{
public:

	// default constructor
	SetupForSymmetryMover();

	SetupForSymmetryMover( core::conformation::symmetry::SymmDataOP symmdata );

	SetupForSymmetryMover( std::string const & );

	~SetupForSymmetryMover();

	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new SetupForSymmetryMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose ) override;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	// XRW TEMP  virtual std::string get_name() const;

	// setter
	void slide_into_contact(bool val) { slide_ = val; }

	/// @brief Sets whether or not the input asymmetric pose's datacache should be copied into
	///        the new symmetric pose.
	/// @param[in] preserve_cache If true, input pose's datacache is copied into new symmetric pose
	///                           If false, input pose's datacache is cleared (default = false)
	void
	set_preserve_datacache( bool const preserve_cache );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void process_symmdef_file(std::string tag);

	/// @brief   constructs a symmetric pose with a symmetric conformation and energies object.
	/// @details Calls core::pose::make_symmetric_pose().  If preserve_datacache is set, this
	///          also copies the datacache into the new symmetric pose.
	/// @param[in,out] pose Input asymmetric pose, output symmetric pose
	void
	make_symmetric_pose( core::pose::Pose & pose ) const;

private:
	bool slide_;
	bool cryst1_; //use cryst1 line
	bool preserve_datacache_;
	core::conformation::symmetry::SymmDataOP symmdef_;
};

///////////////

class ExtractAsymmetricUnitMover : public protocols::moves::Mover
{
public:

	// default constructor
	ExtractAsymmetricUnitMover();

	~ExtractAsymmetricUnitMover();

	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new ExtractAsymmetricUnitMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose ) override;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	// XRW TEMP  virtual std::string get_name() const;

public:
	/// @brief if keep_virtual_residues is true, virtual residues will remain in the pose, otherwise
	///        they will be removed
	void
	set_keep_virtual_residues( bool const keep_virt );

	/// @brief If true, residues with aa() == aa_unk will be kept in the asymmetric unit. If false,
	///        residues of aa type aa_unk will be ignored in the conversion and left out of the
	///        asymmetric unit.
	/// @param[in] keep_unk Desired value for keep_unknown_aas (default=false)
	/// @details If there are NCAAs in the pose, this must be set to false, or the NCAAs will be
	///          ignored.  The keep_unknown_aas defaults to false for historical reasons.
	void
	set_keep_unknown_aas( bool const keep_unk );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// @brief if keep_virtual_residues is true, virtual residues will remain in the pose, otherwise
	///        they will be removed
	bool keep_virtual_residues_;

	/// @brief If false, unknown residue types will be ignored
	///        Must be true if there are NCAAs in the pose, otherwise they will be ignored.
	bool keep_unknown_aas_;
};

class ExtractAsymmetricPoseMover : public protocols::moves::Mover
{
public:

	// default constructor
	ExtractAsymmetricPoseMover();

	~ExtractAsymmetricPoseMover();

	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new ExtractAsymmetricPoseMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose ) override;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	// XRW TEMP  virtual std::string get_name() const;

public:
	/// @brief if true, clears symmetry_definition option from pose.
	void
	clear_sym_def( bool const clear_sym_def );
	
	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	
private:
	/// @brief if true, clears symmetry_definition option from pose.
	bool clear_sym_def_;
};


}
} // symmetric_docking
} // rosetta
#endif
