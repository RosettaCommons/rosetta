#ifndef INCLUDED_protocols_denovo_design_constraints_CoordinateConstraintGenerator_hh
#define INCLUDED_protocols_denovo_design_constraints_CoordinateConstraintGenerator_hh

// unit headers
#include <protocols/denovo_design/constraints/CoordinateConstraintGenerator.fwd.hh>
#include <protocols/moves/ConstraintGenerator.hh>

// protocol headers
#include <protocols/rosetta_scripts/SavePoseMover.hh>

// core headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace constraints {

// coordinate cst rcg class
class CoordinateConstraintGenerator : public protocols::moves::ConstraintGenerator {
public:
	CoordinateConstraintGenerator();
	virtual ~CoordinateConstraintGenerator();
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	virtual core::scoring::constraints::ConstraintCOPs
	generate_constraints( core::pose::Pose const & pose );

	void set_selector( core::select::residue_selector::ResidueSelectorCOP selector );

private:
	void add_constraints( core::scoring::constraints::ConstraintCOPs & csts, core::pose::Pose const & pose ) const;

	/// @brief creates a Ca coordinate constraint for residue resi
	core::scoring::constraints::ConstraintOP
	create_coordinate_cst(
		core::pose::Pose const & pose,
		core::Size const resi ) const;

private:
	core::pose::PoseCOP reference_pose_;
	core::select::residue_selector::ResidueSelectorCOP selector_;
	core::Real weight_;
};

}
}
}

#endif

