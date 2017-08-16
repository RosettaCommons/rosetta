/// @file
/// @brief
/// @author Zhe Zhang
#ifndef INCLUDED_devel_replica_docking_AddEncounterConstraintMover_hh
#define INCLUDED_devel_replica_docking_AddEncounterConstraintMover_hh

#include <devel/replica_docking/AddEncounterConstraintMover.fwd.hh>
#include <protocols/docking/types.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/scoring/InterfaceInfo.fwd.hh>


namespace devel {
namespace replica_docking {

class AddEncounterConstraintMover : public protocols::moves::Mover {

public:

	AddEncounterConstraintMover();

	~AddEncounterConstraintMover() override;

	protocols::moves::MoverOP clone() const override;

	void apply( core::pose::Pose & pose ) override;

	core::scoring::constraints::AtomPairConstraintCOP generate_encounter_cst( core::pose::Pose & pose );

	//  void get_interface_info( core::pose::Pose & pose ) const;
	protocols::scoring::InterfaceInfo const & interface_from_pose( core::pose::Pose const & ) const;


	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Real gap_;
	core::Size interface_jump_;
	core::scoring::constraints::AtomPairConstraintCOP cst_;

};

}
}

#endif
