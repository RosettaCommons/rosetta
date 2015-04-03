// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/kinmatch/FunGroupTK.hh
/// @brief  ToolKit for rapid IK based functional group placement
/// @author Will Sheffler will@sheffler.me


#ifndef INCLUDED_protocols_kinmatch_FunGroupTK_hh
#define INCLUDED_protocols_kinmatch_FunGroupTK_hh

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <protocols/scoring/ImplicitFastClashCheck.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/xyzVector.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <map>


namespace protocols {
namespace kinmatch {

enum KRSQueryType {
    CEN,
    CEN_ANG,
    CEN_AXS,
    CEN_AXS_ORI
};

struct KRSQuery {
    KRSQueryType type;
    numeric::xyzVector<core::Real> cen;
    numeric::xyzVector<core::Real> axs;
    numeric::xyzVector<core::Real> ori;
    core::Real disth;
    core::Real angth;
    core::Real clash;
    KRSQuery(
        KRSQueryType typ
        );
    KRSQuery(
        KRSQueryType typ,
        numeric::xyzVector<core::Real> c,
        numeric::xyzVector<core::Real> a,
        numeric::xyzVector<core::Real> o,
        core::Real dt=0.5,
        core::Real at=0.175,
        core::Real clsh = 2.8*2.8
        );
    KRSQuery(
        KRSQueryType typ,
        numeric::xyzVector<core::Real> c,
        numeric::xyzVector<core::Real> a,
        core::Real dt=0.5,
        core::Real at=0.175,
        core::Real clsh = 2.8*2.8
        );
};

class FunGroupTK : public utility::pointer::ReferenceCount {
protected:
    core::pose::PoseCOP pose_;
    utility::vector1<Size> const pos_;
    protocols::scoring::ImplicitFastClashCheckCOP ifc_;
    std::map<std::string,utility::vector1<core::kinematics::Stub> > stb_;
    std::map<std::string,utility::vector1<core::conformation::ResidueOP> > rsd_;
    core::chemical::ResidueTypeSetCAP frs_;
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FunGroupTK();
    FunGroupTK(
        core::pose::Pose & p_in,
        utility::vector1<Size> & pos
        );

    protocols::scoring::ImplicitFastClashCheck const &
    ifc() const;

    virtual
    void
    place_c(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const = 0;

    virtual
    void
    place_h(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const = 0;
    virtual
    void
    place_d(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const = 0;
};

class BruteFunGroupTK : public FunGroupTK {
public:
    BruteFunGroupTK(
        core::pose::Pose & p_in,
        utility::vector1<Size> pos=utility::vector1<Size>()
        );

    virtual
    void
    place_c(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const;

    virtual
    void
    place_h(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const;

    virtual
    void
    place_d(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const;

};

class KinFunGroupTK : public FunGroupTK {
public:
    KinFunGroupTK(
        core::pose::Pose & p_in,
        utility::vector1<Size> pos=utility::vector1<Size>()
        );


    virtual
    void
    place_c(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits_out
        ) const;

    virtual
    void
    place_h(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const;

    virtual
    void
    place_d(
        KRSQuery const & q,
        core::conformation::Residue const & qrsd,
        utility::vector1<core::conformation::ResidueOP> & hits
        ) const;

};

} // namespace kinmatch
} // namespace protocols

#endif
