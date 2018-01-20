// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project.fwd.h
/// @brief  Project class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#pragma once

#include <memory>

namespace ui {
namespace task {

class Project;

using ProjectSP  = std::shared_ptr< Project >;
using ProjectCSP = std::shared_ptr< Project const >;

using ProjectUP  = std::unique_ptr< Project >;

} // namespace task
} // namespace ui
