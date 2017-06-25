// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project.h
/// @brief  Project class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#ifndef UI_TASK_PROJECT_H
#define UI_TASK_PROJECT_H

#include <ui/task/node.h>


namespace ui {
namespace task {

class Project : public Node
{
public:
	explicit Project(QUuid _project_id);

	std::string type() const override;


private:
	QUuid const project_id;
};

} // namespace task
} // namespace ui

#endif // UI_TASK_PROJECT_H
