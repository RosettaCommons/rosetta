// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/apps/pilot/awatkins/rna_denovo/main.cpp
/// @brief  Main function for the stepwise interface project.
/// @author Andrew Watkins (amw579@stanford.edu)

#include "mainwindow.h"
#include <QApplication>

#include <protocols/init/init.hh>
#include <libxml/parser.h> // dummy include to make sure external excludes path is set

#include <QDebug>

#if defined(MAC) || defined(__APPLE__) || defined (__OSX__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int main(int argc, char *argv[])
{
    qDebug() << "bundle_gui::main: Calling Rosetta init()...";

	glutInit(&argc, argv);
    protocols::init::init(argc, argv);

    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
