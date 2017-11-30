#include "mainwindow.h"
#include <QApplication>

#include <protocols/init/init.hh>
#include <libxml/parser.h> // dummy include to make sure external excludes path is set

#include <ui/ui_lib_test.h>
#include <ui/task/project_view.h>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

#include <QSettings>

#include <QDebug>

int main(int argc, char *argv[])
{
    //qDebug() << "ui::main: Calling Rosetta init()...";

    protocols::init::init(argc, argv);

	{ // creating dummy pose object so later we can create Pose in async mode
		core::pose::Pose p;
		core::import_pose::pose_from_pdbstring(p, "ATOM     17  N   ILE A   1      16.327  47.509  23.466  1.00  0.00\n");
	}

    //ui::ui_lib_test::test();

    QApplication a(argc, argv);

	QCoreApplication::setOrganizationName("RosettaCommons");
    QCoreApplication::setOrganizationDomain("RosettaCommons.org");
    QCoreApplication::setApplicationName("Workbench");

	QSettings::setDefaultFormat(QSettings::IniFormat); // comment this out later, when config editor functionality is in place

    //MainWindow w;
    //w.show();

    ui::task::ProjectView pv;
    pv.show();

    return a.exec();
}
