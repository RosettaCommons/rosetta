#include "mainwindow.h"
#include <QApplication>

#include <protocols/init/init.hh>
#include <libxml/parser.h> // dummy include to make sure external excludes path is set

#include <ui/ui_lib_test.h>
#include <ui/task/project_view.h>

#include <QDebug>

int main(int argc, char *argv[])
{
    qDebug() << "ui::main: Calling Rosetta init()...";

    //protocols::init::init(argc, argv);
    //ui::ui_lib_test::test();

    QApplication a(argc, argv);

    //MainWindow w;
    //w.show();

    ui::task::ProjectView pv;
    pv.show();

    return a.exec();
}
