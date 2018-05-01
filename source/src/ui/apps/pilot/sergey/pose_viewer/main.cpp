#include "mainwindow.h"
#include <QApplication>

#include <protocols/init/init.hh>
#include <libxml/parser.h> // dummy include to make sure external excludes path is set
#include <ui/ui_lib_test.h>


#include <QDebug>

int main(int argc, char *argv[])
{
    qDebug() << "test::main: Calling Rosetta init()...";

    protocols::init::init(argc, argv);
    ui::ui_lib_test::test();

    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
