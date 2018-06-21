#include <ui/widgets/pose_editor.h>
#include <QApplication>

#include <protocols/init/init.hh>
#include <libxml/parser.h> // dummy include to make sure external excludes path is set
#include <ui/ui_lib_test.h>


#include <QDebug>

int main(int argc, char *argv[])
{
	qRegisterMetaType<JSON_SP> ("JSON_SP");
	qRegisterMetaType<JSON_CSP>("JSON_CSP");

	qRegisterMetaType<std::string>("std::string");

    qDebug() << "test::main: Calling Rosetta init()...";

    protocols::init::init(argc, argv);
    ui::ui_lib_test::test();

    QApplication a(argc, argv);
	ui::widgets::PoseEditor w;
    w.show();

    return a.exec();
}
