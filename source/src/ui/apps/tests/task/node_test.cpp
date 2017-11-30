#include <QtTest>

#include <ui/task/node.h>
#include <ui/task/task.h>

#include <QBuffer>
#include <QThread>

using namespace ui::task;


std::pair<QByteArray, int> download_url(QString const & url, QNetworkAccessManager::Operation operation = QNetworkAccessManager::Operation::GetOperation)
{
	QNetworkReply *reply = nullptr;
	if      ( operation == QNetworkAccessManager::Operation::GetOperation )    reply = network_access_manager()->get           ( QNetworkRequest(url) );
	else if ( operation == QNetworkAccessManager::Operation::DeleteOperation ) reply = network_access_manager()->deleteResource( QNetworkRequest(url) );
	else Q_ASSERT_X(false, "download_url", "Unknown operation!");

    QEventLoop loop;
    QObject::connect(reply, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();

	std::pair<QByteArray, int> r = std::make_pair(reply->readAll(), 0);

	QVariant qv_status_code = reply->attribute( QNetworkRequest::HttpStatusCodeAttribute );
	if( qv_status_code.isValid() ) {
		r.second = qv_status_code.toInt();
	}
    delete reply;

	return r;
}


class NodeTest : public QObject
{
    Q_OBJECT

public:
    NodeTest();
    ~NodeTest();

private slots:
    void initTestCase();
    void cleanupTestCase();

	void test_node_serialization();
	void test_task_serialization();
	void test_project_serialization();

	void test_node_upload_and_download();
	void test_task_upload_and_download();

private:

};

NodeTest::NodeTest()
{

}

NodeTest::~NodeTest()
{

}


void NodeTest::initTestCase()
{
	QCoreApplication::setOrganizationName("RosettaCommons");
    QCoreApplication::setOrganizationDomain("RosettaCommons.org");
    QCoreApplication::setApplicationName("Workbench");

	QSettings::setDefaultFormat(QSettings::IniFormat); // comment this out later, when config editor functionality is in place
}

void NodeTest::cleanupTestCase()
{

}

void NodeTest::test_node_serialization()
{
	QUuid const root_id     (0, 0, 0, 0x10, 0x00, 0, 0, 0, 0, 0, 0);
	QUuid const leaf_1_id   (0, 0, 0, 0x11, 0x00, 0, 0, 0, 0, 0, 0);
	QUuid const leaf_1_1_id (0, 0, 0, 0x11, 0x10, 0, 0, 0, 0, 0, 0);
	QUuid const leaf_1_2_id (0, 0, 0, 0x11, 0x20, 0, 0, 0, 0, 0, 0);

	NodeSP root = std::make_shared<Node>(Node::Flags::all, root_id);

	auto leaf_1 = std::make_shared<Node>(Node::Flags::data_in, leaf_1_id);
	root->add("leaf_1", leaf_1);

	auto leaf_1_1 = std::make_shared<Node>(Node::Flags::data_out, leaf_1_1_id);
	leaf_1->add("leaf_1_1", leaf_1_1);

	auto leaf_1_2 = std::make_shared<Node>(Node::Flags::topology_in, leaf_1_2_id);
	leaf_1->add("leaf_1_2", leaf_1_2);


	auto n_root = std::make_shared<Node>();

	QVERIFY( *root != *n_root );


	QBuffer buffer;
    buffer.open(QBuffer::ReadWrite);

	QDataStream out(&buffer);
	out << *root;

	buffer.seek(0);
	QDataStream in(&buffer);

	in >> *n_root;

	QCOMPARE( *root, *n_root );
}


void NodeTest::test_task_serialization()
{
	TaskSP task = std::make_shared<Task>("_t_task_");
	TaskSP n_task = std::make_shared<Task>();

	task->input( File("input test file name", "some test input data") );
	task->flags( File("flags test file name", "some test flags data") );
	task->script( File("script test file name", "some test script data") );

	QVERIFY( *task != *n_task );

	QBuffer buffer;
    buffer.open(QBuffer::ReadWrite);

	QDataStream out(&buffer);
	out << *task;

	buffer.seek(0);
	QDataStream in(&buffer);

	in >> *n_task;

	QCOMPARE( *task, *n_task );
}

void NodeTest::test_project_serialization()
{
	ProjectSP project = std::make_shared<Project>();
	ProjectSP n_project = std::make_shared<Project>();

	project->add("Task N1", std::make_shared<Task>("test task") );

	QVERIFY( *project != *n_project );

	QBuffer buffer;
    buffer.open(QBuffer::ReadWrite);

	QDataStream out(&buffer);
	out << *project;

	buffer.seek(0);
	QDataStream in(&buffer);

	in >> *n_project;

	QCOMPARE( *project, *n_project );
}

void NodeTest::test_node_upload_and_download()
{
	QUuid const root_id     (0xa000000a, 0, 0, 0x10, 0x00, 0, 0, 0, 0, 0, 9);
	QUuid const leaf_1_id   (0xb000000b, 0, 0, 0x11, 0x00, 0, 0, 0, 0, 0, 8);
	QUuid const leaf_1_1_id (0xc000000c, 0, 0, 0x11, 0x10, 0, 0, 0, 0, 0, 7);

	NodeSP root = std::make_shared<Node>(Node::Flags::all, root_id);

	auto leaf_1 = std::make_shared<Node>(Node::Flags::all, leaf_1_id);
	root->add("leaf_1", leaf_1);

	auto leaf_1_1 = std::make_shared<Node>(Node::Flags::all, leaf_1_1_id);
	leaf_1->add("leaf_1_1", leaf_1_1);

	QSignalSpy spy_root(root.get(),         SIGNAL( synced() ) );
	QSignalSpy spy_leaf_1(leaf_1.get(),     SIGNAL( synced() ) );
	QSignalSpy spy_leaf_1_1(leaf_1_1.get(), SIGNAL( synced() ) );


	// Upload
	root->data_is_fresh(true);

	QCOMPARE(spy_root.wait(1000), true);
	QCOMPARE(spy_root.count(), 1);

	QCOMPARE(spy_leaf_1.wait(1000), true);
	QCOMPARE(spy_leaf_1.count(), 1);

	QCOMPARE(spy_leaf_1_1.wait(1000), true);
	QCOMPARE(spy_leaf_1_1.count(), 1);


	// Download
	NodeSP root_dwn = std::make_shared<Node>(Node::Flags::all, root_id);

	QSignalSpy spy_root_dwn(root_dwn.get(), SIGNAL( tree_synced() ) );

	root_dwn->data_is_outdated();

	QCOMPARE(spy_root_dwn.wait(1000), true);
	QCOMPARE(spy_root_dwn.count(), 1);


	// Now checking if uploaded and downloaded tree's are equal...
	QCOMPARE( *root, *root_dwn );


	// Insert new leaf into orignal version
	QUuid const leaf_1_2_id(0xd000000d, 0, 0, 0x11, 0x20, 0, 0, 0, 0, 0, 6);
	auto leaf_1_2 = std::make_shared<Node>(Node::Flags::all, leaf_1_2_id);
	leaf_1->add("leaf_1_2", leaf_1_2);

	QVERIFY( *root != *root_dwn );

	QSignalSpy spy_root2(root.get(), SIGNAL( tree_synced() ) );

	QThread::sleep(1); // to make sure timestamp is changed

	leaf_1->data_is_fresh(true);

	QCOMPARE(spy_root2.wait(2000), true);
	QCOMPARE(spy_root2.count(), 1);

	// start listening to updates and check if tree became synced
	Updater::get().subscribe( root_dwn.get() );  //root_dwn->listen_to_updates();

	QCOMPARE(spy_root_dwn.wait(16000), true);
	QCOMPARE(spy_root_dwn.count(), 2);

	QCOMPARE( *root, *root_dwn );


	// Deleting
	auto r = download_url( server_url() + "/node/" + root_id.toString().mid(1, 36), QNetworkAccessManager::Operation::DeleteOperation);
	QCOMPARE( r.second, 200 );
}


void NodeTest::test_task_upload_and_download()
{
	TaskSP task = std::make_shared<Task>("_t_task_");

	QSignalSpy spy(task.get(), SIGNAL( submitted() ) );

	// Uploading
    task->submit("test");

	QCOMPARE(spy.wait(1000), true);
	QCOMPARE(spy.count(), 1);


	// Deleting
	qDebug() << "Deleting Task " << task->task_id() << "...";
	auto r = download_url( server_url() + "/node/" + task->task_id().toString().mid(1, 36), QNetworkAccessManager::Operation::DeleteOperation);
	QCOMPARE( r.second, 200 );
}


// void NodeTest::test_case1()
// {
// 	QString str = "Hello";
//     QVERIFY(str.toUpper() == "HELLO");
// }



QTEST_MAIN(NodeTest)

#include "node_test.moc"
