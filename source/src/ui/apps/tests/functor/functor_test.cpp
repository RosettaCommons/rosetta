#include <ui/task/functor.h>
#include <ui/task/task.h>

#include <QtEndian>
#include <QThread>
#include <QtTest>

using namespace ui::task;

using CustomHeader = std::pair<QByteArray, QByteArray>;

std::pair<QByteArray, int> download_url(QString const & address, QNetworkAccessManager::Operation operation = QNetworkAccessManager::Operation::GetOperation, QByteArray const &post_data=QByteArray(),
										QString const &content_type="", std::vector<CustomHeader> const & headers=std::vector<CustomHeader>() )
{
	QUrl url(address);
	QNetworkRequest request(url);
	if( not content_type.isEmpty() ) request.setHeader( QNetworkRequest::ContentTypeHeader, content_type);

	for( auto const & h : headers ) request.setRawHeader(h.first, h.second);

	QNetworkReply *reply = nullptr;

	if      ( operation == QNetworkAccessManager::Operation::GetOperation )    reply = network_access_manager()->get           (request);
	else if ( operation == QNetworkAccessManager::Operation::PostOperation )   reply = network_access_manager()->post          (request, post_data);
	else if ( operation == QNetworkAccessManager::Operation::DeleteOperation ) reply = network_access_manager()->deleteResource(request);

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

class FunctorTest : public QObject
{
    Q_OBJECT

public:
    //FunctorTest();
    //~FunctorTest();

private slots:
    void initTestCase();
    void cleanupTestCase();

	void test_basic();

	void test_functor_sequence();

	void test_functor_async_sequence();

	void test_task_upload();

	void test_compression();

private:

};

void FunctorTest::initTestCase()
{
	QCoreApplication::setOrganizationName("RosettaCommons");
    QCoreApplication::setOrganizationDomain("RosettaCommons.org");
    QCoreApplication::setApplicationName("Workbench");

	QSettings::setDefaultFormat(QSettings::IniFormat); // comment this out later, when config editor functionality is in place
}

void FunctorTest::cleanupTestCase()
{
}

void FunctorTest::test_basic()
{
	Functor f("FunctorTest::test_basic");

	QSignalSpy spy_started  ( &f, &Functor::started );
	QSignalSpy spy_tick     ( &f, &Functor::tick );
	QSignalSpy spy_finished ( &f, &Functor::finished );
	QSignalSpy spy_final    ( &f, &Functor::final );

	f.execute();

	QCOMPARE(spy_started . count(), 1);
	QCOMPARE(spy_tick     .count(), 1);
	QCOMPARE(spy_finished .count(), 1);
	QCOMPARE(spy_final    .count(), 1);
}


void FunctorTest::test_functor_sequence()
{
	FunctorSequence s("FunctorTest::test_functor_sequence");

	for(int i=0; i<3; i++) {
		s.add_functor( std::make_shared<Functor>( QString("f-%1").arg(i) ) );
	}

	QSignalSpy spy_started  ( &s, &Functor::started );
	QSignalSpy spy_tick     ( &s, &Functor::tick );
	QSignalSpy spy_finished ( &s, &Functor::finished );
	QSignalSpy spy_final    ( &s, &Functor::final );

	s.execute();

	QCOMPARE(spy_started . count(), 1);
	QCOMPARE(spy_tick     .count(), 5);
	QCOMPARE(spy_finished .count(), 1);
	QCOMPARE(spy_final    .count(), 1);
}


void FunctorTest::test_functor_async_sequence()
{
	FunctorASyncSequence s("FunctorTest::test_functor_async_sequence");

	for(int i=0; i<5; i++) {
		s.add_functor( std::make_shared<Functor>( QString("f-%1").arg(i) ) );
	}

	QSignalSpy spy_started  ( &s, &Functor::started );
	QSignalSpy spy_tick     ( &s, &Functor::tick );
	QSignalSpy spy_finished ( &s, &Functor::finished );
	QSignalSpy spy_final    ( &s, &Functor::final );

	s.execute();

	QCOMPARE(spy_started . count(), 1);
	QCOMPARE(spy_tick     .count(), 6);
	QCOMPARE(spy_finished .count(), 1);
	QCOMPARE(spy_final    .count(), 1);
}


void FunctorTest::test_task_upload()
{
	auto task = std::make_shared<Task>();

	for(int i=0; i<4; i++) {
		task->add_file( QString("file-%1").arg(i), std::make_shared<File>( QByteArray("abc") + QByteArray( std::to_string(i).c_str() ) ) );
	}

	for(auto const & it : task->files() ) {
		QVERIFY( it.second->hash() == "" );
	}

	QSignalSpy spy_task(task.get(), SIGNAL( submitted() ) );
	QSignalSpy spy_task_final(task.get(), SIGNAL( final() ) );

	QVERIFY( task->task_id().isEmpty() );


	task->submit("test");

	QCOMPARE(spy_task.wait(5000), true);
	QCOMPARE(spy_task.count(), 1);

	QVERIFY( not task->task_id().isEmpty() );

	for(auto const & it : task->files() ) {
		QVERIFY( it.second->hash() != "" );
	}

	// post new files
	auto r = download_url( task_api_url() + "/file/" + task->task_id() + "/extra_file_1", QNetworkAccessManager::Operation::PostOperation, "qwerty", "application/octet-stream");
	QCOMPARE( r.second, 200 );

	r = download_url( task_api_url() + "/file/" + task->task_id() + "/extra_file_2", QNetworkAccessManager::Operation::PostOperation, "asdfgh", "application/octet-stream");
	QCOMPARE( r.second, 200 );

	char const * task_json1 = "{\"state\":\"running\",\"app\":\"docking_protocol\",\"name\":\"\",\"queue\":\"test\",\"version\":\"\"}";
	r = download_url( task_api_url() + "/task/" + task->task_id(), QNetworkAccessManager::Operation::PostOperation, task_json1, "application/json");
	QCOMPARE( r.second, 200 );

	char const * task_json2 = "{\"state\":\"finished\",\"app\":\"docking_protocol\",\"name\":\"\",\"queue\":\"test\",\"version\":\"\"}";
	r = download_url( task_api_url() + "/task/" + task->task_id(), QNetworkAccessManager::Operation::PostOperation, task_json2, "application/json");
	QCOMPARE( r.second, 200 );

	// technically subscribe willbe called as result of submit, however it will happened _before_ we posted additional files, so we initate one more subscribe to avoid waiting
	task->subscribe();

	QCOMPARE(spy_task_final.wait(9000), true);
	QCOMPARE(spy_task_final.count(), 1);

	//QThread::sleep(10); // to make sure timestamp is changed

	// Deleting
	qDebug() << "Deleting Task " << task->task_id() << "...";
    r = download_url( task_api_url() + "/task/" + task->task_id(), QNetworkAccessManager::Operation::DeleteOperation);
	QCOMPARE( r.second, 200 );
}

void FunctorTest::test_compression()
{
    // auto r = download_url( task_api_url() + "/_compression_test", QNetworkAccessManager::Operation::GetOperation);
	// qDebug() << "_compression_test: " << r;

	// auto r = download_url( task_api_url() + "/_compression_test_json", QNetworkAccessManager::Operation::GetOperation);
	// qDebug() << "_compression_test_json: " << r;

	// QByteArray p("Qt compression Test"); p+=p; p+=p; p+=p; p+=p; p+=p; p+=p;
	// auto payload = qCompress(p).mid(4); // removing extra 32bit int (size of uncrompressed data from beginning of the data to make result compatible with raw zlib encoding

	//quint32 size = p.size();
	//qToBigEndian(p.size(), &size);
	//qint32_be size(p.size());
	//payload = payload.insert(0, QByteArray(reinterpret_cast<char const*>(&size), sizeof(size) ) );

	// qDebug() << "_compression_test:" << p.size() << " --> " << payload.size();
	// download_url("http://ui.graylab.jhu.edu/api/_compression_test", QNetworkAccessManager::Operation::PostOperation, payload, "application/octet-stream", { {"Content-Encoding", "deflate"} });
	// qDebug() << "_compression_test: payload:" << payload;

	// auto json_string = "{\"state\":\"running\",\"app\":\"docking_protocol\",\"name\":\"\",\"queue\":\"test\",\"version\":\"\"}";
	// QByteArray payload = qCompress( QByteArray(json_string) ).mid(4); // removing extra 32bit int (size of uncrompressed data from beginning of the data to make result compatible with raw zlib encoding
	// download_url("http://ui.graylab.jhu.edu/api/_compression_test_json", QNetworkAccessManager::Operation::PostOperation, payload, "application/json", { {"Content-Encoding", "deflate"} });
}



QTEST_MAIN(FunctorTest)

#include "functor_test.moc"
