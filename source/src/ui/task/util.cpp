#include <ui/task/util.h>
#include <ui/task/node.h>

#include <QSettings>
#include <QTimer>

#include <ctime>

namespace ui {
namespace task {

TimeStamp get_utc_timestamp()
{
	return /*std::time_t t =*/ std::time(nullptr);
}


/// calculate appropriate delay in sec until next network retry should be made
int get_retry_interval()
{
	const double desired_interval = 5;
	const double series_length = 5;

	static /*thread_local*/ double current_interval = desired_interval;
	static /*thread_local*/ double average_interval = desired_interval;
	static /*thread_local*/ TimeStamp last_retry = get_utc_timestamp();

	TimeStamp now = get_utc_timestamp();

	double interval = now - last_retry;

	average_interval = (average_interval * series_length + interval) / (series_length+1);
	if( average_interval < desired_interval ) current_interval = (current_interval+1.0) * 1.1;
	else current_interval = current_interval / 1.1;
	qDebug() << "get_retry_interval() average_interval:" << average_interval << "  current_interval:" << current_interval;

	// if( interval < desired_interval ) current_interval = (current_interval+1.0) * 2.0;
	// else current_interval = current_interval / 2.0;
	// qDebug() << "get_retry_interval() interval:" << interval << "  current_interval:" << current_interval;

	current_interval = std::min(current_interval, 60.0);

	last_retry = now;
	return current_interval;
}

NetworkCall::NetworkCall(QObject * parent) : QObject(parent)
{}

void NetworkCall::abort_network_operation()
{
	if(reply_) {
		//qDebug() << "NetworkCall::abort_network_operation ABORTING previous request!";
		disconnect(reply_, SIGNAL(finished()), 0, 0); // we have to disconnet 'finished' signal because call to abort() will issue 'finished'
		Q_EMIT reply_->abort();
		reply_->deleteLater();
		reply_ = nullptr;
	}
}

void NetworkCall::call(QString const & address, QNetworkAccessManager::Operation operation, QJsonDocument const &post_json_data)
{
	call(address, operation, post_json_data.toJson(QJsonDocument::Compact), "application/json; charset=utf-8");
}


void NetworkCall::call(QString const & address, QNetworkAccessManager::Operation operation, QByteArray const &post_data, QString const &content_type)
{
	abort_network_operation();
	result_.clear();

	address_ = address;  operation_ = operation;  post_data_ = post_data;  content_type_ = content_type;

    QUrl url( server_url() + '/' + address);
	QNetworkRequest request(url);
	request.setHeader( QNetworkRequest::ContentTypeHeader, content_type);

	if      ( operation == QNetworkAccessManager::Operation::GetOperation )    reply_ = network_access_manager()->get           (request);
	else if ( operation == QNetworkAccessManager::Operation::PostOperation )   reply_ = network_access_manager()->post          (request, post_data);
	else if ( operation == QNetworkAccessManager::Operation::DeleteOperation ) reply_ = network_access_manager()->deleteResource(request);
	else Q_ASSERT_X(false, "NetworkCall::call", "Unknown operation!");

	//connect(reply_, &QNetworkReply::finished, this, &NetworkCall::network_operation_is_finished); // does not work on old GCC compilers
	connect(reply_, SIGNAL(QNetworkReply::finished()), this, SLOT(NetworkCall::network_operation_is_finished()));
    if( reply_->isFinished() ) { network_operation_is_finished(); }
}

void NetworkCall::network_operation_is_finished()
{
	if(reply_) {
 		reply_->deleteLater();

		QVariant qv_status_code = reply_->attribute( QNetworkRequest::HttpStatusCodeAttribute );
		if( qv_status_code.isValid() ) {
			// if( int status = qv_status_code.toInt()) {
			// 	//qDebug() << "Node::upload_finished for " << node_id_ << " status:" << status;
			// }
		}

        if( reply_->error() == QNetworkReply::NoError ) {

			result_ = reply_->readAll();
			reply_ = nullptr;
			//qDebug() << "NetworkCall::network_operation_is_finished data: " << data_;

			Q_EMIT finished();

		} else {
            qDebug() << "NetworkCall::call with error: " << reply_->error() << " reply:" << reply_->readAll() << " aborting...";

            //qDebug() << "NetworkCall::call with error: " << reply_->error() << " reply:" << reply_->readAll() << " initiating retry...";
			reply_ = nullptr;

			QTimer::singleShot(get_retry_interval()*1000.0, [this]() { call(address_, operation_, post_data_); });
		}
 	}
}

QByteArray NetworkCall::result() const
{
	return result_;
}

QJsonDocument NetworkCall::json() const
{
	//QJsonDocument jd = QJsonDocument::fromJson(data);
	//QJsonObject root = jd.object();
	//return root;

	return QJsonDocument::fromJson(result_);
}

} // namespace task
} // namespace ui
