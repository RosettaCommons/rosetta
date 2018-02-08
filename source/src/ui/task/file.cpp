
#include <ui/task/file.h>

#include <ui/task/task.h>

#include <ui/ui_core/pose_draw/SimplePoseDrawOpenGLWidget.h>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

#include <QPlainTextEdit>
#include <QDir>

namespace ui {
namespace task {

QWidget * create_viewer_for_file(QWidget *parent, FileSP const &f)
{
	QStringList list = f->file_name().split('.');

	if( not list.isEmpty() ) {
		QString ext = list.last();

		if( false  and  ext == "pdb" ) {
			ui_core::pose_draw::SimplePoseDrawOpenGLWidget *pv = new ui_core::pose_draw::SimplePoseDrawOpenGLWidget(parent);

			core::pose::PoseOP pose = std::make_shared<core::pose::Pose>();

			qDebug() << "pose_from_pdbstring...";
			core::import_pose::pose_from_pdbstring(*pose, std::string( f->data().constData(), f->data().length() ), f->file_name().toUtf8().constData() );
			qDebug() << "pose_from_pdbstring... OK";

			pv->set_pose(pose);

			// QProcess * p = new QProcess(this);
			// p->start("MacPyMOL");

			return pv;
		}

		if( ext == "log"  or  ext == "sf"  or true ) { // default viewer
			QPlainTextEdit *e = new QPlainTextEdit(f->data(), parent);
			e->setReadOnly(true);

			e->setLineWrapMode(QPlainTextEdit::NoWrap);

			//QFont font("Courier");
			e->setFont( QFontDatabase::systemFont(QFontDatabase::FixedFont) );

			//QFont font = e->defaultFont();
			//font.setFamily("Courier New");
			//e->setFont(font);

			return e;
		}

	}

	return nullptr;
}

Node::Key const _input_key_ {"input"};
Node::Key const _output_key_{"output"};

Node::Key const _script_key_ {"script"};

File::File(QString const &file_name)
{
	init_from_file(file_name);
}

File::File(QByteArray const & file_data)
{
	data(file_data);
}


File::File(QString const & file_name, QByteArray const & file_data) : file_name_(file_name)
{
	data(file_data);
}


File& File::operator=(File&& other) noexcept
{
    if(this != &other) { // no-op on self-move-assignment
		file_name_.swap(other.file_name_);
		file_data_.swap(other.file_data_);
    }
    return *this;
}


void File::init_from_file(QString const & file_name)
{
	file_name_ = file_name;
	QFile file(file_name_);
    if( file.open(QIODevice::ReadOnly) ) {
		data( file.readAll() );
	}
}

QByteArray File::data() const
{
	return qUncompress(file_data_);
}

void File::data(QByteArray const &file_data)
{
	file_data_ = qCompress(file_data);
}


bool File::operator ==(File const &rhs) const
{
	return (file_name_ == rhs.file_name_) and (file_data_ == rhs.file_data_);
}

QDataStream &operator<<(QDataStream &out, File const&f)
{
	out << f.file_name_;
	out << f.file_data_;
	return out;
}

QDataStream &operator>>(QDataStream &in, File &f)
{
	in >> f.file_name_;
	in >> f.file_data_;
	return in;
}






/*
NodeTaskSyncer::NodeTaskSyncer(Task *task) : task_(task)
{
}

void NodeTaskSyncer::submit(QString const & queue)
{
	queue_ = queue;
	state_ = State::_draft_;

	create_sync_tree();
	qDebug() << "Task[" << task_id() << "]: Task::submit() Name:" << (project_ ? *project_->find(this) : "''") << " queue:" << queue_;

	connect(root_.get(), SIGNAL(tree_synced()), this, SLOT(draft_to_queued()));

	root_->data_is_fresh(true);

	Q_EMIT task_->changed();
}



void NodeTaskSyncer::create_sync_tree()
{
	//QUuid task_id = QUuid::createUuid();  /// each distinct submitted Task _must_ have unique id

	root_ = std::make_shared<Node>("task", Node::Flags::data_in | Node::Flags::data_out | Node::Flags::topology_out);

	auto files_node = std::make_shared<Node>(Node::Flags::all);
	root_->add("files", files_node);
	//assign_input_node();




	auto script_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	root_->add("script", script_node);
	//assign_script_node();  //script_node_ = script_node;

	auto flags_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	root_->add("flags", flags_node);
	//assign_flags_node();

	auto output_node = std::make_shared<Node>(Node::Flags::topology_in);
	root_->add("output", output_node);
	//assign_output_node();

	connect_task_and_nodes();
}
*/



QVariant FileTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role == Qt::DisplayRole) {
		if (orientation == Qt::Horizontal) {
			switch (section) {
				case 0: return QString("name");
				case 1: return QString("local path");
			}
		} else if (orientation == Qt::Vertical) {
			return QString::number(section+1);
		}
	}
    return QVariant();
}


int FileTableModel::rowCount(const QModelIndex &/*parent*/) const
{
	return rows_.size();
}

int	FileTableModel::columnCount(const QModelIndex &/*parent*/) const
{
	return 2;
}

Qt::ItemFlags FileTableModel::flags(const QModelIndex &index) const
{
    if (!index.isValid()) return Qt::NoItemFlags;

	if( index.column() == 0  and  editable_ ) return Qt::ItemIsEditable | QAbstractItemModel::flags(index);

	return QAbstractItemModel::flags(index);
}

QVariant FileTableModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid()) return QVariant();
    if (role != Qt::DisplayRole && role != Qt::EditRole) return QVariant();

	//if( index.column() == 0 ) {
	//return QStringLiteral("r%1 c%2").arg(index.row()).arg(index.column());

	if( index.row() < static_cast<int>( rows_.size() ) ) {

		if( index.column() == 0 ) return rows_[index.row()].name;
		else if ( index.column() == 1 ) return rows_[index.row()].path;
	}
	return QVariant();
}

bool FileTableModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if (role != Qt::EditRole) return false;

	QString s = value.toString();
	if( not s.isEmpty()  and  static_cast<std::size_t>(index.row()) < rows_.size() ) {
		QString previous_value = rows_[index.row()].name;
		//QTimer::singleShot(0, this, [=]() { Q_EMIT this->rename_file(previous_value, s); } );  //SLOT( rename(previous_value, s) ) );
		Q_EMIT this->rename_file(previous_value, s);
	}

    return false;
}


void FileTableModel::update_from_task(Task const &task)
{
	beginResetModel();

	editable_ = task.state() == Task::State::_none_;

	QDir project_path;
	if( auto project = task.project() ) project_path = QDir( QFileInfo( project->file_name() ).dir() );

	auto const & files = task.files();

	rows_.clear();
	rows_.reserve( files.size() );
	for(auto const & it : files) rows_.push_back( Row{it.first, project_path.relativeFilePath( it.second->file_name() ) } );

	endResetModel();
}


} // namespace task
} // namespace ui
