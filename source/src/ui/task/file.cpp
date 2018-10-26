
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
	QStringList list = f->name().split('.');

	if( not list.isEmpty() ) {
		QString ext = list.last();

		if( false  and  ext == "pdb" ) {
			ui_core::pose_draw::SimplePoseDrawOpenGLWidget *pv = new ui_core::pose_draw::SimplePoseDrawOpenGLWidget(parent);

			core::pose::PoseOP pose = std::make_shared<core::pose::Pose>();

			qDebug() << "pose_from_pdbstring...";
			core::import_pose::pose_from_pdbstring(*pose, std::string( f->data().constData(), f->data().length() ), f->name().toUtf8().constData() );
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

// File::File()
// {
// 	data( QByteArray() );
// }

// File::File(QString const &file_name)
// {
// 	init_from_file(file_name);
// }

// File::File(QByteArray const & file_data)
// {
// 	data(file_data);
// }


File::File(Kind kind, QString const & name, QByteArray const & file_data) : FileID(kind, name)
{
	data(file_data);
}

FileSP File::init_from_local_file(QString const & file_name)
{
	auto fsp = std::make_shared<File>(Kind::input);

	fsp->local_file_name_  = file_name;
	fsp->name_ = QFileInfo( fsp->local_file_name_ ).fileName();

	QFile file( fsp->local_file_name_ );
    if( file.open(QIODevice::ReadOnly) ) {
		fsp->data( file.readAll() );
	}

	return fsp;
}


File& File::operator=(File&& other) noexcept
{
    if(this != &other) { // no-op on self-move-assignment
		kind_ = other.kind_;
		hash_.swap(other.hash_);
		name_.swap(other.name_);
		file_data_.swap(other.file_data_);
		local_file_name_.swap(other.local_file_name_);
    }
    return *this;
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
	return (kind_ == rhs.kind_) and (name_ == rhs.name_) and (file_data_ == rhs.file_data_);
}

QDataStream &operator<<(QDataStream &out, File::Kind k)
{
	qint8 byte = static_cast<qint8>(k);
	out << byte;
	return out;
}

QDataStream &operator>>(QDataStream &in, File::Kind &k)
{
	qint8 byte;
	in >> byte;
	k = static_cast<File::Kind>(byte);
	return in;
}


QDataStream &operator<<(QDataStream &out, File const&f)
{
	out << f.kind_;
	out << f.name_;
	out << f.hash_;
	out << f.file_data_;
	return out;
}

QDataStream &operator>>(QDataStream &in, File &f)
{
	in >> f.kind_;
	in >> f.name_;
	in >> f.hash_;
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
	state_ = State::draft;

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
				case 2: return QString("kind");
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
	return 3;
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
		else if ( index.column() == 2 ) return rows_[index.row()].kind;
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

	editable_ = task.state() == Task::State::none;

	QDir project_path;
	if( auto project = task.project() ) project_path = QDir( QFileInfo( project->file_name() ).dir() );

	auto const & files = task.files();

	rows_.clear();
	rows_.reserve( files.size() );
	for(auto const & it : files) rows_.push_back( Row{it->name(), project_path.relativeFilePath( it->local_file_name() ), it->kind() == File::Kind::input ? "input" : "output" } );

	endResetModel();
}


} // namespace task
} // namespace ui
