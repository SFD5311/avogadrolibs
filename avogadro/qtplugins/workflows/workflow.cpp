/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2016 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#include "workflow.h"

#include <avogadro/qtgui/avogadropython.h>
#include <avogadro/qtgui/filebrowsewidget.h>
#include <avogadro/qtgui/fileformatdialog.h>
#include <avogadro/qtgui/molecule.h>
#include <avogadro/qtgui/utilities.h>
#include <avogadro/qtgui/interfacescript.h>
#include <avogadro/qtgui/interfacewidget.h>

#include <QtWidgets/QAction>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QVBoxLayout>

#include <QtCore/QCoreApplication>
#include <QtCore/QDebug>
#include <QtCore/QDir>
#include <QtCore/QtPlugin>
#include <QtCore/QSettings>
#include <QtCore/QStringList>
#include <QtCore/QStandardPaths>

namespace Avogadro {
namespace QtPlugins {

using Avogadro::QtGui::InterfaceScript;
using Avogadro::QtGui::InterfaceWidget;

Workflow::Workflow(QObject *parent_) :
  ExtensionPlugin(parent_),
  m_molecule(NULL),
  m_outputFormat(NULL)
{
  refreshScripts();
}

Workflow::~Workflow()
{
  qDeleteAll(m_dialogs.values());
  m_dialogs.clear();
}

QList<QAction *> Workflow::actions() const
{
  return m_actions;
}

QStringList Workflow::menuPath(QAction *) const
{
  QStringList path;
  // TODO: Need to query the scripts for these
  path << tr("&Extensions") << tr("Scripts");
  return path;
}

void Workflow::setMolecule(QtGui::Molecule *mol)
{
  if (m_molecule == mol)
    return;

  m_molecule = mol;

  foreach (InterfaceWidget *dlg, m_dialogs.values())
    dlg->setMolecule(mol);
}

bool Workflow::readMolecule(QtGui::Molecule &mol)
{
  Io::FileFormat *reader = m_outputFormat->newInstance();
  bool success = reader->readFile(m_outputFileName.toStdString(), mol);
  if (!success) {
    QMessageBox::information(qobject_cast<QWidget*>(parent()),
                             tr("Error"),
                             tr("Error reading output file '%1':\n%2")
                             .arg(m_outputFileName)
                             .arg(QString::fromStdString(reader->error())));
  }

  m_outputFormat = NULL;
  m_outputFileName.clear();

  return success;
}

void Workflow::refreshScripts()
{
  updateScripts();
  updateActions();
}

void Workflow::menuActivated()
{
  QAction *theSender = qobject_cast<QAction*>(sender());
  if (!theSender)
    return;

  QString scriptFileName = theSender->data().toString();
  QWidget *theParent = qobject_cast<QWidget*>(parent());
  InterfaceWidget *dlg = m_dialogs.value(scriptFileName, NULL);

  if (!dlg) {
    dlg = new InterfaceWidget(theParent);
    m_dialogs.insert(scriptFileName, dlg);
  }
  dlg->setMolecule(m_molecule);
  dlg->show();
  dlg->raise();
}

void Workflow::configurePython()
{
  // Create objects
  QSettings settings;
  QDialog dlg(qobject_cast<QWidget*>(parent()));
  QLabel *label = new QLabel;
  QVBoxLayout *layout = new QVBoxLayout;
  QtGui::FileBrowseWidget *browser = new QtGui::FileBrowseWidget;
  QDialogButtonBox *buttonBox = new QDialogButtonBox;

  // Configure objects
  // Check for python interpreter in env var
  QString pythonInterp = QString::fromLocal8Bit(
        qgetenv("AVO_PYTHON_INTERPRETER"));
  if (pythonInterp.isEmpty()) {
    // Check settings
    pythonInterp = settings.value("interpreters/python",
                                  QString()).toString();
  }
  // Use compile-time default if still not found.
  if (pythonInterp.isEmpty())
    pythonInterp = QString(pythonInterpreterPath);
  browser->setMode(QtGui::FileBrowseWidget::ExecutableFile);
  browser->setFileName(pythonInterp);

  buttonBox->setStandardButtons(QDialogButtonBox::Ok
                                | QDialogButtonBox::Cancel);

  dlg.setWindowTitle(tr("Set path to Python interpreter:"));
  label->setText(tr("Select the python interpreter to run external scripts.\n"
                    "Avogadro must be restarted for any changes to take effect."));

  // Build layout
  layout->addWidget(label);
  layout->addWidget(browser);
  layout->addWidget(buttonBox);
  dlg.setLayout(layout);

  // Connect
  connect(buttonBox, SIGNAL(accepted()), &dlg, SLOT(accept()));
  connect(buttonBox, SIGNAL(rejected()), &dlg, SLOT(reject()));

  // Show dialog
  QDialog::DialogCode response = static_cast<QDialog::DialogCode>(dlg.exec());
  if (response != QDialog::Accepted)
    return;

  // Handle response
  settings.setValue("interpreters/python", browser->fileName());
}

void Workflow::updateScripts()
{
  m_workflowScripts.clear();

  // List of directories to check.
  /// @todo Custom script locations
  QStringList dirs;

  // add the default paths
  QStringList stdPaths
    = QStandardPaths::standardLocations(QStandardPaths::AppLocalDataLocation);
  foreach (const QString &dirStr, stdPaths) {
    QString path = dirStr + "/scripts/workflows";
    QDir dir( path );
    qDebug() << "Checking for generator scripts in" << path;
    if (dir.exists() && dir.isReadable())
       dirs << path;
  }

  dirs << QCoreApplication::applicationDirPath() + "/../"
          + QtGui::Utilities::libraryDirectory()
          + "/avogadro2/scripts/workflows";

  foreach (const QString &dirStr, dirs) {
    qDebug() << "Checking for generator scripts in" << dirStr;
    QDir dir(dirStr);
    if (dir.exists() && dir.isReadable()) {
      foreach (const QFileInfo &file, dir.entryInfoList(QDir::Files |
                                                        QDir::NoDotAndDotDot)) {
        QString filePath = file.absoluteFilePath();
        QString displayName;
        if (queryProgramName(filePath, displayName))
          m_workflowScripts.insert(displayName, filePath);
      }
    }
  }
}

void Workflow::updateActions()
{
  m_actions.clear();

  QAction *action = new QAction(tr("Set Python Path..."), this);
  connect(action, SIGNAL(triggered()), SLOT(configurePython()));
  m_actions << action;

  foreach (const QString &programName, m_workflowScripts.uniqueKeys()) {
    QStringList scripts = m_workflowScripts.values(programName);
    // Include the full path if there are multiple generators with the same name.
    if (scripts.size() == 1) {
      addAction(programName, scripts.first());
    }
    else {
      foreach (const QString &filePath, scripts) {
        addAction(QString("%1 (%2)").arg(programName, filePath), filePath);
      }
    }
  }
}

void Workflow::addAction(const QString &label,
                             const QString &scriptFilePath)
{
  QAction *action = new QAction(label, this);
  action->setData(scriptFilePath);
  action->setEnabled(true);
  connect(action, SIGNAL(triggered()), SLOT(menuActivated()));
  m_actions << action;
}

bool Workflow::queryProgramName(const QString &scriptFilePath,
                                    QString &displayName)
{
  InterfaceScript gen(scriptFilePath);
  displayName = gen.displayName();
  if (gen.hasErrors()) {
    displayName.clear();
    qWarning() << "Workflow::queryProgramName: Unable to retrieve program "
                  "name for" << scriptFilePath << ";"
               << gen.errorList().join("\n\n");
    return false;
  }
  return true;
}

}
}
