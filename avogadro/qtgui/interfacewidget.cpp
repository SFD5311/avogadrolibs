/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2013 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#include "interfacewidget.h"

#include <avogadro/qtgui/filebrowsewidget.h>
#include <avogadro/qtgui/molecule.h>

#include <QtWidgets/QComboBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QTextBrowser>

#include <QtCore/QJsonDocument>
#include <QtCore/QJsonArray>
#include <QtCore/QDebug>
#include <QtCore/QPointer>
#include <QtCore/QSettings>
#include <QtCore/QTimer>

namespace Avogadro {

InterfaceWidget::InterfaceWidget(QWidget *parent_) :
  QWidget(parent_),
  m_molecule(NULL),
  m_interfaceScript(QString())
{
}

InterfaceWidget::~InterfaceWidget()
{
}

void InterfaceWidget::setInterfaceScript(const QString &scriptFile)
{
  m_interfaceScript.setScriptFilePath(scriptFile);
  updateOptions();
}

void InterfaceWidget::setMolecule(QtGui::Molecule *mol)
{
  if (mol == m_molecule)
    return;

  if (m_molecule)
    m_molecule->disconnect(this);

  m_molecule = mol;
}

void InterfaceWidget::defaultsClicked()
{
  setOptionDefaults();
}

void InterfaceWidget::setWarningText(const QString &warn)
{
  qWarning() << tr("Script returns warnings:\n") << warn;
}

QString InterfaceWidget::warningText() const
{
  return QString();
}

void InterfaceWidget::showError(const QString &err)
{
  qWarning() << err;

  QWidget *theParent = this->isVisible() ? this
                                         : qobject_cast<QWidget*>(parent());
  QDialog dlg(theParent);
  QVBoxLayout *vbox = new QVBoxLayout();
  QLabel *label = new QLabel(tr("An error has occurred:"));
  vbox->addWidget(label);
  QTextBrowser *textBrowser = new QTextBrowser();

  // adjust the size of the text browser to ~80 char wide, ~20 lines high
  QSize theSize = textBrowser->sizeHint();
  QFontMetrics metrics(textBrowser->currentFont());
  int charWidth = metrics.width("i7OPlmWn9/") / 10;
  int charHeight = metrics.lineSpacing();
  theSize.setWidth(80 * charWidth);
  theSize.setHeight(20 * charHeight);
  textBrowser->setMinimumSize(theSize);
  textBrowser->setText(err);
  vbox->addWidget(textBrowser);
  dlg.setLayout(vbox);

  dlg.exec();
}

QString InterfaceWidget::settingsKey(const QString &identifier) const
{
  return QString("scriptPlugin/%1/%2").arg(m_interfaceScript.displayName(),
                                           identifier);
}

QString InterfaceWidget::lookupOptionType(const QString &name) const
{
  if (!m_options.contains("userOptions") ||
      !m_options["userOptions"].isObject()) {
    qWarning() << tr("'userOptions' missing, or not an object.");
    return QString();
  }

  QJsonObject userOptions = m_options["userOptions"].toObject();

  if (!userOptions.contains(name)) {
    qWarning() << tr("Option '%1' not found in userOptions.").arg(name);
    return QString();
  }

  if (!userOptions.value(name).isObject()) {
    qWarning() << tr("Option '%1' does not refer to an object.");
    return QString();
  }

  QJsonObject obj = userOptions[name].toObject();

  if (!obj.contains("type") ||
      !obj.value("type").isString()) {
    qWarning() << tr("'type' is not a string for option '%1'.").arg(name);
    return QString();
  }

  return obj["type"].toString();
}

void InterfaceWidget::updateOptions()
{
  // Create the widgets, etc for the gui
  buildOptionGui();
  setOptionDefaults();
}

void InterfaceWidget::buildOptionGui()
{
  // Clear old widgets from the layout
  m_widgets.clear();
  delete layout(); // kill my layout
  QFormLayout *form = new QFormLayout;
  setLayout(form);

  if (!m_options.contains("userOptions") ||
      !m_options["userOptions"].isObject()) {
    showError(tr("'userOptions' missing, or not an object:\n%1")
              .arg(QString(QJsonDocument(m_options).toJson())));
    return;
  }

  QJsonObject userOptions = m_options.value("userOptions").toObject();

  // Title first
  if (userOptions.contains("Title"))
    addOptionRow(tr("Title"), userOptions.take("Title"));

  // File basename next:
  if (userOptions.contains("Filename Base"))
    addOptionRow(tr("Filename Base"), userOptions.take("Filename Base"));

  // Number of cores next:
  if (userOptions.contains("Processor Cores"))
    addOptionRow(tr("Processor Cores"), userOptions.take("Processor Cores"));

  // Calculation Type next:
  if (userOptions.contains("Calculation Type"))
    addOptionRow(tr("Calculation Type"), userOptions.take("Calculation Type"));

  // Theory/basis next. Combine into one row if both present.
  bool hasTheory = userOptions.contains("Theory");
  bool hasBasis = userOptions.contains("Basis");
  if (hasTheory && hasBasis) {
    QWidget *theoryWidget = createOptionWidget(userOptions.take("Theory"));
    QWidget *basisWidget = createOptionWidget(userOptions.take("Basis"));
    QHBoxLayout *hbox = new QHBoxLayout;
    if (theoryWidget) {
      theoryWidget->setObjectName("Theory");
      hbox->addWidget(theoryWidget);
      m_widgets.insert("Theory", theoryWidget);
    }
    if (basisWidget) {
      basisWidget->setObjectName("Basis");
      hbox->addWidget(basisWidget);
      m_widgets.insert("Basis", basisWidget);
    }
    hbox->addStretch();

    form->addRow(tr("Theory:"), hbox);
  }
  else {
    if (hasTheory)
      addOptionRow(tr("Theory"), userOptions.take("Theory"));
    if (hasBasis)
      addOptionRow(tr("Basis"), userOptions.take("Basis"));
  }

  // Other special cases:
  if (userOptions.contains("Charge"))
    addOptionRow(tr("Charge"), userOptions.take("Charge"));
  if (userOptions.contains("Multiplicity"))
    addOptionRow(tr("Multiplicity"), userOptions.take("Multiplicity"));

  // Add remaining keys at bottom.
  for (QJsonObject::const_iterator it = userOptions.constBegin(),
       itEnd = userOptions.constEnd(); it != itEnd; ++it) {
    addOptionRow(it.key(), it.value());
  }

  // Make connections for standard options:
  if (QComboBox *combo = qobject_cast<QComboBox*>(
        m_widgets.value("Calculation Type", NULL))) {
    connect(combo, SIGNAL(currentIndexChanged(int)),
            SLOT(updateTitlePlaceholder()));
  }
  if (QComboBox *combo = qobject_cast<QComboBox*>(
        m_widgets.value("Theory", NULL))) {
    connect(combo, SIGNAL(currentIndexChanged(int)),
            SLOT(updateTitlePlaceholder()));
  }
  if (QComboBox *combo = qobject_cast<QComboBox*>(
        m_widgets.value("Basis", NULL))) {
    connect(combo, SIGNAL(currentIndexChanged(int)),
            SLOT(updateTitlePlaceholder()));
  }
}

void InterfaceWidget::addOptionRow(const QString &label,
                                        const QJsonValue &option)
{
  QWidget *widget = createOptionWidget(option);
  if (!widget)
    return;

  QFormLayout *form = qobject_cast<QFormLayout*>(this->layout());
  if (!form) {
    qWarning() << "Cannot add option" << label
               << "to GUI -- layout is not a form.";
    widget->deleteLater();
    return;
  }

  // For lookups during unit testing:
  widget->setObjectName(label);

  form->addRow(label + ":", widget);
  m_widgets.insert(label, widget);
}

QWidget *InterfaceWidget::createOptionWidget(const QJsonValue &option)
{
  if (!option.isObject())
    return NULL;

  QJsonObject obj = option.toObject();

  if (!obj.contains("type") ||
      !obj.value("type").isString())
    return NULL;

  QString type = obj["type"].toString();

  if (type == "stringList")
    return createStringListWidget(obj);
  else if (type == "string")
    return createStringWidget(obj);
  else if (type == "filePath")
    return createFilePathWidget(obj);
  else if (type == "integer")
    return createIntegerWidget(obj);
  else if (type == "float")
    return createFloatWidget(obj);
  else if (type == "boolean")
    return createBooleanWidget(obj);

  qDebug() << "Unrecognized option type:" << type;
  return NULL;
}

QWidget *InterfaceWidget::createStringListWidget(const QJsonObject &obj)
{
  if (!obj.contains("values") || !obj["values"].isArray()) {
    qDebug() << "QuantumInputDialog::createStringListWidget()"
                "values missing, or not array!";
    return NULL;
  }

  QJsonArray valueArray = obj["values"].toArray();

  QComboBox *combo = new QComboBox(this);

  for (QJsonArray::const_iterator vit = valueArray.constBegin(),
       vitEnd = valueArray.constEnd(); vit != vitEnd; ++vit) {
    if ((*vit).isString())
      combo->addItem((*vit).toString());
    else
      qDebug() << "Cannot convert value to string for stringList:" << *vit;
  }
  connect(combo, SIGNAL(currentIndexChanged(int)), SLOT(updatePreviewText()));

  return combo;
}

QWidget *InterfaceWidget::createStringWidget(const QJsonObject &obj)
{
  Q_UNUSED(obj);
  QLineEdit *edit = new QLineEdit(this);
  connect(edit, SIGNAL(textChanged(QString)), SLOT(updatePreviewText()));
  return edit;
}

QWidget *InterfaceWidget::createFilePathWidget(const QJsonObject &obj)
{
  Q_UNUSED(obj);
  QtGui::FileBrowseWidget *fileBrowse = new QtGui::FileBrowseWidget(this);
  connect(fileBrowse, SIGNAL(fileNameChanged(QString)),
          SLOT(updatePreviewText()));
  return fileBrowse;
}

QWidget *InterfaceWidget::createIntegerWidget(const QJsonObject &obj)
{
  QSpinBox *spin = new QSpinBox(this);
  if (obj.contains("minimum") &&
      obj.value("minimum").isDouble()) {
    spin->setMinimum(static_cast<int>(obj["minimum"].toDouble() + 0.5));
  }
  if (obj.contains("maximum") &&
      obj.value("maximum").isDouble()) {
    spin->setMaximum(static_cast<int>(obj["maximum"].toDouble() + 0.5));
  }
  if (obj.contains("prefix") &&
      obj.value("prefix").isString()) {
    spin->setPrefix(obj["prefix"].toString());
  }
  if (obj.contains("suffix") &&
      obj.value("suffix").isString()) {
    spin->setSuffix(obj["suffix"].toString());
  }
  connect(spin, SIGNAL(valueChanged(int)), SLOT(updatePreviewText()));
  return spin;
}

QWidget *InterfaceWidget::createFloatWidget(const QJsonObject &obj)
{
  QDoubleSpinBox *spin = new QDoubleSpinBox(this);
  if (obj.contains("minimum") &&
      obj.value("minimum").isDouble()) {
    spin->setMinimum(obj["minimum"].toDouble());
  }
  if (obj.contains("maximum") &&
      obj.value("maximum").isDouble()) {
    spin->setMaximum(obj["maximum"].toDouble());
  }
  if (obj.contains("precision") &&
      obj.value("precision").isDouble()) {
    spin->setDecimals(static_cast<int>(obj["precision"].toDouble()));
  }
  if (obj.contains("prefix") &&
      obj.value("prefix").isString()) {
    spin->setPrefix(obj["prefix"].toString());
  }
  if (obj.contains("suffix") &&
      obj.value("suffix").isString()) {
    spin->setSuffix(obj["suffix"].toString());
  }
  connect(spin, SIGNAL(valueChanged(double)), SLOT(updatePreviewText()));
  return spin;
}

QWidget *InterfaceWidget::createBooleanWidget(const QJsonObject &obj)
{
  Q_UNUSED(obj);
  QCheckBox *checkBox = new QCheckBox(this);
  connect(checkBox, SIGNAL(toggled(bool)), SLOT(updatePreviewText()));
  return checkBox;
}

void InterfaceWidget::setOptionDefaults()
{
  if (!m_options.contains("userOptions") ||
      !m_options["userOptions"].isObject()) {
    showError(tr("'userOptions' missing, or not an object:\n%1")
              .arg(QString(QJsonDocument(m_options).toJson())));
    return;
  }

  QJsonObject userOptions = m_options["userOptions"].toObject();

  for (QJsonObject::ConstIterator it = userOptions.constBegin(),
       itEnd = userOptions.constEnd(); it != itEnd; ++it) {
    QString label = it.key();
    QJsonValue val = it.value();

    if (!val.isObject()) {
      qWarning() << tr("Error: value must be object for key '%1'.")
                    .arg(label);
      continue;
    }

    QJsonObject obj = val.toObject();
    if (obj.contains("default"))
      setOption(label, obj["default"]);
    else if (m_interfaceScript.debug())
      qWarning() << tr("Default value missing for option '%1'.").arg(label);
  }
}

void InterfaceWidget::setOption(const QString &name,
                                     const QJsonValue &defaultValue)
{
  QString type = lookupOptionType(name);

  if (type == "stringList")
    return setStringListOption(name, defaultValue);
  else if (type == "string")
    return setStringOption(name, defaultValue);
  else if (type == "filePath")
    return setFilePathOption(name, defaultValue);
  else if (type == "integer")
    return setIntegerOption(name, defaultValue);
  else if (type == "boolean")
    return setBooleanOption(name, defaultValue);

  qWarning() << tr("Unrecognized option type '%1' for option '%2'.")
                .arg(type).arg(name);
  return;
}

void InterfaceWidget::setStringListOption(const QString &name,
                                               const QJsonValue &value)
{
  QComboBox *combo = qobject_cast<QComboBox*>(m_widgets.value(name, NULL));
  if (!combo) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad widget type.")
                  .arg(name);
    return;
  }

  if (!value.isDouble() && !value.isString()) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad default value:")
                  .arg(name)
               << value;
    return;
  }

  int index = -1;
  if (value.isDouble())
    index = static_cast<int>(value.toDouble() + 0.5);
  else if (value.isString())
    index = combo->findText(value.toString());

  if (index < 0) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Could not find valid combo entry index from value:")
                  .arg(name)
               << value;
    return;
  }

  combo->setCurrentIndex(index);
}

void InterfaceWidget::setStringOption(const QString &name,
                                           const QJsonValue &value)
{
  QLineEdit *lineEdit = qobject_cast<QLineEdit*>(m_widgets.value(name, NULL));
  if (!lineEdit) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad widget type.")
                  .arg(name);
    return;
  }

  if (!value.isString()) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad default value:")
                  .arg(name)
               << value;
    return;
  }

  lineEdit->setText(value.toString());
}

void InterfaceWidget::setFilePathOption(const QString &name,
                                             const QJsonValue &value)
{
  QtGui::FileBrowseWidget *fileBrowse =
      qobject_cast<QtGui::FileBrowseWidget*>(m_widgets.value(name, NULL));
  if (!fileBrowse) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad widget type.")
                  .arg(name);
    return;
  }

  if (!value.isString()) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad default value:")
                  .arg(name)
               << value;
    return;
  }

  fileBrowse->setFileName(value.toString());
}

void InterfaceWidget::setIntegerOption(const QString &name,
                                            const QJsonValue &value)
{
  QSpinBox *spin = qobject_cast<QSpinBox*>(m_widgets.value(name, NULL));
  if (!spin) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad widget type.")
                  .arg(name);
    return;
  }

  if (!value.isDouble()) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad default value:")
                  .arg(name)
               << value;
    return;
  }

  int intVal = static_cast<int>(value.toDouble() + 0.5);
  spin->setValue(intVal);
}

void InterfaceWidget::setBooleanOption(const QString &name,
                                            const QJsonValue &value)
{
  QCheckBox *checkBox = qobject_cast<QCheckBox*>(m_widgets.value(name, NULL));
  if (!checkBox) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad widget type.")
                  .arg(name);
    return;
  }

  if (!value.isBool()) {
    qWarning() << tr("Error setting default for option '%1'. "
                     "Bad default value:")
                  .arg(name)
               << value;
    return;
  }

  checkBox->setChecked(value.toBool());
}

bool InterfaceWidget::optionString(const QString &option,
                                        QString &value) const
{
  QWidget *widget = m_widgets.value(option, NULL);
  bool retval = false;
  value.clear();

  if (QLineEdit *edit = qobject_cast<QLineEdit*>(widget)) {
    retval = true;
    value = edit->text();
  }
  else if (QComboBox *combo = qobject_cast<QComboBox*>(widget)) {
    retval = true;
    value = combo->currentText();
  }
  else if (QSpinBox *spinbox = qobject_cast<QSpinBox*>(widget)) {
    retval = true;
    value = QString::number(spinbox->value());
  }
  else if (QDoubleSpinBox *dspinbox = qobject_cast<QDoubleSpinBox*>(widget)) {
    retval = true;
    value = QString::number(dspinbox->value());
  }
  else if (QtGui::FileBrowseWidget *fileBrowse
           = qobject_cast<QtGui::FileBrowseWidget*>(widget)) {
    retval = true;
    value = fileBrowse->fileName();
  }

  return retval;
}

QJsonObject InterfaceWidget::collectOptions() const
{
  QJsonObject ret;

  foreach (QString label, m_widgets.keys()) {
    QWidget *widget = m_widgets.value(label, NULL);
    if (QComboBox *combo = qobject_cast<QComboBox*>(widget)) {
      ret.insert(label, combo->currentText());
    }
    else if (QLineEdit *lineEdit = qobject_cast<QLineEdit*>(widget)) {
      QString value(lineEdit->text());
      if (value.isEmpty() && label == "Title")
        value = generateJobTitle();
      ret.insert(label, value);
    }
    else if (QSpinBox *spinBox = qobject_cast<QSpinBox*>(widget)) {
      ret.insert(label, spinBox->value());
    }
    else if (QCheckBox *checkBox = qobject_cast<QCheckBox*>(widget)) {
      ret.insert(label, checkBox->isChecked());
    }
    else if (QtGui::FileBrowseWidget *fileBrowse
             = qobject_cast<QtGui::FileBrowseWidget*>(widget)) {
      ret.insert(label, fileBrowse->fileName());
    }
    else {
      qWarning() << tr("Unhandled widget in collectOptions for option '%1'.")
                    .arg(label);
    }
  }

  return ret;
}

void InterfaceWidget::applyOptions(const QJsonObject &opts)
{
  foreach (const QString &label, opts.keys())
    setOption(label, opts[label]);
}

QString InterfaceWidget::generateJobTitle() const
{
  QString calculation;
  bool haveCalculation(optionString("Calculation Type", calculation));

  QString theory;
  bool haveTheory(optionString("Theory", theory));

  QString basis;
  bool haveBasis(optionString("Basis", basis));

  // Merge theory/basis into theory
  if (haveBasis) {
    if (haveTheory)
      theory += "/";
    theory += basis;
    theory.replace(QRegExp("\\s+"), "");
    haveTheory = true;
  }

  QString formula(m_molecule ? QString::fromStdString(m_molecule->formula())
                             : tr("[no molecule]"));

  return QString("%1%2%3").arg(formula)
      .arg(haveCalculation ? " | " + calculation : QString())
      .arg(haveTheory      ? " | " + theory      : QString());
  }

} // namespace Avogadro
