/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2013-2016 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#ifndef AVOGADRO_QTPLUGINS_TEMPLATETOOLWIDGET_H
#define AVOGADRO_QTPLUGINS_TEMPLATETOOLWIDGET_H

#include <QtWidgets/QWidget>

namespace Avogadro {
namespace QtGui {
class PeriodicTableView;
}

namespace QtPlugins {

namespace Ui {
class TemplateToolWidget;
}

class TemplateToolWidget : public QWidget
{
  Q_OBJECT

public:
  explicit TemplateToolWidget(QWidget *parent_ = 0);
  ~TemplateToolWidget();

  void setAtomicNumber(unsigned char atomicNum);
  unsigned char atomicNumber() const;

  void setCoordination(unsigned char order);
  unsigned char coordination() const;

private slots:
  void elementChanged(int index);
  void updateElementCombo();
  void addUserElement(unsigned char element);
  void elementSelectedFromTable(int element);
  void selectElement(unsigned char element);

private:
  void buildElements();
  void buildBondOrders();
  void saveElements();

  Ui::TemplateToolWidget *m_ui;
  QtGui::PeriodicTableView *m_elementSelector;
  QList<unsigned char> m_defaultElements;
  QList<unsigned char> m_userElements;
  unsigned char m_currentElement;
};

} // namespace QtPlugins
} // namespace Avogadro
#endif // AVOGADRO_QTPLUGINS_TEMPLATETOOLWIDGET_H
