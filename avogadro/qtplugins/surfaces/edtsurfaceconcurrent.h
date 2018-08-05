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

#ifndef EDTSURFACECONCURRENT_H
#define EDTSURFACECONCURRENT_H

#include <QtCore/QFuture>
#include <QtCore/QFutureWatcher>
#include <QtCore/QObject>

typedef struct volumePixel
{
  int atomId;
  float distance;
  bool inOut;
  bool isBound;
  bool isDone;
}volumePixel;

typedef struct dataStruct{
  Vector3 pTran;
  int boxLength;
  double probeRadius;
  double fixSf;
  double scaleFactor;
  Vector3 pMin, pMax;
  int pHeight, pWidth, pLength;
  int widXz[13];
  int* deptY[13];
  double cutRadius;
  int positIn, positOut, eliminate;
  int certificate;
  int totalSurfaceVox;
  int totalInnerVox;
  Vector3i *inArray, *outArray
}dataStruct;//End struct dataStruct

typedef struct atomStruct{
  Core::Atom *atom;
  volumePixel*** volumePixels;
  int index;
  bool atomType;
}atomStruct;

typedef struct subCube{
  Core::Cube *cube;
  volumePixel** volumePixelsRow;
  int pWidth;
  int pHeight;
  int index;
}subCube;

namespace Avogadro {

namespace Core {
class Cube;
class Molecule;
class Atom;
class EDTSurface;
}

namespace QtPlugins {

class EDTSurfaceConcurrent{
public:
  //Constructor
  EDTSurface();

  //Destructor
  virtual ~EDTSurface();

  /*@brief Populates a cube with values generated by doing a Euclidean distance
  *transform on the molecule provided
  *@param mol A pointer to a the molecule from which the cube is to be generated
  *@param surfType an enum class representing the type of surface (VdW, SES, SAS)
  *@returns a pointer to the cube
  */

  Core::Cube *EDTCube(Core::Molecule *mol, Surfaces::Type surfType);

  //The copying over from array to Cube can and should be done in parallel

  /*@brief Populates a cube with values generated by doing a Euclidean distance
  *transform on the molecule provided
  *@param mol A pointer to a the molecule from which the cube is to be generated
  *@param surfType an enum class representing the type of surface (VdW, SES, SAS)
  *@param probeRadius a double representing the molecular radius of the solvent
  *@returns a pointer to the cube
  */

  Core::Cube *EDTCube(Core::Molecule *mol, Surfaces::Type surfType, double probeRadius);
  // Takes a molecule, a surface type and a probeRadius and

  /*@brief Sets a pointer to the desired molecule
  *@param mol a pointer to the molecule to be set
  */

  void setMolecule(Core::Molecule *mol);

  /*@brief Sets the probe radius to a desired value (default is 1.4 - water)
  *@param probeRadius The molecular radius of the solvent
  */

  void setProbeRadius(double probeRadius);

  /*@brief Copies cube from volumePixel array into Cube object
  */

  QFutureWatcher<void> & watcher() { return m_watcher; }

  private Q_SLOTS:
      /**
       * Slot to set the cube data once Qt Concurrent is done
       */
     void calculationComplete();

private:

  void initPara(bool atomType, bool bType, int surfaceType);
  //This can be done concurrently, but maybe doesn't need to be
  void fillVoxels(bool atomType);
  //This can (and should be done concurrently)
  void fillVoxelsConcurrent(subCube* someVolumePixels);
  //part 1 will take an atomStruct
  //part 2 will take a subCube

  void fillAtom(int indx);
  //This cannot be done concurrently, but there's nor eason for it to be
  void fillAtomWaals(int indx);
  //This cannot be done concurrently, but there's no reason for it to be
  void fillVoxelsWaals(bool atomType);
  //This can and should be done concurrently
  void fillVoxelsWaalsConcurrent(bool atomType);

  void fastOneShell(int* inNum, int* allocOut, Vector3i*** boundPoint,
                    int* outNum, int* elimi);
  //This cannot be done concurrently, we run into issues breaking up the cube
  void fastDistanceMap();
  //This can be done concurrently except for the call to fastOneShell
  void buildBoundary();
  //This cannot be done concurrently (boundary issues),
  //But that's fine because it only happens once
  void boundBox(bool atomType);
  //This can be done concurrently, but maybe doesn't need to be

  void boundingAtom(bool bType);
  //This can be done concurrently but maybe doesn't need to be

  //So mostly, we need a way to make fillVoxels and fillVoxelsWaals run concurrently

  Vector3i vectorFromArray(int* array);
    // Takes an array of integers and returns a vector3i

  int detail(unsigned char atomicNumber);
    // Takes an atomic number and returns an index for rasRad

  void copyCube();

  void copyCubeConcurrent(subCube *someVolumePixels);

  Molecule* m_mol;

  Cube* m_cube;

  volumePixel*** volumePixels;

  dataStruct *data;

  Q_SIGNALS:

    QFuture<void> m_future;
    QFutureWatcher<void> m_watcher;
    Cube *m_cube; // Cube to put the results into
    QVector<atomStruct> m_atomVector;
    QVector<subCube> m_subCubeVector;

  //so the concurrent versions of things that run on a cube should rewritten to run on a subCube
  //and the concurrent versions of things that run on a molecule should be rewritten to run on an atomStruct

};
}
}

#endif // EDTSURFACECONCURRENT_H