/******************************************************************************
  This source file is part of the Avogadro project.

  This source code is released under the New BSD License, (the "License").
******************************************************************************/

#ifndef AVOGADRO_QTPLUGINS_EDTSURFACE_H
#define AVOGADRO_QTPLUGINS_EDTSURFACE_H

#include "surfaces.h"
#include "bitvector.h"
#include "boolcube.h"
#include <vector>
#include <avogadro/core/avogadrocore.h>
#include <avogadro/core/vector.h>
// for the enum

namespace Avogadro {
namespace Core {
class Cube;
class Molecule;
class Atom;
}
namespace QtGui {
class Molecule;
}

namespace QtPlugins {

typedef struct dataStruct
{
  int pHeight, pWidth, pLength, boxLength;
  Vector3 pMin, pMax, pTran;
  double probeRadius, cutRadius;
  double scaleFactor, fixSf;
  bool ignoreHydrogens;//if we want tis feature, write a way of setting it
} dataStruct; // End struct dataStruct

class EDTSurface
{
public:
  // Constructor
  EDTSurface();

  // Destructor
  virtual ~EDTSurface();

  /*@brief Populates a cube with values generated by doing a Euclidean distance
   *transform on the molecule provided
   *@param mol A pointer to a the molecule from which the cube is to be
   *generated
   *@param surfType an enum class representing the type of surface (VdW, SES,
   *SAS)
   *@returns a pointer to the cube
   */

  Core::Cube* EDTCube(QtGui::Molecule* mol, Core::Cube* cube,
                      Surfaces::Type surfType);

  // The copying over from array to Cube can and should be done in parallel

  /*@brief Populates a cube with values generated by doing a Euclidean distance
   *transform on the molecule provided
   *@param mol A pointer to a the molecule from which the cube is to be
   *generated
   *@param surfType an enum class representing the type of surface (VdW, SES,
   *SAS)
   *@param probeRadius a double representing the molecular radius of the solvent
   *@returns a pointer to the cube
   */

  Core::Cube* EDTCube(QtGui::Molecule* mol, Core::Cube* cube,
                      Surfaces::Type surfType, double probeRadius);
  // Takes a molecule, a surface type and a probeRadius and

  /*@brief Sets a pointer to the desired molecule
   *@param mol a pointer to the molecule to be set
   */

  void setMolecule(QtGui::Molecule* mol);

  /*@brief Sets the probe radius to a desired value (default is 1.4 - water)
   *@param probeRadius The molecular radius of the solvent
   */

  void setProbeRadius(double probeRadius);

  void setCube(Core::Cube* cube);

  double getScaleFactor();

  Vector3 getPTran();

private:
  /*
   *@brief Initializes the data members of the class
   *@param atomType
   *@param surfaceType
   */

  void initPara(bool atomType);

  /*
   *@brief For each atom in the molecule, fills the appropriate voxels
   *@param atomType
   */

  void fillVoxels(bool atomType);

  void fillAtom(int indx);

  void fillVoxelsWaals();

  void fastDistanceMap();

  void buildBoundary();

  void boundBox(bool atomType);

  void computeSphere(unsigned char atomicNumber);

  /*
   *@brief Takes a vector and tells whether or not it's within the bounds of the
   *box
   */
  bool inBounds(Vector3i vec);

  // Can I inline these?

  Vector3i round(Vector3 vec);

  Vector3 promote(Vector3i vec);

  QtGui::Molecule* m_mol;

  Core::Cube* m_cube;

  // These bool arrays should probably be converted into BoolCubes
  BoolCube* inSolid;
  BoolCube* onSurface;

  //We can do things as BitVectors, too
//  BitVector* inSolid;
//  BitVector* onSurface;

  Vector3i* neighbors;

  Vector3i** spheres;//An array of pointers to arrays of vectors representing all the points each sphere
  int* numbersOfVectors;//The number of vectors in each sphere
  bool* computed;//An array of bools that tells us if we've already computed the sphere for this element

  dataStruct* data;

  int numberOfSurfaceVoxels;
  Vector3i* surfaceVoxels;

  int numberOfInnerVoxels;//this is a debugging value
  int numberOfInBoundsVoxels;//this is a debugging value
  int alreadyInSolid;//this is a debugging value
}; // End class EDTSurface

} // End namespace QtPlugins
} // End namespace Avogadro

#endif
