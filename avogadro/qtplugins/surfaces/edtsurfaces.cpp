/******************************************************************************
  This source file is part of the Avogadro project.

  This source code is released under the New BSD License, (the "License").
******************************************************************************/

#include "edtsurface.h"

#include <Eigen/Dense>

#include <avogadro/core/cube.h>
#include <avogadro/core/elementdata.h>
#include <avogadro/core/molecule.h>

#define VWS 0
#define MS 1
#define SAS 2
#define SES 3

#define X 0
#define Y 1
#define Z 2

#define I 0
#define J 1
#define K 2

static int nb[26][3] = { 1,  0,  0,  -1, 0,  0,  0, 1,  0,  0,  -1, 0,  0,
                         0,  1,  0,  0,  -1, 1,  1, 0,  1,  -1, 0,  -1, 1,
                         0,  -1, -1, 0,  1,  0,  1, 1,  0,  -1, -1, 0,  1,
                         -1, 0,  -1, 0,  1,  1,  0, 1,  -1, 0,  -1, 1,  0,
                         -1, -1, 1,  1,  1,  1,  1, -1, 1,  -1, 1,  -1, 1,
                         1,  1,  -1, -1, -1, -1, 1, -1, 1,  -1, -1, -1, -1 };

static bool bTypes[4] = { false, true, true, true };

static bool atomTypes[4] = { false, false, true, false };

// static double rasrad[
// ]={1.872,1.872,1.507,1.4,1.848,1.1,1.88,1.872,1.507,1.948,1.5,
// 1.4, 1.1};//lsms
static double rasRad[] = { 1.90, 1.88, 1.63, 1.48, 1.78, 1.2, 1.87,
                           1.96, 1.63, 0.74, 1.8,  1.48, 1.2 }; // liang
//                         ca   c    n    o    s    h   p   cb    ne   fe  other
//                         ox  hx

using namespace Avogadro::Core;

namespace Avogadro {
namespace QtPlugins {

// Constructor
EDTSurfaces::EDTSurfaces()
{
  int i;

  //	data->pTran(0.0,0.0,0.0);
  data->pTran(X) = 0.0;
  data->pTran(Y) = 0.0;
  data->pTran(Z) = 0.0;

  data->boxLength = 0;
  data->probeRadius = 0;
  data->fixSf = 0;
  data->scaleFactor = 0;
  //  data->pMin(0.0,0.0,0.0);
  data->pMin(X) = 0.0;
  data->pMin(Y) = 0.0;
  data->pMin(Z) = 0.0;
  //	data->pMax(0.0,0.0,0.0);
  data->pMax(X) = 0.0;
  data->pMax(Y) = 0.0;
  data->pMax(Z) = 0.0;

  data->pHeight = 0;
  data->pWidth = 0;
  data->pLength = 0;

  for (i = 0; i < 13; i++) {
    data->widXz[i] = 0;
    data->deptY[i] = NULL;
  }
  data->cutRadius = 0;
  volumePixels = NULL;
}

// Destructor
EDTSurfaces::~EDTSurfaces()
{
  int i, j;

  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      free(volumePixels[i][j]);
      volumePixels[i][j] = NULL;
    }
    free(volumePixels[i]);
    volumePixels[i] = NULL;
  }
  free(volumePixels);
  volumePixels = NULL;

  for (i = 0; i < 13; i++) {
    free(data->deptY[i]);
    data->deptY[i] = NULL;
  }
  free(data->deptY);
  //	data->deptY = NULL;

  free(data->widXz);
  //	data->widXz = NULL;
}

// Takes a molecule and a surface type and returns a cube

Core::Cube* EDTSurfaces::EDTCube()
{
  return m_Cube;
}

void EDTSurfaces::fastDistanceMap()
{
  int i, j, k;
  Vector3i ijk;
  data->totalSurfaceVox = 0;
  data->totalInnerVox = 0;

  Vector3i*** boundPoint;
  boundPoint = new Vector3i**[data->pLength];

  for (i = 0; i < data->pLength; i++) {
    boundPoint[i] = new Vector3i*[data->pWidth];
    for (j = 0; j < data->pWidth; j++) {
      boundPoint[i][j] = new Vector3i[data->pHeight];
      for (k = 0; k < data->pHeight; k++) {
        volumePixels[i][j][k].isDone = false;
        if (volumePixels[i][j][k].inOut) {
          if (volumePixels[i][j][k].isBound) {
            data->totalSurfaceVox++;
            ijk(I) = i;
            ijk(J) = j;
            ijk(K) = k;
            boundPoint[i][j][k] = ijk;
            volumePixels[i][j][k].distance = 0;
            volumePixels[i][j][k].isDone = true;
          } else {
            data->totalInnerVox++;
          }
        }
      }
    }
  }

  int allocIn = int(1.2 * data->totalSurfaceVox);
  int allocOut = int(1.2 * data->totalSurfaceVox);
  if (allocIn > data->totalInnerVox)
    allocIn = data->totalInnerVox;
  if (allocIn < data->totalSurfaceVox)
    allocIn = data->totalSurfaceVox;
  if (allocOut > data->totalInnerVox)
    allocOut = data->totalInnerVox;
  data->inArray = new Vector3i[allocIn];
  data->outArray = new Vector3i[allocOut];
  data->positIn = 0;
  data->positOut = 0;

  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      for (k = 0; k < data->pHeight; k++) {
        if (volumePixels[i][j][k].isBound) {
          ijk(I) = i;
          ijk(J) = j;
          ijk(K) = k;
          data->inArray[data->positIn] = ijk;
          data->positIn++;
          volumePixels[i][j][k].isBound = false; // as flag of data->outArray
        }
      }
    }
  }
  data->certificate = data->totalInnerVox;
  ///////////////////////////////////////////////////
  // if(type==0)//do part
  //{
  // type == 0 when we're not doing depth
  do {
    fastOneShell(&data->positIn, &allocOut, boundPoint, &data->positOut, &data->eliminate);
    //	printf("%d %d %d %d
    //%d\n",data->positIn,allocOut,data->positOut,data->totalSurfaceVox,data->totalInnerVox);
    data->certificate -= data->eliminate;

    /*		for(i=0;i<data->positOut;i++)
                            {
                              data->inArray[i](X)=data->outArray[i](X);
                              data->inArray[i](Y)=data->outArray[i](Y);
                              data->inArray[i](Z)=data->outArray[i](Z);
                            }
                            data->positIn=data->positOut;
                    //new code only less dist
            */
    data->positIn = 0;
    for (i = 0; i < data->positOut; i++) {
      volumePixels[data->outArray[i](X)][data->outArray[i](Y)][data->outArray[i](Z)].isBound =
        false;
      if (volumePixels[data->outArray[i](X)][data->outArray[i](Y)][data->outArray[i](Z)]
            .distance <= 1.02 * data->cutRadius) {
        data->inArray[data->positIn] = data->outArray[i];
        data->positIn++;
      }
      if (data->positIn >= allocIn) {
        allocIn *= 2;
        if (allocIn > data->totalInnerVox)
          allocIn = data->totalInnerVox;
        data->inArray = (Vector3i*)realloc(data->inArray, allocIn * sizeof(Vector3i));
      }
    }
  } while (data->positIn != 0);

  //}
  /*else if(type==1)//do all
  {
          Vector3i *tpoint;
          do {

                  fastOneShell( &data->positIn, &allocOut, boundPoint,
  &data->positOut,&data->eliminate);//data->inArray, data->outArray, data->certificate-=data->eliminate;
  //	/*
  //			for(i=0;i<data->positOut;i++)
  //			{
  //
  volumePixels[data->outArray[i](X)][data->outArray[i](Y)][data->outArray[i](Z)].isBound=false;
  //			  data->inArray[i](X)=data->outArray[i](X);
  //			  data->inArray[i](Y)=data->outArray[i](Y);
  //			  data->inArray[i](Z)=data->outArray[i](Z);
  //			}
                          tpoint=data->inArray;
                          data->inArray=data->outArray;
                          data->outArray=tpoint;
                          data->positIn=data->positOut;
                          int alloctmp;
                          alloctmp=allocIn;
                          allocIn=allocOut;
                          allocOut=alloctmp;
                          for(i=0;i<data->positIn;i++)
                                  volumePixels[data->inArray[i](X)][data->inArray[i](Y)][data->inArray[i](Z)].isBound=false;
  //			*/
  // new code only less dist
  /*
  data->positIn=0;
  for(i=0;i<data->positOut;i++)
  {
          volumePixels[data->outArray[i](X)][data->outArray[i](Y)][data->outArray[i](Z)].isBound=false;
          if(volumePixels[data->outArray[i](X)][data->outArray[i](Y)][data->outArray[i](Z)].distance<=1.0*data->cutRadius)
          {
                  data->inArray[data->positIn](X)=data->outArray[i](X);
                  data->inArray[data->positIn](Y)=data->outArray[i](Y);
                  data->inArray[data->positIn](Z)=data->outArray[i](Z);
                  data->positIn++;
          }
  }
  */

  //	}
  //	while(data->positIn!=0);
  //	while(data->positOut!=0);
  //}
  // while(data->positOut!=0);
  if (data->certificate != 0) {
    //	printf("wrong number\n");
  }

  free(data->inArray);
  free(data->outArray);

  double cutsf = data->scaleFactor - 0.5;
  if (cutsf < 0)
    cutsf = 0;
  //	 cutsf=100000000;
  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      for (k = 0; k < data->pHeight; k++) {
        volumePixels[i][j][k].isBound = false;
        // ses solid
        if (volumePixels[i][j][k].inOut) {
          if (!volumePixels[i][j][k].isDone ||
              (volumePixels[i][j][k].isDone &&
               volumePixels[i][j][k].distance >=
                 data->cutRadius - 0.50 / (0.1 + cutsf)) // 0.33  0.75/data->scaleFactor
          ) {
            volumePixels[i][j][k].isBound = true;
            // new add
            if (volumePixels[i][j][k].isDone)
              volumePixels[i][j][k].atomId =
                volumePixels[boundPoint[i][j][k](X)][boundPoint[i][j][k](Y)]
                            [boundPoint[i][j][k](Z)]
                              .atomId;
          }
        }
      }
    }
  }

  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      delete[] boundPoint[i][j];
    }
    delete[] boundPoint[i];
  }
  delete[] boundPoint;
}

void EDTSurface::fastOneShell(int* inNum, int* allocOut, Vector3i*** boundPoint,
                              int* outNum, int* elimi)
{
  int i, number, data->positOut;
  int tx, ty, tz;
  Vector3i dxyz;
  Vector3i txyz;
  int dx, dy, dz;
  int data->eliminate = 0;
  float squre;
  data->positOut = 0;
  number = *inNum;
  if (number == 0)
    return;
  // new code
  int j;
  Vector3i tnv;
  //	Vector3i tnv;
  for (i = 0; i < number; i++) {
    if (data->positOut >= (*allocOut) - 6) {
      (*allocOut) = int(1.2 * (*allocOut));
      if (*allocOut > data->totalInnerVox)
        *allocOut = data->totalInnerVox;
      data->outArray = (Vector3i*)realloc(data->outArray, (*allocOut) * sizeof(Vector3i));
    }
    txyz = data->inArray[i];

    for (j = 0; j < 6; j++) {

      tnv = txyz + vectorFromArray(nb[j]);

      if (tnv(X) < data->pLength && tnv(X) > -1 && tnv(Y) < data->pWidth && tnv(Y) > -1 &&
          tnv(Z) < data->pHeight && tnv(Z) > -1 &&
          volumePixels[tnv(X)][tnv(Y)][tnv(Z)].inOut &&
          !volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone) {
        boundPoint[tnv(X)][tnv(Y)][tnv(Z)] =
          boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
        dxyz = tnv - boundPoint[txyz(X)][txyz(Y)][txyz(Z)];

        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance = dxyz.norm();
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone = true;
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound = true;

        data->outArray[data->positOut] = tnv;

        data->positOut++;
        data->eliminate++;
      } else if (tnv(X) < data->pLength && tnv(X) > -1 && tnv(Y) < data->pWidth &&
                 tnv(Y) > -1 && tnv(Z) < data->pHeight && tnv(Z) > -1 &&
                 volumePixels[tnv(X)][tnv(Y)][tnv(Z)].inOut &&
                 volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone) {
        dxyz = tnv - boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
        if (squre < volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance) {
          boundPoint[tnv(X)][tnv(Y)][tnv(Z)] =
            boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
          volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance = dxyz.norm();
          if (!volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound) {
            volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound = true;
            data->outArray[data->positOut] = tnv;
            data->positOut++;
          }
        }
      }
    }
  }

  for (i = 0; i < number; i++) {
    if (data->positOut >= (*allocOut) - 12) {
      (*allocOut) = int(1.2 * (*allocOut));
      if (*allocOut > data->totalInnerVox)
        *allocOut = data->totalInnerVox;
      data->outArray = (Vector3i*)realloc(data->outArray, (*allocOut) * sizeof(Vector3i));
    }
    txyz = data->inArray[i];
    for (j = 6; j < 18; j++) {
      tnv = txyz + vectorFromArray(nb[j]);
      ;
      if (tnv(X) < data->pLength && tnv(X) > -1 && tnv(Y) < data->pWidth && tnv(Y) > -1 &&
          tnv(Z) < data->pHeight && tnv(Z) > -1 &&
          volumePixels[tnv(X)][tnv(Y)][tnv(Z)].inOut &&
          !volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone) {
        boundPoint[tnv(X)][tnv(Y)][tnv(Z)] =
          boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
        dxyz = tnv - boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance = dxyz.norm();
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone = true;
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound = true;
        data->outArray[data->positOut] = tnv;
        data->positOut++;
        data->eliminate++;
      } else if (tnv(X) < data->pLength && tnv(X) > -1 && tnv(Y) < data->pWidth &&
                 tnv(Y) > -1 && tnv(Z) < data->pHeight && tnv(Z) > -1 &&
                 volumePixels[tnv(X)][tnv(Y)][tnv(Z)].inOut &&
                 volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone) {
        dxyz = tnv - boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
        squre = dxyz.norm();
        if (squre < volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance) {
          boundPoint[tnv(X)][tnv(Y)][tnv(Z)] =
            boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
          volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance = float(squre);
          if (!volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound) {
            volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound = true;
            data->outArray[data->positOut] = tnv;
            data->positOut++;
          }
        }
      }
    }
  }
  for (i = 0; i < number; i++) {
    if (data->positOut >= (*allocOut) - 9) {
      (*allocOut) = int(1.2 * (*allocOut));
      if (*allocOut > data->totalInnerVox)
        *allocOut = data->totalInnerVox;
      data->outArray = (Vector3i*)realloc(data->outArray, (*allocOut) * sizeof(Vector3i));
    }
    txyz = data->inArray[i];

    for (j = 18; j < 26; j++) {
      tnv = txyz + vectorFromArray(nb[j]);
      ;
      if (tnv(X) < data->pLength && tnv(X) > -1 && tnv(Y) < data->pWidth && tnv(Y) > -1 &&
          tnv(Z) < data->pHeight && tnv(Z) > -1 &&
          volumePixels[tnv(X)][tnv(Y)][tnv(Z)].inOut &&
          !volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone) {
        boundPoint[tnv(X)][tnv(Y)][tnv(Z)] =
          boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
        dxyz = tnv - boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance = dxyz.norm();
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone = true;
        volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound = true;
        data->outArray[data->positOut] = tnv;
        data->positOut++;
        data->eliminate++;
      } else if (tnv(X) < data->pLength && tnv(X) > -1 && tnv(Y) < data->pWidth &&
                 tnv(Y) > -1 && tnv(Z) < data->pHeight && tnv(Z) > -1 &&
                 volumePixels[tnv(X)][tnv(Y)][tnv(Z)].inOut &&
                 volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isDone) {

        dxyz = tnv - boundPoint[txyz(X)][txyz(Y)][txyz(Z)];

        if (squre < volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance) {
          boundPoint[tnv(X)][tnv(Y)][tnv(Z)] =
            boundPoint[txyz(X)][txyz(Y)][txyz(Z)];
          volumePixels[tnv(X)][tnv(Y)][tnv(Z)].distance = dxyz.norm();
          if (!volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound) {
            volumePixels[tnv(X)][tnv(Y)][tnv(Z)].isBound = true;
            data->outArray[data->positOut] = tnv;
            data->positOut++;
          }
        }
      }
    }
  }

  *outNum = data->positOut;
  *elimi = data->eliminate;
}

void EDTSurfaces::fillAtom(int indx)
{
  int cx, cy, cz;
  Vector3i cxyz;
  int ox, oy, oz;
  Vector3i oxyz;
  Vector3 cp;

  Array<Vector3> positions = m_mol->atomPositions3d();

  cp = (positions[indx] + data->pTran) * data->scaleFactor;
  cx = int(cp(X) + 0.5);
  cy = int(cp(Y) + 0.5);
  cz = int(cp(Z) + 0.5);

  cxyz(X) = cx;
  cxyz(Y) = cy;
  cxyz(Z) = cz;
  /*Translate, scale, and round positions*/

  Index index = indx;
  Atom current = m_mol->atom(index);
  int at = detail(current.atomicNumber());
  int i, j, k;
  int ii, jj, kk;
  int mi, mj, mk;
  Vector3i mijk;
  int si, sj, sk;
  Vector3i sijk;
  int tIndex;
  int nIndex = 0;
  for (i = 0; i < data->widXz[at]; i++) {
    for (j = 0; j < data->widXz[at]; j++) {
      if (data->deptY[at][nIndex] != -1) {

        for (ii = -1; ii < 2; ii++) {
          for (jj = -1; jj < 2; jj++) {
            for (kk = -1; kk < 2; kk++) {
              if (ii != 0 && jj != 0 && kk != 0) {
                mi = ii * i;
                mk = kk * j;
                for (k = 0; k <= data->deptY[at][nIndex]; k++) {
                  mj = k * jj;
                  mijk(I) = mi;
                  mijk(J) = mj;
                  mijk(K) = mk;
                  sijk = cxyz + mijk;
                  if (sijk(I) < 0 || sijk(J) < 0 || sijk(K) < 0 ||
                      sijk(I) >= data->pLength || sijk(J) >= data->pWidth ||
                      sijk(K) >= data->pHeight) {
                    continue;
                  }

                  else {
                    if (volumePixels[sijk(I)][sijk(J)][sijk(K)].inOut ==
                        false) {
                      volumePixels[sijk(I)][sijk(J)][sijk(K)].inOut = true;
                      volumePixels[sijk(I)][sijk(J)][sijk(K)].atomId = indx;
                    }
                    // no atomic info to each Vector3i change above line
                    //*
                    else if (volumePixels[sijk(I)][sijk(J)][sijk(K)].inOut) {
                      tIndex = volumePixels[sijk(I)][sijk(J)][sijk(K)].atomId;

                      cp = (positions[tIndex] + data->pTran) * data->scaleFactor;
                      // Translating and scaling

                      ox = int(cp(X) + 0.5) - sijk(I);
                      oy = int(cp(Y) + 0.5) - sijk(J);
                      oz = int(cp(Z) + 0.5) - sijk(K);
                      oxyz(X) = ox;
                      oxyz(Y) = oy;
                      oxyz(Z) = oz;
                      // Rounding to the nearest integer

                      if (mijk.squaredNorm() < oxyz.squaredNorm())
                        volumePixels[sijk(I)][sijk(J)][sijk(K)].atomId = indx;
                    }
                    //	*/
                  } // k
                }   // else
              }     // if
            }       // kk
          }         // jj
        }           // ii

      } // if
      nIndex++;
    } // j
  }   // i
}
// sas use inOut
void EDTSurfaces::fillVoxels(bool atomType)
{

  int i, j, k;

  int numberOfAtoms = m_mol->atomCount();

  for (i = 0; i < numberOfAtoms; i++) {
    Index index = i;
    Atom current = m_mol->atom(index);
    if (!atomType || current.atomicNumber() != 1)
      fillAtom(i, mol);
    //			totalNumber++;
  }
  //	printf("%d\n",totalNumber);
  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      for (k = 0; k < data->pHeight; k++) {
        if (volumePixels[i][j][k].inOut) {
          volumePixels[i][j][k].isDone = true;
        }
      }
    }
  }
}
// use isDone
void EDTSurfaces::fillVoxelsWaals(bool atomType)
{
  int i, j, k;

  int numberOfAtoms = m_mol->atomCount();

  for (i = 0; i < numberOfAtoms; i++) {
    Index index = i;
    Atom current = m_mol->atom(index);
    if (!atomType || current.atomicNumber() != 1) {
      fillAtomWaals(i, mol);
    }
  }
}

void EDTSurfaces::fillAtomWaals(int indx)
{
  int cx, cy, cz;
  Vector3i cxyz;
  int ox, oy, oz;
  Vector3i oxyz;
  Vector3 cp;

  int numberOfAtoms = m_mol->atomCount();
  Array<Vector3> positions = m_mol->atomPositions3d();
  Atom current = m_mol->atom(indx);

  cp = (positions[indx] + data->pTran) * data->scaleFactor;
  // Translating and scaling

  cx = int(cp(X) + 0.5);
  cy = int(cp(Y) + 0.5);
  cz = int(cp(Z) + 0.5);

  cxyz(X) = cx;
  cxyz(Y) = cy;
  cxyz(Z) = cz;
  // Rounding to the nearest integer

  int at = detail(current.atomicNumber());
  int i, j, k;
  int ii, jj, kk;
  Vector3i mijk;
  Vector3i sijk;
  int tIndex;
  int nIndex = 0;

  for (i = 0; i < data->widXz[at]; i++) {
    for (j = 0; j < data->widXz[at]; j++) {
      if (data->deptY[at][nIndex] != -1) {
        for (ii = -1; ii < 2; ii++) {
          for (jj = -1; jj < 2; jj++) {
            for (kk = -1; kk < 2; kk++) {
              if (ii != 0 && jj != 0 && kk != 0) {
                mijk(I) = ii * i;
                mijk(K) = kk * j;
                for (k = 0; k <= data->deptY[at][nIndex]; k++) {
                  mijk(J) = jj * k;
                  sijk = cxyz + mijk;
                  if (sijk(I) < 0 || sijk(J) < 0 || sijk(K) < 0) {
                    continue;
                  }

                  else {
                    if (volumePixels[sijk(I)][sijk(J)][sijk(K)].isDone ==
                        false) {
                      volumePixels[sijk(I)][sijk(J)][sijk(K)].isDone = true;
                      volumePixels[sijk(I)][sijk(J)][sijk(K)].atomId = indx;
                    }
                    // with atomic info change above line
                    //*
                    else if (volumePixels[sijk(I)][sijk(J)][sijk(K)].isDone) {
                      tIndex = volumePixels[sijk(I)][sijk(J)][sijk(K)].atomId;
                      cp = (positions[tIndex] + data->pTran) * data->scaleFactor;
                      // Translating and scaling
                      oxyz = cxyz - sijk;
                      if (mijk.squaredNorm() < oxyz.squaredNorm())
                        volumePixels[sijk(I)][sijk(J)][sijk(K)].atomId = indx;
                    }
                    //	 */
                  } // else
                }   // k

              } // if
            }   // kk
          }     // jj
        }       // ii

      } // if
      nIndex++;
    } // j
  }   // i
}

void EDTSurfaces::buildBoundary()
{
  int i, j, k;
  int ii;
  bool flagBound;

  //! Thread this!

  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pHeight; j++) {
      for (k = 0; k < data->pWidth; k++) {
        if (volumePixels[i][k][j].inOut) {
          // 6 neighbors
          //					if(( k-1>-1 && !volumePixels[i][k-1][j].inOut)
          //|| ( k+1<data->pWidth &&!volumePixels[i][k+1][j].inOut)
          //					|| ( j-1>-1 && !volumePixels[i][k][j-1].inOut)
          //|| ( j+1<data->pHeight &&!volumePixels[i][k][j+1].inOut)
          //					|| ( i-1>-1 && !volumePixels[i-1][k][j].inOut)
          //|| ( i+1<data->pLength &&!volumePixels[i+1][k][j].inOut))
          //						volumePixels[i][k][j].isBound=true;
          //	/*
          // 26 neighbors
          flagBound = false;
          ii = 0;
          while (!flagBound && ii < 26) {
            if (i + vectorFromArray(nb[ii])(X) > -1 &&
                i + vectorFromArray(nb[ii])(X) < data->pLength &&
                k + vectorFromArray(nb[ii])(Y) > -1 &&
                k + vectorFromArray(nb[ii])(Y) < data->pWidth &&
                j + vectorFromArray(nb[ii])(Z) > -1 &&
                j + vectorFromArray(nb[ii])(Z) < data->pHeight &&
                !volumePixels[i + vectorFromArray(nb[ii])(X)]
                             [k + vectorFromArray(nb[ii])(Y)]
                             [j + vectorFromArray(nb[ii])(Z)]
                               .inOut) {
              volumePixels[i][k][j].isBound = true;
              flagBound = true;
            } else
              ii++;
          }
          //		*/
        }
      }
    }
  }
}

void EDTSurfaces::boundBox(bool atomType)
{
  /**
   *Finds the bound box of the sequence of atoms
   *The smallest rectangular prism that contains all the atoms
   *@param minPoint A pointer to a vector representing the minimum point
   *@param maxPoint A pointer to a vector representing the maximum point
   **/

  int i;

  int numberOfAtoms = m_mol->atomCount();
  Array<Vector3> positions = m_mol->atomPositions3d();

  data->pMin(X) = 100000;
  data->pMin(Y) = 100000;
  data->pMin(Z) = 100000;
  data->pMax(X) = -100000;
  data->pMax(Y) = -100000;
  data->pMax(Z) = -100000;

  for (i = 0; i < numberOfAtoms; i++) {
    Atom current = m_mol->atom(i);
    if (!atomType || current.atomicNumber() != 1) {
      if (positions[i](X) < data->pMin(X))
        data->pMin(X) = positions[i](X);
      if (positions[i](Y) < data->pMin(Y))
        data->pMin(Y) = positions[i](Y);
      if (positions[i](Z) < data->pMin(Z))
        data->pMin(Z) = positions[i](Z);
      if (positions[i](X) > data->pMax(X))
        data->pMax(X) = positions[i](X);
      if (positions[i](Y) > data->pMax(Y))
        data->pMax(Y) = positions[i](Y);
      if (positions[i](Z) > data->pMax(Z))
        data->pMax(Z) = positions[i](Z);
    }
  }
}

void EDTSurfaces::initPara(bool atomType, bool bType)
{
  int i, j;
  int data->fixSf = 4;
  double fMargin = 2.5;
  int data->pLength, data->pWidth;

  boundBox(atomType, mol);
  if (bType == false) {
    data->pMin(X) -= fMargin;
    data->pMin(Y) -= fMargin;
    data->pMin(Z) -= fMargin;
    data->pMax(X) += fMargin;
    data->pMax(Y) += fMargin;
    data->pMax(Z) += fMargin;
  } else {
    data->pMin(X) -= data->probeRadius + fMargin;
    data->pMin(Y) -= data->probeRadius + fMargin;
    data->pMin(Z) -= data->probeRadius + fMargin;
    data->pMax(X) += data->probeRadius + fMargin;
    data->pMax(Y) += data->probeRadius + fMargin;
    data->pMax(Z) += data->probeRadius + fMargin;
  }

  data->pTran = -data->pMin;

  // data->pTran is the vector to get us to our minimum x, minimum y, and minimum z
  // points

  data->scaleFactor = data->pMax(X) - data->pMin(X);
  if ((data->pMax(Y) - data->pMin(Y)) > data->scaleFactor)
    data->scaleFactor = data->pMax(Y) - data->pMin(Y);
  if ((data->pMax(Z) - data->pMin(Z)) > data->scaleFactor)
    data->scaleFactor = data->pMax(Z) - data->pMin(Z);

  // data->scaleFactor is the maximum distance between our mins and maxes

  data->scaleFactor = (data->boxLength - 1.0) / double(data->scaleFactor);
  ///////////////////////////add this automatically first fix sf then fix
  ///data->boxLength
  //	/*
  data->boxLength = int(data->boxLength * data->fixSf / data->scaleFactor);
  data->scaleFactor = data->fixSf;
  double threshBox = 300;
  if (data->boxLength > threshBox) {
    double sfThresh = threshBox / double(data->boxLength);
    data->boxLength = int(threshBox);
    data->scaleFactor = data->scaleFactor * sfThresh;
  }
  //	*/

  data->pLength = int(ceil(data->scaleFactor * (data->pMax(X) - data->pMin(X))) + 1);
  data->pWidth = int(ceil(data->scaleFactor * (data->pMax(Y) - data->pMin(Y))) + 1);
  data->pHeight = int(ceil(data->scaleFactor * (data->pMax(Z) - data->pMin(Z))) + 1);

  if (data->pLength > data->boxLength)
    data->pLength = data->boxLength;
  if (data->pWidth > data->boxLength)
    data->pWidth = data->boxLength;
  if (data->pHeight > data->boxLength)
    data->pHeight = data->boxLength;

  boundingAtom(bType);
  data->cutRadius = data->probeRadius * data->scaleFactor;
}

void EDTSurfaces::boundingAtom(bool bType)
{
  int i, j, k;
  double tRadius[13];
  double tXz, tDepth, sRadius;
  int indx;
  for (i = 0; i < 13; i++) {
    if (data->deptY[i] != NULL)
      free(data->deptY[i]);
  }

  for (i = 0; i < 13; i++) {
    if (bType == false)
      tRadius[i] = rasRad[i] * data->scaleFactor + 0.5;
    else
      tRadius[i] = (rasRad[i] + data->probeRadius) * data->scaleFactor + 0.5;

    // Multiply by data->scaleFactor
    // Maybe add data->probeRadius first

    sRadius = tRadius[i] * tRadius[i];
    // Square that
    data->widXz[i] = int(tRadius[i]) + 1;
    data->deptY[i] = new int[data->widXz[i] * data->widXz[i]];
    indx = 0;
    for (j = 0; j < data->widXz[i]; j++) {
      for (k = 0; k < data->widXz[i]; k++) {
        tXz = j * j + k * k;
        if (tXz > sRadius) {
          data->deptY[i][indx] = -1;
        } else {
          tDepth = sqrt(sRadius - tXz);
          data->deptY[i][indx] = int(tDepth + 0.0);
        }
        indx++;
      }
    }
  }
}

Vector3i EDTSurfaces::vectorFromArray(int* array)
{
  Vector3i vec;
  vec(0) = array[0];
  vec(1) = array[1];
  vec(2) = array[2];
  return vec;
}

int EDTSurfaces::detail(unsigned char atomicNumber)
{
  // Takes an atomic number and returns the index for rasRad
  switch (atomicNumber) {
    // Hydrogen
    case 1:
      return 5;
    // Carbon
    case 6:
      return 1;
    // Nitrogen
    case 7:
      return 2;
    // Oxygen
    case 8:
      return 3;
    // Phosphorous
    case 15:
      return 6;
    // Sulfur
    case 16:
      return 4;
    // Other
    default:
      return 10;
  }
}

int EDTSurfaces::setMolecule(Molecule *mol){
  m_mol = mol;
  return;
}

bool EDTSurfaces::VanDerWaalsSurface(){

  int i, j, k;

  this->setMolecule(mol);
  //Set molecule

  this->initPara(atomTypes[VWS], bTypes[VWS]);
  // Initialize everything

  this->fillVoxels(atomTypes[VWS]);

  this->buildBoundary();

  m_Cube = new Cube();
  // Initialize cube

  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      for (k = 0; k < data->pHeight; k++) {
        m_Cube->setValue(i, j, k, volumePixels[i][j][k].distance);
      }
      delete[] volumePixels[i][j];
    }
    delete[] volumePixels[i];
  }
  delete[] volumePixels;
  // Copy VolumePixels array into Cube and delete volumePixels
  return true;
}//end VanDerWaalsSurface

bool EDTSurfaces::SolventAccessibleSurface(){

  int i, j, k;

  this->setMolecule(mol);
  //Set molecule

  this->initPara(atomTypes[SAS], bTypes[SAS]);
  // Initialize everything

  this->fastDistanceMap();

  // EDT

  this->fillVoxels(atomTypes[SAS]);

  this->buildBoundary();

  m_Cube = new Cube();
  // Initialize cube

  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      for (k = 0; k < data->pHeight; k++) {
        m_Cube->setValue(i, j, k, volumePixels[i][j][k].distance);
      }
      delete[] volumePixels[i][j];
    }
    delete[] volumePixels[i];
  }
  delete[] volumePixels;
  // Copy VolumePixels array into Cube and delete volumePixels
  return true;
}//end SolventAccessibleSurface

bool EDTSurfaces::SolventExcludedSurface(){

  int i, j, k;

  this->setMolecule(mol);
  //Set molecule

  this->initPara(atomTypes[SES], bTypes[SES]);
  // Initialize everything

  this->fastDistanceMap();
  // EDT

  this->fillVoxels(atomTypes[SES]);

  this->buildBoundary();

  this->fillVoxelsWaals(atomTypes[SES]);
  this->boundingAtom(false);

  m_Cube = new Cube();
  // Initialize cube

  for (i = 0; i < data->pLength; i++) {
    for (j = 0; j < data->pWidth; j++) {
      for (k = 0; k < data->pHeight; k++) {
        m_Cube->setValue(i, j, k, volumePixels[i][j][k].distance);
      }
      delete[] volumePixels[i][j];
    }
    delete[] volumePixels[i];
  }
  delete[] volumePixels;
  // Copy VolumePixels
  return true;
}//end SolventExcludedSurface

} // End namespace Core

} // End namespace Avogadro
