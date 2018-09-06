/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2008 Marcus D. Hanwell
  Copyright 2012 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#include "mesh.h"

#include "mutex.h"

using std::vector;

namespace Avogadro {
namespace Core {

Mesh::Mesh() : m_stable(true), m_other(0), m_cube(0), m_lock(new Mutex)
{
  m_vertices.reserve(100);
  m_normals.reserve(100);
  m_colors.reserve(1);
}

Mesh::Mesh(const Mesh& other)
  : m_vertices(other.m_vertices), m_normals(other.m_normals),
    m_colors(other.m_colors), m_name(other.m_name), m_stable(true),
    m_isoValue(other.m_isoValue), m_other(other.m_other), m_cube(other.m_cube),
    m_lock(new Mutex)
{
}

Mesh::~Mesh()
{
  delete m_lock;
  m_lock = 0;
}

bool Mesh::reserve(unsigned int size, bool useColors)
{
  m_vertices.reserve(size);
  m_normals.reserve(size);
  if (useColors)
    m_colors.reserve(size);
  return true;
}

void Mesh::setStable(bool isStable)
{
  m_stable = isStable;
}

bool Mesh::stable()
{
  return m_stable;
}

const Core::Array<Vector3f>& Mesh::vertices() const
{
  return m_vertices;
}

const Vector3f* Mesh::vertex(int n) const
{
  return &(m_vertices[n * 3]);
}

bool Mesh::setVertices(const Core::Array<Vector3f>& values)
{
  m_vertices.clear();
  m_vertices = values;
  return true;
}

bool Mesh::addVertices(const Core::Array<Vector3f>& values)
{
  if (m_vertices.capacity() < m_vertices.size() + values.size())
    m_vertices.reserve(m_vertices.capacity() * 2);
  if (values.size() % 3 == 0) {
    for (unsigned int i = 0; i < values.size(); ++i)
      m_vertices.push_back(values.at(i));
      //reset volume and surfaceArea so they can be recalculated
      m_volume = 0;
      m_surfaceArea = 0;
    return true;
  } else {
    return false;
  }
}

const Core::Array<Vector3f>& Mesh::normals() const
{
  return m_normals;
}

const Vector3f* Mesh::normal(int n) const
{
  return &(m_normals[n * 3]);
}

bool Mesh::setNormals(const Core::Array<Vector3f>& values)
{
  m_normals.clear();
  m_normals = values;
  return true;
}

bool Mesh::addNormals(const Core::Array<Vector3f>& values)
{
  if (m_normals.capacity() < m_normals.size() + values.size())
    m_normals.reserve(m_normals.capacity() * 2);
  if (values.size() % 3 == 0) {
    for (unsigned int i = 0; i < values.size(); ++i)
      m_normals.push_back(values.at(i));
    return true;
  } else {
    return false;
  }
}

const Core::Array<Color3f>& Mesh::colors() const
{
  return m_colors;
}

const Color3f* Mesh::color(int n) const
{
  // If there is only one color return that, otherwise colored by vertex.
  if (m_colors.size() == 1)
    return &(m_colors[0]);
  else
    return &(m_colors[n * 3]);
}

bool Mesh::setColors(const Core::Array<Color3f>& values)
{
  m_colors.clear();
  m_colors = values;
  return true;
}

bool Mesh::addColors(const Core::Array<Color3f>& values)
{
  if (m_colors.capacity() < m_colors.size() + values.size())
    m_colors.reserve(m_colors.capacity() * 2);
  if (values.size() % 3 == 0) {
    for (unsigned int i = 0; i < values.size(); ++i)
      m_colors.push_back(values.at(i));
    return true;
  } else {
    return false;
  }
}

bool Mesh::valid() const
{
  if (m_vertices.size() == m_normals.size()) {
    if (m_colors.size() == 1 || m_colors.size() == m_vertices.size())
      return true;
    else
      return false;
  } else {
    return false;
  }
}

bool Mesh::clear()
{
  m_vertices.clear();
  m_normals.clear();
  m_colors.clear();
  m_surfaceArea = 0;
  m_volume = 0;
  return true;
}

double Mesh::surfaceArea(){
  if(m_surfaceArea != 0){
    //If we've already calculated this, don't do it again
    return m_surfaceArea;
  }

  else{
    for(int i = 0; i < m_vertices.size() / 3; i++){
      //For each face, find the three vertices
      Vector3f vertexOne = m_vertices[i * 3];
      Vector3f vertexTwo = m_vertices[i * 3 + 1];
      Vector3f vertexThree = m_vertices[i * 3 + 2];

      //Make vectors connecting the vertices
      Vector3f oneTwo = vertexTwo - vertexOne;
      Vector3f twoThree = vertexThree - vertexTwo;
      //The cross product will be the normal vector to the plane containing
      //the vertices
      Vector3f normal = oneTwo.cross(twoThree);
      //And its magnitude is equal to the area of the parallelogram formed by
      //Them, so the triangle is half of that
      m_surfaceArea += (normal.norm() / 2);
    }
  }
  return m_surfaceArea;
}

double Mesh::volume(){

  if(m_volume != 0){
    //If we've already calculated this, don't do it again
    return m_volume;
  }

  else{
    for(int i = 0; i < m_vertices.size() / 3; i++){
      //For each face, find the three vertices
      Vector3f vertexOne = m_vertices[i * 3];
      Vector3f vertexTwo = m_vertices[i * 3 + 1];
      Vector3f vertexThree = m_vertices[i * 3 + 2];
      //These three places plus the origin make a tetrahedron
      //Find the volume of this tetrahedron

      double unsignedVolume = tetrahedronVolume(vertexOne, vertexTwo, vertexThree);
      double signedVolume;
      //Make vectors connecting the vertices
      Vector3f oneTwo = vertexTwo - vertexOne;
      Vector3f twoThree = vertexThree - vertexTwo;
      //The cross product will be the normal vector to the plane containing
      //the vertices
      Vector3f normal = oneTwo.cross(twoThree);
      //If the normal vector points in the direction of the origin
      //Then sign this volume negatively
      if(normal.dot(vertexOne) < 0){
        signedVolume = - unsignedVolume;
      }
      else{
        signedVolume = unsignedVolume;
      }
      //Increment total volume by signed volume
      m_volume += signedVolume;
    }
  }
  return m_volume;

}

//calculates the volume of a tetrahedron, assuming one point is the origin
double Mesh::tetrahedronVolume(Vector3f a, Vector3f b, Vector3f c){
  return( - (a(2) * b(1) * c(0) + a(1) * b(2) * c(0) + a(2) * b(0) * c(1) -
    a(0) * b(2) * c(1) - a(1) * b(0) * c(2) + a(0) * b(1) * c(2)) / 6);
}


Mesh& Mesh::operator=(const Mesh& other)
{
  m_vertices = other.m_vertices;
  m_normals = other.m_vertices;
  m_colors = other.m_colors;
  m_name = other.m_name;
  m_isoValue = other.m_isoValue;

  return *this;
}

} // End namespace QtGui
} // End namespace Avogadro
