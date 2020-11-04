#pragma once

#include "precision.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <OpenMesh/Core/Mesh/SmartHandles.hh>

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::VectorT<real, 3> Point; // use double-values points
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TriMesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> PolyMesh;