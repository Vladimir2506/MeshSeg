#pragma once

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <queue>

#include <iostream>
#include <string>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

using MeshType = OpenMesh::TriMesh_ArrayKernelT<>;
using NormalType = OpenMesh::DefaultTraits::Normal;
using PointType = OpenMesh::DefaultTraits::Point;
using ColorType = OpenMesh::DefaultTraits::Color;

constexpr float UBND_FLT = 1.e7f;
constexpr float LBND_FLT = 1.0e-7f;

class PaletteHSV
{
private:
	void HSVtoRGB(int& R, int& G, int& B, float H, float S = 100.0, float V = 100.0);
public:
	std::vector<ColorType> colors;
	PaletteHSV(uint cls);
};

bool IsConvex(const NormalType& n1, const PointType& p1, const PointType& p2);
float GetGeodesicDistance(const NormalType& n1, const NormalType& n2, const PointType& p1, const PointType& p2);
float GetAngleDistance(const NormalType& n1, const NormalType& n2, const PointType& p1, const PointType& p2);