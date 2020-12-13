#pragma once

#include "common.h"

struct NeighborTraits
{
	float angDist;
	float geoDist;
	float dist;
};

class GeoGraph
{
private:
	std::vector<std::unordered_map<int, NeighborTraits>> dualVerts;
	std::vector<std::vector<float>> minDistMatrix;
	size_t numVerts;

public:
	GeoGraph(const MeshType &mesh, float delta);
	std::vector<std::vector<float>>  RunShortestPath();
};

