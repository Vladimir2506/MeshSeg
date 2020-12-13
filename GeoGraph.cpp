#include "GeoGraph.h"
#include "APSPCuda.cuh"

GeoGraph::GeoGraph(const MeshType &mesh, float delta)
{
	// Initialization of variables.
	this->numVerts = mesh.n_faces();
	this->dualVerts = std::vector<std::unordered_map<int, NeighborTraits>>(this->numVerts);
	this->minDistMatrix = std::vector<std::vector<float>>(this->numVerts, std::vector<float>(this->numVerts, UBND_FLT));

	// Construct the basic graph.
	float totalGeod = 0, totalAngDist = 0;
	int denCounter = 0;

	for (const auto &faceI : mesh.faces())
	{
		int indexI = faceI.idx();
		auto centroidI = mesh.calc_face_centroid(faceI);
		auto normalVecI = mesh.calc_normal(faceI);

		const auto &edgesISet = faceI.edges().to_set();

		for (const auto& faceJ : faceI.faces())
		{
			int indexJ = faceJ.idx();
			const auto &edgesJ = faceJ.edges();
			OpenMesh::SmartEdgeHandle e;

			for (const auto& ej : edgesJ)
			{
				if (edgesISet.contains(ej))
				{
					e = ej;
					break;
				}
			}

			auto centroidJ = mesh.calc_face_centroid(faceJ);
			auto normalVecJ = mesh.calc_normal(faceJ);

			float geodIJ = GetGeodesicDistance(centroidI, centroidJ, mesh.point(e.v0()), mesh.point(e.v1()));
			float angDistIJ = GetAngleDistance(normalVecI, normalVecJ, centroidI, centroidJ);

			this->dualVerts[indexI][indexJ] = NeighborTraits {angDistIJ, geodIJ, UBND_FLT};

			totalGeod += geodIJ;
			totalAngDist += angDistIJ;
			denCounter += 1;
		}
	}

	float avgGeod = totalGeod / denCounter, avgAngd = totalAngDist / denCounter; 

	for (int i = 0; i < this->numVerts; ++i)
	{
		for (auto &kv : this->dualVerts[i])
		{
			int j = kv.first;
			kv.second.dist = kv.second.geoDist * delta / avgGeod + kv.second.angDist * (1.0f - delta) / avgAngd;
			this->minDistMatrix[i][j] = kv.second.dist;
			this->minDistMatrix[j][i] = kv.second.dist;
		}
	}

}


std::vector<std::vector<float>> GeoGraph::RunShortestPath()
{
	// Wrapper for the CUDA APSP.
	std::unique_ptr<graphAPSPTopology> proxy = std::unique_ptr<graphAPSPTopology>(new graphAPSPTopology(this->numVerts));

	for (size_t i = 0; i < this->numVerts; ++i)
	{
		for (size_t j = i; j < this->numVerts; ++j)
		{
			proxy->graph[i * this->numVerts + j] = this->minDistMatrix[i][j];
			proxy->graph[j * this->numVerts + i] = proxy->graph[i * this->numVerts + j];
		}
	}

	cudaBlockedFW(proxy);

	for (size_t i = 0; i < this->numVerts; ++i)
	{
		for (size_t j = i; j < this->numVerts; ++j)
		{
			this->minDistMatrix[i][j] = proxy->graph[i * this->numVerts + j];
			this->minDistMatrix[j][i] = this->minDistMatrix[i][j];
		}
	}
	return this->minDistMatrix;
}
