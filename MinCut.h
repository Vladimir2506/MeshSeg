#pragma once

#include "common.h"

struct UnsurePoint
{
	int index;
	int cls1Idx;
	int cls2Idx;
	float cls1Prob;
	float cls2Prob;
};

class MinCut
{
private:
	const MeshType& mesh;
	const std::vector<std::vector<float>>& prob;
	int totalVertices;
	int totalCategories;
	float epsilon;
	float eta;
	std::vector<std::vector<UnsurePoint>> allRefines;
	std::vector<std::set<int>> fuzzyLabels;
	std::vector<int> labels;

	void GetRefinementArea();
	void FFFlowGraph(const std::vector<UnsurePoint>& pts, int cls1, int cls2);
	
public:
	MinCut(const std::vector<std::vector<float>>& p, float eps, const MeshType& m, float eta);
	void Execute();
	std::vector<int> GetLabels();
	std::unordered_set<int> vcx, vcax, vcbx;
};
