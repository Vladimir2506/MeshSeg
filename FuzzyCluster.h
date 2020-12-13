#pragma once

#include "common.h"

class FuzzyCluster
{
private:
	int maxIters;
	int samples;
	int categories;
	const std::vector<std::vector<float>>& distances;
	std::vector<std::vector<float>> probs;
	std::vector<int> centroids;
	std::vector<std::unordered_set<int>> assignments;
	float gVal;
	void initCentroids();

public:
	FuzzyCluster(int M, int C, const std::vector<std::vector<float>> &dist);
	void RunCluster();
	std::vector<std::vector<float>> GetProbabilities();
	float GetGValue();
};

