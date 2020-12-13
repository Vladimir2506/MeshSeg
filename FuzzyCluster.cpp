#include "FuzzyCluster.h"

void FuzzyCluster::initCentroids()
{
	const auto C = this->categories;
	const auto N = this->distances.size();

	/*float minSumRow = UBND_FLT;
	int primeCentroid = -1;
	for (int i = 0; i < N; ++i)
	{
		float sumRow = 0.0;
		for (int j = 0; j < N; ++j)
		{
			if (i != j)
			{
				sumRow += this->distances[i][j];
			}
		}
		if (sumRow < minSumRow)
		{
			minSumRow = sumRow;
			primeCentroid = i;
		}
	}*/
	this->gVal = -1.0;
	float maxPairDist = LBND_FLT;
	int p1 = -1, p2 = -1;
	for (int i = 0; i < N; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			if (distances[i][j] > maxPairDist)
			{
				maxPairDist = distances[i][j];
				p1 = i;
				p2 = j;
			}
		}
	}

	this->centroids.push_back(p1);
	this->centroids.push_back(p2);

	std::vector<float> temp;

	for (int c1 = 2; c1 < C; ++c1)
	{
		float maxDist = LBND_FLT;
		int maxIdx = -1;

		for (int i = 0; i < N; ++i)
		{
			float centerSetDist = UBND_FLT;
			for (auto source : this->centroids)
			{
				centerSetDist = this->distances[i][source] < centerSetDist ? this->distances[i][source] : centerSetDist;
			}
			if (centerSetDist > maxDist)
			{
				maxDist = centerSetDist;
				maxIdx = i;
			}
		}

		this->centroids.push_back(maxIdx);
		this->gVal = maxDist;
	}

	if (this->gVal < 0)
	{
		this->gVal = maxPairDist;
	}
}

FuzzyCluster::FuzzyCluster(int M, int C, const std::vector<std::vector<float>>& dist)
	: distances(dist)
{
	this->maxIters = M;
	this->categories = C;
	this->samples = dist.size();

	const auto N = this->samples;

	this->probs = std::vector<std::vector<float>>(N, std::vector<float>(C, 0));

	this->initCentroids();

	this->assignments = std::vector<std::unordered_set<int>>(C);
	
	for (int k = 0; k < N; ++k)
	{
		float sumK = 0.0;
		float maxProb = LBND_FLT;
		int maxIdx = 0;
		for (int c = 0; c < C; ++c)
		{
			auto& cntr = this->centroids[c];
			this->probs[k][c] = 1.0 / this->distances[cntr][k];
			sumK += 1.0 / this->distances[cntr][k];
			if (this->probs[k][c] > maxProb)
			{
				maxProb = this->probs[k][c];
				maxIdx = c;
			}
		}
		for (auto& p : this->probs[k])
		{
			p /= sumK;
		}
		this->assignments[maxIdx].insert(k);
	}
}

void FuzzyCluster::RunCluster()
{
	const auto M = this->maxIters;
	const auto C = this->categories;
	const auto N = this->samples;

	for (int iter = 0; iter < M; ++iter)
	{
		for (int k = 0; k < N; ++k)
		{
			float sumK = 0.0;
			for (int c = 0; c < C; ++c)
			{
				//auto cntr = this->centroids[c];
				float meanDist = 0.0;
				for (const auto& ass : this->assignments[c])
				{
					meanDist += this->distances[k][ass];
				}
				meanDist /= this->assignments[c].size();

				/*this->probs[k][c] = 1.0 / this->distances[cntr][k];
				sumK += 1.0 / this->distances[cntr][k];*/
				this->probs[k][c] = 1.0 / meanDist;
				sumK += this->probs[k][c];
			}
			for (auto& p : this->probs[k])
			{
				p /= sumK;
			}
		}

		std::vector<int> nextCntrs(C, 0);
		std::vector<float> minSums(C, UBND_FLT);
		this->assignments.assign(C, std::unordered_set<int>());
		for (int i = 0; i < N; ++i)
		{
			std::vector<float> sums(C, 0.0);
			int maxC = -1;
			float maxP = LBND_FLT;
			for (int j = 0; j < N; ++j)
			{
				for (int c = 0; c < C; ++c)
				{
					sums[c] += this->probs[j][c] * this->distances[i][j];
				}
			}
			for (int c = 0; c < C; ++c)
			{
				if (sums[c] < minSums[c])
				{
					nextCntrs[c] = i;
					minSums[c] = sums[c];
				}
				
				if (this->probs[i][c] > maxP)
				{
					maxP = this->probs[i][c];
					maxC = c;
				}
			}
			this->assignments[maxC].insert(i);
		}

		if (nextCntrs == this->centroids)
		{
			//std::cout << "Converged." << std::endl;
			break;
		}

		this->centroids = nextCntrs;
	}

}

std::vector<std::vector<float>> FuzzyCluster::GetProbabilities()
{
	return this->probs;
}

float FuzzyCluster::GetGValue()
{
	return this->gVal;
}
