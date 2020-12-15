#include "MinCut.h"

void MinCut::GetRefinementArea()
{
    auto& partitions = this->fuzzyLabels;
    auto& unsure = this->allRefines;
    auto N = this->totalVertices, C = this->totalCategories;
    this->labels = std::vector<int>(N, -1);
    for (int k = 0; k < N; ++k)
    {
        int cls1 = -1, cls2 = -1;
        float p1 = LBND_FLT, p2 = LBND_FLT;
        // First largest and second largest
        // This can be optimized
        for (int c = 0; c < C; ++c)
        {
            if (p1 < this->prob[k][c])
            {
                p1 = this->prob[k][c];
                cls1 = c;
            }
        }
        for (int c = 0; c < C; ++c)
        {
            if (p2 < this->prob[k][c] && p1 > this->prob[k][c])
            {
                p2 = this->prob[k][c];
                cls2 = c;
            }
        }
        // If the top 2 probs are too close ...
        if (p1 - p2 < this->epsilon)
        {
            this->labels[k] = -1;
            UnsurePoint point{ k, cls1, cls2, p1, p2 };
            auto pp = std::find(partitions.begin(), partitions.end(), std::set<int>{ cls1, cls2 });
            if (pp != partitions.end())
            {
                int resultIdx = pp - partitions.begin();
                unsure[resultIdx].emplace_back(point);
                
            }
            else
            {
                partitions.emplace_back(std::set<int>{ cls1, cls2 });
                unsure.emplace_back(std::vector<UnsurePoint>{point});
            }
        }
        else
        {
            this->labels[k] = cls1;
        }
    }
}

void MinCut::FFFlowGraph(const std::vector<UnsurePoint>& pts, const int cls1, const int cls2)
{
    int N0 = pts.size();
    std::unordered_set<int> unsuredPool;
    for (const auto& p : pts)
    {
        unsuredPool.insert(p.index);
    }

    /*std::unordered_map<int, int> connects;
    int cnt = 0;
    std::unordered_set<int> unsuredPool;

    for (const auto& p : pts)
    {
        unsuredPool.insert(p.index);
        connects[p.index] = -1;
    }*/

    /*for(int i = 0; i < N0; ++i)
    {
        int j = pts[i].index;

        if (connects[j] >= 0)
        {
            continue;
        }
        std::stack<int> s;
        s.push(j);
        
        while (!s.empty())
        {
            int cur = s.top();
            s.pop();
            connects[cur] = cnt;

            for (const auto& faceNeighbor : (this->mesh.faces_sbegin() + cur)->faces())
            {
                int neighbor = faceNeighbor.idx();
                if (unsuredPool.contains(neighbor) && connects.at(neighbor) == -1)
                {
                    s.push(neighbor);
                }
            }
        }
        ++cnt;
    }*/
    
    std::unordered_set<int> vc, vca, vcb;

    for (int i = 0; i < N0; ++i)
    {
        int j = pts[i].index;
        vc.insert(j);
        
        for (const auto& neighborJ : (this->mesh.faces_begin() + j)->faces())
        {
            int idxNeighbor = neighborJ.idx();
            if (!unsuredPool.contains(idxNeighbor))
            {
                if (this->labels[idxNeighbor] == cls1)
                {
                    vca.insert(idxNeighbor);
                }
                if (this->labels[idxNeighbor] == cls2)
                {
                    vcb.insert(idxNeighbor);
                }
            }
        }
    }
    this->vcax = vca;
    this->vcbx = vcb;
    this->vcx = vc;

    int V = vc.size() + vca.size() + vcb.size();
    std::vector<int> verts(V);
    int vcnt = 0;
    for (auto v : vc)
    {
        verts[vcnt++] = v;
    }
    for (auto v : vca)
    {
        verts[vcnt++] = v;
    }
    for (auto v : vcb)
    {
        verts[vcnt++] = v;
    }
    // Construct Flow Graph
    std::vector<std::vector<float>> flows(V + 2, std::vector<float>(V + 2, LBND_FLT));
    float adSum = 0.0, adDen = 0.0;
    for (int i = 0; i < V; ++i)
    {
        const auto& faceI = *(this->mesh.faces_begin() + verts[i]);
        const auto& n1 = mesh.calc_face_normal(faceI);
        const auto& p1 = mesh.calc_face_centroid(faceI);
        const auto& neighborI = faceI.faces().to_set();
        for (int j = i + 1; j < V; ++j)
        {
            const auto& faceJ = *(this->mesh.faces_begin() + verts[j]);
            if (neighborI.contains(faceJ))
            {
                const auto& n2 = mesh.calc_face_normal(faceJ);
                const auto& p2 = mesh.calc_face_centroid(faceJ);
                float ad = GetAngleDistance(n1, n2, p1, p2, this->eta);
                adSum += ad;
                adDen += 1;
                flows[i][j] = flows[j][i] = ad;
            }
        }
    }
    adSum /= adDen;
    for (int i = 0; i < V; ++i)
    {
        const auto& faceI = *(this->mesh.faces_begin() + verts[i]);
        const auto& neighborI = faceI.faces().to_set();
        for (int j = i + 1; j < V; ++j)
        {
            const auto& faceJ = *(this->mesh.faces_begin() + verts[j]);
            if (neighborI.contains(faceJ))
            {
                flows[i][j] = 1.0 / (1.0 + flows[i][j] / adSum);
                flows[j][i] = flows[i][j];
            }
        }
    }
    for (int i = 0; i < V; ++i)
    {
        if (vca.contains(verts[i]))
        {
            flows[V][i] = UBND_FLT;
            flows[i][V] = UBND_FLT;
        }
        if (vcb.contains(verts[i]))
        {
            flows[V + 1][i] = UBND_FLT;
            flows[i][V + 1] = UBND_FLT;
        }
    }
    std::vector<int> pre(V + 2, -1);
    std::vector<float> flow(V + 2, LBND_FLT);
    auto _BFS = [V, &pre, &flow, &flows](int S, int T) -> float
    {
        std::queue<int> bfsQ;
        pre.assign(V + 2, -1);
        pre[S] = -1;
        flow[S] = UBND_FLT;
        bfsQ.push(S);
        while (!bfsQ.empty())
        {
            int index = bfsQ.front();
            bfsQ.pop();
            if (index == T)
            {
                break;
            }
            for (int i = 0; i < V + 2; ++i)
            {
                if (i != S && pre[i] == -1 && flows[index][i] > LBND_FLT)
                {
                    pre[i] = index;
                    flow[i] = flows[index][i] < flow[index] ? flows[index][i] : flow[index];
                    bfsQ.push(i);
                }
            }
        }
        if (pre[T] == -1)
            return -1.0;
        else
            return flow[T];
    };
    // BFS find max flow
    float flowIncrement = 0;
    while ((flowIncrement = _BFS(V, V + 1)) > LBND_FLT)
    {
        int k = V + 1;
        while (k != V)
        {
            int last = pre[k];
            if (flows[last][k] < UBND_FLT)
            {
                flows[last][k] -= flowIncrement;
            }
            if (flows[k][last] < UBND_FLT)
            {
                flows[k][last] += flowIncrement;
            }
            k = last;
        }
    }
    // Last BFS indicates min cut, now pre contains the points in S-cut with non -1
    for (int j = 0; j < V; ++j)
    {
        int clsJ = pre[j] != -1 ? cls1 : cls2;
        if (vc.contains(verts[j]))
        {
            this->labels[verts[j]] = clsJ;
        }

    }
}

MinCut::MinCut(const std::vector<std::vector<float>>& p, float eps, const MeshType& m, float eta)
    : mesh(m), prob(p)
{
    this->epsilon = eps;
    this->eta = eta;
    this->totalVertices = prob.size();
    this->totalCategories = prob[0].size();

    this->GetRefinementArea();
}

void MinCut::Execute()
{
    int kkPairs = fuzzyLabels.size();
    for (int k = 0; k < kkPairs; ++k)
    {
        auto fuzzyLabel = fuzzyLabels[k].begin();
        auto cls1 = *fuzzyLabel, cls2 = *(++fuzzyLabel);
        const auto& kkRefines = this->allRefines[k];
        this->FFFlowGraph(kkRefines, cls1, cls2);
    }
}

std::vector<int> MinCut::GetLabels()
{
    return this->labels;
}
