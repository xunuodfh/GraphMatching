
#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <chrono>
#include <omp.h>
#include <numeric>

using namespace std;

int nodeCount = 0;
//size_t totalCount = 0;
vector<size_t> totalCount(64,0);

class Timer {
public:
    Timer() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    ~Timer() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        std::cout << "Elapsed time: " << (float)duration/1000/1000/60 << " minutes" << std::endl;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

void updateNodeCount(const string& filename){
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "无法打开文件: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    int u = 0, v = 0;
    while (file >> u >> v) {
        int bigger = max(u,v);
        if(bigger>nodeCount){
            nodeCount = bigger;
        }
    }
    nodeCount++;
}
// CSR数据结构
struct CSR {
    //uint32_t
    //&colIndices[rowPointers[nodeID]]
    vector<int> colIndices;    // 每条边的目标节点
    //int64_t
    vector<int> rowPointers;   // 每个节点的邻居列表在colIndices中的起始位置，下一个为终止位置
};

// 从TXT文件中读取边列表并转换为邻接表
vector<vector<int>> readEdgesFromFile(const string& filename) {
    //第二个vector要排序，还要去重。
//    unordered_map<int, vector<int>> adjList(548411);
    vector<vector<int>>adjList(nodeCount);
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "无法打开文件: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    int u = 0, v = 0;
    while (file >> u >> v) {
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }
    for(auto& vec:adjList){
        sort(vec.begin(),vec.end());
        vec.erase(unique(vec.begin(),vec.end()),vec.end());
    }
    file.close();
//    cout<<adjList.size()<<endl;

    return adjList;

}

// 将邻接表转换为CSR格式
//find部分也要改。
CSR convertToCSR(const vector<vector<int>>& adjList) {
    CSR csr;
    csr.rowPointers.push_back(0); // 第一个元素永远是0

    for (int i = 0; i < nodeCount; ++i) {
//        if (adjList.find(i) != adjList.end()) {
        if(!adjList[i].empty()){
            for (int neighbor : adjList.at(i)) {
                csr.colIndices.push_back(neighbor);
            }
        }
        csr.rowPointers.push_back(csr.colIndices.size());
    }

    return csr;
}

// 输出CSR格式内容
void printCSR(const CSR& csr) {
    cout << "Column Indices: ";
    for (int c : csr.colIndices) {
        cout << c << " ";
    }
    cout << endl;

    cout << "Row Pointers: ";
    for (int r : csr.rowPointers) {
        cout << r << " ";
    }
    cout << endl;
}
// 获取CSR中某个节点的邻居节点集合
//这个要改成返回指针 uint32_t 首地址，和size
//unordered_set<int> getNeighbors(const CSR& csr, int vertex) {
//    unordered_set<int> neighbors;
//    int start = csr.rowPointers[vertex];
//    int end = csr.rowPointers[vertex + 1];
//    for (int i = start; i < end; ++i) {
//        neighbors.insert(csr.colIndices[i]);
//    }
//    return neighbors;
//}
pair<const int*, int> getNeighbors(const CSR& csr, int vertex) {
    // Find the start and end positions of the neighbors for the given vertex
    int start = csr.rowPointers[vertex];
    int end = csr.rowPointers[vertex + 1];

    // Return the pointer to the first neighbor and the number of neighbors
    return make_pair(&csr.colIndices[start], end - start);
}

// 计算两个集合的交集
//改成stl库，指针版本
//unordered_set<int> intersect(const unordered_set<int>& set1, const unordered_set<int>& set2) {
//    unordered_set<int> result;
//    for (int elem : set1) {
//        if (set2.find(elem) != set2.end()) {
//            result.insert(elem);
//        }
//    }
//    return result;
//}
vector<int> intersect(const int* set1_begin, const int* set1_end, const int* set2_begin, const int* set2_end) {
    vector<int> result;
    // 预分配结果向量的大小为两个集合大小中较小的那个，避免动态扩展
//    result.reserve(min(set1_end - set1_begin, set2_end - set2_begin));

    // 计算两个集合的交集，结果存储在 result 中
    set_intersection(set1_begin, set1_end, set2_begin, set2_end, back_inserter(result));
//    sort(result.begin(),result.end());
    return result;
}

// 输出同构子图，去除重复子图
void output(int v0, int v1, int v2, int v3, set<tuple<int, int, int, int>>& uniqueSubgraphs) {
    // 将子图节点排序以避免重复
    vector<int> vertices = {v0, v1, v2, v3};
    sort(vertices.begin(), vertices.end());
    tuple<int, int, int, int> subgraph = make_tuple(vertices[0], vertices[1], vertices[2], vertices[3]);

    // 如果该子图已经输出过，跳过
    if (uniqueSubgraphs.find(subgraph) == uniqueSubgraphs.end()) {
        uniqueSubgraphs.insert(subgraph);
        cout << "Subgraph found: (" << vertices[0] << ", " << vertices[1] << ", " << vertices[2] << ", " << vertices[3] << ")" << endl;
    }
}

// 查找所有与查询图同构的子图
//intersec数组的声明放到循环外面.
void findIsomorphicSubgraphs(const CSR& csr) {
//    int nodeCount = csr.rowPointers.size() - 1;  // 节点数
//    vector<unordered_set<int>> C(4); // 存储候选节点集
    set<tuple<int, int, int, int>> uniqueSubgraphs;  // 存储已经输出的子图


    vector<vector<int>> C2(64,vector<int>(40960,0));
    vector<vector<int>> C3(64,vector<int>(40960,0));
    // C[0] = 所有节点
    #pragma omp parallel for schedule(dynamic, 8) num_threads(64)
    for (int v0 = 0; v0 < nodeCount; ++v0) {  //v0
        int tid0 = omp_get_thread_num();
        // C[1] = N(v0)
//        C[1] = getNeighbors(csr, v0);
        auto [neighbors_ptr_c0, size_c0] = getNeighbors(csr, v0);
        for (int i = 0; i<size_c0; i++) {  //v1
            // C[2] = N(v0) ∩ N(v1)
            auto [neighbors_ptr_c1, size_c1] = getNeighbors(csr, neighbors_ptr_c0[i]);  //把地址变量在循环外存储，然后在getNeighbors里面修改地址变量的内容
            C2[tid0] = intersect(neighbors_ptr_c0,neighbors_ptr_c0+size_c0*sizeof(int),neighbors_ptr_c1,neighbors_ptr_c1+size_c1*sizeof(int));
            for (int v2 : C2[tid0]) {   //v2
                // C[3] = N(v1) ∩ N(v2) − {v0}
//                auto [neighbors_ptr_c2, size_c2] = getNeighbors(csr, neighbors_ptr_c1[i]);
                int v1 = neighbors_ptr_c0[i];
                if(v2>=v1)
                    continue;

                auto [neighbors_ptr_c2, size_c2] = getNeighbors(csr, v2);
                C3[tid0] = intersect(neighbors_ptr_c1,neighbors_ptr_c1+size_c1*sizeof(int), neighbors_ptr_c2,neighbors_ptr_c2+size_c2*sizeof(int));
//                C3.erase(v0);
//                C3.erase(remove(C3.begin(), C3.end(), v0), C3.end());
                for (int v3 : C3[tid0]) {
                    if(v3!=v0 && v3<v0)
                        totalCount[tid0]++;


//                    output(v0, neighbors_ptr_c0[i], v2, v3, uniqueSubgraphs);

                }
            }
        }
    }
}

//void findIsomorphicSubgraphs2(const CSR& csr) {
//    int nodeCount = csr.rowPointers.size() - 1;  // 节点数
//    vector<unordered_set<int>> C(4); // 存储候选节点集
//    set<tuple<int, int, int, int>> uniqueSubgraphs;  // 存储已经输出的子图
//
//    // C[0] = 所有节点
//    for (int v0 = 0; v0 < nodeCount; ++v0) {
//        // C[1] = N(v0)
//        C[1] = getNeighbors(csr, v0);
//        for (int v1 : C[1]) {
//            // C[2] = N(v0) ∩ N(v1)
//            C[2] = intersect(getNeighbors(csr, v0), getNeighbors(csr, v1));
//            for (int v2 : C[2]) {
//                // C[3] = N(v1) ∩ N(v2) − {v0}
//                C[3] = intersect(getNeighbors(csr, v1), getNeighbors(csr, v2));
//                C[3].erase(v0);
//                for (int v3 : C[3]) {
//                    output(v0, v1, v2, v3, uniqueSubgraphs);
//                }
//            }
//        }
//    }
//}
int main() {
    Timer t;
    string filename = "C:\\CLionProject\\GraphMatching\\cmake-build-debug\\amazon.txt";
    updateNodeCount(filename);
    // 从文件中读取边列表并转换为邻接表
    vector<vector<int>> adjList = readEdgesFromFile(filename);

    // 假设节点编号是连续的，从0开始
//    int nodeCount = adjList.size();

    // 转换为CSR格式
    CSR csr = convertToCSR(adjList);

    // 输出CSR格式
//    printCSR(csr);

//    CSR csr = {
//            // colIndices数组，存储所有边的目标节点
//            {1, 2, 0, 2, 3, 0, 1, 3, 1, 2},
//            // rowPointers数组，存储每个节点的邻居列表起始位置
//            {0, 2, 5, 8, 10}
//    };

    findIsomorphicSubgraphs(csr);
    size_t total = accumulate(totalCount.begin(),totalCount.end(),0);
    cout << total<<endl;

    return 0;
}
