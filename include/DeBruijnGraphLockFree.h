#ifndef GRAPH_DEBRUIJNGRAPHLOCKFREE_H
#define GRAPH_DEBRUIJNGRAPHLOCKFREE_H

#include <atomic>
#include <boost/pool/object_pool.hpp>
#include <cmath>
#include <deque>
#include <iostream>
#include <queue>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

static const int A = 0;
static const int C = 1;
static const int T = 2;
static const int G = 3;
static const int N = 4;

class DeBruijnGraphLockFree
{
  private:
    struct Node
    {
        std::string_view                         k_mer;
        std::array<std::atomic<std::uint8_t>, N> counters{};
        std::array<std::atomic<Node *>, N>       pointers{};
        Node                                    *next;
        std::atomic<bool>                        first;
        Node(std::string_view mer, Node *n, bool f) : k_mer(mer), next(n), first(f) {}
        int onlyOnePath() const;
    };
    const unsigned long long table_size;
    const unsigned int       k_mer;
    const double             redundancy;
    const int                threads_for_assembly;

    std::deque<std::atomic<Node *>> hash_table;

    std::deque<std::string>                contigs;
    std::deque<std::deque<std::string>>    contigs_pools;
    std::deque<boost::object_pool<Node> *> object_pools;
    bool                                   already_searched = false;

    void makeGraph(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end, int index);
    void normalizeGraph(unsigned int start, unsigned int end);
    void findContig(Node *node, const unsigned int &index);
    void makePath(Node *node, const unsigned int &index);
    void findContigs(unsigned int start, unsigned int end, long unsigned int i);
    void delPool(long unsigned int i);

  public:
    DeBruijnGraphLockFree(unsigned int              k,
                          std::vector<std::string> &reads,
                          unsigned int              read_length,
                          unsigned int              genome_size,
                          int                       threads_for_assembly,
                          unsigned long long        table_size,
                          bool                      normalize = true);
    ~DeBruijnGraphLockFree();
    std::deque<std::string>  getContigs();
    std::vector<std::string> getGraph(bool withEdges);
};
#endif // GRAPH_DEBRUIJNGRAPHLOCKFREE_H
