// This library with source codes is available under MIT license.

// Copyright (c) 2022  Daniel Gorniak, Robert Nowak, Warsaw University of Technology

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef LFdBG_H
#define LFdBG_H

#include <atomic>
#include <boost/pool/object_pool.hpp>
#include <cmath>
#include <deque>
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

class LFdBG
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
    LFdBG(unsigned int              k,
          std::vector<std::string> &reads,
          unsigned int              read_length,
          unsigned int              genome_size,
          int                       threads_for_assembly,
          unsigned long long        table_size,
          bool                      normalize = true);
    ~LFdBG();
    std::deque<std::string>  getContigs();
    std::vector<std::string> getGraph(bool withEdges);
};
#endif // LFdBG_H
