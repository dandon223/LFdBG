#include "../include/LFdBG.h"

// This method converts 0 to 'A', 1 to 'C', 2 to 'T' and 3 to 'G'.
inline char getChar(int i) {
    switch (i) {
        case A:
            return 'A';
        case C:
            return 'C';
        case T:
            return 'T';
        default:
            return 'G';
    }
}

// This method converts 'A'to 0, 'C' to 1, 'T' to 2 and 'G' to 3.
inline int findIndex(char c) {
    return (c >> 1) & 0x3;
}

static const int NO_NEXT_NODES   = -2; // used in onlyOnePath
static const int MANY_NEXT_NODES = -1; // used in onlyOnePath

// This method checks if node has only one outgoing path.
int LFdBG::Node::onlyOnePath() const {
    if (counters[A] == 0 && counters[C] == 0 && counters[T] == 0 && counters[G] == 0)
        return NO_NEXT_NODES;
    if (counters[A] > 0 && counters[C] == 0 && counters[T] == 0 && counters[G] == 0)
        return A;
    if (counters[C] > 0 && counters[A] == 0 && counters[T] == 0 && counters[G] == 0)
        return C;
    if (counters[T] > 0 && counters[C] == 0 && counters[A] == 0 && counters[G] == 0)
        return T;
    if (counters[G] > 0 && counters[C] == 0 && counters[T] == 0 && counters[A] == 0)
        return G;
    return MANY_NEXT_NODES;
}

LFdBG::~LFdBG() {
    std::vector<std::thread> threads;
    for (long unsigned int i = 0; i < object_pools.size(); ++i)
        threads.emplace_back(&LFdBG::delPool, this, i);
    for (long unsigned int i = 0; i < object_pools.size(); ++i)
        threads[i].join();
}
void LFdBG::delPool(long unsigned int i) {
    delete object_pools[i];
}
// The constructor also makes graph.
LFdBG::LFdBG(unsigned int              k,
                                             std::vector<std::string> &reads,
                                             unsigned int              read_length,
                                             unsigned int              genome_size,
                                             int                       t_f_a,
                                             unsigned long long        t_s,
                                             bool                      normalize) :
    table_size(t_s),
    k_mer(k), redundancy(double(reads.size()) * double(read_length - k_mer + 1) / genome_size),
    threads_for_assembly(t_f_a) {

    for (unsigned int i = 0; i <= table_size; ++i)
        hash_table.emplace_back(nullptr);

    std::vector<std::thread> threads;

    for (int i = 0; i < threads_for_assembly; ++i)
        object_pools.push_back(new boost::object_pool<Node>{32, 1048576});

    for (int i = 0; i < threads_for_assembly; ++i)
        contigs_pools.emplace_back();

    unsigned int for_one_thread = reads.size() / threads_for_assembly;

    // Every thread gets its share of reads to work on.
    for (int i = 0; i < threads_for_assembly; ++i) {
        if (i == threads_for_assembly - 1) {
            threads.emplace_back(
                &LFdBG::makeGraph, this, reads.begin() + i * for_one_thread, reads.end(), i);
        } else {
            threads.emplace_back(&LFdBG::makeGraph,
                                 this,
                                 reads.begin() + i * for_one_thread,
                                 reads.begin() + (i + 1) * for_one_thread,
                                 i);
        }
    }

    for (int i = 0; i < threads_for_assembly; ++i)
        threads[i].join();
    if(!normalize)
        return;

    for_one_thread = table_size / threads_for_assembly;
    std::vector<std::thread> threads2;

    for (int i = 0; i < threads_for_assembly; ++i) {
        if (i == threads_for_assembly - 1) {
            threads2.emplace_back(&LFdBG::normalizeGraph, this, i * for_one_thread, table_size);
        } else {
            threads2.emplace_back(
                &LFdBG::normalizeGraph, this, i * for_one_thread, (i + 1) * for_one_thread);
        }
    }

    for (int i = 0; i < threads_for_assembly; ++i)
        threads2[i].join();
}

// This method is main method to make graph.
void LFdBG::makeGraph(std::vector<std::string>::iterator begin,
                                      std::vector<std::string>::iterator end,
                                      int                                index) {

    std::hash<std::string_view> hasher;
    Node                       *left_new_node  = nullptr;
    Node                       *right_new_node = nullptr;

    // We loop through every created mer.
    for (; begin != end; ++begin) {
        if (begin->length() < k_mer)
            continue;
        std::string_view read{*begin};
        for (long unsigned int i = 0; i <= begin->length() - k_mer; ++i) {
            auto mer = read.substr(i, k_mer);
            // We make left node from our mer.
            auto         left       = mer.substr(0, k_mer - 1);
            int          left_index = findIndex(mer[mer.size() - 1]);
            unsigned int left_hash  = hasher(left) % table_size;
            bool         left_bool  = false;
            // Node* left_new_node = nullptr;
            // We try as long as we do not add left node to the structure or (if it's already present) do not increment
            // proper counter.
            Node *find;
            for (find = hash_table[left_hash].load(); find != nullptr; find = find->next) {
                if (find->k_mer == left) {
                    ++find->counters[left_index];
                    left_new_node = find;
                    left_bool     = true;
                    break;
                }
            }
            if (!left_bool) {
                left_new_node = object_pools[index]->construct(left, hash_table[left_hash].load(), true);
                ++left_new_node->counters[left_index];
                do {

                    // We loop through already present nodes. If it's present we increment proper counter and break the
                    // loop.
                    for (find = left_new_node->next; find != nullptr; find = find->next) {
                        if (find->k_mer == left) {
                            ++find->counters[left_index];
                            object_pools[index]->free(left_new_node);
                            left_new_node = find;
                            left_bool     = true;
                            break;
                        }
                    }
                } while (!left_bool && !std::atomic_compare_exchange_strong(
                                           &hash_table[left_hash], &left_new_node->next, left_new_node));
            }
            // Node* right_new_node = nullptr;
            bool         right_bool = false;
            auto         right      = mer.substr(1, k_mer - 1);
            unsigned int right_hash = hasher(right) % table_size;
            for (find = hash_table[right_hash].load(); find != nullptr; find = find->next) {
                if (find->k_mer == right) {
                    right_new_node        = find;
                    right_new_node->first = false;
                    right_bool            = true;
                    break;
                }
            }
            if (!right_bool) {
                right_new_node = object_pools[index]->construct(right, hash_table[right_hash].load(), false);
                do {

                    // We loop through already present nodes. If it's present we increment proper counter and break the
                    // loop.
                    for (find = right_new_node->next; find != nullptr; find = find->next) {
                        if (find->k_mer == right) {
                            object_pools[index]->free(right_new_node);
                            right_new_node        = find;
                            right_new_node->first = false;
                            right_bool            = true;
                            break;
                        }
                    }
                } while (!right_bool && !std::atomic_compare_exchange_strong(
                                            &hash_table[right_hash], &right_new_node->next, right_new_node));
            }
            left_new_node->pointers[left_index].store(right_new_node, std::memory_order_relaxed);
        }
    }
}
// This method makes contigs from graph that was created.
void LFdBG::findContigs(unsigned int start, unsigned int end, long unsigned int index) {
    for (unsigned int i = start; i != end; i++) {
        if (hash_table[i] == nullptr)
            continue;
        Node *node = hash_table[i];
        do {
            if (node->first || node->onlyOnePath() == MANY_NEXT_NODES)
                findContig(node, index);
            node = node->next;
        } while (node != nullptr);
    }
}
// This method creates contig starting from node.
void LFdBG::findContig(Node *node, const unsigned int &index) {
    if (node->onlyOnePath() == MANY_NEXT_NODES) {
        if (node->counters[A] > 0) {
            makePath(node->pointers[A], index);
            node->counters[A] = 0;
        }
        if (node->counters[C] > 0) {
            makePath(node->pointers[C], index);
            node->counters[C] = 0;
        }
        if (node->counters[T] > 0) {
            makePath(node->pointers[T], index);
            node->counters[T] = 0;
        }
        if (node->counters[G] > 0) {
            makePath(node->pointers[G], index);
            node->counters[G] = 0;
        }
    } else if (node->first) {
        makePath(node, index);
    }
}
void LFdBG::makePath(Node *n, const unsigned int &index) {
    Node       *node = n;
    std::string contig;
    contig = node->k_mer;
    int next_index;
    while (true) {
        next_index = node->onlyOnePath();
        if (next_index < 0)
            break;

        --node->counters[next_index];
        contig.append(1, getChar(next_index));
        node = node->pointers[next_index];
    }
    contigs_pools[index].push_back(contig);
}

// This method normalizes our created graph.
// In our implementation, counter is set to one when dividing counter by weight gives as number in range from 0.1
// exclusive to 1.5 exclusive.
void LFdBG::normalizeGraph(unsigned int start, unsigned int end) {
    if (redundancy < 1)
        return;

    for (unsigned int i = start; i < end; i++) {
        if (hash_table[i] == nullptr)
            continue;
        Node *node = hash_table[i];
        do {

            if (node->counters[A] > 0 && std::round(double(node->counters[A]) / redundancy) == 0 &&
                double(node->counters[A]) / redundancy > 0.1)
                node->counters[A] = 1;
            else
                node->counters[A] = std::round(double(node->counters[A]) / redundancy);
            if (node->counters[C] > 0 && std::round(double(node->counters[C]) / redundancy) == 0 &&
                double(node->counters[C]) / redundancy > 0.1)
                node->counters[C] = 1;
            else
                node->counters[C] = std::round(double(node->counters[C]) / redundancy);
            if (node->counters[T] > 0 && std::round(double(node->counters[T]) / redundancy) == 0 &&
                double(node->counters[T]) / redundancy > 0.1)
                node->counters[T] = 1;
            else
                node->counters[T] = std::round(double(node->counters[T]) / redundancy);
            if (node->counters[G] > 0 && std::round(double(node->counters[G]) / redundancy) == 0 &&
                double(node->counters[G]) / redundancy > 0.1)
                node->counters[G] = 1;
            else
                node->counters[G] = std::round(double(node->counters[G]) / redundancy);

            node = node->next;
        } while (node != nullptr);
    }
}

// This method returns found contigs.
std::deque<std::string> LFdBG::getContigs() {
    if (already_searched)
        return contigs;

    int                      for_one_thread = table_size / threads_for_assembly;
    std::vector<std::thread> threads2;

    for (int i = 0; i < threads_for_assembly; ++i) {
        if (i == threads_for_assembly - 1) {
            threads2.emplace_back(&LFdBG::findContigs, this, i * for_one_thread, table_size, i);
        } else {
            threads2.emplace_back(
                &LFdBG::findContigs, this, i * for_one_thread, (i + 1) * for_one_thread, i);
        }
    }

    for (int i = 0; i < threads_for_assembly; ++i)
        threads2[i].join();

    already_searched = true;
    for (long unsigned int i = 0; i < contigs_pools.size(); i++) {
        long unsigned int size = contigs_pools[i].size();
        for (long unsigned int j = 0; j < size; j++) {
            contigs.push_back(contigs_pools[i].front());
            contigs_pools[i].pop_front();
        }
    }
    return contigs;
}
std::vector<std::string> LFdBG::getGraph(bool withEdges) {
    std::vector<std::string> graph;
    for (unsigned int i = 0; i < table_size; i++) {
        if (hash_table[i] == nullptr)
            continue;
        Node *node = hash_table[i];
        while (node != nullptr) {
            if (withEdges)
                graph.push_back(std::string(node->k_mer) + " " + std::to_string(node->counters[A]) + " " +
                                std::to_string(node->counters[C]) + " " + std::to_string(node->counters[T]) + " " +
                                std::to_string(node->counters[G]));
            else
                graph.push_back(std::string(node->k_mer));

            node = node->next;
        }
    }
    return graph;
}
