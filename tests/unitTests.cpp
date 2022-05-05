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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE UnitTests
#include "../include/LFdBG.h"
#include "../src/LFdBG.cpp"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <thread>
#include <vector>
bool findMer(std::vector<std::string> &kmers, std::string mer) {
    for (auto kmer : kmers)
        if (kmer == mer)
            return true;
    return false;
}
BOOST_AUTO_TEST_CASE(kmers) {
    std::vector<std::string> reads;
    reads.push_back("ACTGACCTGG");
    LFdBG d        = LFdBG(6, reads, 10, 10, 1, 100, false);
    auto  vertices = d.getGraph(false);
    BOOST_ASSERT(findMer(vertices, "ACTGA"));
    BOOST_ASSERT(findMer(vertices, "CTGAC"));
    BOOST_ASSERT(findMer(vertices, "TGACC"));
    BOOST_ASSERT(findMer(vertices, "GACCT"));
    BOOST_ASSERT(findMer(vertices, "ACCTG"));
    BOOST_ASSERT(findMer(vertices, "CCTGG"));
    BOOST_ASSERT(vertices.size() == 6);
}
BOOST_AUTO_TEST_CASE(weights) {
    std::vector<std::string> reads;
    reads.push_back("ACTGACCTGG");
    LFdBG d        = LFdBG(6, reads, 10, 10, 1, 100, false);
    auto  vertices = d.getGraph(true);
    BOOST_ASSERT(findMer(vertices, "GACCT 0 0 0 1"));
    BOOST_ASSERT(findMer(vertices, "CTGAC 0 1 0 0"));
    BOOST_ASSERT(findMer(vertices, "ACCTG 0 0 0 1"));
    BOOST_ASSERT(findMer(vertices, "ACTGA 0 1 0 0"));
    BOOST_ASSERT(findMer(vertices, "TGACC 0 0 1 0"));
    BOOST_ASSERT(findMer(vertices, "CCTGG 0 0 0 0"));
    BOOST_ASSERT(vertices.size() == 6);
}
BOOST_AUTO_TEST_CASE(weights_2) {
    std::vector<std::string> reads;
    reads.push_back("ACTGACCTGG");
    reads.push_back("ACTGACCTGG");
    LFdBG d        = LFdBG(6, reads, 10, 10, 1, 100, false);
    auto  vertices = d.getGraph(true);
    BOOST_ASSERT(findMer(vertices, "GACCT 0 0 0 2"));
    BOOST_ASSERT(findMer(vertices, "CTGAC 0 2 0 0"));
    BOOST_ASSERT(findMer(vertices, "ACCTG 0 0 0 2"));
    BOOST_ASSERT(findMer(vertices, "ACTGA 0 2 0 0"));
    BOOST_ASSERT(findMer(vertices, "TGACC 0 0 2 0"));
    BOOST_ASSERT(findMer(vertices, "CCTGG 0 0 0 0"));
    BOOST_ASSERT(vertices.size() == 6);
}
BOOST_AUTO_TEST_CASE(weights_many_threads) {
    std::vector<std::string> reads;
    reads.push_back("ACTGACCTGG");
    reads.push_back("ACTGACCTGG");
    reads.push_back("ACTGACCTGG");
    reads.push_back("ACTGACCTGG");
    LFdBG d        = LFdBG(6, reads, 10, 10, 4, 1, false);
    auto  vertices = d.getGraph(true);
    BOOST_ASSERT(findMer(vertices, "GACCT 0 0 0 4"));
    BOOST_ASSERT(findMer(vertices, "CTGAC 0 4 0 0"));
    BOOST_ASSERT(findMer(vertices, "ACCTG 0 0 0 4"));
    BOOST_ASSERT(findMer(vertices, "ACTGA 0 4 0 0"));
    BOOST_ASSERT(findMer(vertices, "TGACC 0 0 4 0"));
    BOOST_ASSERT(findMer(vertices, "CCTGG 0 0 0 0"));
    BOOST_ASSERT(vertices.size() == 6);
}
BOOST_AUTO_TEST_CASE(contigs) {
    std::vector<std::string> reads;
    reads.push_back("ACTGACCTGG");
    LFdBG d       = LFdBG(6, reads, 10, 10, 1, 100);
    auto  contigs = d.getContigs();
    for (auto contig : contigs)
        BOOST_ASSERT(reads[0].find(contig, 0) != std::string::npos);
}