# LFdBG (Lock-Free-de-Bruijn-Graph)

LFdBG is fast and easy to work with C++ library that uses Lock-Free structures for de-Novo genome assembly. 

## Table of contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)


## Installation

For the library to work properly it needs Boost library to be installed on the computer. On linux it can be isntalled with this command.

    sudo apt-get install libboost-all-dev
This library needs version 11 of C++, because of std::atomic.

Other than that using LFdBG library does not need anything else installed.

## Usage
When creating the LFdBG object, the user must provide a few fundamental parameters: the graph dimension (lenght of k-mers), the number of threads to be used for assembly, the size of the hash table, average read length and approximate genome size.
There is also type bool parameter called normalize. If we do not want to try to normalize our graph we can simply set it to false.

The LFdBG class provides two public methods:

* Using the getContigs method we can get a copy of the created contigs. When used for the first time, the algorithm will search the graph and remember the paths found. On following calls, the previously saved contigs will be returned.

* The getGraph method returns the graph as a set of all created nodes with their mers as a series of strings. If we also want to obtain weights, we must pass the boolean value true as a parameter of the method.

## Licence
This library with source codes is available under MIT license.

Copyright (c) 2022  Daniel Gorniak, Robert Nowak, Warsaw University of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
