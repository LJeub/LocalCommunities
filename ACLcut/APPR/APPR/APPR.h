// Approximate Personalised PageRank computation using the algorithm from:
//
//      Andersen, R., Chung, F. R. K., & Lang, K. J. (2006).
//      Local Graph Partitioning using PageRank Vectors (pp. 475–486).
//      in Proceedings of the 47th Annual Symposium on Foundations of Computer
//      Science, IEEE. http://doi.org/10.1109/FOCS.2006.44
//
// Version: 2.0.2
// Date: Wed 20 Jun 2018 16:01:01 CEST
// Author: Lucas Jeub
// Email: lucasjeub@gmail.com

#ifndef __APPR__main__
#define __APPR__main__

#include <stdio.h>
#include "matlab_matrix.h"
#include "mex.h"
#include "matrix.h"
#include <queue>
#include <vector>

void push(full & p, full & r, const full & d, const sparse & W, std::queue<mwIndex> & update, double alpha, double epsilon);

#endif /* defined(__APPR__main__) */
