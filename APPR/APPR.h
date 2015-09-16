//
//  main.h
//  APPR
//
//  Created by Lucas Jeub on 29/01/2015.
//
//

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
