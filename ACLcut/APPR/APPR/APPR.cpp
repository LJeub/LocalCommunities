// Approximate Personalised PageRank computation using the algorithm from:
//
//      Andersen, R., Chung, F. R. K., & Lang, K. J. (2006).
//      Local Graph Partitioning using PageRank Vectors (pp. 475–486).
//      in Proceedings of the 47th Annual Symposium on Foundations of Computer
//      Science, IEEE. http://doi.org/10.1109/FOCS.2006.44
//
// Version: 2.0
// Date: Mon 25 Jul 2016 17:06:57 BST
// Author: Lucas Jeub
// Email: jeub@maths.ox.ac.uk

#include "APPR.h"

#define MAXITER 1000000

using namespace std;

// Provides APPR function:
// [p, flag, r] = APPR(alpha, epsilon, seed, W, d)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    //input check
    if (nrhs!=5) {
        mexErrMsgIdAndTxt("APPR:input", "APPR needs 5 input arguments");
    }
    
    //input:
    double alpha= * mxGetPr(prhs[0]);
    double epsilon= * mxGetPr(prhs[1]);
    full seed(prhs[2]); //make sure to copy as this will be modified
    full r=seed;

    sparse W(prhs[3]);
    full d(prhs[4]);
    
    
    //number of nodes
    mwIndex N=W.m;
    
    
    //initialise APPR vector
    full p(N,1);
    for (mwIndex i=0; i<N; ++i) {
        p.get(i)=0;
    }
    
    
    //initialise queue
    queue<mwIndex> update;
    for (mwIndex i=0; i<N; ++i) {
        if (r.get(i)>epsilon*d.get(i)) {
            update.push(i);
        }
    }
    
    
    //push until convergence
    mwIndex iter=0;
    while (!update.empty() && iter < MAXITER ) {
        push(p, r, d, W, update, alpha, epsilon);
        ++iter;
    }
    
    //output:
    p.export_matlab(plhs[0]);
    if (nlhs>1) {
       plhs[1]=mxCreateLogicalScalar(iter==MAXITER);
    }
    if (nlhs>2) {
        r.export_matlab(plhs[2]);
    }

}


//push algorithm
void push(full & p, full & r, const full & d, const sparse & W, queue<mwIndex> & update, double alpha, double epsilon) {
    mwIndex node=update.front();
    update.pop();
    if (r.get(node)>epsilon*d.get(node)) {
        double r_node=r.get(node);
        p.get(node)+= alpha*r_node; //push to p
        r.get(node)=(1-alpha)*r_node/2;
        if (r.get(node)>epsilon*d.get(node)) {
            update.push(node);
        }
        
        for (mwIndex i=W.col[node]; i<W.col[node+1]; ++i) {
            mwIndex target = W.row[i];
            double delta = (1-alpha)*r_node*W.val[i]/(2*d.get(node));
            r.get(target)+=delta;
            if (0<r.get(target)-epsilon*d.get(target) && r.get(target)-epsilon*d.get(target)<=delta) {
                update.push(target);
            }
        }
    }
}
