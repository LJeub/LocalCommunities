//
//  main.cpp
//  APPR
//
//  Created by Lucas Jeub on 29/01/2015.
//
//

#include "main.h"

#define MAXITER 1000000

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nrhs!=5) {
        mexErrMsgIdAndTxt("APPR:input", "APPR needs 5 input arguments");
    }
    //input:
    
    double alpha= * mxGetPr(prhs[0]);
    double epsilon= * mxGetPr(prhs[1]);
    
    full seed(prhs[2]); //make sure to copy as this will be modified
    sparse W(prhs[3]);
    full d(prhs[4]);
    
    full r=seed;
    
    //number of nodes
    mwIndex N=W.m;
    
    
    //initialise queue
    queue<mwIndex> update;
    
    full p(N,1);
    for (mwIndex i=0; i<N; ++i) {
        p.get(i)=0;
    }
    
    
    for (mwIndex i=0; i<N; ++i) {
        if (r.get(i)>epsilon*d.get(i)) {
            update.push(i);
        }
    }
    
    mwIndex iter=0;
    
    while (!update.empty() && iter < MAXITER ) {
        push(p, r, d, W, update, alpha, epsilon);
        ++iter;
    }
    
    p.export_matlab(plhs[0]);
    if (nlhs>1) {
       plhs[1]=mxCreateLogicalScalar(iter==MAXITER);
    }
    if (nlhs>2) {
        r.export_matlab(plhs[2]);
    }

}

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
