/*
 * arrayProduct.c - example in MATLAB External mwSizeerfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "matrix.h"
#include <omp.h>
#include "vema.h"
#include "eig3.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

Vector closestPointTriangle(Vector&, Vector&, Vector&, Vector&, double&, double&, double&);

// Generates pomwSize-triangle proximity lists using the linked cell algorithm
vector<std::vector<int>>& createNNLtriangle(vector<std::vector<int>>& NNLt, vector<std::vector<double>> V, vector<std::vector<int>> Fb, vector<int> SN, mwSize nsn, mwSize nf, double hs, double bw, double mw) {

    Vector* Ut = new Vector[V.size()]();
    Vector* faces = new Vector[Fb.size()]();
    for (int i=0; i<V.size(); i++){
            Ut[i].x = V[i][0];
            Ut[i].y = V[i][1];
            Ut[i].z = V[i][2];
    }
    for (int i=0; i<Fb.size(); i++){
            faces[i].x = Fb[i][0];
            faces[i].y = Fb[i][1];
            faces[i].z = Fb[i][2];
    }
      
	int mx = max(1, (int)(bw/mw)); // ** = 40 cells bw=3.2, mw=0.08
    vector<int> head(mx*mx*mx, -1);
	vector<int> list(nf);
// 	std::vector<int> head(mx*mx*mx, -1); //****** mx*mx*mx cells nomber, size mx*mx*mx vector with all values are -1, 40*40*40 = 64000
// 	std::vector<int> list(nf); // **** nf = 101882
	int xa, ya, za, xi, yi, zi;
	double ub, vb, wb;
	int pt, tri;
	Vector cog;
	for (int i = 0; i < nf; i++) { // Divide triangle faces mwSizeo cells, i index of face
		//Vector cog = (Ut[faces[i].n1] + Ut[faces[i].n2] + Ut[faces[i].n3])/3.0;
        cog = (Ut[(int)faces[i].x] + Ut[(int)faces[i].y] + Ut[(int)faces[i].z])/3.0;
		int xa = (int)((cog.x + 0.5*bw)/bw*mx);
		int ya = (int)((cog.y + 0.5*bw)/bw*mx);
		int za = (int)((cog.z + 0.5*bw)/bw*mx);
        int tmp =  mx*mx*za + mx*ya + xa; // *** 1641838 > 64000

		list[i]=head[mx*mx*za + mx*ya + xa];
		head[mx*mx*za + mx*ya + xa] = i;
	}
	#pragma omp parallel for
	for (int i = 0; i < nsn; i++) { // Search cells around each pomwSize and build proximity list
		int pt = SN[i];
		NNLt[i].clear();
		int xa = (int)((Ut[pt].x + 0.5*bw)/bw*mx);
		int ya = (int)((Ut[pt].y + 0.5*bw)/bw*mx);
		int za = (int)((Ut[pt].z + 0.5*bw)/bw*mx);

		for (int xi = max(0, xa-1); xi <= min(mx-1, xa+1); xi++)// *** Browse head list
		for (int yi = max(0, ya-1); yi <= min(mx-1, ya+1); yi++)
		for (int zi = max(0, za-1); zi <= min(mx-1, za+1); zi++) {
			int tri = head[mx*mx*zi + mx*yi + xi];
			while (tri != -1) {
				if ( pt != (int)faces[tri].x && pt != (int)faces[tri].y && pt != (int)faces[tri].z ) {				
					if ( (closestPointTriangle(Ut[pt], Ut[(int)faces[tri].x], Ut[(int)faces[tri].y], Ut[(int)faces[tri].z], ub, vb, wb) - Ut[pt]).length() < hs) {
						NNLt[i].push_back(tri);
					}
				}
				tri = list[tri];
			}
		}
	}
    return NNLt;
}


// Returns the closest pomwSize of triangle abc to pomwSize p ***** a or b or c, if not, pt projection through the barycenter inside the triangle 
Vector closestPointTriangle(Vector& p, Vector& a, Vector& b, Vector& c, double& u, double& v, double& w) {
	
	Vector ab = b - a;
	Vector ac = c - a;
	Vector ap = p - a;
	double d1 = ab.dot(ap);
	double d2 = ac.dot(ap);
	if (d1 <= 0.0 && d2 <= 0.0) {
		u = 1.0;
		v = 0.0;
		w = 0.0;
		return a;
	}
	Vector bp = p - b;
	double d3 = ab.dot(bp);
	double d4 = ac.dot(bp);
	if (d3 >= 0.0 && d4 <= d3) {
		u = 0.0;
		v = 1.0;
		w = 0.0;
		return b;
	}
	double vc = d1*d4 - d3*d2;
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
		v = d1 / (d1 - d3);
		u = 1.0 - v;
		w = 0.0;
		return a + ab * v;
	}
	Vector cp = p - c;
	double d5 = ab.dot(cp);
	double d6 = ac.dot(cp);
	if (d6 >= 0.0 && d5 <= d6) {
		u = 0.0;
		v = 0.0;
		w = 1.0;
		return c;
	}
	double vb = d5*d2 - d1*d6;
	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
		w = d2 / (d2 - d6);
		u = 1.0 - w;
		v = 0.0;	
		return a + ac * w;
	}
	double va = d3*d6 - d5*d4;
	if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
		w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		u = 0.0;
		v = 1.0 - w;
		return b + (c - b) * w;
	}
	double denom = 1.0 / (va + vb + vc);
	v = vb * denom;
	w = vc * denom;
	u = 1.0 - v - w;
	return a + ab * v + ac * w;
}


/* The gateway function */
void mexFunction(mwSize nlhs, mxArray *plhs[],
                 mwSize nrhs, const mxArray *prhs[])
{   
//     double multiplier;      /* input scalar */
//     double *inMatrix;       /* 1xN input matrix */
//     mwSize ncols;           /* size of matrix */
//     double *outMatrix;      /* output matrix */
//     vector<int>* NNLt;
//     Vector* Ut;
//     vector<Face>& faces;
    
//     vector<int>* NNLt[prhs[4]];
//     Vector *V = new Vector[mxGetN(prhs[1]),]();
    vector<std::vector<double>> V(mxGetN(prhs[1]), std::vector<double>(3));
    vector<std::vector<int>> Fb(mxGetN(prhs[2]), std::vector<int>(3));
//     Vector *Fb = new Vector[mxGetM(prhs[2])]();
//     int *sn;
    mwSize nsn; 
    mwSize nf; 
    double hs;
    double bw;
    double mw; 
//     mwSize ncols; 
    mwSize i;
    mwSize j;
           
    /* check for proper number of arguments */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Nine inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }    
  
    
    /* get the value of the scalar input  */
    nsn = mxGetScalar(prhs[4]);
    nf = mxGetScalar(prhs[5]);
    hs = mxGetScalar(prhs[6]);  
    bw = mxGetScalar(prhs[7]);
    mw = mxGetScalar(prhs[8]);
    /* create a pomwSizeer to the real data in the input matrix  */
//     NNLt = (vector<int> *)mxGetData(prhs[0]);
//     V = (Vector *)mxGetData(prhs[1]);
//     V = static_cast<double *>(mxGetData(prhs[1]));
//     Fb = static_cast<double *>(mxGetData(prhs[2]));
//     Fb = (Vector *)mxGetData(prhs[2]);
//     sn = (int *)mxGetData(prhs[3]);
    vector<std::vector<int>> NNLt(nsn, std::vector<int>(1));
    
//     for(i=0; i<nsn; i++) {
//         mxArray *tmp = mxGetCell(prhs[0],i);
// //         std::size_t M = mxGetM(tmp); // number of rows        
// //         std::size_t N = mxGetN(tmp); // number of columns
// //         double* ptr = mxGetPr(tmp);
// //         copy(ptr, ptr, NNLt[i].begin());         
//    }    
   double* ptr = mxGetPr(prhs[1]);
   for(i=0; i<mxGetN(prhs[1]); i++) {
        copy(ptr,ptr+3, V[i].begin()); 
        ptr += 3;
   }


   double* ptrr = mxGetPr(prhs[2]);
   for(i=0; i<mxGetN(prhs[2]); i++) {
        copy(ptrr,ptrr+3, Fb[i].begin()); 
        ptrr += 3;
   }
    
   std::vector<int> sn(nsn);
   copy(mxGetPr(prhs[3]),mxGetPr(prhs[3])+nsn, sn.begin());

//     double* ptr = mxGetPr(prhs[1]);
//     for (i=0; i<mxGetM(prhs[1]); i++){
//     copy(ptr, ptr+3, V.begin());
//     for(i=0; i<nsn; i++) {
//         mxArray *tmp = mxGetCell(prhV[i].begin()s[0],i);
//         std::size_t M = mxGetM(tmp); // number of rows
//         vector<std::vector<int>> NNLt(nsn, std::vector<int>(M));
//         double* ptr = mxGetPr(tmp);
//         copy(ptr, ptr+M, NNLt[i].begin());  
//     };    

    /* get dimensions of the input matrix */
//     ncols = mxGetN(prhs[0]); 

    /* call the computational routine */    
    NNLt = createNNLtriangle(NNLt, V, Fb, sn, nsn, nf, hs, bw, mw); 
//     for (i=0;i<nsn;i++){
// //         for (j=0;j<NNLt[i].size();j++){
//             mexPrintf("%d",NNLt[i][0]);       
// //         }
//     }
    /* create the output matrix */
//     plhs[0] = mxCreateCellMatrix(nsn,50);
//     plhs[0] = (mxArray *)mxGetData(NNLt);
    plhs[0] = mxCreateCellMatrix(1,nsn);
//        NNLtout = mxGetPr(plhs[0]);u
    /* get a pomwSizeer to the real data in the output matrix */
//     NNLtout = plhs[0];
    for(i=0;i<nsn;i++){
//             copy(NNLt[i].begin(),NNLt[i].end(),NNLtout[i*50;i*50+NNLt[i].size()]);
//            NNLtoutt=mxCreatStrucMatrix(1,50,1,fnom);
        
         mxArray* tmp = mxCreateDoubleMatrix(1, NNLt[i].size(), mxREAL);
         copy(NNLt[i].begin(), NNLt[i].end(), mxGetPr(tmp));     
         mxSetCell(plhs[0], i, tmp);
    }
//         NNLtout = (double *)mxGetData(NNLt);
//     mxSetCell(NNLtout,mxDuplicateArray(mxGetPr(NNLt[i])));

//         memcpy(mxGetData(plhs[0]),NNLt[0],500);

}
