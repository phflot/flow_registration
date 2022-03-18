// Author   : Philipp Flotho
// Copyright 2021 by Philipp Flotho, All rights reserved.

#include <mex.h>
#include <math.h>
#include <stdlib.h>

/* input arguments */
#define J11 0
#define J22 1
#define J33 2
#define J12 3
#define J13 4
#define J23 5
#define WEIGHT 6
#define U 7
#define V 8
#define ALPHA 9
#define ITERATIONS 10
#define UPDATE_LAG 11
#define VERBOSE 12
#define A_DATA 13
#define A_SMOOTH 14
#define HX 15
#define HY 16

/* output arguments */
#define DU 0
#define DV 1

/* constants */
#define OMEGA 1.95

void nonlinearity_smoothness(double* psi_smooth, mxArray *du, const mxArray *u, 
	mxArray *dv, const mxArray *v, double a, double hx, double hy);
void nonlinearity(double *psi,
	const double *j11, const double *j22, const double *j33, const double *j12,
	const double *j13, const double *j23, double *du, double *dv, int n,
	int n_channels, const double *a);
mxArray* p_sum(const mxArray* m1, const mxArray* m2);

void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {

	/* input arguments */
	double hx;
	double hy;
	int verbose; /* if 0 no terminal output */

	/* checking number of inputs */
	if (nrhs < 15)
		mexErrMsgIdAndTxt("OF:input_parameters",
			"at least 14 inputs required.");

	/* setting inputs */
	const double* j11 = mxGetPr(prhs[J11]);
	const double* j22 = mxGetPr(prhs[J22]);
	const double* j12 = mxGetPr(prhs[J12]);
	const double* j13 = mxGetPr(prhs[J13]);
	const double* j23 = mxGetPr(prhs[J23]);
	const double* j33 = mxGetPr(prhs[J33]);
	const double* data_weight = mxGetPr(prhs[WEIGHT]);

	const mxArray *u_mat = prhs[U];
	const mxArray *v_mat = prhs[V];

	const double* u = mxGetPr(prhs[U]);
	const double* v = mxGetPr(prhs[V]);
	mwSize ny = mxGetM(prhs[U]);
	mwSize nx = mxGetN(prhs[U]);

	mwSize weight_ny = mxGetM(prhs[WEIGHT]);
	mwSize weight_nx = mxGetN(prhs[WEIGHT]);
	mwSize weight_size = weight_ny * weight_nx;

	mwSize n_channels;
	if (mxGetNumberOfDimensions(prhs[J11]) < 3)
		n_channels = 1;
	else
		n_channels = mxGetDimensions(prhs[J11])[2];

	const double* a_data = mxGetPr(prhs[A_DATA]);
	const double* alpha = mxGetPr(prhs[ALPHA]);
    double a_smooth = mxGetScalar(prhs[A_SMOOTH]);
	int iterations = (int)mxGetScalar(prhs[ITERATIONS]);
	int update_lag = (int)mxGetScalar(prhs[UPDATE_LAG]);
	if (nrhs == 17) {
		hx = mxGetScalar(prhs[HX]);
		hy = mxGetScalar(prhs[HY]);
	}
	else {
		hx = 1.0;
		hy = 1.0;
	}

	verbose = (int) mxGetScalar(prhs[VERBOSE]);

	/* checking number of outputs */
	if (nlhs != 2) 
		mexErrMsgIdAndTxt("OF:output_parameters",
			"Two outputs required.");

	mxArray* du_mat = mxCreateDoubleMatrix(ny, nx, mxREAL);
	mxArray* dv_mat = mxCreateDoubleMatrix(ny, nx, mxREAL);

	plhs[DU] = du_mat;
	plhs[DV] = dv_mat;

	double* du = mxGetPr(du_mat);
	double* dv = mxGetPr(dv_mat);

	if (verbose) {
		printf("Starting OF calculation with alpha = (%f, %f) and %i iterations and update lag of %i. /n",
			alpha[0], alpha[1], iterations, update_lag);
	}

	mxArray* psi_mat = mxCreateNumericArray(
		mxGetNumberOfDimensions(prhs[J11]), 
		mxGetDimensions(prhs[J11]), mxDOUBLE_CLASS, mxREAL); //mxCreateDoubleMatrix(ny, nx, n_channels, mxREAL);
	double *psi = mxGetPr(psi_mat);
	mxArray* psi_smooth_mat = mxCreateDoubleMatrix(ny, nx, mxREAL);
	double* psi_smooth = mxGetPr(psi_smooth_mat);

	for (int i = 0; i < nx * ny; i++) {
		du[i] = 0;
		dv[i] = 0;
		psi_smooth[i] = 1;
		psi[i] = 1;
	}

	/* printf("hx = %f, hy = %f \n", hx, hy); */

	double w_right, w_left, w_down, w_up, w_sum;

	double smooth_u, smooth_v;

	int idx, idx_left, idx_right, idx_up, idx_down;
	
	double denom_u, denom_v, num_u, num_v, tmp, du_kp1, dv_kp1;
	const double weight = 1.0; // / (double)n_channels;

	int nd_idx;
	int iteration_counter = 0;

	double alpha_stencil[4] = {alpha[0] / (hx * hx), alpha[0] / (hx * hx), alpha[1] / (hy * hy), alpha[1] / (hy * hy)};
	int s_idx[4];

	int idx_inner;

	if (weight_size < nx * ny) {

        while (iteration_counter++ < iterations) {

            // updating the non-linearities, compare Bruhn dissertation chapter 4 and 5:
            // we treat a_smooth == 1 separately for computational efficiency, as we know psi_smooth = 1 in that case
            if (iteration_counter % update_lag == 0) {

                for (int k = 0; k < n_channels; k++) {
                    for (int i = 0; i < nx * ny; i++) {
                        idx = i + k * nx * ny;
                        tmp = j11[idx] * du[i] * du[i] + j22[idx] * dv[i] * dv[i] + j23[idx] * dv[i] +
                              +2 * j12[idx] * du[i] * dv[i] + 2 * j13[idx] * du[i] + j23[idx] * dv[i] + j33[idx];
                        tmp = tmp < 0 ? 0 : tmp;
                        psi[idx] = a_data[k] * pow(tmp +  0.00001, a_data[k] - 1);
                    }
                }

                if (a_smooth != 1) {
                    nonlinearity_smoothness(psi_smooth, du_mat, u_mat, dv_mat, v_mat, a_smooth, hx, hy);
                }
            }

            /// setting virtual nodes at the boundary
            for (int i = 0; i < nx; i++) {
                idx = 0 + i * ny;
                idx_inner = 1 + i * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];

                idx = ny - 1 + i * ny;
                idx_inner = ny - 2 + i * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];
            }

            for (int j = 0; j < ny; j++) {
                idx = j + 0 * ny;
                idx_inner = j + 1 * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];

                idx = j + (nx - 1) * ny;
                idx_inner = j + (nx - 2) * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];
            }

            for (int i = 1; i < nx - 1; i++) {
                for (int j = 1; j < ny - 1; j++) {

                    denom_u = 0;
                    denom_v = 0;
                    num_u = 0;
                    num_v = 0;

                    idx = j + i * ny;
                    s_idx[0] = j + (i - 1) * ny;
                    s_idx[1] = j + (i + 1) * ny;
                    s_idx[2] = (j + 1) + i * ny;
                    s_idx[3] = (j - 1) + i * ny;

                    smooth_u = 0;
                    smooth_v = 0;

                    // convolution for the smoothness term, compare Bruhn dissertation chapter 4 and 5,
                    // we treat a_smooth == 1 separately for computational efficiency
                    if (a_smooth != 1) {
                        for (int d = 0; d < 4; d++) {
                            tmp = 0.5 * (psi_smooth[idx] + psi_smooth[s_idx[d]]) * alpha_stencil[d];
                            num_u += tmp * (u[s_idx[d]] + du[s_idx[d]] - u[idx]);
                            num_v += tmp * (v[s_idx[d]] + dv[s_idx[d]] - v[idx]);
                            denom_u += tmp;
                            denom_v += tmp;
                        }
                    }
                    else {
                        for (int d = 0; d < 4; d++) {
                            num_u += alpha_stencil[d] * (u[s_idx[d]] + du[s_idx[d]] - u[idx]);
                            num_v += alpha_stencil[d] * (v[s_idx[d]] + dv[s_idx[d]] - v[idx]);
                            denom_u += alpha_stencil[d];
                            denom_v += alpha_stencil[d];
                        }
                    }

                    // summing up the data terms for each channel
                    for (int k = 0; k < n_channels; k++) {
                        nd_idx = j + i * ny + k * (nx * ny);

                        num_u -= data_weight[k] * psi[nd_idx] * (j13[nd_idx] + j12[nd_idx] * dv[idx]);

                        denom_u += data_weight[k] * psi[nd_idx] * j11[nd_idx];
                        denom_v += data_weight[k] * psi[nd_idx] * j22[nd_idx];
                    }
                    du_kp1 = num_u / denom_u;

                    // SOR interpolation step:
                    du[idx] = (1 - OMEGA) * du[idx] + OMEGA * du_kp1;

                    for (int k = 0; k < n_channels; k++) {
                        nd_idx = j + i * ny + k * (nx * ny);
                        num_v -= data_weight[k] * psi[nd_idx] * (j23[nd_idx] + j12[nd_idx] * du[idx]);
                    }
                    dv_kp1 = num_v / denom_v;

                    // SOR interpolation step:
                    dv[idx] = (1 - OMEGA) * dv[idx] + OMEGA *  dv_kp1;
                }
            }
        }
	} else {
                
        while (iteration_counter++ < iterations) {

            // updating the non-linearities, compare Bruhn dissertation chapter 4 and 5:
            // we treat a_smooth == 1 separately for computational efficiency, as we know psi_smooth = 1 in that case
            if (iteration_counter % update_lag == 0) {

                for (int k = 0; k < n_channels; k++) {
                    for (int i = 0; i < nx * ny; i++) {
                        idx = i + k * nx * ny;
                        tmp = j11[idx] * du[i] * du[i] + j22[idx] * dv[i] * dv[i] + j23[idx] * dv[i] +
                              +2 * j12[idx] * du[i] * dv[i] + 2 * j13[idx] * du[i] + j23[idx] * dv[i] + j33[idx];
                        tmp = tmp < 0 ? 0 : tmp;
                        psi[idx] = a_data[k] * pow(tmp +  0.00001, a_data[k] - 1);
                    }
                }

                if (a_smooth != 1) {
                    nonlinearity_smoothness(psi_smooth, du_mat, u_mat, dv_mat, v_mat, a_smooth, hx, hy);
                }
            }

            /// setting virtual nodes at the boundary
            for (int i = 0; i < nx; i++) {
                idx = 0 + i * ny;
                idx_inner = 1 + i * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];

                idx = ny - 1 + i * ny;
                idx_inner = ny - 2 + i * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];
            }

            for (int j = 0; j < ny; j++) {
                idx = j + 0 * ny;
                idx_inner = j + 1 * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];

                idx = j + (nx - 1) * ny;
                idx_inner = j + (nx - 2) * ny;
                du[idx] = du[idx_inner];
                dv[idx] = dv[idx_inner];
            }

            for (int i = 1; i < nx - 1; i++) {
                for (int j = 1; j < ny - 1; j++) {

                    denom_u = 0;
                    denom_v = 0;
                    num_u = 0;
                    num_v = 0;

                    idx = j + i * ny;
                    s_idx[0] = j + (i - 1) * ny;
                    s_idx[1] = j + (i + 1) * ny;
                    s_idx[2] = (j + 1) + i * ny;
                    s_idx[3] = (j - 1) + i * ny;

                    smooth_u = 0;
                    smooth_v = 0;

                    // convolution for the smoothness term, compare Bruhn dissertation chapter 4 and 5,
                    // we treat a_smooth == 1 separately for computational efficiency
                    if (a_smooth != 1) {
                        for (int d = 0; d < 4; d++) {
                            tmp = 0.5 * (psi_smooth[idx] + psi_smooth[s_idx[d]]) * alpha_stencil[d];
                            num_u += tmp * (u[s_idx[d]] + du[s_idx[d]] - u[idx]);
                            num_v += tmp * (v[s_idx[d]] + dv[s_idx[d]] - v[idx]);
                            denom_u += tmp;
                            denom_v += tmp;
                        }
                    }
                    else {
                        for (int d = 0; d < 4; d++) {
                            num_u += alpha_stencil[d] * (u[s_idx[d]] + du[s_idx[d]] - u[idx]);
                            num_v += alpha_stencil[d] * (v[s_idx[d]] + dv[s_idx[d]] - v[idx]);
                            denom_u += alpha_stencil[d];
                            denom_v += alpha_stencil[d];
                        }
                    }

                    // summing up the data terms for each channel
                    for (int k = 0; k < n_channels; k++) {
                        nd_idx = j + i * ny + k * (nx * ny);

                        num_u -= data_weight[nd_idx] * psi[nd_idx] * (j13[nd_idx] + j12[nd_idx] * dv[idx]);

                        denom_u += data_weight[nd_idx] * psi[nd_idx] * j11[nd_idx];
                        denom_v += data_weight[nd_idx] * psi[nd_idx] * j22[nd_idx];
                    }
                    du_kp1 = num_u / denom_u;

                    // SOR interpolation step:
                    du[idx] = (1 - OMEGA) * du[idx] + OMEGA * du_kp1;

                    for (int k = 0; k < n_channels; k++) {
                        nd_idx = j + i * ny + k * (nx * ny);
                        num_v -= data_weight[nd_idx] * psi[nd_idx] * (j23[nd_idx] + j12[nd_idx] * du[idx]);
                    }
                    dv_kp1 = num_v / denom_v;

                    // SOR interpolation step:
                    dv[idx] = (1 - OMEGA) * dv[idx] + OMEGA *  dv_kp1;
                }
            }
        }
    }
	mxDestroyArray(psi_smooth_mat);
	mxDestroyArray(psi_mat);
}

void nonlinearity_smoothness(double* psi_smooth, mxArray *du, const mxArray *u, mxArray *dv, 
	const mxArray *v, double a, double hx, double hy) {

	double eps = 0.00001;

	int m = mxGetM(du);
	int n = mxGetN(du);

	mxArray* u_full = p_sum(u, du);
	mxArray* v_full = p_sum(v, dv);

    mxArray* G_u[2];
    mxArray* mx_hx = mxCreateDoubleScalar(hx);
    mxArray* mx_hy = mxCreateDoubleScalar(hy);
    mxArray* rhs_u[3] = { u_full, mx_hx, mx_hy };
    mexCallMATLAB(2, G_u, 3, rhs_u, "gradient");
    mxDestroyArray(mx_hx);
    mxDestroyArray(mx_hy);
    
    mxArray* G_v[2];
    mx_hx = mxCreateDoubleScalar(hx);
    mx_hy = mxCreateDoubleScalar(hy);
    mxArray* rhs_v[3] = { v_full, mx_hx, mx_hy };
    mexCallMATLAB(2, G_v, 3, rhs_v, "gradient");
    mxDestroyArray(mx_hx);
    mxDestroyArray(mx_hy);
    
//     mexPutVariable("base", "test", G_v[0]);

	double *ux = mxGetPr(G_u[0]);
	double *uy = mxGetPr(G_u[1]);
	double *vx = mxGetPr(G_v[0]);
	double *vy = mxGetPr(G_v[1]);
	
	double tmp;
	for (int i = 0; i < m * n; i++) {
		tmp = ux[i] * ux[i] + uy[i] * uy[i] + vx[i] * vx[i] + vy[i] * vy[i];
		tmp = tmp < 0 ? 0 : tmp;
		psi_smooth[i] = a * pow(tmp + eps, a - 1);
	}

	
	mxDestroyArray(u_full);
	mxDestroyArray(v_full);
	for (int i = 0; i < 2; i++) {
		mxDestroyArray(G_u[i]);
		mxDestroyArray(G_v[i]);
	}
}

mxArray* p_sum(const mxArray* m1, const mxArray* m2) {
    mxArray* out[1];
    mxArray* tmp1 = mxDuplicateArray(m1);
    mxArray* tmp2 = mxDuplicateArray(m2);
    mxArray* rhs[2] = { tmp1, tmp2 };
    mexCallMATLAB(1, out, 2, rhs, "+");
    mxDestroyArray(tmp1);
    mxDestroyArray(tmp2);
    return out[0];
}