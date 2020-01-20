/**
 * Copyright (c) 2019 d-gfx
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */
#ifndef dgfx_vex_utils_h
#define dgfx_vex_utils_h

#define mod(a, b)		((a) % (b))
#define mix(a, b, t)	lerp((a), (b), (t))
#define clamp01(v)		clamp(v, 0.0, 1.0)

function vector dgfx_Calc_HeatMap_Color(float value_01)
{
	vector colors[] = array({0,0,1}, {0,1,0}, {1,1,0}, {1,0,0}); // enable to add color table.
	int num_colors = len(colors);
	float fract = 0;
	float value = clamp(value_01, 0, 1) * (num_colors - 1);
	int idx1 = floor(value);
	int idx2 = clamp(idx1+1, 0, num_colors-1);
	fract = frac(value);
	return lerp(colors[idx1], colors[idx2], fract);
}

/**
 *	Quadratic Bezier Curve	
 */
function vector dgfx_BezierQuadratic_Position(float t_01; vector p0, p1, p2)
{
	float t = clamp(t_01, 0.0, 1.0);
	return (1-t)*(1-t) * p0 + 2*(1-t)*t * p1 + t*t * p2;
}

function vector dgfx_BezierQuadratic_Derivative_1st(float t_01; vector p0, p1, p2)
{
	float t = clamp(t_01, 0.0, 1.0);
	return 2*(1-t) * (p1 - p0) + 2*t*(p2 - p1);
}

function vector dgfx_BezierQuadratic_Derivative_2nd(float t_01; vector p0, p1, p2)
{
	float t = clamp(t_01, 0.0, 1.0);
	return 2 * (p2 - 2*p1 + p0);
}

/**
 *	Cubic Bezier Curve
 */
function vector dgfx_BezierCubic_Position(float t_01; vector p0, p1, p2, p3)
{
	float t = clamp(t_01, 0.0, 1.0);
	#if 0
	return (1 - t) * bezierQuadratic_Position(t, p0, p1, p2) + t * bezierQuadratic_Position(t, p1, p2, p3);
	#else
	// explicit form
	float rev_t = clamp(1.0 - t, 0.0, 1.0);
	return rev_t*rev_t*rev_t * p0 + 3*rev_t*rev_t*t * p1 + 3*rev_t*t*t * p2 + t*t*t * p3;
	#endif
}

function vector dgfx_BezierCubic_Derivative_1st(float t_01; vector p0, p1, p2, p3)
{
	float t = clamp(t_01, 0.0, 1.0);
	float rev_t = clamp(1.0 - t, 0.0, 1.0);
	return 3*rev_t*rev_t * (p1 - p0) + 6*rev_t*t * (p2 - p1) + 3*t*t * (p3-p2);
}

function vector dgfx_BezierCubic_Derivative_2nd(float t_01; vector p0, p1, p2, p3)
{
	float t = clamp(t_01, 0.0, 1.0);
	float rev_t = clamp(1.0 - t, 0.0, 1.0);
	return 6*rev_t * (p2 - 2*p1 + p0) + 6*t * (p3 - 2*p2 + p1);
}

/**
 *	is Primitive Vertex Looped ?
 */
function int dgfx_IsLooped_PrimVertex(int geometry, prim_num)
{
	int num_vtx = primvertexcount(geometry, prim_num);
	int src_vtx_1st = primvertex(geometry, prim_num, 0);
	int src_vtx_end = primvertex(geometry, prim_num, num_vtx-1);
	int src_pt_1st = vertexpoint(geometry, src_vtx_1st);
	int src_pt_end = vertexpoint(geometry, src_vtx_end);
	return (src_pt_1st == src_pt_end);
}

/**
 *	append polyline from array of point number
 */
function void dgfx_Append_PolyLine_from_Points(vector pts[]; int is_prim_looped)
{
	int num_pt = len(pts);
	int prim = addprim(0, "polyline");
	int first_pt = -1;
	for (int i=0; i<num_pt; ++i)
	{
		int pt = addpoint(0, pts[i]);
		addvertex(0, prim, pt);
		if (i == 0)
		{
			first_pt = pt;
		}
	}
	if (is_prim_looped == 1)
	{
		addvertex(0, prim, first_pt);
	}
}

/**
 *	append mid point to edges to create rounded corner
 */
function void dgfx_Append_Mid_Point_Edge_Array(vector pts[]; int src_geo, src_prim, num_vtx, is_prim_looped, order; float round_rate)
{
	int num_pt = (is_prim_looped) ? num_vtx-1 : num_vtx;
	int num_corner = (is_prim_looped) ? num_pt : num_pt-2;
	float round_01 = (order == 3) ? 0.0 : clamp(round_rate, 0, 1);
	float round = round_01 * 0.5;
	int div_num = order-1;
	int is_exist_corner_vtx = ((order % 2) == 1);
	float branch_index = float(div_num)/2.0;
	int end_j_1st = int(floor(branch_index));
	int start_j_2nd = int(ceil(branch_index));
	if (is_exist_corner_vtx)
	{
		start_j_2nd++;
		end_j_1st--;
	}
//	printf("num_vtx = %d, order = %d, num_corner = %d, round_01 = %f\n", num_vtx, order, num_corner, round_01);
//	printf("branch_index = %d, end_j_1st = %d, start_j_2nd = %d, div_num = %d, is_prim_looped = %d\n", branch_index, end_j_1st, start_j_2nd, div_num, is_prim_looped);
	for (int i=0; i<num_corner; ++i)
	{
		int vtx_0 = primvertex(src_geo, src_prim, (i)%num_vtx);
		int vtx_1 = primvertex(src_geo, src_prim, (i+1)%num_vtx);
		int vtx_2 = primvertex(src_geo, src_prim, (i+2)%num_vtx);
		int pt_0  = vertexpoint(src_geo, vtx_0);
		int pt_1  = vertexpoint(src_geo, vtx_1);
		int pt_2  = vertexpoint(src_geo, vtx_2);
		if (pt_1 == pt_2)
		{
			vtx_2 = (vtx_2 + 1) % num_vtx;
			pt_2 = vertexpoint(src_geo, vtx_2);
		}
//		printf("pts = [%d, %d, %d]\n", pt_0, pt_1, pt_2);
		vector pos_0 = point(src_geo, "P", pt_0);
		vector pos_1 = point(src_geo, "P", pt_1);
		vector pos_2 = point(src_geo, "P", pt_2);

		// if prim is not looped, add first straight line control points
		if ((i==0) && (is_prim_looped == 0))
		{
			int append_1st_num = order-1; // number of control points
			for (int s=0; s<append_1st_num; ++s)
			{
				float rate = fit(clamp(float(s)/(append_1st_num), 0, 1), 0, 1, 0, 0.5);
				vector pos = lerp(pos_0, pos_1, rate);
				append(pts, pos);
//				printf("01 : [%d - %d]rate = %f\n", pt_0, pt_1, rate);
			}
		}

		// Resample for Curve
		for (int j=0; j<div_num; ++j)
		{
			vector pos_s, pos_e;
			float rate = 0.0;
			if (j < branch_index)
			{
				pos_s = pos_0; pos_e = pos_1;
				rate = clamp(float(j) / end_j_1st, 0.0, 1.0);
				rate = fit(rate, 0.0, 1.0, 0.5, 0.5 + 0.5 * (1.0-round_01));
			//	printf("02 : [%d - %d]rate = %f j = %d\n", pt_0, pt_1, rate, j);
			}
			else if (j == branch_index) // just corner vertex
			{
				pos_s = pos_0; pos_e = pos_1;
				rate = 1.0;
			//	printf("03 : [%d - %d]rate = %f j = %d\n", pt_0, pt_1, rate, j);
			}
			else
			{
				pos_s = pos_1; pos_e = pos_2;
				rate = clamp((float(j - start_j_2nd) / (div_num - start_j_2nd)), 0.0, 1.0);
				rate = fit(rate, 0.0, 1.0, 0.5*round_01, 0.5);
			//	printf("04 : [%d - %d]rate = %f j = %d\n", pt_1, pt_2, rate, j);
			}
			vector append_pos = lerp(pos_s, pos_e, rate);
			append(pts, append_pos);
		}

		// if prim is not looped, add last straight line control points
		if ((i==(num_corner-1)) && (is_prim_looped == 0))
		{
			int append_last_num = order; // number of control points
			for (int s=0; s<append_last_num; ++s)
			{
				float rate = fit(clamp(float(s)/(append_last_num-1), 0, 1), 0, 1, 0.5, 1.0);
				vector pos = lerp(pos_1, pos_2, rate);
				append(pts, pos);
//				printf("04 : [%d - %d]rate = %f\n", pt_1, pt_2, rate);
			}
		}
	}
}

/**
 *	Controlable Smooth Step
 */
function float dgfx_SmoothStep(const float x, edge, ofs)
{
	float edge_ = clamp01(edge - 0.5);
	float rate  = clamp01(edge * 2.0);
	// smooth step
	float edge0 = clamp01(edge_+ofs), edge1 = clamp01(1 - edge_+ofs);
	float smoothstep = smooth(edge0, edge1, x);
	return lerp(x, smoothstep, rate);
}
#endif // dgfx_vex_utils_h
