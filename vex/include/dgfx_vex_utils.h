/**
 * Copyright (c) 2019 d-gfx
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */
#ifndef dgfx_vex_utils_h
#define dgfx_vex_utils_h

#define mod(a, b)			((a) % (b))
#define mix(a, b, t)		lerp((a), (b), (t))
#define clamp01(v)			clamp(v, 0.0, 1.0)
#define vector3_ctor(v, f)	(set((v).x, (v).y, f))
#define vector4_ctor(v, f)	(set((v).x, (v).y, (v).z, f))
#define swap(type, a, b)			{ type __tmp = a; a = b; b = __tmp; }

// repeat(4, 0, 4) => 0
#define repeat(value, min, max) (((max-min) <= 0) ? min : mod(value-min, max-min) + min)
#define log_base(value, base)	(log(value)/log(base))
#define log2(value)	(log_base(value, 2.0))
#define exp2(value) (pow(2.0, value))
#define factorial(n)	(dgfx_Calc_Factorial(n))

/**
 *	calc factorial
 *	ex: 3! = 3*2*1
 */
function int dgfx_Calc_Factorial(const int n)
{
    int ret = 1;
    for (int i=2; i<=n; ++i) { ret *= i; }
    return ret;
}

/**
 *	Array Intersection
 */
#define ArrayIntersection(out_array, array_a, array_b)	\
{\
	for (int i=0; i<len(array_a); ++i)\
	{\
		if (0 <= find(array_b, array_a[i]))\
			append(out_array, array_a[i]);\
	}\
}

function int[] dgfx_Array_Intersection(const int a[], b[]) { int ret[]; ArrayIntersection(ret, a, b); return ret; }
function float[] dgfx_Array_Intersection(const float a[], b[]) { float ret[]; ArrayIntersection(ret, a, b); return ret; }
function vector[] dgfx_Array_Intersection(const vector a[], b[]) { vector ret[]; ArrayIntersection(ret, a, b); return ret; }

/**
 *	Shuffle Array (in place)
 */
#define ArrayShuffle(type, inout_array, seed_01)	\
{\
	int size = len(inout_array);\
	for (int i=size-1; 1<=i; --i)\
	{\
		float rnd = rand(seed_01);\
		int j = int(rint(rnd * i));\
		swap(type, inout_array[i], inout_array[j]);\
	}\
}

function void dgfx_Array_Shuffle(int arr[]; float seed_01) { ArrayShuffle(int, arr, seed_01); }
function void dgfx_Array_Shuffle(float arr[]; float seed_01) { ArrayShuffle(float, arr, seed_01); }
function void dgfx_Array_Shuffle(vector arr[]; float seed_01) { ArrayShuffle(vector, arr, seed_01); }

/**
 *	Calc Vector Angle (radian)
 */
function float dgfx_Calc_Vector_Radian(const vector vec1, vec2)
{
	vector dir1 = normalize(vec1);
	vector dir2 = normalize(vec2);
	return qdistance(dihedral(dir1, dir2), {0,0,0,1});
}

/**
 *	Vector Projection
 */
function float dgfx_Vector_Extract_Vertical(vector tgt; const vector dir)
{
	vector nrm = normalize(dir);
	float d = dot(tgt, nrm);
	tgt = tgt - d*nrm;
	return d/length(dir);
}
function float dgfx_Vector_Extract_Parallel(vector tgt; const vector dir)
{
	vector nrm = normalize(dir);
	float d = dot(tgt, nrm);
	tgt = d*nrm;
	return d/length(dir);
}

function float dgfx_Calc_Perpendicular_Foot(vector ret; const vector pos, edge_begin, edge_end)
{
	vector edge_vec = edge_end - edge_begin;
	vector begin_to_pos = pos - edge_begin;
	float d = dgfx_Vector_Extract_Parallel(begin_to_pos, edge_vec);
	ret = edge_begin + begin_to_pos;
	return d;
}

/**
 *	Calc orient (Quaternion) from Two Vectors
 */
function vector4 dgfx_Calc_Orient_Front_Side(const vector front, side)
{
	vector _up = normalize(cross(front, side));
	vector _side = normalize(cross(_up, front));
	matrix3 mtx = set(_side, _up, front);
	return normalize(quaternion(mtx));
}

/**
 *	Create Axis Geometry (for Visualize)
 */
function void dgfx_Create_Axis_Geomety(const float axis_len)
{
	int pt_orig = addpoint(0, {0,0,0});
	int pt_x = addpoint(0, set(axis_len,0,0));
	int pt_y = addpoint(0, set(0,axis_len,0));
	int pt_z = addpoint(0, set(0,0,axis_len));
	int prim_x = addprim(0, "polyline", pt_orig, pt_x);
	int prim_y = addprim(0, "polyline", pt_orig, pt_y);
	int prim_z = addprim(0, "polyline", pt_orig, pt_z);
	setprimattrib(0, "Cd", prim_x, {1,0,0});
	setprimattrib(0, "Cd", prim_y, {0,1,0});
	setprimattrib(0, "Cd", prim_z, {0,0,1});
}

/**
 *	Normalized Value to Color
 */
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
 *	Calc Circle Center Position
 */
function void dgfx_Calc_Circumscribed_Circle(float ret_r; vector ret_center, ret_plane_nrm; const vector a, b, c)
{
	vector ac = c - a;
	vector ab = b - a;
	vector ab_x_ac = cross(ab, ac);
	vector a_center = (cross(ab_x_ac, ab) * length2(ac) + cross(ac, ab_x_ac) * length2(ab)) / (2.0 * length2(ab_x_ac));

	float radius = length(a_center);
	vector center = a + a_center;

	ret_r = radius;
	ret_center = center;
	ret_plane_nrm = normalize(ab_x_ac);
}

/**
 *	is Primitive Vertex Looped ?
 */
function int dgfx_IsLooped_PrimVertex(int geometry, prim_num)
{
//	int is_closed = primintrinsic(geometry, "closed", prim_num); // this is whether open or closed prim
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

/**
 *	Count Edges (not half edges)
 *	from https://www.sidefx.com/docs/houdini/vex/functions/pointhedgenext.html
 */
function int dgfx_CountEdges(int geo, pt)
{
    int edge_count = 0;
    int hout = pointhedge(geo, pt);
    while ( hout != -1 )
    {
        if (hedge_isprimary(geo, hout))	{ edge_count++; }
        int hin = hedge_prev(geo, hout);
        if (hedge_isprimary(geo, hin))	{ edge_count++; }
        hout = pointhedgenext(geo, hout);
    }
    return edge_count;
}

/**
 *	Make 3D Shear Matrix
 *	ref https://www.sidefx.com/docs/houdini/vex/functions/maketransform.html
 */
function matrix dgfx_MakeShearMatrix(const vector shear)
{
	vector zero = set(0,0,0);
	vector ones = set(1,1,1);
	vector t = zero;
	vector r = zero;
	vector s = ones;
	vector p = zero;
	vector pr = zero;
	matrix shear_mtx = maketransform(XFORM_SRT, XFORM_XYZ, t, r, s, p, pr, shear);
	return shear_mtx;
}

/**
 *	Calc Color Palette
 *	from http://iquilezles.org/www/articles/palettes/palettes.htm
 */
function vector dgfx_ColorPalette(float t; vector a, b, c, d)
{
	vector color = a + b * cos(6.28318 * (c * t + d));
	return color;
}

/**
 *	Rainbow Color Palette
 */
function vector dgfx_ColorPaletteRainbow(float t)
{
	vector a = set(0.5, 0.5, 0.5);
	vector b = set(0.5, 0.5, 0.5);
	vector c = set(1, 1, 1);
	vector d = set(0.0, 0.33, 0.67);
	return dgfx_ColorPalette(t, a, b, c, d);
}

/**
 *	Break Matrix Component
 */
function void dgfx_BreakMatrix(vector A_col_0, A_col_1, A_col_2; const matrix3 A)
{
	assign(A_col_0.x, A_col_0.y, A_col_0.z
		 , A_col_1.x, A_col_1.y, A_col_1.z
		 , A_col_2.x, A_col_2.y, A_col_2.z
		 , A);
}

/**
 *	Project vector to plane
 */
function vector dgfx_Proj_to_Plane(const vector src, plane_nrm)
{
    vector delete_dir_nrm = normalize(plane_nrm);
    float elem_dot = dot(src, delete_dir_nrm);
    return src - elem_dot*delete_dir_nrm;
}

/**
 *	OpenGL Projection Matrix
 */
function matrix dgfx_Calc_Projection_Mtx_GL(const float fovy, near, far, h_div_w)
{
	float e = 1.0/tan(fovy/2.0);
	float a = h_div_w;
	float inv_fn = 1.0/(far - near);
	matrix proj_mtx = set(e,   0,				 0,					  0,
						  0, e/a,				 0,					  0,
						  0,   0,-(far+near)*inv_fn, -2*far*near*inv_fn,
						  0,   0,				-1,					  0);
	proj_mtx = transpose(proj_mtx);
	return proj_mtx;
}

/**
 *	Get Local Box Points
 */
function void dgfx_Calc_Local_Box_Points(vector points[])
{
	append(points, set(-1, -1, -1));
	append(points, set( 1, -1, -1));
	append(points, set( 1,  1, -1));
	append(points, set(-1,  1, -1));
	append(points, set(-1, -1,  1));
	append(points, set( 1, -1,  1));
	append(points, set( 1,  1,  1));
	append(points, set(-1,  1,  1));
}

/**
 *	Apply Transform using Perspective Divide
 *	World <--> NDC
 */
function vector dgfx_Transform_Projection(const vector P; const matrix A, B)
{
    vector4 vec4_P = vector4_ctor(P, 1);
    vector4 proj_P = vec4_P * A * B;
    float w = proj_P.w;
    return set(proj_P.x/w, proj_P.y/w, proj_P.z/w);
}
#endif // dgfx_vex_utils_h
