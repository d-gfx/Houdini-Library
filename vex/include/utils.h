/**
 * Copyright (c) 2019 d-gfx
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

function vector calcHeatMapColor(float value_01)
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
function vector bezierQuadratic_Position(float t_01; vector p0, p1, p2)
{
	float t = clamp(t_01, 0.0, 1.0);
	return (1-t)*(1-t) * p0 + 2*(1-t)*t * p1 + t*t * p2;
}

function vector bezierQuadratic_Derivative_1st(float t_01; vector p0, p1, p2)
{
	float t = clamp(t_01, 0.0, 1.0);
	return 2*(1-t) * (p1 - p0) + 2*t*(p2 - p1);
}

function vector bezierQuadratic_Derivative_2nd(float t_01; vector p0, p1, p2)
{
	float t = clamp(t_01, 0.0, 1.0);
	return 2 * (p2 - 2*p1 + p0);
}

/**
 *	Cubic Bezier Curve
 */
function vector bezierCubic_Position(float t_01; vector p0, p1, p2, p3)
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

function vector bezierCubic_Derivative_1st(float t_01; vector p0, p1, p2, p3)
{
	float t = clamp(t_01, 0.0, 1.0);
	float rev_t = clamp(1.0 - t, 0.0, 1.0);
	return 3*rev_t*rev_t * (p1 - p0) + 6*rev_t*t * (p2 - p1) + 3*t*t * (p3-p2);
}

function vector bezierCubic_Derivative_2nd(float t_01; vector p0, p1, p2, p3)
{
	float t = clamp(t_01, 0.0, 1.0);
	float rev_t = clamp(1.0 - t, 0.0, 1.0);
	return 6*rev_t * (p2 - 2*p1 + p0) + 6*t * (p3 - 2*p2 + p1);
}
