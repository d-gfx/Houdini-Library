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

function void calcBezierPointQuadratic(vector interp_pos, derivative_1st, derivative_2nd; float t; vector p0, p1, p2)
{
	float t_ = clamp(t, 0.0, 1.0);
	interp_pos = (1-t_)*(1-t_) * p0 + 2*(1-t_)*t_ * p1 + t_*t_ * p2;
	derivative_1st = 2*(1-t) * (p1 - p0) + 2*t*(p2 - p1);
	derivative_2nd = 2 * (p2 - 2*p1 + p0);
}
