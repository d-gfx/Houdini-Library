/**
 * Copyright (c) 2019 d-gfx
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Dual Number
 */
struct dual
{
	float x, dx;
}

dual dvar(float x)   { return dual(x, 1); } // variable
dual dconst(float x) { return dual(x, 0); } // constant
float dreal(const dual c) { return c.x;  }
float deps(const dual c)  { return c.dx; }

/**
 * addition
 */
dual dadd(const dual a; const dual b)
{
	dual ret;
	ret.x  = a.x + b.x;
	ret.dx = a.dx + b.dx;
	return ret;
}

/**
 * subtraction
 */
dual dsub(const dual a; const dual b)
{
	dual ret;
	ret.x  = a.x - b.x;
	ret.dx = a.dx - b.dx;
	return ret;
}

/**
 * multiplication
 */
dual dmul(const dual a; const dual b)
{
	dual ret;
	ret.x  = a.x * b.x;
	ret.dx = (a.dx * b.x + a.x * b.dx);
	return ret;
}

/**
 * division
 */
dual ddiv(const dual a; const dual b)
{
	dual ret;
	ret.x  = a.x / b.x;
	ret.dx = (a.dx * b.x - a.x * b.dx) / (b.x*b.x);
	return ret;
}

/**
 * sin
 */
dual dsin(const dual a)
{
	dual ret;
	ret.x = sin(a.x);
	ret.dx = a.dx*cos(a.x);
	return ret;
}

/**
 * cos
 */
dual dcos(const dual a)
{
	dual ret;
	ret.x = cos(a.x);
	ret.dx = -a.dx*sin(a.x);
	return ret;
}

/**
 * exponential
 */
dual dexp(const dual a)
{
	dual ret;
	ret.x = exp(a.x);
	ret.dx = a.dx * exp(a.x);
	return ret;
}

/**
 * logarithm
 */
dual dlog(const dual a)
{
	dual ret;
	ret.x = log(a.x);
	ret.dx = a.dx / a.x;
	return ret;
}

/**
 * power
 */
dual dpow(const dual a; float k)
{
	dual ret;
	ret.x  = pow(a.x, k);
	ret.dx = k * pow(a.x, k-1) * a.dx;
	return ret;
}

/**
 * absolute
 */
dual dabs(const dual a)
{
	dual ret;
	ret.x = abs(a.x);
	ret.dx = a.dx * sign(a.x);
	return ret;
}
