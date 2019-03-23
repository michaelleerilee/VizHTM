/*
 * Edges.C
 *
 *  Created on: Feb 28, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "tests.h"

bool Edges1( VizHTM *viz ) {

	SpatialVector origin(0,0,0);
	SpatialVector xhat(1,0,0);
	SpatialVector yhat(0,1,0);
	SpatialVector zhat(0,0,1);

	float64 scale = 1.0;
	float64 a     = -1;

	if(true) {

		viz->addEdge(origin, xhat, 1.0, 0.0, 0.0, a, scale);
		viz->addEdge(origin, yhat, 0.0, 1.0, 0.0, a, scale);
		viz->addEdge(origin, zhat, 0.0, 0.0, 1.0, a, scale);

		viz->addSphere(origin,1.0,1.0,1.0,0.125);
		viz->addSphere(xhat,1.0,0,0,0.125);
		viz->addSphere(yhat,0,1.0,0,0.125);
		viz->addSphere(zhat,0,0,1.0,0.125);

	}


	SpatialVector x0(0.1,0.1,0.1);
	SpatialVector x1(0.1,-0.1,0.1);
	for( int i = 0; i < 10; ++i ) {
		viz->addEdge(x0,x1,1.0,0.0,0.0,a,scale); scale += 0.1;
		viz->addEdge(x0,x1,0.0,1.0,0.0,a,scale); scale += 0.1;
		viz->addEdge(x0,x1,0.0,0.0,1.0,a,scale); scale += 0.1;
	}

	SpatialVector *xplus = new SpatialVector(xhat*1.2);
	viz->addAnnotation(xplus, "x\n", 0.25, 1.0, 0.0, 0.0);

	SpatialVector *zplus = new SpatialVector(zhat*1.2);
	viz->addAnnotation(zplus, "z\n", 0.25, 0.0, 0.0, 1.0);

	return true;
}


// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Edges2Barycentric0(
		SpatialVector p, SpatialVector a, SpatialVector b, SpatialVector c,
		float64 &u, float64 &v, float64 &w, bool verbose)
{
	SpatialVector v0 = b - a, v1 = c - a, v2 = p - a;
    float64 d00 = v0*v0;
    float64 d01 = v0*v1;
    float64 d11 = v1*v1;
    float64 d20 = v2*v0;
    float64 d21 = v2*v1;
    float64 idenom = 1.0/(d00 * d11 - d01 * d01);
    v = (d11 * d20 - d01 * d21) * idenom;
    w = (d00 * d21 - d01 * d20) * idenom;
    u = 1.0 - v - w;

    if(verbose) {
    	cout
		<< setprecision(17)
		<< setw(20)
		<< scientific
    	<< endl;
    	cout << " - " <<  endl << flush;
//    	print_vector("a",a);
//    	print_vector("b",b);
//    	print_vector("c",c);
//    	print_vector("p",p);
    	cout << " - " <<  endl << flush;
    	cout << " - " <<  endl << flush;
    	cout << "a " << a << endl << flush;
    	cout << "b " << b << endl << flush;
    	cout << "c " << c << endl << flush;
    	cout << endl << flush;
    	cout << "v0 = " << v0 << ", l = " << v0.length() << endl << flush;
    	cout << "v1 = " << v1 << ", l = " << v1.length() << endl << flush;
    	cout << "v2 = " << v2 << ", l = " << v2.length() << endl << flush;
    	cout << endl << flush;
    	cout << "d00    " << d00 << endl << flush;
    	cout << "d01    " << d01 << endl << flush;
    	cout << "d11    " << d11 << endl << flush;
    	cout << "d20    " << d20 << endl << flush;
    	cout << "d21    " << d21 << endl << flush;
    	cout << "idenom " << idenom << endl << flush;
    	cout << "v      " << v << endl << flush;
    	cout << "w      " << w << endl << flush;
    	cout << "u      " << u << endl << flush;
    	cout << endl << flush;
    }
}

void Edges2Barycentric1(
		SpatialVector p, SpatialVector a, SpatialVector b, SpatialVector c,
		float64 &u, float64 &v, float64 &w, bool verbose)
{
	cout << "200" << endl << flush;
	SpatialVector normal = (b-a)^(c-a); normal.normalize();
	cout << "201" << endl << flush;
	float64 area_abc = normal * ((b-a)^(c-a));
	cout << "202" << endl << flush;
	float64 area_pbc = normal * ((b-p)^(c-p));
	cout << "203" << endl << flush;
	float64 area_pca = normal * ((c-p)^(a-p));
	cout << "204" << endl << flush;
	u = area_pbc/area_abc;
	v = area_pca/area_abc;
	w = 1 - u - v;
}

bool
isInside0(const SpatialVector & v, const SpatialVector & v0,
	       const SpatialVector & v1, const SpatialVector & v2)
{
	float64 gEpsilon = 1.0e-16;
	float64 d01 = (v0 ^ v1) * v;
	float64 d12 = (v1 ^ v2) * v;
	float64 d20 = (v2 ^ v0) * v;
	cout << "ii d01 = " << d01 << endl << flush;
	cout << "ii d12 = " << d12 << endl << flush;
	cout << "ii d20 = " << d20 << endl << flush;
	if( (v0 ^ v1) * v < -gEpsilon) return false;
	if( (v1 ^ v2) * v < -gEpsilon) return false;
	if( (v2 ^ v0) * v < -gEpsilon) return false;
	return true;
}

float64 triple_product0(const SpatialVector &a, const SpatialVector &b, const SpatialVector &c)
{
	return (a ^ b) * c;
}
float64 min_triangle_quality0(
		const SpatialVector & v, const SpatialVector & v0,
		const SpatialVector & v1, const SpatialVector & v2)
{
	float64 d01 = triple_product0(v0,v1,v);
	float64 d12 = triple_product0(v1,v2,v);
	float64 d20 = triple_product0(v2,v0,v);
	return min(d01,min(d12,d20));
}

bool Edges2( VizHTM * viz ) {
	float64 scale          = 1.0;
	float64 a_transparency = 0.0;

    SpatialVector a( -1.81695612430945508e-01, -8.18712639798629960e-01, 5.44698373283143189e-01 );
    SpatialVector b( -1.90494010835855154e-01, -8.26148716846798692e-01, 5.30273825007474109e-01 );
    SpatialVector c( -1.72486271871561264e-01, -8.29116447465894857e-01, 5.31802973437892623e-01 );
    SpatialVector q = a + b + c; q = q * (1.0/3.0);
    SpatialVector p( -1.87801794532830940e-01, -8.23897295103509841e-01, 5.34718368013824996e-01 );
    SpatialVector dab = a - b;

	// SpatialVector normal = (b-a)^(c-a); normal.normalize();
    SpatialVector normal = q; normal.normalize();
	SpatialVector p1     =  p - a;
	SpatialVector dproj  = normal * (normal * p1);
	SpatialVector proj   = p - dproj;

    cout
	<< setprecision(17)
	<< setw(20)
	<< scientific;

    cout << "a: " << a << endl << flush;
    cout << "b: " << b << endl << flush;
    cout << "c: " << c << endl << flush;
    cout << "p: " << p << endl << flush;

    // cout << "ab" << endl << flush;
	viz->addEdge(
			a,b,
			1.0,0.0,0.0,
			a_transparency,scale);

	// cout << "bc" << endl << flush;
	viz->addEdge(
			b,c,
			0.0,1.0,0.0,
			a_transparency,scale);

	// cout << "ca" << endl << flush;
	viz->addEdge(
			c,a,
			0.0,0.0,1.0,
			a_transparency,scale);

	// cout << "radius" << endl << flush;
	float64 radius = dab.length()*0.0001;
	// cout << "sphere" << endl << flush;
	// viz->addSphere(p,0.75,1.0,0.75,radius);

	float64 eps = 1.0e-3;
	viz->addEdge(
			p*(1-eps),p*(1+eps),
			0.75,1.0,0.75,
			a_transparency,scale
			);
	viz->addEdge(
			p,q,
			1.0,1.0,1.0,
			a_transparency,scale);
	viz->addEdge(
			proj,q,
			1.0,1.0,0.0,
			a_transparency,scale);
	// viz->addSphere(zhat,0,0,1.0,0.125);

	// cout << "TriangleCheck1 done" << endl << flush;

	float64 u,v,w; bool verbose = true;
	cout << 100 << endl << flush;
	Edges2Barycentric0( p,  a,  b,  c, u, v, w, verbose );
	cout << "u,v,w = " << u << "," << v << "," << w << endl << flush;

	cout << 110 << endl << flush;
	Edges2Barycentric0( p, c, a, b, u, v, w, verbose );

	cout << 120 << endl << flush;
	Edges2Barycentric1( p,  a,  b,  c, u, v, w, verbose );
	cout << "u,v,w = " << u << "," << v << "," << w << endl << flush;

	cout << 130 << endl << flush;
	Edges2Barycentric1( proj,  a,  b,  c, u, v, w, verbose );
	cout << "u,v,w = " << u << "," << v << "," << w << endl << flush;

	cout << 140 << endl << flush;
	Edges2Barycentric0( proj,  a,  b,  c, u, v, w, verbose );
	cout << "u,v,w = " << u << "," << v << "," << w << endl << flush;

	cout << 150 << endl << flush;
	bool inside = isInside0(p,a,b,c);
	cout << "isInside:    " << inside  << endl << flush;
	cout << "min quality: " << min_triangle_quality0(p,a,b,c) << endl << flush;

	return true;
}
