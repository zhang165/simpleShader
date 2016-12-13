#include "object.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>
#include <cfloat>
#include <iostream>

// 2016 Version

bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // Assert the correct values of the W coords of the origin and direction.
    // You can comment this out if you take great care to construct the rays
    // properly.
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    Ray localRay(i_transform * ray.origin, i_transform * ray.direction);
	//!!! USEFUL NOTES: to calculate depth in localIntersect(), if the intersection happens at
	//ray.origin + ray.direction * t, then t is just the depth
	//!!! USEFUL NOTES: Here direction might be scaled, so you must not renormalize it in
	//localIntersect(), or you will get a depth in local coordinate system,
	//which can't be compared with intersections with other objects
    if (localIntersect(localRay, hit))	{
        // Assert correct values of W.
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
		// Transform intersection coordinates into global coordinates.
        hit.position = transform * hit.position;
        hit.normal = (n_transform * hit.normal).normalized();
        
		return true;
    }

    return false;
}

//INTERSECTIONS 

// helper function from Scratchapixel
bool solveQuadratic(const double &a, const double &b, const double &c, double &x0, double &x1){
	double discr = b * b - 4 * a * c; // b^2 - 4*ac
	if (discr < 0) return false;
	else if (discr == 0) x0 = x1 = -0.5 * b / a;
	else {
		double q = (b > 0) ? -0.5 * (b + sqrt(discr)) :	-0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}
	if (x0 > x1) std::swap(x0, x1);
	return true;
}

bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const {
	//////////////////
    // YOUR CODE HERE 
	// For this part you are required to calculate if a ray has intersected your sphere.
	// Be sure to cover all the possible intersection scenarios(zero, one, and two points of intersection).
	// Test your result by comparing the output depth image of your algorithm with the provided example solution's results.

	// Here in local coordinate system, the sphere is centered on (0, 0, 0)
	//with radius 1.0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
	double radius = this->radius;
	radius = radius * radius;
	Vector p0 = Vector(0,0,0); 

	double t0, t1;

	Vector L = ray.origin - p0;
	double a = ray.direction.dot(ray.direction);
	double b = 2 * ray.direction.dot(L);
	double c = L.dot(L) - radius;
	if (!solveQuadratic(a, b, c, t0, t1)) return false;

    if (t0 > t1) std::swap(t0, t1); // make t0 the smallest
         if (t0 < 0) { 
            t0 = t1; // if t0 is negative, let's use t1 instead 
            if (t0 < 0) return false; // both t0 and t1 are negative 
        } 
 	
 	if(hit.depth < t0) return false; // take the closest depth

 	hit.depth = t0;
 	hit.position = ray.origin + ray.direction * hit.depth;
 	hit.normal = hit.position;
	return true;
}


bool Plane::localIntersect(Ray const &ray, Intersection &hit) const{
	//////////////////
	// YOUR CODE HERE 
	// The implementation of this part is similar to the previous part in that you are 
	// calculating if a line has intersected your plane.Test this function in a similar 
	// way to the previous part(note that as you do new objects will appear).

	// Do not accept intersection while ray is inside the plane
	// Here in local coordinate system, the plane is at z = 0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
	Vector n = Vector(0, 0, 1);
	Vector p0 = Vector(0, 0, 0);

	double ln = n.dot(ray.direction);
	if (abs(ln) <= 1e-6) return false; // parallel

	double op = (ray.origin - p0).dot(n);
	if (abs(op) <= 1e-6) return false; // inside

	double t = -(op) / ln;

	if ((t < 0 ) || (hit.depth < t)) return false; // take the closest valid ray

	hit.depth = t;
	hit.normal = n;
	hit.position = ray.origin + ray.direction * t;
	return true;
}

bool Mesh::intersectTriangle(Ray const &ray, Triangle const &tri, Intersection &hit) const
{
	// Extract vertex positions from the mesh data.
	Vector const &p0 = positions[tri[0].pi];
	Vector const &p1 = positions[tri[1].pi];
	Vector const &p2 = positions[tri[2].pi];

	//////////////////
	// YOUR CODE HERE 
	// Think back through the course and try to decide what equations might help you decide on which 
	// side of the bounding lines of the triangle the ray intersects.
	// Test this part just like the above two parts; when triangle intersection is working properly, 
	// you should be able to see full meshes appear in your scenes.
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
	//!!! USEFUL NOTES: for the intersection point, its normal should satisfy hit.normal.dot(ray.direction) < 0
	Vector v0v1 = p1 - p0;
	Vector v0v2 = p2 - p0;

	// PLANE TEST
	Vector n = v0v1.cross(v0v2);
	
	double ln = n.dot(ray.direction);
	if (abs(ln) <= 1e-6) return false; // parallel

	double op = (ray.origin - p0).dot(n);
	if (abs(op) <= 1e-6) return false; // inside

	double t = -(op) / ln;

	Vector P = ray.origin + t * ray.direction; // the intersection point

	// INSIDE TEST
	Vector C; // vector perpendicular to triangle's plane 
	Vector edge0 = p1 - p0;
	Vector vp0 = P - p0;
	C = edge0.cross(vp0);
	if (n.dot(C) < 0) return false; // P is on the right side 

	Vector edge1 = p2 - p1;
	Vector vp1 = P - p1;
	C = edge1.cross(vp1);
	if (n.dot(C) < 0)  return false; // P is on the right side 

	Vector edge2 = p0 - p2;
	Vector vp2 = P - p2;
	C = edge2.cross(vp2);
	if (n.dot(C) < 0) return false; // P is on the right side; 

	if (hit.depth < t || t < 0) return false;

	hit.depth = t;
	hit.normal = n;
	hit.position = P;

	return true;
}

bool Conic::localIntersect(Ray const &ray, Intersection &hit) const {
	//////////////////
	// YOUR CODE HERE (creative license)
    return false;
}


// Intersections!
bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
	// Bounding box check
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++) {
		if (ray.direction[i] == 0.0) {
			if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) {
				// Ray parallel to bounding box plane and outside of box!
				return false;
			}
			// Ray parallel to bounding box plane and inside box: continue;
		}
		else {
			double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Ensure t1 <= t2

			if (t1 > tNear) tNear = t1; // We want the furthest tNear
			if (t2 < tFar) tFar = t2; // We want the closest tFar

			if (tNear > tFar) return false; // Ray misses the bounding box.
			if (tFar < 0) return false; // Bounding box is behind the ray.
		}
	}
	// If we made it this far, the ray does intersect the bounding box.

	// The ray hits the bounding box, so check each triangle.
	bool isHit = false;
	for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
		Triangle const &tri = triangles[tri_i];

		if (intersectTriangle(ray, tri, hit)) {
			isHit = true;
		}
	}
	return isHit;
}

double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}


