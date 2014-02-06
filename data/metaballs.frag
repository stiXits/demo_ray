#version 150

const float INFINITY = 1e+4 ;

const int SIZE = 16;
const float THRESHOLD = 0.66;

struct Sphere
{
	vec3 position;
	float radius;
};

struct Material
{
	vec4 sr; // vec3 specular, float reflectance
	vec4 dr; // vec3 diffuse, float roughness
};

struct Ray
{
	vec3 origin;
	vec3 direction;
};


Sphere blobs[16];
Material materials[16];

void cache(
	sampler2D positions
,	sampler2D materials0
,	sampler2D materials1)
{
	for(int i = 0; i < SIZE; ++i)
	{
		ivec2 uv = ivec2(i % 4, i / 4);

		vec4 v = texelFetch(positions, uv, 0);

		Sphere blob;
		blob.position = v.xyz;
		blob.radius = v.w;

		blobs[i] = blob;

		Material mat;
		mat.sr = texelFetch(materials0, uv, 0);
		mat.dr = texelFetch(materials1, uv, 0);

		materials[i] = mat;
	}
}

// sphere intersection
bool intersect(
    const in Sphere blob
,   const in Ray ray // needs to be normalized!
,   out float t0
,   out float t1)
{
	// Task_5_1 - ToDo Begin
	// implement a ray sphere intersection

	// hint: the case of single intersection with sphere can be neglected, so
	// only t > 0.0 intersections are relevant (if you like you can also
	// account for t == 0.0 intersections, which should have a rare occurrence).

	// hint: use ray.origin, ray.direction, blob.position, blob.radius
	// hint: make it fast (if you like)! ;)

    vec3 normRay = ray.origin - blob.position;

    float a = dot(ray.direction, ray.direction);
    float b = 2 * dot(ray.direction, normRay);
    float c = dot(normRay, normRay) - (blob.radius * blob.radius);

    float disc = b * b - 4 * a * c;
    if (disc < 0)
        return false;


    float distSqrt = sqrt(disc);

    t0 = (-b - distSqrt)/(2.0*a);
    t1 = (-b + distSqrt)/(2.0*a);

    if (t0 > t1)
    {
        float temp = t0;
        t0 = t1;
        t1 = temp;
    }

	return true;
}

bool rcast(in Ray ray, out vec3 normal, out Material material, out float t)
{
    // Task_5_1 - ToDo Begin

	// return normal for the nearest intersected sphere, as well as
	// the intersection parameter t (that should allow to retrieve the
	// itnersection point I = ray.origin + ray.direction * t).

	// the function should return true, if at least one sphere was hit,
	// and false if no when no sphere was hit.

	// (Task_5_2 - ToDo return material of nearest intersected sphere)
	t =  INFINITY;

	for(int i = 0; i < SIZE; ++i)
	{
		float t0; // = ?
		float t1; // = ?

		// ...

		if(intersect(blobs[i], ray, t0, t1) && t0 < t)
		{
            t = t0;
            normal = normalize(ray.origin + ray.direction * t0 - blobs[i].position);
			material = materials[i];
		}
	}
    if(t != INFINITY)
        return true;
    else
    {
        normal = vec3(0.0f);
        return false;
    }
	// ToDo: End Task 5_1
}


// Task_5_3 - ToDo Begin

// ... your helper functions

// ... more ...
bool nearestIntersection(in Ray ray, out float t)
{
    t = INFINITY;
    float t0;
    float t1;

	for(int i = 0; i < SIZE; ++i)
	{
        if(intersect(blobs[i], ray, t0, t1) && t0 < t)
        {
            t = t0;
        }
	}

	return (t < INFINITY);
}

bool farthestIntersection(in Ray ray, out float t)
{
    t = 0;
    float t0;
    float t1;

	for(int i = 0; i < SIZE; ++i)
	{
        if(intersect(blobs[i], ray, t0, t1) && t0 > t)
        {
            t = t0;
        }
	}

	return (t > 0);
}

bool getDistance(in vec3 position, in int sphere, out float t)
{
    t = sqrt(    pow(blobs[sphere].position.x - position.x, 2) +
                    pow(blobs[sphere].position.y - position.y, 2) +
                    pow(blobs[sphere].position.z - position.z, 2)) - blobs[sphere].radius;
    return true;
}

bool nearestDistance(in vec3 position, out float t)
{
    t = INFINITY;
    float t0;

	for(int i = 0; i < SIZE; ++i)
	{
        if(getDistance(position, i, t0) && t0 < t)
        {
            t = t0;
        }
	}

	return (t < INFINITY);
}

bool hitSphere(in vec3 position, out int sphere)
{
    float t;

    for(int i = 0; i < SIZE; ++i)
	{
        if(getDistance(position, i, t) && (t <= 0))
        {
            sphere = i;
            return true;
        }
	}
    return false;
}

bool trace(in Ray ray, out vec3 normal, out Material material, out float t)
{
	// Task_5_3 - ToDo Begin

	// find nearest and farthest intersection for all metaballs
	// hint: use for loop, INFINITE, SIZE, intersect, min and max...



	float tmin;
	float tmax;

	if(!farthestIntersection(ray, tmax) || !nearestIntersection(ray, tmin))
    {
        return false;
    }

	// implement raymarching within your tmin and tmax
	// hint: e.g., use while loop, THRESHOLD, and implment yourself
	// an attribute interpolation function (e.g., interp(pos, normal, material, actives?))
	// as well as a summation function (e.g., sum(pos))

	// your shader should terminate!

	// return true if iso surface was hit, fals if not
    float marchedDistance = tmin;
    float currentStep;
    bool hit = false;
    int sphere;
    Ray marchingRay = ray;

	while(marchedDistance <= tmax && !hit)
    {
        nearestDistance(marchingRay.origin, currentStep);
        marchingRay.origin += currentStep * marchingRay.direction;
        marchedDistance += currentStep;
        hit = hitSphere(marchingRay.origin, sphere);
    }

    material = materials[sphere];
    normal = normalize(marchingRay.origin - blobs[sphere].position);
    t = marchedDistance;

	return hit;
}

// Task_5_3 - ToDo End
