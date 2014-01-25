#version 150

const float INFINITY = 1e+4;

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

	vec3 rayNormal;
	float distanceToNormal;


    {
        vec3 r = ray.direction;
        vec3 s = ray.origin;
        vec3 c = blob.position;

        distanceToNormal = - (r.x * (s.x - c.x) + r.y * (s.y - c.y) + r.z * (s.z - c.z))/(pow(2, r.x) + pow(2, r.x) + pow(2, r.z));
    }

	rayNormal = blob.position - ray.origin - distanceToNormal * normalize(ray.direction);

	if(length(rayNormal) < blob.radius)
    {
        t0 = length((blob.position - rayNormal) - (normalize(ray.direction) * sqrt(pow(blob.radius, 2.0) - pow(length(rayNormal), 2.0))));
        t1 = length((blob.position - rayNormal) + (normalize(ray.direction) * sqrt(pow(blob.radius, 2.0) - pow(length(rayNormal), 2.0))));
        return true;
    }
	// Task_5_1 - ToDo End
	return false;
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

		if(intersect(blobs[i], ray, t0, t1) /* todo, more? */)
		{
		    if(t0 < t)
            {
                t = t0;
    			t /= length(ray.direction);
                normal = normalize(blobs[i].position - (ray.origin + t * ray.direction));
                return true;
            }
			// Task 5_2: material = ?;
		}
	}
	normal = vec3(0.0f);
	return false; // ?

	// ToDo: End Task 5_1
}


// Task_5_3 - ToDo Begin

// ... your helper functions

// ... more ...

bool trace(in Ray ray, out vec3 normal, out Material material, out float t)
{
	// Task_5_3 - ToDo Begin

	// find nearest and farthest intersection for all metaballs
	// hint: use for loop, INFINITE, SIZE, intersect, min and max...

	// ...

	// implement raymarching within your tmin and tmax
	// hint: e.g., use while loop, THRESHOLD, and implment yourself
	// an attribute interpolation function (e.g., interp(pos, normal, material, actives?))
	// as well as a summation function (e.g., sum(pos))

	// your shader should terminate!

	// return true if iso surface was hit, fals if not

	return false;
}

// Task_5_3 - ToDo End
