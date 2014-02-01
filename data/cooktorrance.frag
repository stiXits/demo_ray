#version 150

const float EPSILON  = 1e-6;
const float c_e = 2.7182;

struct Material
{
	vec4 sr; // vec3 specular, float reflectance
	vec4 dr; // vec3 diffuse, float roughness
};

// Schlick's Approximation of the Fresnel Factor
float fresnel(in float VdotH, in float r)
{
	// r: reflectance
	// Task_5_2 - ToDo Begin

	return r + (1.0f - r) * pow((1 - VdotH), 5.0f);

	// Task_5_2 - ToDo End
}

// Beckmann's distribution for roughness
float roughness(in float NdotH, in float r)
{
	// r: roughness
	// Task_5_2 - ToDo Begin

	float exponent = (pow(NdotH, 2) - 1)/(pow(r, 2) * pow(NdotH, 2));
	float base = 1/(pow(r, 2) * pow(NdotH, 4));

	return base*exp(exponent);

	// Task_5_2 - ToDo End
}

// Geometric attenuation accounts for the shadowing and
// self masking of one microfacet by another.
float geom(in float NdotH, in float NdotV, in float VdotH, in float NdotL)
{
	// Task_5_2 - ToDo Begin

    float term1 = (2 * NdotH * NdotV)/VdotH;
    float term2 = (2 * NdotH * NdotL)/VdotH;


    return min(min(1, term1), term2);

	// Task_5_2 - ToDo End
}

vec3 CookTorrance(in vec3 V, in vec3 N, in vec3 L, in Material m, in vec3 R, in vec3 ambient)
{
	vec3 H = normalize(L + V);

	float VdotH = clamp(dot(V, H), 0.0, 1.0);
	float NdotV = clamp(dot(N, V), 0.0, 1.0);
	float NdotH = clamp(dot(N, H), 0.0, 1.0);
	float NdotL = clamp(dot(N, L), 0.0, 1.0);

	// Task_5_2 - ToDo Begin

	// hint: R is reflection (e.g., ray in envmap)

    float Rs = (fresnel(VdotH, m.sr.a) * roughness(NdotH, m.dr.a)  * geom(NdotH, NdotV, VdotH, NdotL))/(NdotV * NdotL);
    vec3 final = max(0, NdotL) * (m.sr.rgb * Rs  + m.dr.rgb * ambient);
    final += R * max(0, NdotL);
	return final;

	// Task_5_2 - ToDo End
}
