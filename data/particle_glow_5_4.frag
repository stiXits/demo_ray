#version 150

uniform sampler2D color;
uniform sampler2D glow;

in vec2 v_uv;

out vec4 fragColor;

void main()
{
	vec4 c = texelFetch(color, ivec2(gl_FragCoord.xy), 0);
	vec4 g = texture(glow, v_uv/40, 0);

	// Task_5_4 - ToDo Begin

//	vec4 color = mix(c, g, vec4(0.5f));

	fragColor = vec4(g.xyz, 1.0);

	// Task_5_4 - ToDo End
}
