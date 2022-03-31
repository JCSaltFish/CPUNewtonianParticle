#version 410

layout (location=0) out vec4 oColor;

uniform int particlesMode;

void main()
{
	if (particlesMode == 0)
		oColor = vec4(0.3);
	else
		oColor = vec4(1.0);
}
