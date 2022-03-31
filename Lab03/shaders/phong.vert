#version 410

layout (location=0) in vec4 iPosition;
layout (location=1) in vec3 iNormal;

uniform mat4 view;
uniform mat4 proj;

void main()
{
	gl_Position = proj*view*iPosition;  //standard vertex out
}
