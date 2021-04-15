#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexcoord;

out vec2 fTexcoord;

uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;

void main()
{
	fTexcoord = vec2(aTexcoord.x, 1 - aTexcoord.y);
	gl_Position = projection * view * model * vec4(aPos, 1.0);
}