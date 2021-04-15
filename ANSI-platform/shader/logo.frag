#version 330 core

in vec2 fTexcoord;

layout (location = 0) out vec4 color;

uniform sampler2D Sampler;

void main()
{
    vec3 texture = vec3(texture(Sampler, fTexcoord));
    // color = vec4(1.0f, 1.0f, 1.0f, 1.0f);
    color = vec4(texture, 0.2f);
}
