// rayterm: ray-tracer for the terminal
// Copyright (C) 2026 0xdeadf1.sh
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

///////////////////////////////////////////////////////////////////////////
typedef struct {
    float x, y, z;
} rt_vec3_t;

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_create(float x, float y, float z)
{
    rt_vec3_t p = { x, y, z };
    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_add(rt_vec3_t p, rt_vec3_t q)
{
    p.x += q.x;
    p.y += q.y;
    p.z += q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_add_scalar(rt_vec3_t p, float k)
{
    p.x += k;
    p.y += k;
    p.z += k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_sub(rt_vec3_t p, rt_vec3_t q)
{
    p.x -= q.x;
    p.y -= q.y;
    p.z -= q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_sub_scalar(rt_vec3_t p, float k)
{
    p.x -= k;
    p.y -= k;
    p.z -= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_mul(rt_vec3_t p, rt_vec3_t q)
{
    p.x *= q.x;
    p.y *= q.y;
    p.z *= q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_mul_scalar(rt_vec3_t p, float k)
{
    p.x *= k;
    p.y *= k;
    p.z *= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_div(rt_vec3_t p, rt_vec3_t q)
{
    p.x /= q.x;
    p.y /= q.y;
    p.z /= q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_div_scalar(rt_vec3_t p, float k)
{
    p.x /= k;
    p.y /= k;
    p.z /= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline float rt_vec3_len(rt_vec3_t p)
{
    return sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
}

///////////////////////////////////////////////////////////////////////////
static inline float rt_vec3_sqrlen(rt_vec3_t p)
{
    return p.x * p.x + p.y * p.y + p.z * p.z;
}

///////////////////////////////////////////////////////////////////////////
static inline float rt_vec3_dot(rt_vec3_t p, rt_vec3_t q)
{
    return p.x * q.x + p.y * q.y + p.z * q.z;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_cross(rt_vec3_t p, rt_vec3_t q)
{
    rt_vec3_t r = {
        .x = p.y * q.z - p.z * q.y,
        .y = p.z * q.x - p.x * q.z,
        .z = p.x * q.y - p.y * q.x,
    };

    return r;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_norm(rt_vec3_t p)
{
    float k = p.x * p.x + p.y * p.y + p.z * p.z;
    if (k > 0.0f) {
        k = 1.0f / sqrtf(k);
        p.x *= k;
        p.y *= k;
        p.z *= k;
    }
    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline uint32_t rt_vec3_to_uint32(rt_vec3_t p)
{
    uint32_t x = (uint32_t)(p.x * 255.99f);
    uint32_t y = (uint32_t)(p.y * 255.99f);
    uint32_t z = (uint32_t)(p.z * 255.99f);

    return (z << 16) | (y << 8) | (x << 0);
}

///////////////////////////////////////////////////////////////////////////
static inline uint32_t rt_vec3_to_uint32_alpha(rt_vec3_t p, float a)
{
    uint32_t x = (uint32_t)(p.x * 255.99f);
    uint32_t y = (uint32_t)(p.y * 255.99f);
    uint32_t z = (uint32_t)(p.z * 255.99f);
    uint32_t w = (uint32_t)(  a * 255.99f);

    return (w << 24) | (z << 16) | (y << 8) | (x << 0);
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_apply_0(rt_vec3_t p, float (*fun)(void))
{
    p.x = fun();
    p.y = fun();
    p.z = fun();

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_apply_1(rt_vec3_t p, float (*fun)(float))
{
    p.x = fun(p.x);
    p.y = fun(p.y);
    p.z = fun(p.z);

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_apply_2(rt_vec3_t p, float (*fun)(float, float), float k)
{
    p.x = fun(p.x, k);
    p.y = fun(p.y, k);
    p.z = fun(p.z, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_lerp(rt_vec3_t p, rt_vec3_t q, float k)
{
    p.x += k * (q.x - p.x);
    p.y += k * (q.y - p.y);
    p.z += k * (q.z - p.z);

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_rand(uint32_t seed)
{
    srand(seed);
    rt_vec3_t p = {
        .x = (float)(rand() % RAND_MAX),
        .y = (float)(rand() % RAND_MAX),
        .z = (float)(rand() % RAND_MAX),
    };

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_max(rt_vec3_t p, float k)
{
    p.x = fmaxf(p.x, k);
    p.y = fmaxf(p.y, k);
    p.z = fmaxf(p.z, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_min(rt_vec3_t p, float k)
{
    p.x = fminf(p.x, k);
    p.y = fminf(p.y, k);
    p.z = fminf(p.z, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_clamp(rt_vec3_t p, float min, float max)
{
    p.x = fminf(max, fmaxf(p.x, min));
    p.y = fminf(max, fmaxf(p.y, min));
    p.z = fminf(max, fmaxf(p.z, min));

    return p;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_vec3_reflect(rt_vec3_t p, rt_vec3_t n)
{
    float d = 2.0f * (p.x * n.x + p.y * n.y + p.z * n.z);
    p.x -= n.x * d;
    p.y -= n.y * d;
    p.z -= n.z * d;

    return p;
}
