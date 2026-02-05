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
//////////////////////////////// DEFINES //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define RT_API                          static inline
#define RT_GAMMA                        2.2f
#define RT_GAMMA_INVERSE                0.454545f
#define RT_PI                           3.1415926f

///////////////////////////////////////////////////////////////////////////
//////////////////////////// MACRO FUNCTIONS //////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define RT_BUFFER_LEN(BUFFER)           (sizeof(BUFFER) / sizeof((BUFFER)[0]))
#define RT_TO_RADIANS(DEGREES)          (DEGREES * 0.0174533f)
#define RT_TO_DEGREES(RADIANS)          (RADIANS * 57.295779f)

///////////////////////////////////////////////////////////////////////////
/////////////////////////////// UTILITIES /////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
RT_API float rt_apply_inverse_gamma(float x)
{
    return powf(x, RT_GAMMA_INVERSE);
}

///////////////////////////////////////////////////////////////////////////
RT_API float rt_apply_gamma(float x)
{
    return powf(x, RT_GAMMA);
}

///////////////////////////////////////////////////////////////////////////
RT_API float rt_apply_gamma_custom(float x, float gamma)
{
    return powf(x, gamma);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// MATH ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    float x, y;
}
rt_vec2_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    float x, y, z;
}
rt_vec3_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_add(rt_vec3_t p, rt_vec3_t q)
{
    p.x += q.x;
    p.y += q.y;
    p.z += q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_add_scalar(rt_vec3_t p, float k)
{
    p.x += k;
    p.y += k;
    p.z += k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_sub(rt_vec3_t p, rt_vec3_t q)
{
    p.x -= q.x;
    p.y -= q.y;
    p.z -= q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_sub_scalar(rt_vec3_t p, float k)
{
    p.x -= k;
    p.y -= k;
    p.z -= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_mul(rt_vec3_t p, rt_vec3_t q)
{
    p.x *= q.x;
    p.y *= q.y;
    p.z *= q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_mul_scalar(rt_vec3_t p, float k)
{
    p.x *= k;
    p.y *= k;
    p.z *= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_div(rt_vec3_t p, rt_vec3_t q)
{
    p.x /= q.x;
    p.y /= q.y;
    p.z /= q.z;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_div_scalar(rt_vec3_t p, float k)
{
    p.x /= k;
    p.y /= k;
    p.z /= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_negate(rt_vec3_t p)
{
    p.x *= -1.0f;
    p.y *= -1.0f;
    p.z *= -1.0f;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API float rt_vec3_len(rt_vec3_t p)
{
    return sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
}

///////////////////////////////////////////////////////////////////////////
RT_API float rt_vec3_sqrlen(rt_vec3_t p)
{
    return p.x * p.x + p.y * p.y + p.z * p.z;
}

///////////////////////////////////////////////////////////////////////////
RT_API float rt_vec3_dot(rt_vec3_t p, rt_vec3_t q)
{
    return p.x * q.x + p.y * q.y + p.z * q.z;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_cross(rt_vec3_t p, rt_vec3_t q)
{
    rt_vec3_t r = {
        .x = p.y * q.z - p.z * q.y,
        .y = p.z * q.x - p.x * q.z,
        .z = p.x * q.y - p.y * q.x,
    };

    return r;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_norm(rt_vec3_t p)
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
RT_API uint32_t rt_vec3_to_uint32(rt_vec3_t p)
{
    uint32_t x = (uint32_t)(p.x * 255.99f);
    uint32_t y = (uint32_t)(p.y * 255.99f);
    uint32_t z = (uint32_t)(p.z * 255.99f);

    return (z << 16) | (y << 8) | (x << 0);
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_vec3_to_uint32_alpha(rt_vec3_t p, float a)
{
    uint32_t x = (uint32_t)(p.x * 255.99f);
    uint32_t y = (uint32_t)(p.y * 255.99f);
    uint32_t z = (uint32_t)(p.z * 255.99f);
    uint32_t w = (uint32_t)(  a * 255.99f);

    return (w << 24) | (z << 16) | (y << 8) | (x << 0);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_apply_0(rt_vec3_t p, float (*fun)(void))
{
    p.x = fun();
    p.y = fun();
    p.z = fun();

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_apply_1(rt_vec3_t p, float (*fun)(float))
{
    p.x = fun(p.x);
    p.y = fun(p.y);
    p.z = fun(p.z);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_apply_2(rt_vec3_t p, float (*fun)(float, float), float k)
{
    p.x = fun(p.x, k);
    p.y = fun(p.y, k);
    p.z = fun(p.z, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_lerp(rt_vec3_t p, rt_vec3_t q, float k)
{
    p.x += k * (q.x - p.x);
    p.y += k * (q.y - p.y);
    p.z += k * (q.z - p.z);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_rand(uint32_t seed)
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
RT_API rt_vec3_t rt_vec3_max(rt_vec3_t p, float k)
{
    p.x = fmaxf(p.x, k);
    p.y = fmaxf(p.y, k);
    p.z = fmaxf(p.z, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_min(rt_vec3_t p, float k)
{
    p.x = fminf(p.x, k);
    p.y = fminf(p.y, k);
    p.z = fminf(p.z, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_clamp(rt_vec3_t p, float min, float max)
{
    p.x = fminf(max, fmaxf(p.x, min));
    p.y = fminf(max, fmaxf(p.y, min));
    p.z = fminf(max, fmaxf(p.z, min));

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_reflect(rt_vec3_t p, rt_vec3_t n)
{
    n       = rt_vec3_norm(n);
    float d = 2.0f * rt_vec3_dot(p, n);

    p.x -= n.x * d;
    p.y -= n.y * d;
    p.z -= n.z * d;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_refract(rt_vec3_t p, rt_vec3_t n, float k)
{
    rt_vec3_t pnegated      = rt_vec3_negate(p);
    float cos_theta         = fminf(1.0f, rt_vec3_dot(pnegated, n));
    rt_vec3_t perpendicular = rt_vec3_mul_scalar(rt_vec3_add(p, rt_vec3_mul_scalar(n, cos_theta)), k);
    float a                 = -sqrtf(fabsf(1.0f - rt_vec3_sqrlen(perpendicular)));
    rt_vec3_t parallel      = rt_vec3_mul_scalar(n, a);
    return rt_vec3_add(perpendicular, parallel);
}

///////////////////////////////////////////////////////////////////////////
RT_API float rt_vec3_dist(rt_vec3_t p, rt_vec3_t q)
{
    float x = p.x - q.x;
    float y = p.y - q.y;
    float z = p.z - q.z;

    return sqrtf(x * x + y * y + z * z);
}

///////////////////////////////////////////////////////////////////////////
RT_API float rt_vec3_sqrdist(rt_vec3_t p, rt_vec3_t q)
{
    float x = p.x - q.x;
    float y = p.y - q.y;
    float z = p.z - q.z;

    return x * x + y * y + z * z;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_rotate_x(rt_vec3_t p, float r)
{
    float py = p.y;
    float pz = p.z;

    float cos_r = cosf(r);
    float sin_r = sinf(r);

    p.y = py * cos_r - pz * sin_r;
    p.z = py * sin_r + pz * cos_r;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_vec3_rotate_y(rt_vec3_t p, float r)
{
    float px = p.x;
    float pz = p.z;

    float cos_r = cosf(r);
    float sin_r = sinf(r);

    p.x =  px * cos_r + pz * sin_r;
    p.z = -px * sin_r + pz * cos_r;

    return p;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// RAYS ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t org;
    rt_vec3_t dir;
}
rt_ray_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t position;
    rt_vec3_t normal;
    float t;
    bool is_front_facing;
}
rt_hit_info;

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_ray_at(rt_ray_t ray, float t)
{
    rt_vec3_t dir_scaled = rt_vec3_mul_scalar(ray.dir, t);
    return rt_vec3_add(ray.org, dir_scaled);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////// MATERIALS ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t color;
}
rt_emissive_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t color_0;
    rt_vec3_t color_1;
    float shadow_factor;
}
rt_checkerboard_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t ambient;
    rt_vec3_t diffuse;
}
rt_diffuse_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t ambient;
    rt_vec3_t specular;
}
rt_metallic_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t ambient;
    rt_vec3_t specular;
    float refractive_index;
}
rt_dielectric_material_t;

///////////////////////////////////////////////////////////////////////////
//////////////////////////// GEOMETRY /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t center;
    float radius;
}
rt_sphere_t;

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_sphere_hit(rt_sphere_t sphere,
                          rt_ray_t ray,
                          rt_vec2_t range,
                          rt_hit_info* info)
{
    rt_vec3_t oc    = rt_vec3_sub(sphere.center, ray.org);
    float a         = rt_vec3_sqrlen(ray.dir);
    float h         = rt_vec3_dot(ray.dir, oc);
    float c         = rt_vec3_sqrlen(oc) - sphere.radius * sphere.radius;

    float disc      = h * h - a * c;
    if (disc < 0.0f) {
        return NULL;
    }

    float sqrt_disc     = sqrtf(disc);
    float root_nearest  = (h - sqrt_disc) / a;
    float root_farthest = (h + sqrt_disc) / a;

    float range_min     = range.x;
    float range_max     = range.y;

    float root_chosen   = root_nearest;
    if (root_chosen <= range_max || root_chosen >= range_max) {
        root_chosen = root_farthest;
        if (root_chosen <= range_min || root_chosen >= range_max) {
            return false;
        }
    }

    rt_vec3_t n             = rt_vec3_sub(info->position, sphere.center);
    n                       = rt_vec3_div_scalar(n, sphere.radius);

    info->t                 = root_chosen;
    info->position          = rt_ray_at(ray, root_chosen);
    info->normal            = rt_vec3_norm(n);
    info->is_front_facing   = rt_vec3_dot(info->normal, ray.dir) < 0.0f;

    return true;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////// LIGHTS //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t color;
    rt_vec3_t position;
    float intensity;
}
rt_point_light_t;

///////////////////////////////////////////////////////////////////////////
//////////////////////// FRAGMENT SHADERS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_fragment_shader_emissive(rt_emissive_material_t* material)
{
    return material->color;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_fragment_shader_diffuse(const rt_hit_info* hit_info,
                                            const rt_diffuse_material_t* material,
                                            const rt_point_light_t* point_lights,
                                            uint32_t point_light_count,
                                            bool is_in_shadow)
{
    if (is_in_shadow) {
        return material->ambient;
    }

    rt_vec3_t final_color = {};
    for (uint32_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light   = &point_lights[i];

        rt_vec3_t light_color           = rt_vec3_mul_scalar(light->color, light->intensity);
        rt_vec3_t light_dir             = rt_vec3_norm(rt_vec3_sub(light->position, hit_info->position));

        rt_vec3_t diffuse_color         = rt_vec3_mul(material->diffuse, light_color);
        float diff                      = fmaxf(rt_vec3_dot(light_dir, hit_info->normal), 0.0f);
        rt_vec3_t diffuse               = rt_vec3_mul_scalar(diffuse_color, diff);

        final_color                     = rt_vec3_add(final_color, diffuse);
    }

    final_color = rt_vec3_add(final_color, material->ambient);
    return rt_vec3_clamp(final_color, 0.0f, 1.0f);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_fragment_shader_checkerboard(const rt_hit_info* hit_info,
                                                 const rt_checkerboard_material_t* material,
                                                 const rt_point_light_t* point_lights,
                                                 uint32_t point_light_count,
                                                 bool is_in_shadow)
{
    uint32_t pos_x_quant    = (uint32_t)(ceilf(hit_info->position.x * 0.5f));
    uint32_t pos_z_quant    = (uint32_t)(ceilf(hit_info->position.z * 0.5f));
    uint32_t pos_quant      = pos_x_quant + pos_z_quant;

    rt_vec3_t material_color = (pos_quant & 1) ? material->color_0
                                               : material->color_1;

    if (is_in_shadow) {
        return rt_vec3_mul_scalar(material_color, material->shadow_factor);
    }

    rt_vec3_t final_color = {};
    for (uint32_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light   = &point_lights[i];

        rt_vec3_t light_color           = rt_vec3_mul_scalar(light->color, light->intensity);
        rt_vec3_t light_dir             = rt_vec3_norm(rt_vec3_sub(light->position, hit_info->position));

        rt_vec3_t diffuse_color         = rt_vec3_mul(material_color, light_color);
        float diff                      = fmaxf(rt_vec3_dot(light_dir, hit_info->normal), 0.0f);
        rt_vec3_t diffuse               = rt_vec3_mul_scalar(diffuse_color, diff);

        final_color                     = rt_vec3_add(final_color, diffuse);
    }

    return rt_vec3_clamp(final_color, material->shadow_factor, 1.0f);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// WORLD ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////
