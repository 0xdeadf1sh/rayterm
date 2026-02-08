// rayterm: ray-tracer for the terminal tritten in modern C (version 23)
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

///////////////////////////////////////////////////////////////////////////
//////////////////////////////// RAYTERM.H ////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifndef RT_RAYTERM_H
#define RT_RAYTERM_H

#include <assert.h>
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
#define RT_INIT_CAP                     8
#define RT_SHADOW_BIAS                  0.001f

///////////////////////////////////////////////////////////////////////////
//////////////////////////// MACRO FUNCTIONS //////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define RT_BUFFER_LEN(BUFFER)           (sizeof(BUFFER) / sizeof((BUFFER)[0]))
#define RT_TO_RADIANS(DEGREES)          (DEGREES * 0.0174533f)
#define RT_TO_DEGREES(RADIANS)          (RADIANS * 57.295779f)
#define RT_ASSERT(EXPR)                 assert(EXPR && "RT_ASSERT: " #EXPR " failed")

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
///////////////////////////// MATERIALS ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
enum rt_material_type
{
    RT_MATERIAL_NONE,
    RT_MATERIAL_EMISSIVE,
    RT_MATERIAL_CHECKERBOARD,
    RT_MATERIAL_DIFFUSE,
    RT_MATERIAL_METALLIC,
    RT_MATERIAL_DIELECTRIC
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   color;
}
rt_emissive_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   color_0;
    rt_vec3_t   color_1;
    float       shadow_factor;
    bool        receives_shadows;
}
rt_checkerboard_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   ambient;
    rt_vec3_t   diffuse;
    bool        receives_shadows;
}
rt_diffuse_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   ambient;
    rt_vec3_t   specular;
    bool        receives_shadows;
}
rt_metallic_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   ambient;
    rt_vec3_t   specular;
    float       refractive_index;
    bool        receives_shadows;
}
rt_dielectric_material_t;

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
enum rt_hit_geometry_type
{
    RT_HIT_NONE,
    RT_HIT_SPHERE,
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   position;
    rt_vec3_t   normal;
    float       t;
    bool        is_front_facing;
}
rt_hit_info_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_ray_at(rt_ray_t ray, float t)
{
    rt_vec3_t dir_scaled = rt_vec3_mul_scalar(ray.dir, t);
    return rt_vec3_add(ray.org, dir_scaled);
}

///////////////////////////////////////////////////////////////////////////
//////////////////////////// GEOMETRY /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
enum rt_geometry_type
{
    RT_GEOMETRY_NONE,
    RT_GEOMETRY_SPHERE
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t               center;
    float                   radius;

    uint32_t                material_index;
    enum rt_material_type   material_type;
}
rt_sphere_t;

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_sphere_hit(rt_sphere_t sphere,
                          rt_ray_t ray,
                          rt_vec2_t range,
                          rt_hit_info_t* info)
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

    info->position          = rt_ray_at(ray, root_chosen);
    info->normal            = rt_vec3_norm(n);
    info->t                 = root_chosen;
    info->is_front_facing   = rt_vec3_dot(info->normal, ray.dir) < 0.0f;

    return true;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////// LIGHTS //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
enum rt_light_type
{
    RT_LIGHT_NONE,
    RT_LIGHT_DIRECTIONAL,
    RT_LIGHT_POINT,
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   color;
    rt_vec3_t   direction;
    float       intensity;
    bool        casts_shadows;
}
rt_directional_light_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec3_t   color;
    rt_vec3_t   position;
    float       intensity;
    bool        casts_shadows;
}
rt_point_light_t;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// WORLD ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
enum rt_face_cull_mode
{
    RT_FACE_CULL_NONE,
    RT_FACE_CULL_FRONT   = 1 << 0,
    RT_FACE_CULL_BACK    = 1 << 1,
    RT_FACE_CULL_BOTH    = RT_FACE_CULL_FRONT | RT_FACE_CULL_BACK,
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    ///////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// SKY ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    rt_vec3_t clear_color;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// LIGHTS //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    rt_directional_light_t* directional_lights;
    uint32_t directional_light_count;
    uint32_t directional_light_capacity;

    ///////////////////////////////////////////////////////////////////////////
    rt_point_light_t* point_lights;
    uint32_t point_light_count;
    uint32_t point_light_capacity;

    ///////////////////////////////////////////////////////////////////////////
    //////////////////////////////// GEOMETRY /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    rt_sphere_t* spheres;
    uint32_t sphere_count;
    uint32_t sphere_capacity;

    ////////////////////////////////////////////////////////////////////RT_INIT_CAP///////
    /////////////////////////////// MATERIALS /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    rt_emissive_material_t* emissive_materials;
    uint32_t emissive_material_count;
    uint32_t emissive_material_capacity;

    ///////////////////////////////////////////////////////////////////////////
    rt_checkerboard_material_t* checkerboard_materials;
    uint32_t checkerboard_material_count;
    uint32_t checkerboard_material_capacity;

    ///////////////////////////////////////////////////////////////////////////
    rt_diffuse_material_t* diffuse_materials;
    uint32_t diffuse_material_count;
    uint32_t diffuse_material_capacity;

    ///////////////////////////////////////////////////////////////////////////
    rt_metallic_material_t* metallic_materials;
    uint32_t metallic_material_count;
    uint32_t metallic_material_capacity;

    ///////////////////////////////////////////////////////////////////////////
    rt_dielectric_material_t* dielectric_materials;
    uint32_t dielectric_material_count;
    uint32_t dielectric_material_capacity;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// STATE ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    enum rt_face_cull_mode face_cull_mode;
}
rt_world_t;

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_PUSH(WORLD_NAME,                                           \
                      BUFFER_TYPE,                                          \
                      BUFFER_NAME,                                          \
                      COUNT_NAME,                                           \
                      CAPACITY_NAME)                                        \
{                                                                           \
    uint32_t count     = WORLD_NAME->COUNT_NAME;                            \
    uint32_t capacity  = WORLD_NAME->CAPACITY_NAME;                         \
                                                                            \
    if (count == UINT32_MAX - 1) {                                          \
        return UINT32_MAX;                                                  \
    }                                                                       \
    else if (count >= capacity) {                                           \
        uint32_t capacity_new   = capacity ? 2 * capacity : RT_INIT_CAP;    \
                                                                            \
        uint32_t bytes_required = capacity_new * sizeof(BUFFER_TYPE);       \
        BUFFER_TYPE* buffer_new = (BUFFER_TYPE*)malloc(bytes_required);     \
                                                                            \
        if (!buffer_new) {                                                  \
            return UINT32_MAX;                                              \
        }                                                                   \
                                                                            \
        for (uint32_t i = 0; i < count; ++i) {                              \
            buffer_new[i] = WORLD_NAME->BUFFER_NAME[i];                     \
        }                                                                   \
                                                                            \
        free(WORLD_NAME->BUFFER_NAME);                                      \
                                                                            \
        WORLD_NAME->BUFFER_NAME     = buffer_new;                           \
        WORLD_NAME->CAPACITY_NAME   = capacity_new;                         \
    }                                                                       \
                                                                            \
    WORLD_NAME->BUFFER_NAME[count]  = (BUFFER_TYPE){};                      \
    WORLD_NAME->COUNT_NAME          += 1;                                   \
    return count;                                                           \
}

///////////////////////////////////////////////////////////////////////////
#define RT_ASSERT_PUSH(COUNT) RT_ASSERT(COUNT != UINT32_MAX)

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_POP(WORLD_NAME,                                            \
                     BUFFER_TYPE,                                           \
                     BUFFER_NAME,                                           \
                     COUNT_NAME)                                            \
{                                                                           \
    uint32_t count          = WORLD_NAME->COUNT_NAME;                       \
                                                                            \
    if (!count) {                                                           \
        return NULL;                                                        \
    }                                                                       \
                                                                            \
    uint32_t last_index     = count - 1;                                    \
                                                                            \
    BUFFER_TYPE* popped     = &WORLD_NAME->BUFFER_NAME[last_index];         \
                                                                            \
    WORLD_NAME->COUNT_NAME -= 1;                                            \
    return popped;                                                          \
}

///////////////////////////////////////////////////////////////////////////
#define RT_ASSERT_POP(POPPED) RT_ASSERT(POPPED)

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_RESERVE(WORLD_NAME,                                        \
                         BUFFER_TYPE,                                       \
                         BUFFER_NAME,                                       \
                         COUNT_NAME,                                        \
                         CAPACITY_NAME,                                     \
                         NEW_CAPACITY)                                      \
{                                                                           \
    uint32_t count          = WORLD_NAME->COUNT_NAME;                       \
    uint32_t capacity_new   = NEW_CAPACITY;                                 \
                                                                            \
    uint32_t bytes_required = capacity_new * sizeof(BUFFER_TYPE);           \
    BUFFER_TYPE* buffer_new = (BUFFER_TYPE*)malloc(bytes_required);         \
                                                                            \
    if (!buffer_new) {                                                      \
        return false;                                                       \
    }                                                                       \
                                                                            \
    for (uint32_t i = 0; i < count; ++i) {                                  \
        buffer_new[i] = WORLD_NAME->BUFFER_NAME[i];                         \
    }                                                                       \
                                                                            \
    free(WORLD_NAME->BUFFER_NAME);                                          \
                                                                            \
    WORLD_NAME->BUFFER_NAME     = buffer_new;                               \
    WORLD_NAME->CAPACITY_NAME   = capacity_new;                             \
    return true;                                                            \
}

///////////////////////////////////////////////////////////////////////////
#define RT_ASSERT_RESERVE(RESERVED) RT_ASSERT(RESERVED)

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_FREE(WORLD_NAME,                                           \
                      BUFFER_NAME,                                          \
                      COUNT_NAME,                                           \
                      CAPACITY_NAME)                                        \
{                                                                           \
    free(WORLD_NAME->BUFFER_NAME);                                          \
    WORLD_NAME->BUFFER_NAME     = NULL;                                     \
    WORLD_NAME->COUNT_NAME      = 0;                                        \
    WORLD_NAME->CAPACITY_NAME   = 0;                                        \
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_directional_light(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_directional_light_t,
                  directional_lights,
                  directional_light_count,
                  directional_light_capacity)
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_point_light(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_point_light_t,
                  point_lights,
                  point_light_count,
                  point_light_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_sphere(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_sphere_t,
                  spheres,
                  sphere_count,
                  sphere_capacity)
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_emissive_material(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_emissive_material_t,
                  emissive_materials,
                  emissive_material_count,
                  emissive_material_capacity)
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_checkerboard_material(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_checkerboard_material_t,
                  checkerboard_materials,
                  checkerboard_material_count,
                  checkerboard_material_capacity)
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_diffuse_material(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_diffuse_material_t,
                  diffuse_materials,
                  diffuse_material_count,
                  diffuse_material_capacity)
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_metallic_material(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_metallic_material_t,
                  metallic_materials,
                  metallic_material_count,
                  metallic_material_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_push_dielectric_material(rt_world_t* world)
{
    RT_WORLD_PUSH(world,
                  rt_dielectric_material_t,
                  dielectric_materials,
                  dielectric_material_count,
                  dielectric_material_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_directional_light_t* rt_world_pop_directional_light(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_directional_light_t,
                 directional_lights,
                 directional_light_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_point_light_t* rt_world_pop_point_light(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_point_light_t,
                 point_lights,
                 point_light_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_sphere_t* rt_world_pop_sphere(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_sphere_t,
                 spheres,
                 sphere_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_emissive_material_t* rt_world_pop_emissive_material(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_emissive_material_t,
                 emissive_materials,
                 emissive_material_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_checkerboard_material_t* rt_world_pop_checkerboard_material(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_checkerboard_material_t,
                 checkerboard_materials,
                 checkerboard_material_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_diffuse_material_t* rt_world_pop_diffuse_material(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_diffuse_material_t,
                 diffuse_materials,
                 diffuse_material_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_metallic_material_t* rt_world_pop_metallic_material(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_metallic_material_t,
                 metallic_materials,
                 metallic_material_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_dielectric_material_t* rt_world_pop_dielectric_material(rt_world_t* world)
{
    RT_WORLD_POP(world,
                 rt_dielectric_material_t,
                 dielectric_materials,
                 dielectric_material_count);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_directional_lights(rt_world_t* world,
                                                uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_directional_light_t,
                     directional_lights,
                     directional_light_count,
                     directional_light_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_point_lights(rt_world_t* world,
                                          uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_point_light_t,
                     point_lights,
                     point_light_count,
                     point_light_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_spheres(rt_world_t* world,
                                     uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_sphere_t,
                     spheres,
                     sphere_count,
                     sphere_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_emissive_materials(rt_world_t* world,
                                                uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_emissive_material_t,
                     emissive_materials,
                     emissive_material_count,
                     emissive_material_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_checkerboard_materials(rt_world_t* world,
                                                    uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_checkerboard_material_t,
                     checkerboard_materials,
                     checkerboard_material_count,
                     checkerboard_material_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_diffuse_materials(rt_world_t* world,
                                               uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_diffuse_material_t,
                     diffuse_materials,
                     diffuse_material_count,
                     diffuse_material_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_metallic_materials(rt_world_t* world,
                                                uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_metallic_material_t,
                     metallic_materials,
                     metallic_material_count,
                     metallic_material_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_reserve_dielectric_materials(rt_world_t* world,
                                                  uint32_t capacity)
{
    RT_WORLD_RESERVE(world,
                     rt_dielectric_material_t,
                     dielectric_materials,
                     dielectric_material_count,
                     dielectric_material_capacity,
                     capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_directional_lights(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  directional_lights,
                  directional_light_count,
                  directional_light_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_point_lights(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  point_lights,
                  point_light_count,
                  point_light_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_spheres(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  spheres,
                  sphere_count,
                  sphere_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_emissive_materials(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  emissive_materials,
                  emissive_material_count,
                  emissive_material_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_checkerboard_materials(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  checkerboard_materials,
                  checkerboard_material_count,
                  checkerboard_material_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_diffuse_materials(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  diffuse_materials,
                  diffuse_material_count,
                  diffuse_material_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_metallic_materials(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  metallic_materials,
                  metallic_material_count,
                  metallic_material_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free_dielectric_materials(rt_world_t* world)
{
    RT_WORLD_FREE(world,
                  dielectric_materials,
                  dielectric_material_count,
                  dielectric_material_capacity);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free(rt_world_t* world)
{
    rt_world_free_directional_lights            (world);
    rt_world_free_point_lights                  (world);
    rt_world_free_spheres                       (world);
    rt_world_free_emissive_materials            (world);
    rt_world_free_checkerboard_materials        (world);
    rt_world_free_diffuse_materials             (world);
    rt_world_free_metallic_materials            (world);
    rt_world_free_dielectric_materials          (world);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_sphere_unset_material(rt_world_t* world,
                                     uint32_t sphere_index)
{
    RT_ASSERT(sphere_index      < world->sphere_count);

    rt_sphere_t* sphere         = &world->spheres[sphere_index];
    sphere->material_type       = RT_MATERIAL_NONE;
    sphere->material_index      = UINT32_MAX;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_sphere_set_emissive_material(rt_world_t* world,
                                            uint32_t sphere_index,
                                            uint32_t material_index)
{
    RT_ASSERT(sphere_index      < world->sphere_count);
    RT_ASSERT(material_index    < world->emissive_material_count);

    rt_sphere_t* sphere         = &world->spheres[sphere_index];
    sphere->material_type       = RT_MATERIAL_EMISSIVE;
    sphere->material_index      = material_index;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_sphere_set_checkerboard_material(rt_world_t* world,
                                                uint32_t sphere_index,
                                                uint32_t material_index)
{
    RT_ASSERT(sphere_index      < world->sphere_count);
    RT_ASSERT(material_index    < world->checkerboard_material_count);

    rt_sphere_t* sphere         = &world->spheres[sphere_index];
    sphere->material_type       = RT_MATERIAL_CHECKERBOARD;
    sphere->material_index      = material_index;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_sphere_set_diffuse_material(rt_world_t* world,
                                           uint32_t sphere_index,
                                           uint32_t material_index)
{
    RT_ASSERT(sphere_index      < world->sphere_count);
    RT_ASSERT(material_index    < world->diffuse_material_count);

    rt_sphere_t* sphere         = &world->spheres[sphere_index];
    sphere->material_type       = RT_MATERIAL_DIFFUSE;
    sphere->material_index      = material_index;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_sphere_set_metallic_material(rt_world_t* world,
                                            uint32_t sphere_index,
                                            uint32_t material_index)
{
    RT_ASSERT(sphere_index      < world->sphere_count);
    RT_ASSERT(material_index    < world->metallic_material_count);

    rt_sphere_t* sphere         = &world->spheres[sphere_index];
    sphere->material_type       = RT_MATERIAL_METALLIC;
    sphere->material_index      = material_index;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_sphere_set_dielectric_material(rt_world_t* world,
                                              uint32_t sphere_index,
                                              uint32_t material_index)
{
    RT_ASSERT(sphere_index      < world->sphere_count);
    RT_ASSERT(material_index    < world->dielectric_material_count);

    rt_sphere_t* sphere         = &world->spheres[sphere_index];
    sphere->material_type       = RT_MATERIAL_DIELECTRIC;
    sphere->material_index      = material_index;
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_world_sphere_closest_hit(const rt_world_t* world,
                                            rt_ray_t ray,
                                            rt_vec2_t range,
                                            rt_hit_info_t* info)
{
    if (RT_FACE_CULL_BOTH == world->face_cull_mode) {
        return UINT32_MAX;
    }

    rt_hit_info_t closest_hit_info    = {};
    uint32_t closest_hit_index      = UINT32_MAX;

    uint32_t sphere_count           = world->sphere_count;
    const rt_sphere_t* spheres      = world->spheres;

    enum rt_face_cull_mode cull_mode = world->face_cull_mode;

    for (uint32_t i = 0; i < sphere_count; ++i) {

        rt_hit_info_t hit_info = {};

        if (rt_sphere_hit(spheres[i], ray, range, &hit_info)) {

            bool is_front_facing = hit_info.is_front_facing;

            if ((RT_FACE_CULL_FRONT == cull_mode &&  is_front_facing) ||
                (RT_FACE_CULL_BACK  == cull_mode && !is_front_facing)) {
                continue;
            }

            if (UINT32_MAX == closest_hit_index ||
                hit_info.t < closest_hit_info.t) {

                closest_hit_info    = hit_info;
                closest_hit_index   = i;
            }
        }
    }

    if (UINT32_MAX != closest_hit_index) {
        *info = closest_hit_info;
        return closest_hit_index;
    }

    return UINT32_MAX;
}

///////////////////////////////////////////////////////////////////////////
#define RT_HAS_HIT(HIT_INDEX) (UINT32_MAX != (HIT_INDEX))

///////////////////////////////////////////////////////////////////////////
RT_API enum rt_hit_geometry_type rt_world_any_closest_hit(const rt_world_t* world,
                                                          rt_ray_t ray,
                                                          rt_vec2_t range,
                                                          rt_hit_info_t* info,
                                                          uint32_t* geometry_hit_index)
{
    uint32_t hit_index = rt_world_sphere_closest_hit(world, ray, range, info);
    if (RT_HAS_HIT(hit_index)) {
        *geometry_hit_index = hit_index;
        return RT_HIT_SPHERE;
    }

    return RT_HIT_NONE;
}

///////////////////////////////////////////////////////////////////////////
//////////////////////// FRAGMENT SHADERS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_fragment_shader_emissive(rt_emissive_material_t* material)
{
    return material->color;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_fragment_shader_diffuse(const rt_hit_info_t* hit_info,
                                            const rt_diffuse_material_t* material,
                                            const rt_world_t* world,
                                            bool is_in_shadow)
{
    if (is_in_shadow && material->receives_shadows) {
        return material->ambient;
    }

    rt_vec3_t final_color = {};

    uint32_t directional_light_count = world->directional_light_count;
    const rt_directional_light_t* directional_lights = world->directional_lights;

    for (uint32_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];

        rt_vec3_t light_color           = rt_vec3_mul_scalar(light->color,
                                                             light->intensity);
        rt_vec3_t light_dir             = rt_vec3_norm(rt_vec3_negate(light->direction));

        rt_vec3_t diffuse_color         = rt_vec3_mul(material->diffuse,
                                                      light_color);
        float diff                      = fmaxf(rt_vec3_dot(light_dir,
                                                            hit_info->normal), 0.0f);
        rt_vec3_t diffuse               = rt_vec3_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec3_add(final_color,
                                                      diffuse);
    }

    uint32_t point_light_count = world->point_light_count;
    const rt_point_light_t* point_lights = world->point_lights;

    for (uint32_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light   = &point_lights[i];

        rt_vec3_t light_color           = rt_vec3_mul_scalar(light->color,
                                                             light->intensity);
        rt_vec3_t light_dir             = rt_vec3_norm(rt_vec3_sub(light->position,
                                                                   hit_info->position));

        rt_vec3_t diffuse_color         = rt_vec3_mul(material->diffuse,
                                                      light_color);
        float diff                      = fmaxf(rt_vec3_dot(light_dir,
                                                            hit_info->normal), 0.0f);
        rt_vec3_t diffuse               = rt_vec3_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec3_add(final_color,
                                                      diffuse);
    }

    final_color = rt_vec3_add(final_color, material->ambient);
    return rt_vec3_clamp(final_color, 0.0f, 1.0f);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_fragment_shader_checkerboard(const rt_hit_info_t* hit_info,
                                                 const rt_checkerboard_material_t* material,
                                                 bool is_in_shadow)
{
    uint32_t pos_x_quant    = (uint32_t)(ceilf(hit_info->position.x * 0.5f));
    uint32_t pos_z_quant    = (uint32_t)(ceilf(hit_info->position.z * 0.5f));
    uint32_t pos_quant      = pos_x_quant + pos_z_quant;

    rt_vec3_t material_color = (pos_quant & 1) ? material->color_0
                                               : material->color_1;

    if (is_in_shadow && material->receives_shadows) {
        return rt_vec3_mul_scalar(material_color, material->shadow_factor);
    }

    return material_color;
}

///////////////////////////////////////////////////////////////////////////
////////////////////////////// FINAL COLOR ////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_is_in_shadow(const rt_world_t* world,
                            const rt_hit_info_t* hit_info,
                            rt_vec2_t range)
{
    uint32_t directional_light_count = world->directional_light_count;
    const rt_directional_light_t* directional_lights = world->directional_lights;

    for (uint32_t i = 0; i < directional_light_count; ++i) {
        const rt_directional_light_t* light = &directional_lights[i];

        rt_ray_t new_ray = {
            .org = rt_vec3_add(hit_info->position,
                               rt_vec3_mul_scalar(hit_info->normal, RT_SHADOW_BIAS)),
            .dir = rt_vec3_norm(rt_vec3_negate(light->direction)),
        };

        rt_hit_info_t new_hit_info = {};

        uint32_t geometry_hit_index = UINT32_MAX;

        enum rt_hit_geometry_type shadow_hit_type = rt_world_any_closest_hit(world,
                                                                             new_ray,
                                                                             range,
                                                                             &new_hit_info,
                                                                             &geometry_hit_index);

        if (RT_HIT_NONE != shadow_hit_type) {
            return true;
        }
    }
    
    uint32_t point_light_count = world->point_light_count;
    const rt_point_light_t* point_lights = world->point_lights;

    for (uint32_t i = 0; i < point_light_count; ++i) {
        const rt_point_light_t* light = &point_lights[i];

        rt_ray_t new_ray = {
            .org = rt_vec3_add(hit_info->position,
                               rt_vec3_mul_scalar(hit_info->normal, RT_SHADOW_BIAS)),
            .dir = rt_vec3_norm(rt_vec3_sub(light->position,
                                            hit_info->position)),
        };

        rt_hit_info_t new_hit_info = {};
        uint32_t geometry_hit_index = UINT32_MAX;

        enum rt_hit_geometry_type shadow_hit_type = rt_world_any_closest_hit(world,
                                                                             new_ray,
                                                                             range,
                                                                             &new_hit_info,
                                                                             &geometry_hit_index);

        if (RT_HIT_NONE != shadow_hit_type) {
            return true;
        }
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec3_t rt_world_compute_color(const rt_world_t* world,
                                        rt_ray_t ray,
                                        rt_vec2_t range,
                                        uint32_t depth)
{
    if (!depth) {
        return world->clear_color;
    }

    rt_hit_info_t hit_info = {};
    uint32_t geometry_hit_index = UINT32_MAX;

    enum rt_hit_geometry_type hit_type = rt_world_any_closest_hit(world,
                                                                  ray,
                                                                  range,
                                                                  &hit_info,
                                                                  &geometry_hit_index);

    uint32_t material_index             = UINT32_MAX;
    enum rt_material_type material_type = RT_MATERIAL_NONE;

    switch (hit_type) {
        case RT_HIT_SPHERE:
            material_index = world->spheres[geometry_hit_index].material_index;
            material_type = world->spheres[geometry_hit_index].material_type;
            break;
        case RT_HIT_NONE:
            return world->clear_color;
        default:
            RT_ASSERT(false && "Unhandled hit type");
            return world->clear_color;
    }

    if (UINT32_MAX == material_index || RT_MATERIAL_NONE == material_type) {
        return world->clear_color;
    }

    bool is_in_shadow = rt_is_in_shadow(world,
                                        &hit_info,
                                        range);

    switch (material_type) {
        case RT_MATERIAL_EMISSIVE: {
            rt_emissive_material_t* emissive_material = &world->emissive_materials[material_index];
            return rt_fragment_shader_emissive(emissive_material);
        }
        case RT_MATERIAL_CHECKERBOARD: {
            rt_checkerboard_material_t* checkerboard_material = &world->checkerboard_materials[material_index];
            return rt_fragment_shader_checkerboard(&hit_info,
                                                   checkerboard_material,
                                                   is_in_shadow);
        }
        case RT_MATERIAL_DIFFUSE: {
            rt_diffuse_material_t* diffuse_material = &world->diffuse_materials[material_index];
            return rt_fragment_shader_diffuse(&hit_info,
                                              diffuse_material,
                                              world,
                                              is_in_shadow);
        }
        /*
        case RT_MATERIAL_METALLIC: {
            rt_metallic_material_t* metallic_material = &world->metallic_materials[material_index];
            return rt_fragment_shader_metallic(&hit_info, metallic_material, world->point_lights, world->point_light_count, false);
        }
        case RT_MATERIAL_DIELECTRIC: {
            rt_dielectric_material_t* dielectric_material = &world->dielectric_materials[material_index];
            return rt_fragment_shader_dielectric(&hit_info, dielectric_material, world->point_lights, world->point_light_count, false);
        }
        */
        default:
            RT_ASSERT(false && "Unhandled material type");
            return world->clear_color;
    }
}

#endif // RT_RAYTERM_H
