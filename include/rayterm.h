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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

///////////////////////////////////////////////////////////////////////////
//////////////////////////////// DEFINES //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define RT_API                          static inline

///////////////////////////////////////////////////////////////////////////
#ifdef RT_USE_FLOAT64
#define RT_FLOAT(X)                     X
#else
#define RT_FLOAT(X)                     X ## f
#endif

///////////////////////////////////////////////////////////////////////////
#define RT_GAMMA                        RT_FLOAT(2.2)
#define RT_GAMMA_INVERSE                RT_FLOAT(0.454545)
#define RT_PI                           RT_FLOAT(3.1415926)
#define RT_EPSILON                      RT_FLOAT(0.000001)
#define RT_SHADOW_BIAS                  RT_FLOAT(0.001)

///////////////////////////////////////////////////////////////////////////
#define RT_INIT_CAP                     8

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// TYPES ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
#ifdef RT_USE_FLOAT64
typedef double rt_float_t;
#else
typedef float rt_float_t;
#endif

///////////////////////////////////////////////////////////////////////////
#ifdef RT_USE_IDX64
typedef int64_t rt_idx_t;

#define RT_IDX_MIN INT64_MIN
#define RT_IDX_MAX INT64_MAX

#else
typedef int32_t rt_idx_t;

#define RT_IDX_MIN INT32_MIN
#define RT_IDX_MAX INT32_MAX

#endif

///////////////////////////////////////////////////////////////////////////
enum rt_status
{
    RT_STATUS_success,
    RT_STATUS_failure,
};
typedef enum rt_status rt_status_t;

///////////////////////////////////////////////////////////////////////////
#define RT_ASSERT(EXPR)                                                     \
    do {                                                                    \
        if (!(EXPR)) {                                                      \
            fprintf(stderr, #EXPR " failed at %s (%s) : %u\n",              \
                            __FILE__,                                       \
                            __PRETTY_FUNCTION__,                            \
                            __LINE__);                                      \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    }                                                                       \
    while (0);

///////////////////////////////////////////////////////////////////////////
#define RT_ASSERT_STATUS(STATUS)        RT_ASSERT(STATUS == RT_STATUS_success)

///////////////////////////////////////////////////////////////////////////
//////////////////////////// MACRO FUNCTIONS //////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define RT_BUFFER_LEN(BUFFER)           (sizeof(BUFFER) / sizeof((BUFFER)[0]))

///////////////////////////////////////////////////////////////////////////
#define RT_TO_RADIANS(DEGREES)          (DEGREES * RT_FLOAT(0.0174533))
#define RT_TO_DEGREES(RADIANS)          (RADIANS * RT_FLOAT(57.295779))

///////////////////////////////////////////////////////////////////////////
#define RT_CONCAT2(A, B)                A ## B
#define RT_CONCAT3(A, B, C)             A ## B ## C
#define RT_CONCAT4(A, B, C, D)          A ## B ## C ## D

///////////////////////////////////////////////////////////////////////////
#define RT_ALLOC                        malloc
#define RT_FREE                         free

///////////////////////////////////////////////////////////////////////////
/////////////////////////////// UTILITIES /////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define pow powf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define sqrt sqrtf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define fmin fminf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define fmax fmaxf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define fabs fabsf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define cos cosf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define sin sinf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define ceil ceilf
#endif

///////////////////////////////////////////////////////////////////////////
#ifndef RT_USE_FLOAT64
#define floor floorf
#endif

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_apply_inverse_gamma(rt_float_t x)
{
    return pow(x, RT_GAMMA_INVERSE);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_apply_gamma(rt_float_t x)
{
    return pow(x, RT_GAMMA);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_apply_gamma_custom(rt_float_t x,
                                        rt_float_t gamma)
{
    return pow(x, gamma);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_clamp(rt_float_t x,
                           rt_float_t minimum,
                           rt_float_t maximum)
{
    return fmin(maximum, fmax(x, minimum));
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// MATH ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_float_t x, y, z, w;
}
rt_vec4_t;

///////////////////////////////////////////////////////////////////////////
RT_API void rt_vec4_zero(rt_vec4_t* p)
{
    RT_ASSERT(p != NULL);

    p->x = RT_FLOAT(0.0);
    p->y = RT_FLOAT(0.0);
    p->z = RT_FLOAT(0.0);
    p->w = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_vec4_set(rt_vec4_t* p,
                        rt_float_t x,
                        rt_float_t y,
                        rt_float_t z,
                        rt_float_t w)
{
    RT_ASSERT(p != NULL);

    p->x = x;
    p->y = y;
    p->z = z;
    p->w = w;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_add(rt_vec4_t p,
                             rt_vec4_t q)
{
    p.x += q.x;
    p.y += q.y;
    p.z += q.z;
    p.w += q.w;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_add_scalar(rt_vec4_t   p,
                                    rt_float_t  k)
{
    p.x += k;
    p.y += k;
    p.z += k;
    p.w += k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_sub(rt_vec4_t p,
                             rt_vec4_t q)
{
    p.x -= q.x;
    p.y -= q.y;
    p.z -= q.z;
    p.w -= q.w;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_sub_scalar(rt_vec4_t   p,
                                    rt_float_t  k)
{
    p.x -= k;
    p.y -= k;
    p.z -= k;
    p.w -= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_mul(rt_vec4_t p,
                             rt_vec4_t q)
{
    p.x *= q.x;
    p.y *= q.y;
    p.z *= q.z;
    p.w *= q.w;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_mul_scalar(rt_vec4_t   p,
                                    rt_float_t  k)
{
    p.x *= k;
    p.y *= k;
    p.z *= k;
    p.w *= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_div(rt_vec4_t p,
                             rt_vec4_t q)
{
    p.x /= q.x;
    p.y /= q.y;
    p.z /= q.z;
    p.w /= q.w;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_div_scalar(rt_vec4_t   p,
                                    rt_float_t  k)
{
    p.x /= k;
    p.y /= k;
    p.z /= k;
    p.w /= k;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_negate(rt_vec4_t p)
{
    p.x *= RT_FLOAT(-1.0);
    p.y *= RT_FLOAT(-1.0);
    p.z *= RT_FLOAT(-1.0);
    p.w *= RT_FLOAT(-1.0);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_vec4_len(rt_vec4_t p)
{
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z + p.w * p.w);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_vec4_sqrlen(rt_vec4_t p)
{
    return p.x * p.x + p.y * p.y + p.z * p.z + p.w * p.w;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_vec4_dot(rt_vec4_t p,
                              rt_vec4_t q)
{
    return p.x * q.x + p.y * q.y + p.z * q.z + p.w * q.w;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_cross(rt_vec4_t p,
                               rt_vec4_t q)
{
    rt_vec4_t r = {
        .x = p.y * q.z - p.z * q.y,
        .y = p.z * q.x - p.x * q.z,
        .z = p.x * q.y - p.y * q.x,
        .w = 0.0f,
    };

    return r;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_norm(rt_vec4_t p)
{
    rt_float_t k = p.x * p.x + p.y * p.y + p.z * p.z + p.w * p.w;
    if (k > RT_FLOAT(0.0)) {
        k = RT_FLOAT(1.0) / sqrt(k);
        p.x *= k;
        p.y *= k;
        p.z *= k;
        p.w *= k;
    }
    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_vec4_to_uint32(rt_vec4_t p)
{
    uint32_t x = (uint32_t)(p.x * RT_FLOAT(255.99));
    uint32_t y = (uint32_t)(p.y * RT_FLOAT(255.99));
    uint32_t z = (uint32_t)(p.z * RT_FLOAT(255.99));
    uint32_t w = (uint32_t)(p.w * RT_FLOAT(255.99));

    return (w << 24) | (z << 16) | (y << 8) | (x << 0);
}

///////////////////////////////////////////////////////////////////////////
RT_API uint32_t rt_vec4_to_uint32_alpha(rt_vec4_t   p,
                                        rt_float_t  a)
{
    uint32_t x = (uint32_t)(p.x * RT_FLOAT(255.99));
    uint32_t y = (uint32_t)(p.y * RT_FLOAT(255.99));
    uint32_t z = (uint32_t)(p.z * RT_FLOAT(255.99));
    uint32_t w = (uint32_t)(  a * RT_FLOAT(255.99));

    return (w << 24) | (z << 16) | (y << 8) | (x << 0);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_apply_0(rt_vec4_t  p,
                                 rt_float_t (*fun)(void))
{
    p.x = fun();
    p.y = fun();
    p.z = fun();
    p.w = fun();

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_apply_1(rt_vec4_t  p,
                                 rt_float_t (*fun)(rt_float_t))
{
    p.x = fun(p.x);
    p.y = fun(p.y);
    p.z = fun(p.z);
    p.w = fun(p.w);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_apply_2(rt_vec4_t  p,
                                 rt_float_t (*fun)(rt_float_t, rt_float_t),
                                 rt_float_t k)
{
    p.x = fun(p.x, k);
    p.y = fun(p.y, k);
    p.z = fun(p.z, k);
    p.w = fun(p.w, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_lerp(rt_vec4_t     p,
                              rt_vec4_t     q,
                              rt_float_t    k)
{
    p.x += k * (q.x - p.x);
    p.y += k * (q.y - p.y);
    p.z += k * (q.z - p.z);
    p.w += k * (q.w - p.w);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_rand(uint32_t seed)
{
    srand(seed);
    rt_vec4_t p = {
        .x = (rt_float_t)(rand() % RAND_MAX),
        .y = (rt_float_t)(rand() % RAND_MAX),
        .z = (rt_float_t)(rand() % RAND_MAX),
        .w = (rt_float_t)(rand() % RAND_MAX),
    };

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_rand_alpha(uint32_t seed,
                                    rt_float_t alpha)
{
    srand(seed);
    rt_vec4_t p = {
        .x = (rt_float_t)(rand() % RAND_MAX),
        .y = (rt_float_t)(rand() % RAND_MAX),
        .z = (rt_float_t)(rand() % RAND_MAX),
        .w = alpha,
    };

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_max(rt_vec4_t  p,
                             rt_float_t k)
{
    p.x = fmax(p.x, k);
    p.y = fmax(p.y, k);
    p.z = fmax(p.z, k);
    p.w = fmax(p.w, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_min(rt_vec4_t  p,
                             rt_float_t k)
{
    p.x = fmin(p.x, k);
    p.y = fmin(p.y, k);
    p.z = fmin(p.z, k);
    p.w = fmin(p.w, k);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_clamp(rt_vec4_t            p,
                               rt_float_t           min,
                               rt_float_t           max)
{
    p.x = rt_clamp(p.x, min, max);
    p.y = rt_clamp(p.y, min, max);
    p.z = rt_clamp(p.z, min, max);
    p.w = rt_clamp(p.w, min, max);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_reflect(rt_vec4_t p,
                                 rt_vec4_t n)
{
    n               = rt_vec4_norm(n);
    rt_float_t d    = RT_FLOAT(2.0) * rt_vec4_dot(p, n);

    p.x -= n.x * d;
    p.y -= n.y * d;
    p.z -= n.z * d;
    p.w  = 0.0f;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_refract(rt_vec4_t      p,
                                 rt_vec4_t      n,
                                 rt_float_t     k)
{
    rt_vec4_t pnegated      = rt_vec4_negate(p);
    rt_float_t cos_theta    = fmin(RT_FLOAT(1.0), rt_vec4_dot(pnegated, n));
    rt_vec4_t perpendicular = rt_vec4_mul_scalar(rt_vec4_add(p,
                                                             rt_vec4_mul_scalar(n,
                                                                                cos_theta)),
                                                             k);
    rt_float_t a            = -sqrt(fabs(RT_FLOAT(1.0) - rt_vec4_sqrlen(perpendicular)));
    rt_vec4_t parallel      = rt_vec4_mul_scalar(n, a);
    return rt_vec4_add(perpendicular, parallel);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_vec4_dist(rt_vec4_t p,
                               rt_vec4_t q)
{
    rt_float_t x = p.x - q.x;
    rt_float_t y = p.y - q.y;
    rt_float_t z = p.z - q.z;
    rt_float_t w = p.w - q.w;

    return sqrt(x * x + y * y + z * z + w * w);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_vec4_sqrdist(rt_vec4_t p,
                                  rt_vec4_t q)
{
    rt_float_t x = p.x - q.x;
    rt_float_t y = p.y - q.y;
    rt_float_t z = p.z - q.z;
    rt_float_t w = p.w - q.w;

    return x * x + y * y + z * z + w * w;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_rotate_x(rt_vec4_t     p,
                                  rt_float_t    r)
{
    rt_float_t py = p.y;
    rt_float_t pz = p.z;

    rt_float_t cos_r = cos(r);
    rt_float_t sin_r = sin(r);

    p.y = py * cos_r - pz * sin_r;
    p.z = py * sin_r + pz * cos_r;

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_rotate_y(rt_vec4_t     p,
                                  rt_float_t    r)
{
    rt_float_t px = p.x;
    rt_float_t pz = p.z;

    rt_float_t cos_r = cos(r);
    rt_float_t sin_r = sin(r);

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
    RT_MATERIAL_TYPE_null_material,
    RT_MATERIAL_TYPE_emissive_material,
    RT_MATERIAL_TYPE_checkerboard_material,
    RT_MATERIAL_TYPE_diffuse_material,
    RT_MATERIAL_TYPE_metallic_material,
    RT_MATERIAL_TYPE_dielectric_material,
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   color;
}
rt_emissive_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   color_0;
    rt_vec4_t   color_1;
    rt_float_t  shadow_factor;
    bool        receives_shadows;
}
rt_checkerboard_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   ambient;
    rt_vec4_t   diffuse;
    bool        receives_shadows;
}
rt_diffuse_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   ambient;
    rt_vec4_t   specular;
    bool        receives_shadows;
}
rt_metallic_material_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   ambient;
    rt_vec4_t   specular;
    rt_float_t  refractive_index;
    bool        receives_shadows;
}
rt_dielectric_material_t;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// RAYS ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t org;
    rt_vec4_t dir;
}
rt_ray_t;

///////////////////////////////////////////////////////////////////////////
enum rt_hit_geometry_type
{
    RT_HIT_null,
    RT_HIT_sphere,
    RT_HIT_plane,
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   position;
    rt_vec4_t   normal;
    rt_float_t  t;
    bool        is_front_facing;
}
rt_hit_info_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_hit_info_t               info;
    enum rt_hit_geometry_type   geometry_type;
    rt_idx_t                    geometry_index;
}
rt_hit_ext_info_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_ray_at(rt_ray_t         ray,
                           rt_float_t       t)
{
    rt_vec4_t dir_scaled = rt_vec4_mul_scalar(ray.dir, t);
    return rt_vec4_add(ray.org, dir_scaled);
}

///////////////////////////////////////////////////////////////////////////
//////////////////////////// GEOMETRY /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t               center;
    rt_float_t              radius;

    rt_idx_t                material_index;
    enum rt_material_type   material_type;
}
rt_sphere_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t               position;
    rt_vec4_t               normal;

    rt_idx_t                material_index;
    enum rt_material_type   material_type;
}
rt_plane_t;

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_sphere_hit(rt_sphere_t               sphere,
                          rt_ray_t                  ray,
                          rt_float_t                nearestZ,
                          rt_float_t                farthestZ,
                          rt_hit_info_t*            info)
{
    RT_ASSERT(info      != NULL);
    RT_ASSERT(nearestZ  <= farthestZ);

    rt_vec4_t oc    = rt_vec4_sub(sphere.center, ray.org);
    rt_float_t a    = rt_vec4_sqrlen(ray.dir);
    rt_float_t h    = rt_vec4_dot(ray.dir, oc);
    rt_float_t c    = rt_vec4_sqrlen(oc) - sphere.radius * sphere.radius;

    rt_float_t disc = h * h - a * c;
    if (disc < RT_FLOAT(0.0)) {
        return false;
    }

    rt_float_t sqrt_disc     = sqrt(disc);
    rt_float_t root_nearest  = (h - sqrt_disc) / a;
    rt_float_t root_farthest = (h + sqrt_disc) / a;

    rt_float_t range_min     = nearestZ;
    rt_float_t range_max     = farthestZ;

    rt_float_t root_chosen   = root_nearest;

    if (root_chosen <= range_min || root_chosen >= range_max) {

        root_chosen = root_farthest;

        if (root_chosen <= range_min || root_chosen >= range_max) {
            return false;
        }
    }

    info->position          = rt_ray_at(ray, root_chosen);

    rt_vec4_t n             = rt_vec4_sub(info->position, sphere.center);
    n                       = rt_vec4_div_scalar(n, sphere.radius);

    info->normal            = rt_vec4_norm(n);
    info->t                 = root_chosen;
    info->is_front_facing   = rt_vec4_dot(info->normal, ray.dir) < RT_FLOAT(0.0);

    return true;
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_plane_hit(rt_plane_t             plane,
                         rt_ray_t               ray,
                         rt_float_t             nearestZ,
                         rt_float_t             farthestZ,
                         rt_hit_info_t*         info)
{
    RT_ASSERT(info      != NULL);
    RT_ASSERT(nearestZ  <= farthestZ);

    rt_float_t denom = rt_vec4_dot(plane.normal, ray.dir);

    if (denom < RT_EPSILON) {
        return false;
    }

    rt_vec4_t dir = rt_vec4_sub(plane.position, ray.org);
    rt_float_t t = rt_vec4_dot(dir, plane.normal) / denom;

    if (t < nearestZ || t > farthestZ) {
        return false;
    }

    info->position          = rt_ray_at(ray, t);
    info->normal            = plane.normal;
    info->t                 = t;
    info->is_front_facing   = true;

    return true;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////// LIGHTS //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
enum rt_light_type
{
    RT_LIGHT_null_light,
    RT_LIGHT_directional_light,
    RT_LIGHT_point_light,
};

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   color;
    rt_vec4_t   direction;
    rt_float_t  intensity;
    bool        casts_shadows;
}
rt_directional_light_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   color;
    rt_vec4_t   position;
    rt_float_t  intensity;
    bool        casts_shadows;
}
rt_point_light_t;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// WORLD ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
enum rt_face_cull_mode
{
    RT_FACE_cull_null,
    RT_FACE_cull_front   = 1 << 0,
    RT_FACE_cull_back    = 1 << 1,
    RT_FACE_cull_both    = RT_FACE_cull_front | RT_FACE_cull_back,
};

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_DEF_BUFFER(NAME)                                       \
                                                                        \
    RT_CONCAT3(rt_, NAME, _t) * RT_CONCAT2(NAME, _buffer);              \
                                                                        \
    rt_idx_t RT_CONCAT2(NAME, _count);                                  \
    rt_idx_t RT_CONCAT2(NAME, _capacity);                               \

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    ///////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// SKY ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    rt_vec4_t clear_color;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// LIGHTS //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    RT_WORLD_DEF_BUFFER(directional_light)
    RT_WORLD_DEF_BUFFER(point_light)

    ///////////////////////////////////////////////////////////////////////////
    //////////////////////////////// GEOMETRY /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    RT_WORLD_DEF_BUFFER(sphere)
    RT_WORLD_DEF_BUFFER(plane);

    ///////////////////////////////////////////////////////////////////////////
    /////////////////////////////// MATERIALS /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    RT_WORLD_DEF_BUFFER(emissive_material)
    RT_WORLD_DEF_BUFFER(checkerboard_material)
    RT_WORLD_DEF_BUFFER(diffuse_material)
    RT_WORLD_DEF_BUFFER(metallic_material)
    RT_WORLD_DEF_BUFFER(dielectric_material)

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// STATE ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    enum rt_face_cull_mode face_cull_mode;
}
rt_world_t;

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_DEFINE_PUSH(OBJECT) [[nodiscard]]                          \
RT_API rt_status_t RT_CONCAT2(rt_world_push_, OBJECT)                       \
(rt_world_t* world, rt_idx_t* idx)                                          \
{                                                                           \
    RT_ASSERT(world != NULL);                                               \
    RT_ASSERT(idx   != NULL);                                               \
                                                                            \
    typedef RT_CONCAT3(rt_, OBJECT, _t) buffer_type;                        \
    buffer_type* buffer = world->RT_CONCAT2(OBJECT, _buffer);               \
                                                                            \
    rt_idx_t count     = world->RT_CONCAT2(OBJECT, _count);                 \
    RT_ASSERT(count >= 0);                                                  \
                                                                            \
    rt_idx_t capacity  = world->RT_CONCAT2(OBJECT, _capacity);              \
    RT_ASSERT(capacity >= 0);                                               \
                                                                            \
    if (count == RT_IDX_MAX) {                                              \
        return RT_STATUS_failure;                                           \
    }                                                                       \
    else if (count >= capacity) {                                           \
                                                                            \
        rt_idx_t capacity_new = capacity ? 2 * capacity : RT_INIT_CAP;      \
        if (capacity_new <= 0) {                                            \
            return RT_STATUS_failure;                                       \
        }                                                                   \
                                                                            \
        size_t bytes_required = (size_t)capacity_new * sizeof(buffer_type); \
        buffer_type* buffer_new = (buffer_type*)RT_ALLOC(bytes_required);   \
                                                                            \
        if (!buffer_new) {                                                  \
            return RT_STATUS_failure;                                       \
        }                                                                   \
                                                                            \
        for (rt_idx_t i = 0; i < count; ++i) {                              \
            buffer_new[i] = buffer[i];                                      \
        }                                                                   \
                                                                            \
        RT_FREE(buffer);                                                    \
                                                                            \
        world->RT_CONCAT2(OBJECT, _buffer)     = buffer_new;                \
        world->RT_CONCAT2(OBJECT, _capacity)   = capacity_new;              \
    }                                                                       \
                                                                            \
    world->RT_CONCAT2(OBJECT, _buffer)[count]  = (buffer_type){};           \
    world->RT_CONCAT2(OBJECT, _count)         += 1;                         \
                                                                            \
    *idx = count;                                                           \
                                                                            \
    return RT_STATUS_success;                                               \
}

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_DEFINE_POP(OBJECT) [[nodiscard]]                           \
RT_API rt_status_t RT_CONCAT2(rt_world_pop_, OBJECT)                        \
(rt_world_t* world, RT_CONCAT3(rt_, OBJECT, _t) *obj)                       \
{                                                                           \
    RT_ASSERT(world     != NULL);                                           \
    RT_ASSERT(obj       != NULL);                                           \
                                                                            \
    rt_idx_t count = world->RT_CONCAT2(OBJECT, _count);                     \
    RT_ASSERT(count >= 0);                                                  \
                                                                            \
    if (!count) {                                                           \
        return RT_STATUS_failure;                                           \
    }                                                                       \
                                                                            \
    world->RT_CONCAT2(OBJECT, _count) -= 1;                                 \
                                                                            \
    typedef RT_CONCAT3(rt_, OBJECT, _t) buffer_type;                        \
    buffer_type* buffer = world->RT_CONCAT2(OBJECT, _buffer);               \
                                                                            \
    RT_ASSERT(buffer != NULL);                                              \
                                                                            \
    memcpy(obj, &buffer[count], sizeof(buffer_type));                       \
                                                                            \
    return RT_STATUS_success;                                               \
}

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_DEFINE_RESERVE(OBJECT) [[nodiscard]]                       \
RT_API rt_status_t RT_CONCAT2(rt_world_reserve_, OBJECT)                    \
(rt_world_t* world, rt_idx_t capacity)                                      \
{                                                                           \
    RT_ASSERT(world     != NULL);                                           \
    RT_ASSERT(capacity  > 0);                                               \
                                                                            \
    rt_idx_t count          = world->RT_CONCAT2(OBJECT, _count);            \
    RT_ASSERT(count >= 0);                                                  \
                                                                            \
    typedef RT_CONCAT3(rt_, OBJECT, _t) buffer_type;                        \
                                                                            \
    size_t bytes_required   = (size_t)capacity * sizeof(buffer_type);       \
    buffer_type* buffer_new = (buffer_type*)RT_ALLOC(bytes_required);       \
                                                                            \
    if (!buffer_new) {                                                      \
        return RT_STATUS_failure;                                           \
    }                                                                       \
                                                                            \
    for (rt_idx_t i = 0; i < count; ++i) {                                  \
        buffer_new[i] = world->RT_CONCAT2(OBJECT, _buffer)[i];              \
    }                                                                       \
                                                                            \
    RT_FREE(world->RT_CONCAT2(OBJECT, _buffer));                            \
                                                                            \
    world->RT_CONCAT2(OBJECT, _buffer)     = buffer_new;                    \
    world->RT_CONCAT2(OBJECT, _capacity)   = capacity;                      \
                                                                            \
    return RT_STATUS_success;                                               \
}

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_DEFINE_FREE(OBJECT)                                        \
RT_API void RT_CONCAT3(rt_world_free_, OBJECT, s)                           \
(rt_world_t* world)                                                         \
{                                                                           \
    rt_idx_t capacity = world->RT_CONCAT2(OBJECT, _capacity);               \
    RT_ASSERT(capacity >= 0);                                               \
                                                                            \
    RT_FREE(world->RT_CONCAT2(OBJECT, _buffer));                            \
    world->RT_CONCAT2(OBJECT, _buffer)     = NULL;                          \
    world->RT_CONCAT2(OBJECT, _capacity)   = 0;                             \
    world->RT_CONCAT2(OBJECT, _count)      = 0;                             \
}

///////////////////////////////////////////////////////////////////////////
RT_WORLD_DEFINE_PUSH(directional_light)
RT_WORLD_DEFINE_PUSH(point_light)
RT_WORLD_DEFINE_PUSH(sphere)
RT_WORLD_DEFINE_PUSH(plane)
RT_WORLD_DEFINE_PUSH(emissive_material)
RT_WORLD_DEFINE_PUSH(checkerboard_material)
RT_WORLD_DEFINE_PUSH(diffuse_material)
RT_WORLD_DEFINE_PUSH(metallic_material)
RT_WORLD_DEFINE_PUSH(dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_WORLD_DEFINE_POP(directional_light)
RT_WORLD_DEFINE_POP(point_light)
RT_WORLD_DEFINE_POP(sphere)
RT_WORLD_DEFINE_POP(plane)
RT_WORLD_DEFINE_POP(emissive_material)
RT_WORLD_DEFINE_POP(checkerboard_material)
RT_WORLD_DEFINE_POP(diffuse_material)
RT_WORLD_DEFINE_POP(metallic_material)
RT_WORLD_DEFINE_POP(dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_WORLD_DEFINE_RESERVE(directional_light)
RT_WORLD_DEFINE_RESERVE(point_light)
RT_WORLD_DEFINE_RESERVE(sphere)
RT_WORLD_DEFINE_RESERVE(plane)
RT_WORLD_DEFINE_RESERVE(emissive_material)
RT_WORLD_DEFINE_RESERVE(checkerboard_material)
RT_WORLD_DEFINE_RESERVE(diffuse_material)
RT_WORLD_DEFINE_RESERVE(metallic_material)
RT_WORLD_DEFINE_RESERVE(dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_WORLD_DEFINE_FREE(directional_light)
RT_WORLD_DEFINE_FREE(point_light)
RT_WORLD_DEFINE_FREE(sphere)
RT_WORLD_DEFINE_FREE(plane)
RT_WORLD_DEFINE_FREE(emissive_material)
RT_WORLD_DEFINE_FREE(checkerboard_material)
RT_WORLD_DEFINE_FREE(diffuse_material)
RT_WORLD_DEFINE_FREE(metallic_material)
RT_WORLD_DEFINE_FREE(dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_free(rt_world_t* world)
{
    rt_world_free_directional_lights            (world);
    rt_world_free_point_lights                  (world);
    rt_world_free_spheres                       (world);
    rt_world_free_planes                        (world);
    rt_world_free_emissive_materials            (world);
    rt_world_free_checkerboard_materials        (world);
    rt_world_free_diffuse_materials             (world);
    rt_world_free_metallic_materials            (world);
    rt_world_free_dielectric_materials          (world);
}

///////////////////////////////////////////////////////////////////////////
#define RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(GEOMETRY, MATERIAL)             \
RT_API void RT_CONCAT4(rt_, GEOMETRY, _link_, MATERIAL)                     \
(rt_world_t* world, rt_idx_t geometry_index, rt_idx_t material_index)       \
{                                                                           \
    RT_ASSERT(world != NULL);                                               \
                                                                            \
    RT_ASSERT(geometry_index    >= 0);                                      \
    RT_ASSERT(material_index    >= 0);                                      \
                                                                            \
    RT_ASSERT(geometry_index    < world->RT_CONCAT2(GEOMETRY, _count));     \
    RT_ASSERT(material_index    < world->RT_CONCAT2(MATERIAL, _count));     \
                                                                            \
    typedef RT_CONCAT3(rt_, GEOMETRY, _t) geometry_type;                    \
                                                                            \
    geometry_type* geometry_buffer  = world->RT_CONCAT2(GEOMETRY, _buffer); \
    RT_ASSERT(geometry_buffer != NULL);                                     \
                                                                            \
    geometry_type* geometry         = &geometry_buffer[geometry_index];     \
    RT_ASSERT(geometry != NULL);                                            \
                                                                            \
    geometry->material_type     = RT_CONCAT2(RT_MATERIAL_TYPE_, MATERIAL);  \
    geometry->material_index    = material_index;                           \
}                                                                           \

///////////////////////////////////////////////////////////////////////////
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(sphere, emissive_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(sphere, checkerboard_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(sphere, diffuse_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(sphere, metallic_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(sphere, dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(plane, emissive_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(plane, checkerboard_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(plane, diffuse_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(plane, metallic_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(plane, dielectric_material)

///////////////////////////////////////////////////////////////////////////
#define RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(GEOMETRY, MATERIAL)           \
RT_API void RT_CONCAT4(rt_, GEOMETRY, _unlink_, MATERIAL)                   \
(rt_world_t* world, rt_idx_t geometry_index)                                \
{                                                                           \
    RT_ASSERT(world != NULL);                                               \
                                                                            \
    RT_ASSERT(geometry_index    >= 0);                                      \
    RT_ASSERT(geometry_index    < world->RT_CONCAT2(GEOMETRY, _count));     \
                                                                            \
    typedef RT_CONCAT3(rt_, GEOMETRY, _t) geometry_type;                    \
                                                                            \
    geometry_type* geometry_buffer  = world->RT_CONCAT2(GEOMETRY, _buffer); \
    RT_ASSERT(geometry_buffer != NULL);                                     \
                                                                            \
    geometry_type* geometry         = &geometry_buffer[geometry_index];     \
    RT_ASSERT(geometry != NULL);                                            \
                                                                            \
    geometry->material_type     = RT_MATERIAL_TYPE_null_material;           \
    geometry->material_index    = 0;                                        \
}                                                                           \

///////////////////////////////////////////////////////////////////////////
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(sphere, emissive_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(sphere, checkerboard_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(sphere, diffuse_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(sphere, metallic_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(sphere, dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(plane, emissive_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(plane, checkerboard_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(plane, diffuse_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(plane, metallic_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(plane, dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_directional_light_params(rt_world_t*   world,
                                                  rt_idx_t      light_index,
                                                  rt_vec4_t     color,
                                                  rt_vec4_t     direction,
                                                  rt_float_t    intensity,
                                                  bool          casts_shadows)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(light_index   >= 0);

    rt_directional_light_t* light_buffer = world->directional_light_buffer;
    RT_ASSERT(light_buffer != NULL);

    rt_directional_light_t* light = &light_buffer[light_index];
    RT_ASSERT(light != NULL);

    light->color            = color;
    light->direction        = direction;
    light->intensity        = intensity;
    light->casts_shadows    = casts_shadows;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_point_light_params(rt_world_t*     world,
                                            rt_idx_t        light_index,
                                            rt_vec4_t       color,
                                            rt_vec4_t       position,
                                            rt_float_t      intensity,
                                            bool            casts_shadows)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(light_index   >= 0);

    rt_point_light_t* light_buffer = world->point_light_buffer;
    RT_ASSERT(light_buffer != NULL);

    rt_point_light_t* light = &light_buffer[light_index];
    RT_ASSERT(light != NULL);

    light->color            = color;
    light->position         = position;
    light->intensity        = intensity;
    light->casts_shadows    = casts_shadows;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_sphere_params(rt_world_t*          world,
                                       rt_idx_t             sphere_index,
                                       rt_vec4_t            sphere_center,
                                       rt_float_t           sphere_radius)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(sphere_index  >= 0);

    rt_sphere_t* sphere_buffer = world->sphere_buffer;
    RT_ASSERT(sphere_buffer != NULL);

    rt_sphere_t* sphere = &sphere_buffer[sphere_index];
    RT_ASSERT(sphere != NULL);

    sphere->center = sphere_center;
    sphere->radius = sphere_radius;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_plane_params(rt_world_t*           world,
                                      rt_idx_t              plane_index,
                                      rt_vec4_t             plane_position,
                                      rt_vec4_t             plane_normal)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(plane_index   >= 0);

    rt_plane_t* plane_buffer = world->plane_buffer;
    RT_ASSERT(plane_buffer != NULL);

    rt_plane_t* plane = &plane_buffer[plane_index];
    RT_ASSERT(plane != NULL);

    plane->position     = plane_position;
    plane->normal       = plane_normal;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_emissive_material_params(rt_world_t*   world,
                                                  rt_idx_t      material_index,
                                                  rt_vec4_t     color)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);

    rt_emissive_material_t* material_buffer = world->emissive_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_emissive_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    material->color = color;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_checkerboard_material_params(rt_world_t*   world,
                                                      rt_idx_t      material_index,
                                                      rt_vec4_t     color_0,
                                                      rt_vec4_t     color_1,
                                                      rt_float_t    shadow_factor,
                                                      bool          receives_shadows)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);

    rt_checkerboard_material_t* material_buffer = world->checkerboard_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_checkerboard_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    material->color_0           = color_0;
    material->color_1           = color_1;
    material->shadow_factor     = shadow_factor;
    material->receives_shadows  = receives_shadows;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_diffuse_material_params(rt_world_t*    world,
                                                 rt_idx_t       material_index,
                                                 rt_vec4_t      ambient,
                                                 rt_vec4_t      diffuse,
                                                 bool           receives_shadows)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);

    rt_diffuse_material_t* material_buffer = world->diffuse_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_diffuse_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    material->ambient           = ambient;
    material->diffuse           = diffuse;
    material->receives_shadows  = receives_shadows;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_metallic_material_params(rt_world_t*   world,
                                                  rt_idx_t      material_index,
                                                  rt_vec4_t     ambient,
                                                  rt_vec4_t     specular,
                                                  bool          receives_shadows)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);

    rt_metallic_material_t* material_buffer = world->metallic_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_metallic_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    material->ambient           = ambient;
    material->specular          = specular;
    material->receives_shadows  = receives_shadows;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_dielectric_material_params(rt_world_t*     world,
                                                    rt_idx_t        material_index,
                                                    rt_vec4_t       ambient,
                                                    rt_vec4_t       specular,
                                                    rt_float_t      refractive_index,
                                                    bool            receives_shadows)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);

    rt_dielectric_material_t* material_buffer = world->dielectric_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_dielectric_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    material->ambient           = ambient;
    material->specular          = specular;
    material->refractive_index  = refractive_index;
    material->receives_shadows  = receives_shadows;
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_sphere_closest_hit(const rt_world_t*       world,
                                        rt_ray_t                ray,
                                        rt_float_t              nearestZ,
                                        rt_float_t              farthestZ,
                                        rt_hit_ext_info_t*      info)
{
    RT_ASSERT(world     != NULL);
    RT_ASSERT(nearestZ  <= farthestZ);
    RT_ASSERT(info      != NULL);

    enum rt_face_cull_mode cull_mode = world->face_cull_mode;

    if (RT_FACE_cull_both == cull_mode) {

        info->geometry_type         = RT_HIT_null;
        info->geometry_index        = 0;

        return false;
    }

    rt_hit_info_t closest_hit_info  = {};
    rt_idx_t closest_hit_index      = -1;

    rt_idx_t sphere_count           = world->sphere_count;
    RT_ASSERT(sphere_count         >= 0);

    const rt_sphere_t* spheres      = world->sphere_buffer;

    for (rt_idx_t i = 0; i < sphere_count; ++i) {

        rt_hit_info_t hit_info = {};

        if (rt_sphere_hit(spheres[i],
                          ray,
                          nearestZ,
                          farthestZ,
                          &hit_info)) {

            bool is_front_facing = hit_info.is_front_facing;

            if ((RT_FACE_cull_front == cull_mode &&  is_front_facing) ||
                (RT_FACE_cull_back  == cull_mode && !is_front_facing)) {

                continue;
            }

            if (-1 == closest_hit_index ||
                hit_info.t < closest_hit_info.t) {

                closest_hit_info    = hit_info;
                closest_hit_index   = i;
            }
        }
    }

    if (-1 != closest_hit_index) {

        info->info              = closest_hit_info;
        info->geometry_type     = RT_HIT_sphere;
        info->geometry_index    = closest_hit_index;

        return true;
    }

    info->geometry_type         = RT_HIT_null;
    info->geometry_index        = 0;

    return false;
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_plane_closest_hit(const rt_world_t*    world,
                                       rt_ray_t             ray,
                                       rt_float_t           nearestZ,
                                       rt_float_t           farthestZ,
                                       rt_hit_ext_info_t*   info)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(nearestZ      <= farthestZ);
    RT_ASSERT(info          != NULL);

    enum rt_face_cull_mode cull_mode = world->face_cull_mode;

    if (RT_FACE_cull_both == cull_mode) {

        info->geometry_type     = RT_HIT_null;
        info->geometry_index    = 0;

        return false;
    }

    rt_hit_info_t closest_hit_info  = {};
    rt_idx_t closest_hit_index      = -1;

    rt_idx_t plane_count            = world->plane_count;
    RT_ASSERT(plane_count          >= 0);

    const rt_plane_t* planes        = world->plane_buffer;

    for (rt_idx_t i = 0; i < plane_count; ++i) {

        rt_hit_info_t hit_info = {};

        if (rt_plane_hit(planes[i],
                         ray,
                         nearestZ,
                         farthestZ,
                         &hit_info)) {

            bool is_front_facing = hit_info.is_front_facing;

            if ((RT_FACE_cull_front == cull_mode &&  is_front_facing) ||
                (RT_FACE_cull_back  == cull_mode && !is_front_facing)) {

                continue;
            }

            if (-1 == closest_hit_index ||
                hit_info.t < closest_hit_info.t) {

                closest_hit_info    = hit_info;
                closest_hit_index   = i;
            }
        }
    }

    if (-1 != closest_hit_index) {

        info->info              = closest_hit_info;
        info->geometry_type     = RT_HIT_plane;
        info->geometry_index    = closest_hit_index;

        return true;
    }

    info->geometry_type         = RT_HIT_null;
    info->geometry_index        = 0;

    return false;
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_any_closest_hit(const rt_world_t*     world,
                                     rt_ray_t              ray,
                                     rt_float_t            nearestZ,
                                     rt_float_t            farthestZ,
                                     rt_hit_ext_info_t*    info)
{
    bool is_hit_bool_array[2]           = {};
    rt_hit_ext_info_t hit_info_array[2] = {};

    RT_ASSERT(RT_BUFFER_LEN(is_hit_bool_array) == RT_BUFFER_LEN(hit_info_array));

    is_hit_bool_array[0] = rt_world_sphere_closest_hit(world,
                                                       ray,
                                                       nearestZ,
                                                       farthestZ,
                                                       &hit_info_array[0]);

    is_hit_bool_array[1] = rt_world_plane_closest_hit(world,
                                                      ray,
                                                      nearestZ,
                                                      farthestZ,
                                                      &hit_info_array[1]);

    bool is_hit_any     = false;
    uint32_t hit_idx    = 0;

    for (uint32_t i = 0; i < RT_BUFFER_LEN(is_hit_bool_array); ++i) {

        if (is_hit_bool_array[i]) {

            is_hit_any  = true;
            hit_idx     = i;
        }
    }

    if (is_hit_any) {
        *info = hit_info_array[hit_idx];
    }

    return is_hit_any;
}

///////////////////////////////////////////////////////////////////////////
//////////////////////// FRAGMENT SHADERS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_fragment_shader_emissive(rt_emissive_material_t* material)
{
    RT_ASSERT(material != NULL);
    return material->color;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_fragment_shader_diffuse(const rt_hit_info_t*            hit_info,
                                            const rt_diffuse_material_t*    material,
                                            const rt_world_t*               world,
                                            bool                            is_in_shadow)
{
    RT_ASSERT(hit_info  != NULL);
    RT_ASSERT(material  != NULL);
    RT_ASSERT(world     != NULL);

    if (is_in_shadow && material->receives_shadows) {
        return material->ambient;
    }

    rt_vec4_t final_color = { 0.0f, 0.0f, 0.0f, 1.0f };

    const rt_idx_t directional_light_count              = world->directional_light_count;
    const rt_directional_light_t* directional_lights    = world->directional_light_buffer;

    for (rt_idx_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];
        RT_ASSERT(light != NULL);

        rt_vec4_t light_color           = rt_vec4_mul_scalar(light->color,
                                                             light->intensity);
        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_negate(light->direction));

        rt_vec4_t diffuse_color         = rt_vec4_mul(material->diffuse,
                                                      light_color);
        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), 0.0f);
        rt_vec4_t diffuse               = rt_vec4_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      diffuse);
    }

    rt_idx_t point_light_count          = world->point_light_count;
    const rt_point_light_t* point_lights= world->point_light_buffer;

    for (rt_idx_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light   = &point_lights[i];

        rt_vec4_t light_color           = rt_vec4_mul_scalar(light->color,
                                                             light->intensity);
        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_sub(light->position,
                                                                   hit_info->position));

        rt_vec4_t diffuse_color         = rt_vec4_mul(material->diffuse,
                                                      light_color);
        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), 0.0f);
        rt_vec4_t diffuse               = rt_vec4_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      diffuse);
    }

    final_color = rt_vec4_add(final_color, material->ambient);
    return rt_vec4_clamp(final_color, 0.0f, 1.0f);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_fragment_shader_checkerboard(const rt_hit_info_t*               hit_info,
                                                 const rt_checkerboard_material_t*  material,
                                                 bool                               is_in_shadow)
{
    RT_ASSERT(hit_info != NULL);
    RT_ASSERT(material != NULL);

    uint32_t pos_x_quant    = (uint32_t)(ceil(hit_info->position.x * 0.5f));
    uint32_t pos_z_quant    = (uint32_t)(ceil(hit_info->position.z * 0.5f));
    uint32_t pos_quant      = pos_x_quant + pos_z_quant;

    rt_vec4_t material_color = (pos_quant & 1) ? material->color_0
                                               : material->color_1;

    if (is_in_shadow && material->receives_shadows) {

        rt_vec4_t final_color = rt_vec4_mul_scalar(material_color,
                                                   material->shadow_factor);

        final_color.w = 1.0f;

        return final_color;
    }

    return material_color;
}

///////////////////////////////////////////////////////////////////////////
////////////////////////////// FINAL COLOR ////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_is_in_shadow(const rt_world_t*           world,
                            rt_hit_ext_info_t*          hit_info,
                            rt_float_t                  nearestZ,
                            rt_float_t                  farthestZ)
{
    RT_ASSERT(world     != NULL);
    RT_ASSERT(hit_info  != NULL);
    RT_ASSERT(nearestZ  <= farthestZ);

    rt_idx_t directional_light_count                    = world->directional_light_count;
    const rt_directional_light_t* directional_lights    = world->directional_light_buffer;

    for (rt_idx_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];
        RT_ASSERT(light != NULL);

        rt_ray_t new_ray = {
            .org = rt_vec4_add(hit_info->info.position,
                               rt_vec4_mul_scalar(hit_info->info.normal, RT_SHADOW_BIAS)),
            .dir = rt_vec4_norm(rt_vec4_negate(light->direction)),
        };

        if (rt_world_any_closest_hit(world,
                                     new_ray,
                                     nearestZ,
                                     farthestZ,
                                     hit_info)) {

            return true;
        }
    }
    
    rt_idx_t point_light_count              = world->point_light_count;
    const rt_point_light_t* point_lights    = world->point_light_buffer;

    for (rt_idx_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light = &point_lights[i];
        RT_ASSERT(light != NULL);

        rt_vec4_t org = rt_vec4_add(hit_info->info.position,
                                    rt_vec4_mul_scalar(hit_info->info.normal, RT_SHADOW_BIAS));

        rt_ray_t new_ray = {
            .org = org,
            .dir = rt_vec4_norm(rt_vec4_sub(light->position, org)),
        };

        rt_float_t newNearestZ  = RT_SHADOW_BIAS;
        rt_float_t newFarthestZ = rt_vec4_dist(org, light->position) - RT_SHADOW_BIAS;

        if (rt_world_any_closest_hit(world,
                                     new_ray,
                                     newNearestZ,
                                     newFarthestZ,
                                     hit_info)) {

            return true;
        }
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_world_compute_color(const rt_world_t*       world,
                                        rt_ray_t                ray,
                                        rt_float_t              nearestZ,
                                        rt_float_t              farthestZ,
                                        rt_idx_t                depth)
{
    RT_ASSERT(world     != NULL);
    RT_ASSERT(nearestZ  <= farthestZ);
    RT_ASSERT(depth     >= 0);

    rt_hit_ext_info_t hit_info = {};

    if (!rt_world_any_closest_hit(world,
                                  ray,
                                  nearestZ,
                                  farthestZ,
                                  &hit_info)) {

        return world->clear_color;;
    }

    enum rt_material_type material_type     = RT_MATERIAL_TYPE_null_material;
    rt_idx_t material_index                 = 0;

    enum rt_hit_geometry_type geometry_type = hit_info.geometry_type;
    RT_ASSERT(geometry_type                != RT_HIT_null);

    rt_idx_t geometry_index                 = hit_info.geometry_index;
    RT_ASSERT(geometry_index               >= 0);

    switch (geometry_type) {
        case RT_HIT_sphere:

            RT_ASSERT(world->sphere_buffer != NULL);

            material_type   = world->sphere_buffer[geometry_index].material_type;
            material_index  = world->sphere_buffer[geometry_index].material_index;

            break;
        case RT_HIT_plane:

            RT_ASSERT(world->plane_buffer != NULL);

            material_type   = world->plane_buffer[geometry_index].material_type;
            material_index  = world->plane_buffer[geometry_index].material_index;

            break;
        default:
            RT_ASSERT(0);
            break;
    }

    if (RT_MATERIAL_TYPE_null_material == material_type) {

        return world->clear_color;
    }

    bool is_in_shadow = rt_is_in_shadow(world,
                                        &hit_info,
                                        nearestZ,
                                        farthestZ);

    switch (material_type) {
        case RT_MATERIAL_TYPE_emissive_material: {
            rt_emissive_material_t* m = &world->emissive_material_buffer[material_index];
            return rt_fragment_shader_emissive(m);
        }
        case RT_MATERIAL_TYPE_checkerboard_material: {
            rt_checkerboard_material_t* m = &world->checkerboard_material_buffer[material_index];
            return rt_fragment_shader_checkerboard(&hit_info.info,
                                                   m,
                                                   is_in_shadow);
        }
        case RT_MATERIAL_TYPE_diffuse_material: {
            rt_diffuse_material_t* m = &world->diffuse_material_buffer[material_index];
            return rt_fragment_shader_diffuse(&hit_info.info,
                                              m,
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

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    uint32_t* rgb_buffer;
    uint32_t width;
    uint32_t height;
}
rt_framebuffer_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_status_t rt_framebuffer_create(uint32_t           width,
                                         uint32_t           height,
                                         rt_framebuffer_t*  framebuffer)
{
    RT_ASSERT(framebuffer != NULL);

    framebuffer->rgb_buffer = (uint32_t*)RT_ALLOC(width *
                                                  height *
                                                  sizeof(uint32_t));
    if (!framebuffer->rgb_buffer) {
        return RT_STATUS_failure;
    }

    framebuffer->width  = width;
    framebuffer->height = height;

    return RT_STATUS_success;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_status_t rt_framebuffer_resize(uint32_t           new_width,
                                         uint32_t           new_height,
                                         rt_framebuffer_t*  framebuffer)
{
    RT_ASSERT(framebuffer != NULL);

    uint32_t* rgb_buffer    = (uint32_t*)RT_ALLOC(new_width *
                                                  new_height *
                                                  sizeof(uint32_t));

    if (!rgb_buffer) {
        return RT_STATUS_failure;
    }

    RT_FREE(framebuffer->rgb_buffer);

    framebuffer->rgb_buffer = rgb_buffer;

    framebuffer->width      = new_width;
    framebuffer->height     = new_height;

    return RT_STATUS_success;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_framebuffer_write(uint32_t               row,
                                 uint32_t               col,
                                 rt_vec4_t              color,
                                 rt_framebuffer_t*      framebuffer)
{
    RT_ASSERT(framebuffer   != NULL);
    RT_ASSERT(row            < framebuffer->height);
    RT_ASSERT(col            < framebuffer->width);

    uint32_t index = row * framebuffer->width + col;
    framebuffer->rgb_buffer[index] = rt_vec4_to_uint32_alpha(color, 1.0f);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_framebuffer_free(rt_framebuffer_t* framebuffer)
{
    RT_ASSERT(framebuffer != NULL);

    if (framebuffer->rgb_buffer) {
        RT_FREE(framebuffer->rgb_buffer);
        framebuffer->rgb_buffer = NULL;
    }
}

#endif // RT_RAYTERM_H
