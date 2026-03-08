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
////////////////////////////////// API ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define RT_API                          static inline

///////////////////////////////////////////////////////////////////////////
///////////////////////////////// TYPES ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
#ifdef RT_USE_FLOAT64

#define RT_FLOAT(X)                     X
typedef double rt_float_t;

#else

#define RT_FLOAT(X)                     X ## f
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
////////////////////////////// CONSTANTS //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define RT_GAMMA                        RT_FLOAT(2.2)
#define RT_GAMMA_INVERSE                RT_FLOAT(0.454545)
#define RT_PI                           RT_FLOAT(3.1415926)
#define RT_EPSILON                      RT_FLOAT(0.000001)
#define RT_SHADOW_BIAS                  RT_FLOAT(0.001)
#define RT_INIT_CAP                     8

///////////////////////////////////////////////////////////////////////////
//////////////////////////////// ASSERTIONS ///////////////////////////////
///////////////////////////////////////////////////////////////////////////

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
}
rt_sphere_params_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_sphere_params_t      geometry_params;
    rt_idx_t                material_index;
    enum rt_material_type   material_type;
}
rt_sphere_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t               position;
    rt_vec4_t               normal;
    rt_float_t              side_length;
}
rt_plane_params_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_plane_params_t       geometry_params;
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

    rt_float_t sphere_radius = sphere.geometry_params.radius;
    rt_vec4_t sphere_center  = sphere.geometry_params.center;

    rt_vec4_t oc    = rt_vec4_sub(sphere.geometry_params.center, ray.org);
    rt_float_t a    = rt_vec4_sqrlen(ray.dir);
    rt_float_t h    = rt_vec4_dot(ray.dir, oc);
    rt_float_t c    = rt_vec4_sqrlen(oc) - sphere_radius * sphere_radius;

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

    rt_vec4_t n             = rt_vec4_sub(info->position, sphere_center);
    n                       = rt_vec4_div_scalar(n, sphere_radius);

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

    rt_vec4_t plane_position        = plane.geometry_params.position;
    rt_vec4_t plane_normal          = plane.geometry_params.normal;
    rt_float_t plane_side_length    = plane.geometry_params.side_length;

    rt_vec4_t negated_plane_normal = rt_vec4_negate(plane_normal);

    rt_float_t denom = rt_vec4_dot(negated_plane_normal,
                                   ray.dir);

    if (denom < RT_EPSILON) {
        return false;
    }

    rt_vec4_t dir = rt_vec4_sub(plane_position,
                                ray.org);

    rt_float_t t = rt_vec4_dot(dir,
                               negated_plane_normal) / denom;

    rt_vec4_t hit_position = rt_ray_at(ray, t);

    rt_float_t diff_x = fabs(hit_position.x - plane_position.x);
    rt_float_t diff_y = fabs(hit_position.y - plane_position.y);
    rt_float_t diff_z = fabs(hit_position.z - plane_position.z);

    if (diff_x > plane_side_length ||
        diff_y > plane_side_length ||
        diff_z > plane_side_length) {

        return false;
    }

    if (t < nearestZ || t > farthestZ) {
        return false;
    }

    info->position          = rt_ray_at(ray, t);
    info->normal            = negated_plane_normal;
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
RT_API void rt_world_set_directional_light_params(rt_world_t*                    world,
                                                  rt_idx_t                       light_index,
                                                  const rt_directional_light_t*  light_params)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(light_index   >= 0);
    RT_ASSERT(light_params  != NULL);

    rt_directional_light_t* light_buffer = world->directional_light_buffer;
    RT_ASSERT(light_buffer != NULL);

    rt_directional_light_t* light = &light_buffer[light_index];
    RT_ASSERT(light != NULL);

    memcpy(light, light_params, sizeof(rt_directional_light_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_point_light_params(rt_world_t*              world,
                                            rt_idx_t                 light_index,
                                            const rt_point_light_t*  light_params)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(light_index   >= 0);
    RT_ASSERT(light_params  != NULL);

    rt_point_light_t* light_buffer = world->point_light_buffer;
    RT_ASSERT(light_buffer != NULL);

    rt_point_light_t* light = &light_buffer[light_index];
    RT_ASSERT(light != NULL);

    memcpy(light, light_params, sizeof(rt_point_light_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_sphere_params(rt_world_t*                  world,
                                       rt_idx_t                     sphere_index,
                                       const rt_sphere_params_t*    sphere_params)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(sphere_index  >= 0);
    RT_ASSERT(sphere_params != NULL);

    rt_sphere_t* sphere_buffer = world->sphere_buffer;
    RT_ASSERT(sphere_buffer != NULL);

    rt_sphere_t* sphere = &sphere_buffer[sphere_index];
    RT_ASSERT(sphere != NULL);

    memcpy(&sphere->geometry_params,
           sphere_params,
           sizeof(rt_sphere_params_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_plane_params(rt_world_t*               world,
                                      rt_idx_t                  plane_index,
                                      const rt_plane_params_t*  plane_params)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(plane_index   >= 0);
    RT_ASSERT(plane_params  != NULL);

    rt_plane_t* plane_buffer = world->plane_buffer;
    RT_ASSERT(plane_buffer != NULL);

    rt_plane_t* plane = &plane_buffer[plane_index];
    RT_ASSERT(plane != NULL);

    memcpy(&plane->geometry_params,
           plane_params,
           sizeof(rt_plane_params_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_emissive_material_params(rt_world_t*                    world,
                                                  rt_idx_t                       material_index,
                                                  const rt_emissive_material_t*  material_params)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);
    RT_ASSERT(material_params   != NULL);

    rt_emissive_material_t* material_buffer = world->emissive_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_emissive_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    memcpy(material, material_params, sizeof(rt_emissive_material_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_checkerboard_material_params(rt_world_t*                        world,
                                                      rt_idx_t                           material_index,
                                                      const rt_checkerboard_material_t*  material_params)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);
    RT_ASSERT(material_params   != NULL);

    rt_checkerboard_material_t* material_buffer = world->checkerboard_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_checkerboard_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    memcpy(material, material_params, sizeof(rt_checkerboard_material_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_diffuse_material_params(rt_world_t*                     world,
                                                 rt_idx_t                        material_index,
                                                 const rt_diffuse_material_t*    material_params)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);
    RT_ASSERT(material_params   != NULL);

    rt_diffuse_material_t* material_buffer = world->diffuse_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_diffuse_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    memcpy(material, material_params, sizeof(rt_diffuse_material_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_metallic_material_params(rt_world_t*                    world,
                                                  rt_idx_t                       material_index,
                                                  const rt_metallic_material_t*  material_params)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);
    RT_ASSERT(material_params   != NULL);

    rt_metallic_material_t* material_buffer = world->metallic_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_metallic_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    memcpy(material, material_params, sizeof(rt_metallic_material_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_dielectric_material_params(rt_world_t*                      world,
                                                    rt_idx_t                         material_index,
                                                    const rt_dielectric_material_t*  material_params)
{
    RT_ASSERT(world             != NULL);
    RT_ASSERT(material_index    >= 0);
    RT_ASSERT(material_params   != NULL);

    rt_dielectric_material_t* material_buffer = world->dielectric_material_buffer;
    RT_ASSERT(material_buffer != NULL);

    rt_dielectric_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material != NULL);

    memcpy(material, material_params, sizeof(rt_dielectric_material_t));
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

    rt_idx_t hit_idx        = -1;
    rt_float_t closest_z    = RT_FLOAT(0.0);

    for (uint32_t i = 0; i < RT_BUFFER_LEN(is_hit_bool_array); ++i) {

        if (is_hit_bool_array[i]) {

            if (-1 == hit_idx || hit_info_array[i].info.t < closest_z) {

                closest_z   = hit_info_array[i].info.t;
                hit_idx     = (rt_idx_t)i;

            }
        }
    }

    if (hit_idx >= 0) {
        memcpy(info, &hit_info_array[hit_idx], sizeof(rt_hit_ext_info_t));
    }

    return hit_idx >= 0;
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

    rt_vec4_t final_color = { RT_FLOAT(0.0),
                              RT_FLOAT(0.0),
                              RT_FLOAT(0.0),
                              RT_FLOAT(1.0) };

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
RT_API rt_vec4_t rt_fragment_shader_metallic(const rt_hit_info_t*           hit_info,
                                             const rt_metallic_material_t*  material,
                                             const rt_world_t*              world,
                                             bool                           is_in_shadow)
{
    RT_ASSERT(hit_info  != NULL);
    RT_ASSERT(material  != NULL);
    RT_ASSERT(world     != NULL);

    if (is_in_shadow && material->receives_shadows) {
        return material->ambient;
    }

    rt_vec4_t final_color = { RT_FLOAT(0.0),
                              RT_FLOAT(0.0),
                              RT_FLOAT(0.0),
                              RT_FLOAT(0.0) };

    const rt_idx_t directional_light_count              = world->directional_light_count;
    const rt_directional_light_t* directional_lights    = world->directional_light_buffer;

    for (rt_idx_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];
        RT_ASSERT(light != NULL);

        rt_vec4_t light_color           = rt_vec4_mul_scalar(light->color,
                                                             light->intensity);

        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_negate(light->direction));

        rt_vec4_t specular_color        = rt_vec4_mul(material->specular,
                                                      light_color);

        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), 0.0f);

        rt_vec4_t specular              = rt_vec4_mul_scalar(specular_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      specular);
    }

    rt_idx_t point_light_count          = world->point_light_count;
    const rt_point_light_t* point_lights= world->point_light_buffer;

    for (rt_idx_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light   = &point_lights[i];
        RT_ASSERT(light != NULL);

        rt_vec4_t light_color           = rt_vec4_mul_scalar(light->color,
                                                             light->intensity);
        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_sub(light->position,
                                                                   hit_info->position));

        rt_vec4_t specular_color        = rt_vec4_mul(material->specular,
                                                      light_color);
        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), 0.0f);
        rt_vec4_t specular              = rt_vec4_mul_scalar(specular_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      specular);
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

        return final_color;
    }

    return material_color;
}

///////////////////////////////////////////////////////////////////////////
////////////////////////////// FINAL COLOR ////////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_should_be_in_shadow(const rt_world_t*         world,
                                   const rt_hit_ext_info_t*  hit_info)
{
        enum rt_material_type material_type = RT_MATERIAL_TYPE_null_material;
        rt_idx_t material_index             = -1;

        switch (hit_info->geometry_type) {
            case RT_HIT_sphere: {

                rt_sphere_t* sphere = &world->sphere_buffer[hit_info->geometry_index];
                RT_ASSERT(sphere != NULL);

                material_type   = sphere->material_type;
                material_index  = sphere->material_index;

                break;
            }
            case RT_HIT_plane: {

                rt_plane_t* plane = &world->plane_buffer[hit_info->geometry_index];
                RT_ASSERT(plane != NULL);

                material_type   = plane->material_type;
                material_index  = plane->material_index;

                break;
            }
            default: {
                RT_ASSERT(0);
                break;
            }
        }

        return !(material_index == -1                                 ||
                 RT_MATERIAL_TYPE_null_material == material_type      ||
                 RT_MATERIAL_TYPE_emissive_material == material_type);
}

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

            return rt_should_be_in_shadow(world, hit_info);
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

            return rt_should_be_in_shadow(world, hit_info);
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

    rt_hit_ext_info_t hit_info_copy = hit_info;;

    bool is_in_shadow = rt_is_in_shadow(world,
                                        &hit_info_copy,
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
        case RT_MATERIAL_TYPE_metallic_material: {

            rt_vec4_t reflected_ray_dir = rt_vec4_reflect(ray.dir,
                                                          hit_info.info.normal);

            rt_vec4_t reflected_ray_pos = rt_vec4_add(hit_info.info.position,
                                                      rt_vec4_mul_scalar(hit_info.info.normal,
                                                                         RT_SHADOW_BIAS));

            rt_ray_t reflected_ray = { .org = reflected_ray_pos,
                                       .dir = reflected_ray_dir };

            rt_vec4_t reflected_color = rt_world_compute_color(world,
                                                               reflected_ray,
                                                               nearestZ,
                                                               farthestZ,
                                                               depth - 1);

            rt_metallic_material_t* m = &world->metallic_material_buffer[material_index];

            rt_vec4_t fragment_shader_output =  rt_fragment_shader_metallic(&hit_info.info,
                                                                            m,
                                                                            world,
                                                                            is_in_shadow);

            return rt_vec4_mul(reflected_color, fragment_shader_output);
        }
        /*
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

///////////////////////////////////////////////////////////////////////////
//////////////////////////////// FPS CAMERA ///////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_vec4_t   y_axis;
    rt_vec4_t   z_axis;

    rt_vec4_t   position;
    rt_vec4_t   rotation;
    rt_vec4_t   rotation_new;

    rt_vec4_t   velocity_x;
    rt_vec4_t   velocity_y;
    rt_vec4_t   velocity_z;

    rt_vec4_t   velocity_new_x;
    rt_vec4_t   velocity_new_y;
    rt_vec4_t   velocity_new_z;

    rt_float_t  fovy;
    rt_float_t  sensitivity;
    rt_float_t  smoothing;
    rt_float_t  movement_speed;
    rt_float_t  rotation_speed;

    rt_float_t  pitch_delta;
    rt_float_t  yaw_delta;

    rt_idx_t    depth;

    rt_float_t  near;
    rt_float_t  far;
}
rt_fps_camera_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_fps_camera_t rt_fps_camera_create()
{
    rt_fps_camera_t camera = {

        .y_axis         = { RT_FLOAT(0.0),
                            RT_FLOAT(1.0),
                            RT_FLOAT(0.0),
                            RT_FLOAT(0.0), },

        .z_axis         = { RT_FLOAT(0.0),
                            RT_FLOAT(0.0),
                            RT_FLOAT(1.0),
                            RT_FLOAT(0.0), },

        .position       = { RT_FLOAT( 0.0),
                            RT_FLOAT( 0.0),
                            RT_FLOAT(-5.0),
                            RT_FLOAT( 1.0), },

        .rotation       = {},
        .rotation_new   = {},

        .velocity_x     = {},
        .velocity_y     = {},
        .velocity_z     = {},

        .velocity_new_x = {},
        .velocity_new_y = {},
        .velocity_new_z = {},

        .fovy           = RT_PI * RT_FLOAT(0.5),

        .sensitivity    = RT_FLOAT(0.1),
        .smoothing      = RT_FLOAT(0.01),
        .movement_speed = RT_FLOAT(0.5),
        .rotation_speed = RT_FLOAT(0.05),

        .pitch_delta    = RT_FLOAT(0.0),
        .yaw_delta      = RT_FLOAT(0.0),

        .depth          = 3,

        .near           = RT_FLOAT(0.3),
        .far            = RT_FLOAT(1000.0),
    };

    return camera;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_forward(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->velocity_new_z = rt_vec4_negate(rt_vec4_mul_scalar(camera->z_axis,
                                                               camera->movement_speed));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_forward(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    rt_vec4_zero(&camera->velocity_new_z);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_backward(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->velocity_new_z = rt_vec4_mul_scalar(camera->z_axis,
                                                camera->movement_speed);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_backward(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    rt_vec4_zero(&camera->velocity_new_z);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_left(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    rt_vec4_t strafe = rt_vec4_mul_scalar(rt_vec4_norm(rt_vec4_cross(camera->y_axis,
                                                                     camera->z_axis)), camera->movement_speed);

    camera->velocity_new_x = rt_vec4_negate(strafe);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_left(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    rt_vec4_zero(&camera->velocity_new_x);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_right(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    rt_vec4_t strafe = rt_vec4_mul_scalar(rt_vec4_norm(rt_vec4_cross(camera->y_axis,
                                                                     camera->z_axis)), camera->movement_speed);

    camera->velocity_new_x = strafe;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_right(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    rt_vec4_zero(&camera->velocity_new_x);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_left(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->yaw_delta = camera->rotation_speed;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_left(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->yaw_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_right(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->yaw_delta = -camera->rotation_speed;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_right(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->yaw_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_down(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->pitch_delta = -camera->rotation_speed;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_down(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->pitch_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_up(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->pitch_delta = camera->rotation_speed;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_up(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->pitch_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_up(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->velocity_new_y.y = camera->movement_speed;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_up(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->velocity_new_y.y = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_down(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->velocity_new_y.y = -camera->movement_speed;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_down(rt_fps_camera_t* camera)
{
    RT_ASSERT(camera != NULL);

    camera->velocity_new_y.y = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_render(rt_fps_camera_t*       camera,
                                 const rt_world_t*      world,
                                 rt_framebuffer_t*      framebuffer,
                                 rt_float_t             delta_time)
{
    RT_ASSERT(camera            != NULL);
    RT_ASSERT(world             != NULL);
    RT_ASSERT(framebuffer       != NULL);

    uint32_t cols               = framebuffer->width;
    uint32_t rows               = framebuffer->height;

    rt_float_t aspect           = (rt_float_t)cols / (rt_float_t)(rows * 2);

    rt_vec4_t camera_dir        = { RT_FLOAT(0.0),
                                    RT_FLOAT(0.0),
                                    RT_FLOAT(1.0),
                                    RT_FLOAT(0.0) };

    camera->rotation_new.x      += camera->pitch_delta;
    camera->rotation_new.x      = rt_clamp(camera->rotation_new.x,
                                           RT_TO_RADIANS(RT_FLOAT(-89.0)),
                                           RT_TO_RADIANS(RT_FLOAT( 89.0)));

    camera->rotation_new.y      += camera->yaw_delta;

    camera->rotation            = rt_vec4_lerp(camera->rotation,
                                               camera->rotation_new,
                                               camera->smoothing * delta_time);

    camera_dir                  = rt_vec4_rotate_x(camera_dir, camera->rotation.x);
    camera_dir                  = rt_vec4_rotate_y(camera_dir, camera->rotation.y);

    camera->z_axis              = camera_dir;

    camera->velocity_z          = rt_vec4_lerp(camera->velocity_z,
                                               camera->velocity_new_z,
                                               camera->smoothing * delta_time);

    camera->velocity_x          = rt_vec4_lerp(camera->velocity_x,
                                               camera->velocity_new_x,
                                               camera->smoothing * delta_time);

    camera->velocity_y          = rt_vec4_lerp(camera->velocity_y,
                                               camera->velocity_new_y,
                                               camera->smoothing * delta_time);

    rt_vec4_t total_velocity    = rt_vec4_add(camera->velocity_z,
                                              camera->velocity_x);

    total_velocity              = rt_vec4_add(total_velocity,
                                              camera->velocity_y);

    camera->position            = rt_vec4_add(camera->position,
                                              total_velocity);

    rt_vec4_t camera_look_at    = rt_vec4_sub(camera->position,
                                              camera_dir);

    rt_float_t focal_length     = rt_vec4_dist(camera->position,
                                               camera_look_at);

    rt_float_t camera_h         = tanf(camera->fovy * RT_FLOAT(0.5));

    rt_float_t viewport_height  = RT_FLOAT(2.0) * camera_h * focal_length;

    rt_float_t viewport_width   = viewport_height * aspect;

    rt_vec4_t camera_w          = rt_vec4_norm(rt_vec4_sub(camera->position,
                                                           camera_look_at));

    rt_vec4_t camera_u          = rt_vec4_norm(rt_vec4_cross(camera->y_axis,
                                                             camera_w));

    rt_vec4_t camera_v          = rt_vec4_norm(rt_vec4_cross(camera_w,
                                                             camera_u));

    rt_vec4_t viewport_u        = rt_vec4_mul_scalar(camera_u,
                                                     viewport_width);

    rt_vec4_t viewport_v        = rt_vec4_mul_scalar(rt_vec4_negate(camera_v),
                                                     viewport_height);

    rt_vec4_t pixel_delta_u     = rt_vec4_div_scalar(viewport_u, (rt_float_t)cols);

    rt_vec4_t pixel_delta_v     = rt_vec4_div_scalar(viewport_v, (rt_float_t)rows);

    rt_vec4_t pixel_delta_diag  = rt_vec4_mul_scalar(rt_vec4_add(pixel_delta_u,
                                                                 pixel_delta_v),
                                                     RT_FLOAT(0.5));

    rt_vec4_t viewport_u_half   = rt_vec4_mul_scalar(viewport_u,
                                                     RT_FLOAT(0.5));

    rt_vec4_t viewport_v_half   = rt_vec4_mul_scalar(viewport_v,
                                                     RT_FLOAT(0.5));

    rt_vec4_t viewport_upper_left = rt_vec4_sub(camera->position,
                                                rt_vec4_mul_scalar(camera_w,
                                                                   focal_length));

    viewport_upper_left         = rt_vec4_sub(viewport_upper_left,
                                              viewport_u_half);

    viewport_upper_left         = rt_vec4_sub(viewport_upper_left,
                                              viewport_v_half);

    rt_vec4_t pixel00_loc       = rt_vec4_add(viewport_upper_left,
                                              pixel_delta_diag);

    for (uint32_t row = 0; row < rows; ++row) {

        rt_vec4_t ver = rt_vec4_mul_scalar(pixel_delta_v, (rt_float_t)row);

        for (uint32_t col = 0; col < cols; ++col) {

            rt_vec4_t hor           = rt_vec4_mul_scalar(pixel_delta_u,
                                                         (rt_float_t)col);

            rt_vec4_t pixel_center  = rt_vec4_add(pixel00_loc,
                                                  rt_vec4_add(hor, ver));

            rt_vec4_t ray_dir       = rt_vec4_sub(pixel_center,
                                                  camera->position);

            rt_ray_t r = {
                .dir = ray_dir,
                .org = camera->position,
            };

            rt_vec4_t pixel_color = rt_world_compute_color(world,
                                                           r,
                                                           camera->near,
                                                           camera->far,
                                                           camera->depth);

            pixel_color = rt_vec4_apply_2(pixel_color,
                                          rt_apply_gamma_custom,
                                          RT_GAMMA_INVERSE);

            rt_framebuffer_write(row, col, pixel_color, framebuffer);
        }
    }
}

///////////////////////////////////////////////////////////////////////////
/////////////////////////// NOTCURSES SURFACE /////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef RT_USE_NOTCURSES

#endif

///////////////////////////////////////////////////////////////////////////
/////////////////////// SDL INPUT (FOR JOYSTICKS) /////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef RT_USE_SDL3

#endif

#endif // RT_RAYTERM_H
