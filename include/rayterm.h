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
#include <inttypes.h>
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
#define RT_SHADOW_BIAS                  RT_FLOAT(0.01)
#define RT_MAX_SHADOW_LIGHTS            8
#define RT_INIT_CAP                     8

///////////////////////////////////////////////////////////////////////////
//////////////////////////// ERROR CALLBACK ///////////////////////////////
///////////////////////////////////////////////////////////////////////////
typedef void (*rt_error_callback_t)(const char*             filename,
                                    uint32_t                line,
                                    const char*             function_name,
                                    const char*             message,
                                    void*                   userParam);

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
static void* rt_default_user_param = NULL;

///////////////////////////////////////////////////////////////////////////
RT_API void rt_default_error_callback(const char*             filename,
                                      uint32_t                line,
                                      const char*             function_name,
                                      const char*             message,
                                      [[maybe_unused]] void*  userParam)
{
    fprintf(stderr, "ERROR: %s:%" PRIu32 " in %s: %s\n", filename,
                                                         line,
                                                         function_name,
                                                         message);
    exit(EXIT_FAILURE);
}

///////////////////////////////////////////////////////////////////////////
rt_error_callback_t rt_error_callback = rt_default_error_callback;

///////////////////////////////////////////////////////////////////////////
RT_API void rt_set_error_callback(rt_error_callback_t   callback,
                                  void*                 userParam)
{
    rt_error_callback       = callback;
    rt_default_user_param   = userParam;
}

///////////////////////////////////////////////////////////////////////////
#define RT_ASSERT(EXPR)                                                     \
    do {                                                                    \
        if (!(EXPR)) {                                                      \
            rt_error_callback(__FILE__,                                     \
                              __LINE__,                                     \
                              __PRETTY_FUNCTION__,                          \
                              #EXPR " failed!",                             \
                              rt_default_user_param);                       \
        }                                                                   \
    }                                                                       \
    while (0);

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
#define tan tanf
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
#ifndef RT_USE_FLOAT64
#define exp expf
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
typedef struct
{
    rt_float_t last_time;
    rt_float_t total_time;
    rt_float_t delay_time;
}
rt_timer_t;

///////////////////////////////////////////////////////////////////////////
typedef void (*rt_timer_callback_t)(rt_timer_t*     timer,
                                    void*           user_param);

// #undef RT_USE_SDL3

///////////////////////////////////////////////////////////////////////////
#ifdef RT_USE_SDL3
#include <SDL3/SDL.h>
#else
#include <time.h>
#endif

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_timer_update(rt_timer_t*              timer,
                                  rt_timer_callback_t      callback,
                                  void*                    user_param)
{
#ifdef RT_USE_SDL3

    rt_float_t current_time     = (rt_float_t)SDL_GetTicks() * RT_FLOAT(0.001);

#else

    rt_float_t current_time     = (rt_float_t)clock() / (rt_float_t)CLOCKS_PER_SEC;

#endif

    rt_float_t delta_time       = (timer->last_time > RT_FLOAT(0.0)) ? (current_time - timer->last_time)
                                                                     : RT_FLOAT(0.0);
    timer->last_time            = current_time;
    timer->total_time           += delta_time;

    if (timer->total_time       >= timer->delay_time) {

        if (callback) {

            callback(timer, user_param);

        }

        timer->total_time       = RT_FLOAT(0.0);

    }

    return delta_time;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_timer_wait(rt_timer_t*   timer,
                          rt_float_t    target_fps)
{
#ifdef RT_USE_SDL3

    rt_float_t elapsed_time     = (rt_float_t)SDL_GetTicks() * RT_FLOAT(0.001) - timer->last_time;

#else

    rt_float_t elapsed_time     = ((rt_float_t)clock() / (rt_float_t)CLOCKS_PER_SEC) - timer->last_time;

#endif

    rt_float_t ms               = RT_FLOAT(1000.0) / target_fps;
    rt_float_t remaining_time   = ms - elapsed_time;

    if (remaining_time > RT_FLOAT(0.0)) {

#ifdef RT_USE_SDL3

        SDL_Delay((uint32_t)remaining_time);

#else

        struct timespec ts = {};

        ts.tv_sec   = (time_t)(remaining_time / RT_FLOAT(1000.0));
        ts.tv_nsec  = (long)((remaining_time - (rt_float_t)ts.tv_sec) * RT_FLOAT(1000.0));

        nanosleep(&ts, NULL);

#endif
    }
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
        .w = RT_FLOAT(0.0),
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
RT_API rt_vec4_t rt_vec4_max_vec4(rt_vec4_t p,
                                  rt_vec4_t q)
{
    p.x = fmax(p.x, q.x);
    p.y = fmax(p.y, q.y);
    p.z = fmax(p.z, q.z);
    p.w = fmax(p.w, q.w);

    return p;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_vec4_min_vec4(rt_vec4_t p,
                                  rt_vec4_t q)
{
    p.x = fmin(p.x, q.x);
    p.y = fmin(p.y, q.y);
    p.z = fmin(p.z, q.z);
    p.w = fmin(p.w, q.w);

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
    p.w  = RT_FLOAT(0.0);

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
    rt_float_t  ambient_factor;
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
    RT_HIT_aabb,
    RT_HIT_count,
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
typedef struct
{
    rt_vec4_t               min_bounding_box;
    rt_vec4_t               max_bounding_box;
}
rt_aabb_params_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_aabb_params_t        geometry_params;
    rt_idx_t                material_index;
    enum rt_material_type   material_type;
}
rt_aabb_t;

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_sphere_hit(rt_sphere_t               sphere,
                          rt_ray_t                  ray,
                          rt_float_t                nearest_z,
                          rt_float_t                farthest_z,
                          rt_hit_info_t*            info)
{
    RT_ASSERT(info          != NULL);
    RT_ASSERT(nearest_z     <= farthest_z);

    rt_float_t sphere_radius = sphere.geometry_params.radius;
    rt_vec4_t sphere_center  = sphere.geometry_params.center;

    rt_vec4_t oc            = rt_vec4_sub(sphere.geometry_params.center, ray.org);
    rt_float_t a            = rt_vec4_sqrlen(ray.dir);
    rt_float_t h            = rt_vec4_dot(ray.dir, oc);
    rt_float_t c            = rt_vec4_sqrlen(oc) - sphere_radius * sphere_radius;

    rt_float_t disc         = h * h - a * c;

    if (disc                < RT_FLOAT(0.0)) {
        return false;
    }

    rt_float_t sqrt_disc     = sqrt(disc);
    rt_float_t root_nearest  = (h - sqrt_disc) / a;
    rt_float_t root_farthest = (h + sqrt_disc) / a;

    rt_float_t range_min     = nearest_z;
    rt_float_t range_max     = farthest_z;

    rt_float_t root_chosen   = root_nearest;

    if (root_chosen         <= range_min ||
        root_chosen         >= range_max) {

        root_chosen = root_farthest;

        if (root_chosen     <= range_min ||
            root_chosen     >= range_max) {

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
                         rt_float_t             nearest_z,
                         rt_float_t             farthest_z,
                         rt_hit_info_t*         info)
{
    RT_ASSERT(info          != NULL);
    RT_ASSERT(nearest_z     <= farthest_z);

    rt_vec4_t plane_position        = plane.geometry_params.position;
    rt_vec4_t plane_normal          = plane.geometry_params.normal;
    rt_float_t plane_side_length    = plane.geometry_params.side_length;

    rt_float_t denom                = rt_vec4_dot(plane_normal,
                                                  ray.dir);

    if (denom < RT_EPSILON) {
        return false;
    }

    rt_vec4_t dir           = rt_vec4_sub(plane_position,
                                          ray.org);

    rt_float_t t            = rt_vec4_dot(dir,
                                          plane_normal) / denom;

    rt_vec4_t hit_position  = rt_ray_at(ray, t);

    rt_float_t diff_x       = fabs(hit_position.x - plane_position.x);
    rt_float_t diff_y       = fabs(hit_position.y - plane_position.y);
    rt_float_t diff_z       = fabs(hit_position.z - plane_position.z);

    if (diff_x              > plane_side_length ||
        diff_y              > plane_side_length ||
        diff_z              > plane_side_length) {

        return false;
    }

    if (t                   < nearest_z ||
        t                   > farthest_z) {
        return false;
    }

    info->position          = rt_ray_at(ray, t);
    info->normal            = rt_vec4_negate(plane_normal);
    info->t                 = t;
    info->is_front_facing   = true;

    return true;
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_aabb_hit(rt_aabb_t              aabb,
                        rt_ray_t               ray,
                        rt_float_t             nearest_z,
                        rt_float_t             farthest_z,
                        rt_hit_info_t*         info)
{
    RT_ASSERT(info      != NULL);
    RT_ASSERT(nearest_z <= farthest_z);

    rt_float_t tmin = (aabb.geometry_params.min_bounding_box.x - ray.org.x) / ray.dir.x;
    rt_float_t tmax = (aabb.geometry_params.max_bounding_box.x - ray.org.x) / ray.dir.x;

    if (tmin > tmax) {
        rt_float_t tmp = tmin;
        tmin = tmax;
        tmax = tmp;
    }

    rt_float_t tymin = (aabb.geometry_params.min_bounding_box.y - ray.org.y) / ray.dir.y;
    rt_float_t tymax = (aabb.geometry_params.max_bounding_box.y - ray.org.y) / ray.dir.y;

    if (tymin > tymax) {
        rt_float_t tmp = tymin;
        tymin = tymax;
        tymax = tmp;
    }

    if ((tmin > tymax) || (tymin > tmax)) {
        return false;
    }

    tmin = fmax(tmin, nearest_z);
    tmax = fmin(tmax, farthest_z);

    rt_float_t tzmin = (aabb.geometry_params.min_bounding_box.z - ray.org.z) / ray.dir.z;
    rt_float_t tzmax = (aabb.geometry_params.max_bounding_box.z - ray.org.z) / ray.dir.z;

    if (tzmin > tzmax) {
        rt_float_t tmp = tzmin;
        tzmin = tzmax;
        tzmax = tmp;
    }

    if ((tmin > tzmax) || (tzmin > tmax)) {
        return false;
    }

    tmin = fmax(tmin, nearest_z);
    tmax = fmin(tmax, farthest_z);

    if (tmin < nearest_z) {
        tmin = tmax;
    }

    if (tmin > farthest_z) {
        return false;
    }

    info->position          = rt_ray_at(ray, tmin);
    info->normal            = (rt_vec4_t){}; // TODO: compute normal
    info->t                 = tmin;
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
    rt_float_t  radius;
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
/////////////////////// ATMOSPHERIC SCATTERING ////////////////////////////
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    rt_float_t  planetary_radius;
    rt_float_t  atmospheric_height;
    rt_vec4_t   rayleigh_scattering;
    rt_float_t  mie_scattering;
    rt_float_t  phase_eccentricity;
    uint32_t    steps;
}
rt_atmosphere_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_atmosphere_t rt_atmosphere_create_default()
{
    rt_atmosphere_t default_atmosphere = {

        .planetary_radius       = RT_FLOAT(6371.0e3),

        .atmospheric_height     = RT_FLOAT(10.0e3),

        .rayleigh_scattering    = { RT_FLOAT(5.5e-6),
                                    RT_FLOAT(13.0e-6),
                                    RT_FLOAT(22.4e-6),
                                    RT_FLOAT(0.0), },

        .mie_scattering         = RT_FLOAT(21e-6),

        .phase_eccentricity     = RT_FLOAT(0.758),

        .steps                  = 2,

    };

    return default_atmosphere;
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_atmosphere_occludes(const rt_atmosphere_t*    atmosphere,
                                   rt_ray_t                  ray)
{
    RT_ASSERT(atmosphere   != NULL);

    rt_float_t radius       = atmosphere->planetary_radius;

    rt_vec4_t center        = { RT_FLOAT(0.0),
                               -radius,
                                RT_FLOAT(0.0),
                                RT_FLOAT(1.0) };

    rt_vec4_t offset        = rt_vec4_sub(center, ray.org);

    rt_float_t t            = rt_vec4_dot(ray.dir, offset);

    rt_vec4_t dir_scaled    = rt_vec4_mul_scalar(ray.dir, t);

    return t >= RT_FLOAT(0.0) &&
           rt_vec4_sqrlen(rt_vec4_sub(offset, dir_scaled)) < (radius * radius);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_float_t rt_atmosphere_trace(const rt_atmosphere_t* atmosphere,
                                      rt_ray_t               ray)
{
    RT_ASSERT(atmosphere   != NULL);

    rt_float_t radius       = atmosphere->planetary_radius;

    rt_float_t r_sky        = radius + atmosphere->atmospheric_height;

    rt_vec4_t center        = { RT_FLOAT(0.0),
                                -radius,
                                RT_FLOAT(0.0),
                                RT_FLOAT(1.0) };

    rt_vec4_t offset        = rt_vec4_sub(center, ray.org);

    rt_float_t t            = rt_vec4_dot(ray.dir, offset);

    rt_float_t discriminant = t * t - rt_vec4_dot(offset, offset) + r_sky * r_sky;

    if (discriminant < RT_FLOAT(0.0)) {

        return RT_FLOAT(0.0);

    }
    return t + sqrt(discriminant);
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_atmosphere_scatter(const rt_atmosphere_t*  atmosphere,
                                       rt_ray_t                ray,
                                       rt_vec4_t               irradiance,
                                       rt_vec4_t               direction)
{
    RT_ASSERT(atmosphere   != NULL);

    direction               = rt_vec4_norm(rt_vec4_negate(direction));

    rt_float_t g            = atmosphere->phase_eccentricity;

    rt_float_t dist_eye     = rt_atmosphere_trace(atmosphere, ray);

    rt_float_t cos_theta    = rt_vec4_dot(ray.dir, direction);

    rt_float_t rayleigh     = RT_FLOAT(3.0) / (RT_FLOAT(16.0) * RT_PI) *
                             (RT_FLOAT(1.0) + cos_theta * cos_theta);

    rt_float_t mie          = RT_FLOAT(1.0) / (RT_FLOAT(4.0) * RT_PI) *
                              pow(RT_FLOAT(1.0) - g, RT_FLOAT(2.0)) /
                              pow(RT_FLOAT(1.0) + g * g - RT_FLOAT(2.0) * g * cos_theta,
                              RT_FLOAT(1.5));

    uint32_t steps          = atmosphere->steps;

    rt_float_t dt           = dist_eye / (rt_float_t)steps;

    rt_vec4_t f             = rt_vec4_apply_1(
                                rt_vec4_mul_scalar(
                                    rt_vec4_negate(
                                        rt_vec4_add_scalar(atmosphere->rayleigh_scattering,
                                                           atmosphere->mie_scattering)), dt), exp);

    rt_vec4_t color         = { RT_FLOAT(0.0),
                                RT_FLOAT(0.0),
                                RT_FLOAT(0.0),
                                RT_FLOAT(1.0), };

    for (uint32_t i = 0; i < steps; ++i) {

        rt_float_t t            = dist_eye * (RT_FLOAT(1.0) - (((rt_float_t)i + RT_FLOAT(0.5)) / (rt_float_t)steps));

        rt_vec4_t orig          = rt_vec4_add(ray.org, rt_vec4_mul_scalar(ray.dir, t));

        rt_float_t sun_dist     = rt_atmosphere_trace(atmosphere, (rt_ray_t){.org = orig, .dir = direction} );

        rt_vec4_t l_in_right    = rt_vec4_apply_1(rt_vec4_mul_scalar(
                                                    rt_vec4_negate(
                                                        rt_vec4_add_scalar(atmosphere->rayleigh_scattering,
                                                                           atmosphere->mie_scattering)),
                                                    sun_dist),
                                                  exp);

        rt_vec4_t l_in_left     = rt_vec4_add_scalar(rt_vec4_mul_scalar(atmosphere->rayleigh_scattering,
                                                                        rayleigh),
                                                 atmosphere->mie_scattering * mie);

        rt_vec4_t l_in          = rt_vec4_mul(irradiance, l_in_left);

        l_in                    = rt_vec4_mul(l_in, l_in_right);

        l_in                    = rt_vec4_mul_scalar(l_in, dt);

        color                   = rt_vec4_mul(color, f);

        color                   = rt_vec4_add(color, l_in);
    }

    return color;
}

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
    rt_atmosphere_t atmosphere;

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
    RT_WORLD_DEF_BUFFER(aabb)

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
        RT_ASSERT(capacity_new > capacity);                                 \
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
RT_WORLD_DEFINE_PUSH(aabb)
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
RT_WORLD_DEFINE_POP(aabb)
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
RT_WORLD_DEFINE_RESERVE(aabb)
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
RT_WORLD_DEFINE_FREE(aabb)
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
    rt_world_free_aabbs                         (world);
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
    RT_ASSERT(geometry_buffer   != NULL);                                   \
                                                                            \
    geometry_type* geometry     = &geometry_buffer[geometry_index];         \
    RT_ASSERT(geometry          != NULL);                                   \
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
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(aabb, emissive_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(aabb, checkerboard_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(aabb, diffuse_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(aabb, metallic_material)
RT_DEFINE_LINK_GEOMETRY_TO_MATERIAL(aabb, dielectric_material)

///////////////////////////////////////////////////////////////////////////
#define RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(GEOMETRY, MATERIAL)           \
RT_API void RT_CONCAT4(rt_, GEOMETRY, _unlink_, MATERIAL)                   \
(rt_world_t* world, rt_idx_t geometry_index)                                \
{                                                                           \
    RT_ASSERT(world             != NULL);                                   \
    RT_ASSERT(geometry_index    >= 0);                                      \
    RT_ASSERT(geometry_index    < world->RT_CONCAT2(GEOMETRY, _count));     \
                                                                            \
    typedef RT_CONCAT3(rt_, GEOMETRY, _t) geometry_type;                    \
                                                                            \
    geometry_type* geometry_buffer  = world->RT_CONCAT2(GEOMETRY, _buffer); \
    RT_ASSERT(geometry_buffer   != NULL);                                   \
                                                                            \
    geometry_type* geometry     = &geometry_buffer[geometry_index];         \
    RT_ASSERT(geometry          != NULL);                                   \
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
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(aabb, emissive_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(aabb, checkerboard_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(aabb, diffuse_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(aabb, metallic_material)
RT_DEFINE_UNLINK_GEOMETRY_TO_MATERIAL(aabb, dielectric_material)

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_directional_light_params(rt_world_t*                    world,
                                                  rt_idx_t                       light_index,
                                                  const rt_directional_light_t*  light_params)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(light_index   >= 0);
    RT_ASSERT(light_params  != NULL);

    rt_directional_light_t* light_buffer = world->directional_light_buffer;
    RT_ASSERT(light_buffer  != NULL);

    rt_directional_light_t* light = &light_buffer[light_index];
    RT_ASSERT(light         != NULL);

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
    RT_ASSERT(light_buffer  != NULL);

    rt_point_light_t* light = &light_buffer[light_index];
    RT_ASSERT(light         != NULL);

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
    RT_ASSERT(sphere        != NULL);

    memcpy(&sphere->geometry_params, sphere_params, sizeof(rt_sphere_params_t));
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
    RT_ASSERT(plane_buffer  != NULL);

    rt_plane_t* plane = &plane_buffer[plane_index];
    RT_ASSERT(plane         != NULL);

    memcpy(&plane->geometry_params, plane_params, sizeof(rt_plane_params_t));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_world_set_aabb_params(rt_world_t*                world,
                                     rt_idx_t                   aabb_index,
                                     const rt_aabb_params_t*    aabb_params)
{
    RT_ASSERT(world         != NULL);
    RT_ASSERT(aabb_index    >= 0);
    RT_ASSERT(aabb_params   != NULL);

    rt_aabb_t* aabb_buffer  = world->aabb_buffer;
    RT_ASSERT(aabb_buffer  != NULL);

    rt_aabb_t* aabb = &aabb_buffer[aabb_index];
    RT_ASSERT(aabb         != NULL);

    memcpy(&aabb->geometry_params, aabb_params, sizeof(rt_aabb_params_t));
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
    RT_ASSERT(material_buffer   != NULL);

    rt_emissive_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material          != NULL);

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
    RT_ASSERT(material_buffer   != NULL);

    rt_checkerboard_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material          != NULL);

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
    RT_ASSERT(material_buffer   != NULL);

    rt_diffuse_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material          != NULL);

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
    RT_ASSERT(material_buffer   != NULL);

    rt_metallic_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material          != NULL);

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
    RT_ASSERT(material_buffer   != NULL);

    rt_dielectric_material_t* material = &material_buffer[material_index];
    RT_ASSERT(material          != NULL);

    memcpy(material, material_params, sizeof(rt_dielectric_material_t));
}

///////////////////////////////////////////////////////////////////////////
#define RT_WORLD_DEF_CLOSEST_HIT(GEOMETRY)                                  \
RT_API bool RT_CONCAT3(rt_world_, GEOMETRY, _closest_hit)                   \
(const rt_world_t*       world,                                             \
 rt_ray_t                ray,                                               \
 rt_float_t              nearest_z,                                         \
 rt_float_t              farthest_z,                                        \
 rt_hit_ext_info_t*      info)                                              \
{                                                                           \
    RT_ASSERT(world     != NULL);                                           \
    RT_ASSERT(nearest_z <= farthest_z);                                     \
    RT_ASSERT(info      != NULL);                                           \
                                                                            \
    enum rt_face_cull_mode cull_mode = world->face_cull_mode;               \
                                                                            \
    if (RT_FACE_cull_both == cull_mode) {                                   \
                                                                            \
        info->geometry_type         = RT_HIT_null;                          \
        info->geometry_index        = 0;                                    \
                                                                            \
        return false;                                                       \
    }                                                                       \
                                                                            \
    rt_hit_info_t closest_hit_info  = {};                                   \
    rt_idx_t closest_hit_index      = -1;                                   \
                                                                            \
    rt_idx_t RT_CONCAT2(GEOMETRY, _count) =                                 \
                                      world->RT_CONCAT2(GEOMETRY, _count);  \
                                                                            \
    RT_ASSERT(RT_CONCAT2(GEOMETRY, _count) >= 0);                           \
                                                                            \
    const RT_CONCAT3(rt_, GEOMETRY, _t)* RT_CONCAT2(GEOMETRY, _buffer) =    \
        world->RT_CONCAT2(GEOMETRY, _buffer);                               \
                                                                            \
    for (rt_idx_t i = 0; i < RT_CONCAT2(GEOMETRY, _count); ++i) {           \
                                                                            \
        rt_hit_info_t hit_info = {};                                        \
                                                                            \
        if (RT_CONCAT3(rt_, GEOMETRY, _hit)(RT_CONCAT2(GEOMETRY, _buffer)[i],\
                                            ray,                            \
                                            nearest_z,                      \
                                            farthest_z,                     \
                                            &hit_info)) {                   \
                                                                            \
            bool is_front_facing = hit_info.is_front_facing;                \
                                                                            \
            if ((RT_FACE_cull_front == cull_mode &&  is_front_facing) ||    \
                (RT_FACE_cull_back  == cull_mode && !is_front_facing)) {    \
                                                                            \
                continue;                                                   \
            }                                                               \
                                                                            \
            if (-1 == closest_hit_index ||                                  \
                hit_info.t < closest_hit_info.t) {                          \
                                                                            \
                closest_hit_info    = hit_info;                             \
                closest_hit_index   = i;                                    \
            }                                                               \
        }                                                                   \
    }                                                                       \
                                                                            \
    if (-1 != closest_hit_index) {                                          \
                                                                            \
        info->info              = closest_hit_info;                         \
        info->geometry_type     = RT_CONCAT2(RT_HIT_, GEOMETRY);            \
        info->geometry_index    = closest_hit_index;                        \
                                                                            \
        return true;                                                        \
    }                                                                       \
                                                                            \
    info->geometry_type         = RT_HIT_null;                              \
    info->geometry_index        = 0;                                        \
                                                                            \
    return false;                                                           \
}                                                                           \

///////////////////////////////////////////////////////////////////////////
RT_WORLD_DEF_CLOSEST_HIT(sphere)
RT_WORLD_DEF_CLOSEST_HIT(plane)
RT_WORLD_DEF_CLOSEST_HIT(aabb)

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_world_any_closest_hit(const rt_world_t*     world,
                                     rt_ray_t              ray,
                                     rt_float_t            nearest_z,
                                     rt_float_t            farthest_z,
                                     rt_hit_ext_info_t*    info)
{
    bool is_hit_bool_array[RT_HIT_count]            = {};
    rt_hit_ext_info_t hit_info_array[RT_HIT_count]  = {};

    RT_ASSERT(RT_BUFFER_LEN(is_hit_bool_array) == RT_BUFFER_LEN(hit_info_array));

    is_hit_bool_array[RT_HIT_sphere] = rt_world_sphere_closest_hit(world,
                                                                   ray,
                                                                   nearest_z,
                                                                   farthest_z,
                                                                   &hit_info_array[RT_HIT_sphere]);

    is_hit_bool_array[RT_HIT_plane] = rt_world_plane_closest_hit(world,
                                                                 ray,
                                                                 nearest_z,
                                                                 farthest_z,
                                                                 &hit_info_array[RT_HIT_plane]);

    is_hit_bool_array[RT_HIT_aabb] = rt_world_aabb_closest_hit(world,
                                                               ray,
                                                               nearest_z,
                                                               farthest_z,
                                                               &hit_info_array[RT_HIT_aabb]);

    rt_idx_t hit_idx        = -1;
    rt_float_t closest_z    = RT_FLOAT(0.0);

    for (uint32_t i = 0; i < RT_HIT_count; ++i) {

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
#define RT_COMPUTE_SHADOW_FOR_LIGHT(LIGHT,                              \
                                    MATERIAL,                           \
                                    LIGHT_INDEX,                        \
                                    LIGHT_TYPE,                         \
                                    SHADOW_LIGHT_INDEX,                 \
                                    SHADOW_LIGHT_TYPE)                  \
{                                                                       \
    if (LIGHT->casts_shadows &&                                         \
        MATERIAL->receives_shadows &&                                   \
        LIGHT_INDEX == SHADOW_LIGHT_INDEX &&                            \
        LIGHT_TYPE == SHADOW_LIGHT_TYPE) {                              \
                                                                        \
        continue;                                                       \
    }                                                                   \
}

///////////////////////////////////////////////////////////////////////////
#define RT_FIND_INDEX_OF_SHADOW_LIGHT(SHADOW_LIGHTS,                    \
                                      SHADOW_LIGHT_TYPES,               \
                                      SHADOW_LIGHT_COUNT,               \
                                      LIGHT_INDEX,                      \
                                      LIGHT_TYPE,                       \
                                      FOUND_INDEX,                      \
                                      FOUND_TYPE)                       \
{                                                                       \
    FOUND_INDEX = -1;                                                   \
    FOUND_TYPE  = RT_LIGHT_null_light;                                  \
                                                                        \
    for (rt_idx_t j = 0; j < SHADOW_LIGHT_COUNT; ++j) {                 \
                                                                        \
        if (SHADOW_LIGHTS[j] == LIGHT_INDEX &&                          \
            SHADOW_LIGHT_TYPES[j] == LIGHT_TYPE) {                      \
                                                                        \
            FOUND_INDEX = SHADOW_LIGHTS[j];                             \
            FOUND_TYPE  = SHADOW_LIGHT_TYPES[j];                        \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}                                                                       \

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_fragment_shader_diffuse(const rt_hit_info_t*            hit_info,
                                            const rt_diffuse_material_t*    material,
                                            const rt_world_t*               world,
                                            const rt_idx_t*                 shadow_lights,
                                            const enum rt_light_type*       shadow_light_types,
                                            rt_idx_t                        shadow_light_count)
{
    RT_ASSERT(hit_info  != NULL);
    RT_ASSERT(material  != NULL);
    RT_ASSERT(world     != NULL);

    rt_vec4_t final_color = material->ambient;

    const rt_idx_t directional_light_count              = world->directional_light_count;
    const rt_directional_light_t* directional_lights    = world->directional_light_buffer;

    for (rt_idx_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];
        RT_ASSERT(light != NULL);

        rt_idx_t shadow_light_index             = -1;
        enum rt_light_type shadow_light_type    = RT_LIGHT_null_light;

        RT_FIND_INDEX_OF_SHADOW_LIGHT(  shadow_lights,
                                        shadow_light_types,
                                        shadow_light_count,
                                        i,
                                        RT_LIGHT_directional_light,
                                        shadow_light_index,
                                        shadow_light_type);

        RT_COMPUTE_SHADOW_FOR_LIGHT(    light,
                                        material,
                                        i,
                                        RT_LIGHT_directional_light,
                                        shadow_light_index,
                                        shadow_light_type);

        rt_ray_t light_ray              = { .org = hit_info->position,
                                            .dir = rt_vec4_negate(light->direction), };

        rt_vec4_t light_color           = rt_atmosphere_scatter(&world->atmosphere,
                                                                light_ray,
                                                                rt_vec4_mul_scalar(light->color,
                                                                                   light->intensity),
                                                                light->direction);

        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_negate(light->direction));

        rt_vec4_t diffuse_color         = rt_vec4_mul(material->diffuse,
                                                      light_color);

        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), RT_FLOAT(0.0));

        rt_vec4_t diffuse               = rt_vec4_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      diffuse);
    }

    rt_idx_t point_light_count          = world->point_light_count;
    const rt_point_light_t* point_lights= world->point_light_buffer;

    for (rt_idx_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light   = &point_lights[i];
        RT_ASSERT(light != NULL);

        rt_idx_t shadow_light_index             = -1;
        enum rt_light_type shadow_light_type    = RT_LIGHT_null_light;

        RT_FIND_INDEX_OF_SHADOW_LIGHT(  shadow_lights,
                                        shadow_light_types,
                                        shadow_light_count,
                                        i,
                                        RT_LIGHT_point_light,
                                        shadow_light_index,
                                        shadow_light_type);

        RT_COMPUTE_SHADOW_FOR_LIGHT(    light,
                                        material,
                                        i,
                                        RT_LIGHT_point_light,
                                        shadow_light_index,
                                        shadow_light_type);

        rt_vec4_t light_color           = rt_vec4_mul_scalar(light->color,
                                                             light->intensity);

        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_sub(light->position,
                                                                   hit_info->position));

        rt_vec4_t diffuse_color         = rt_vec4_mul(material->diffuse,
                                                      light_color);

        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), RT_FLOAT(0.0));

        rt_vec4_t diffuse               = rt_vec4_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      diffuse);
    }

    return final_color;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_fragment_shader_metallic(const rt_hit_info_t*           hit_info,
                                             const rt_metallic_material_t*  material,
                                             const rt_world_t*              world,
                                             const rt_idx_t*                shadow_lights,
                                             const enum rt_light_type*      shadow_light_types,
                                             rt_idx_t                       shadow_light_count)

{
    RT_ASSERT(hit_info  != NULL);
    RT_ASSERT(material  != NULL);
    RT_ASSERT(world     != NULL);

    rt_vec4_t final_color = material->ambient;

    rt_idx_t directional_light_count              = world->directional_light_count;
    const rt_directional_light_t* directional_lights    = world->directional_light_buffer;

    for (rt_idx_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];
        RT_ASSERT(light != NULL);

        rt_idx_t shadow_light_index             = -1;
        enum rt_light_type shadow_light_type    = RT_LIGHT_null_light;

        RT_FIND_INDEX_OF_SHADOW_LIGHT(  shadow_lights,
                                        shadow_light_types,
                                        shadow_light_count,
                                        i,
                                        RT_LIGHT_directional_light,
                                        shadow_light_index,
                                        shadow_light_type);

        RT_COMPUTE_SHADOW_FOR_LIGHT(    light,
                                        material,
                                        i,
                                        RT_LIGHT_directional_light,
                                        shadow_light_index,
                                        shadow_light_type);

        rt_ray_t light_ray             = { .org = hit_info->position,
                                           .dir = rt_vec4_negate(light->direction), };

        rt_vec4_t light_color           = rt_atmosphere_scatter(&world->atmosphere,
                                                                light_ray,
                                                                rt_vec4_mul_scalar(light->color,
                                                                                   light->intensity),
                                                                light->direction);

        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_negate(light->direction));

        rt_vec4_t specular_color        = rt_vec4_mul(material->specular,
                                                      light_color);

        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), RT_FLOAT(0.0));

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

        rt_idx_t shadow_light_index             = -1;
        enum rt_light_type shadow_light_type    = RT_LIGHT_null_light;

        RT_FIND_INDEX_OF_SHADOW_LIGHT(  shadow_lights,
                                        shadow_light_types,
                                        shadow_light_count,
                                        i,
                                        RT_LIGHT_point_light,
                                        shadow_light_index,
                                        shadow_light_type);

        RT_COMPUTE_SHADOW_FOR_LIGHT(    light,
                                        material,
                                        i,
                                        RT_LIGHT_point_light,
                                        shadow_light_index,
                                        shadow_light_type);

        rt_vec4_t light_color           = rt_vec4_mul_scalar(light->color,
                                                             light->intensity);

        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_sub(light->position,
                                                                   hit_info->position));

        rt_vec4_t specular_color        = rt_vec4_mul(material->specular,
                                                      light_color);

        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), RT_FLOAT(0.0));

        rt_vec4_t specular              = rt_vec4_mul_scalar(specular_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      specular);
    }

    return final_color;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_fragment_shader_checkerboard(const rt_hit_info_t*               hit_info,
                                                 const rt_checkerboard_material_t*  material,
                                                 const rt_world_t*                  world,
                                                 const rt_idx_t*                    shadow_lights,
                                                 const enum rt_light_type*          shadow_light_types,
                                                 rt_idx_t                           shadow_light_count)
{
    RT_ASSERT(hit_info     != NULL);
    RT_ASSERT(material     != NULL);

    uint32_t pos_x_quant    = (uint32_t)(ceil(hit_info->position.x * 0.5f));
    uint32_t pos_z_quant    = (uint32_t)(ceil(hit_info->position.z * 0.5f));
    uint32_t pos_quant      = pos_x_quant + pos_z_quant;

    rt_vec4_t material_color = (pos_quant & 1) ? material->color_0
                                               : material->color_1;

    rt_vec4_t final_color   = rt_vec4_mul_scalar(material_color,
                                                 material->ambient_factor);

    const rt_idx_t directional_light_count              = world->directional_light_count;
    const rt_directional_light_t* directional_lights    = world->directional_light_buffer;

    for (rt_idx_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];
        RT_ASSERT(light != NULL);

        rt_idx_t shadow_light_index             = -1;
        enum rt_light_type shadow_light_type    = RT_LIGHT_null_light;

        RT_FIND_INDEX_OF_SHADOW_LIGHT(  shadow_lights,
                                        shadow_light_types,
                                        shadow_light_count,
                                        i,
                                        RT_LIGHT_directional_light,
                                        shadow_light_index,
                                        shadow_light_type);

        RT_COMPUTE_SHADOW_FOR_LIGHT(    light,
                                        material,
                                        i,
                                        RT_LIGHT_directional_light,
                                        shadow_light_index,
                                        shadow_light_type);

        rt_ray_t light_ray             = { .org = hit_info->position,
                                           .dir = rt_vec4_negate(light->direction), };

        rt_vec4_t light_color           = rt_atmosphere_scatter(&world->atmosphere,
                                                                light_ray,
                                                                rt_vec4_mul_scalar(light->color,
                                                                                   light->intensity),
                                                                light->direction);

        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_negate(light->direction));

        rt_vec4_t diffuse_color         = rt_vec4_mul(material_color,
                                                      light_color);

        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), RT_FLOAT(0.0));

        rt_vec4_t diffuse               = rt_vec4_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      diffuse);
    }

    rt_idx_t point_light_count          = world->point_light_count;
    const rt_point_light_t* point_lights= world->point_light_buffer;

    for (rt_idx_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light   = &point_lights[i];
        RT_ASSERT(light != NULL);

        rt_idx_t shadow_light_index             = -1;
        enum rt_light_type shadow_light_type    = RT_LIGHT_null_light;

        RT_FIND_INDEX_OF_SHADOW_LIGHT(  shadow_lights,
                                        shadow_light_types,
                                        shadow_light_count,
                                        i,
                                        RT_LIGHT_point_light,
                                        shadow_light_index,
                                        shadow_light_type);

        RT_COMPUTE_SHADOW_FOR_LIGHT(    light,
                                        material,
                                        i,
                                        RT_LIGHT_point_light,
                                        shadow_light_index,
                                        shadow_light_type);

        rt_vec4_t light_color           = rt_vec4_mul_scalar(light->color,
                                                             light->intensity);

        rt_vec4_t light_dir             = rt_vec4_norm(rt_vec4_sub(light->position,
                                                                   hit_info->position));

        rt_vec4_t diffuse_color         = rt_vec4_mul(material_color,
                                                      light_color);

        rt_float_t diff                 = fmax(rt_vec4_dot(light_dir,
                                                           hit_info->normal), RT_FLOAT(0.0));

        rt_vec4_t diffuse               = rt_vec4_mul_scalar(diffuse_color,
                                                             diff);

        final_color                     = rt_vec4_add(final_color,
                                                      diffuse);
    }

    return final_color;
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
            case RT_HIT_aabb: {
                
                rt_aabb_t* aabb = &world->aabb_buffer[hit_info->geometry_index];
                RT_ASSERT(aabb != NULL);

                material_type   = aabb->material_type;
                material_index  = aabb->material_index;

                break;
            }
            default: {
                RT_ASSERT(0 && "Unspecified geometry type!");
                break;
            }
        }

        return !(material_index == -1                                 ||
                 RT_MATERIAL_TYPE_null_material == material_type      ||
                 RT_MATERIAL_TYPE_emissive_material == material_type);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_is_in_shadow(const rt_world_t*           world,
                            rt_hit_ext_info_t*          hit_info,
                            rt_float_t                  nearest_z,
                            rt_float_t                  farthest_z,
                            rt_idx_t*                   shadow_light_index,
                            enum rt_light_type*         shadow_light_type,
                            rt_idx_t                    max_shadow_lights)
{
    RT_ASSERT(world                 != NULL);
    RT_ASSERT(hit_info              != NULL);
    RT_ASSERT(nearest_z             <= farthest_z);
    RT_ASSERT(shadow_light_index    != NULL);
    RT_ASSERT(shadow_light_type     != NULL);
    RT_ASSERT(max_shadow_lights     >= 0);

    rt_idx_t shadow_light_cnt       = 0;

    rt_idx_t directional_light_count                    = world->directional_light_count;
    const rt_directional_light_t* directional_lights    = world->directional_light_buffer;

    for (rt_idx_t i = 0; i < directional_light_count; ++i) {

        const rt_directional_light_t* light = &directional_lights[i];
        RT_ASSERT(light != NULL);

        if (!light->casts_shadows) {
            continue;
        }

        rt_ray_t new_ray = {

            .org = rt_vec4_add(hit_info->info.position,
                               rt_vec4_mul_scalar(hit_info->info.normal, RT_SHADOW_BIAS)),

            .dir = rt_vec4_norm(rt_vec4_negate(light->direction)),

        };

        if (rt_world_any_closest_hit(world,
                                     new_ray,
                                     nearest_z,
                                     farthest_z,
                                     hit_info)) {

            bool should_be_in_shadow = rt_should_be_in_shadow(world,
                                                              hit_info);

            if (should_be_in_shadow && shadow_light_cnt < max_shadow_lights) {

                shadow_light_index[shadow_light_cnt]   = i;
                shadow_light_type[shadow_light_cnt]    = RT_LIGHT_directional_light;

                ++shadow_light_cnt;
            }
        }
    }
    
    rt_idx_t point_light_count              = world->point_light_count;
    const rt_point_light_t* point_lights    = world->point_light_buffer;

    for (rt_idx_t i = 0; i < point_light_count; ++i) {

        const rt_point_light_t* light = &point_lights[i];
        RT_ASSERT(light != NULL);

        if (!light->casts_shadows) {
            continue;
        }

        rt_vec4_t org = rt_vec4_add(hit_info->info.position,
                                    rt_vec4_mul_scalar(hit_info->info.normal, RT_SHADOW_BIAS));

        rt_ray_t new_ray = {
            .org = org,
            .dir = rt_vec4_norm(rt_vec4_sub(light->position, org)),
        };

        rt_float_t new_nearest_z  = RT_SHADOW_BIAS;
        rt_float_t new_farthest_z = rt_vec4_dist(org, light->position) - RT_SHADOW_BIAS;

        if (rt_world_any_closest_hit(world,
                                     new_ray,
                                     new_nearest_z,
                                     new_farthest_z,
                                     hit_info)) {

            bool should_be_in_shadow = rt_should_be_in_shadow(world,
                                                              hit_info);

            if (should_be_in_shadow && shadow_light_cnt < max_shadow_lights) {

                shadow_light_index[shadow_light_cnt]   = i;
                shadow_light_type[shadow_light_cnt]    = RT_LIGHT_point_light;

                ++shadow_light_cnt;
            }
        }
    }

    for (rt_idx_t i = shadow_light_cnt; i < max_shadow_lights; ++i) {

        shadow_light_index[i]   = -1;
        shadow_light_type[i]    = RT_LIGHT_null_light;

    }
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_world_compute_sky_color(const rt_world_t*   world,
                                            rt_ray_t            ray)
{
    RT_ASSERT(world != NULL);

    rt_idx_t dir_light_count = world->directional_light_count;

    if (!dir_light_count) {

        return world->clear_color;

    }

    rt_vec4_t total_scattered_light = {};

    for (rt_idx_t i = 0; i < dir_light_count; ++i) {

        rt_directional_light_t* light = &world->directional_light_buffer[i];
        RT_ASSERT(light != NULL);

        rt_vec4_t direction = rt_vec4_norm(rt_vec4_negate(light->direction));

        rt_vec4_t scattered_light = rt_atmosphere_scatter(&world->atmosphere,
                                                          ray,
                                                          rt_vec4_mul_scalar(light->color,
                                                                             light->intensity),
                                                          light->direction);

        rt_float_t disc_factor = RT_FLOAT(1.0) - light->radius;

        if (fmax(rt_vec4_dot(ray.dir, direction), RT_FLOAT(0.0)) > disc_factor) {

            total_scattered_light = rt_vec4_mul_scalar(scattered_light,
                                                       light->intensity);
        }
        
        total_scattered_light = rt_vec4_add(total_scattered_light,
                                            scattered_light);
    }

    return total_scattered_light;
}

///////////////////////////////////////////////////////////////////////////
RT_API rt_vec4_t rt_world_compute_color(const rt_world_t*       world,
                                        rt_ray_t                ray,
                                        rt_float_t              nearest_z,
                                        rt_float_t              farthest_z,
                                        rt_idx_t                depth)
{
    RT_ASSERT(world     != NULL);
    RT_ASSERT(nearest_z <= farthest_z);
    RT_ASSERT(depth     >= 0);

    rt_hit_ext_info_t hit_info = {};

    if (!rt_world_any_closest_hit(world,
                                  ray,
                                  nearest_z,
                                  farthest_z,
                                  &hit_info)) {

        return rt_vec4_max_vec4(rt_world_compute_sky_color(world, ray),
                                world->clear_color);
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
        case RT_HIT_aabb:

            RT_ASSERT(world->aabb_buffer != NULL);

            material_type   = world->aabb_buffer[geometry_index].material_type;
            material_index  = world->aabb_buffer[geometry_index].material_index;

            break;
        default:
            RT_ASSERT(0 && "Unhandled geometry type!");
            break;
    }

    if (RT_MATERIAL_TYPE_null_material == material_type) {

        return rt_vec4_max_vec4(rt_world_compute_sky_color(world, ray),
                                world->clear_color);
    }

    rt_hit_ext_info_t hit_info_copy = hit_info;;

    rt_idx_t shadow_lights[RT_MAX_SHADOW_LIGHTS]                    = {};
    enum rt_light_type shadow_light_types[RT_MAX_SHADOW_LIGHTS]     = {};

    rt_is_in_shadow(world,
                    &hit_info_copy,
                    nearest_z,
                    farthest_z,
                    shadow_lights,
                    shadow_light_types,
                    RT_MAX_SHADOW_LIGHTS);

    switch (material_type) {

        case RT_MATERIAL_TYPE_emissive_material: {

            rt_emissive_material_t* m = &world->emissive_material_buffer[material_index];

            return rt_fragment_shader_emissive(m);
        }
        case RT_MATERIAL_TYPE_checkerboard_material: {

            rt_checkerboard_material_t* m = &world->checkerboard_material_buffer[material_index];

            return rt_fragment_shader_checkerboard(&hit_info.info,
                                                   m,
                                                   world,
                                                   shadow_lights,
                                                   shadow_light_types,
                                                   RT_MAX_SHADOW_LIGHTS);
        }
        case RT_MATERIAL_TYPE_diffuse_material: {

            rt_diffuse_material_t* m = &world->diffuse_material_buffer[material_index];

            return rt_fragment_shader_diffuse(&hit_info.info,
                                              m,
                                              world,
                                              shadow_lights,
                                              shadow_light_types,
                                              RT_MAX_SHADOW_LIGHTS);
        }
        case RT_MATERIAL_TYPE_metallic_material: {

            rt_vec4_t reflected_ray_dir = rt_vec4_reflect(ray.dir,
                                                          hit_info.info.normal);

            rt_vec4_t reflected_ray_pos = rt_vec4_add(hit_info.info.position,
                                                      rt_vec4_mul_scalar(hit_info.info.normal,
                                                                         RT_SHADOW_BIAS));

            rt_ray_t reflected_ray = { .org = reflected_ray_pos,
                                       .dir = rt_vec4_norm(reflected_ray_dir) };

            rt_vec4_t reflected_color = rt_world_compute_color(world,
                                                               reflected_ray,
                                                               nearest_z,
                                                               farthest_z,
                                                               depth - 1);

            rt_metallic_material_t* m = &world->metallic_material_buffer[material_index];

            rt_vec4_t fragment_shader_output =  rt_fragment_shader_metallic(&hit_info.info,
                                                                            m,
                                                                            world,
                                                                            shadow_lights,
                                                                            shadow_light_types,
                                                                            RT_MAX_SHADOW_LIGHTS);

            return rt_vec4_mul(reflected_color, fragment_shader_output);
        }
        /*
        case RT_MATERIAL_DIELECTRIC: {
            rt_dielectric_material_t* dielectric_material = &world->dielectric_materials[material_index];
            return rt_fragment_shader_dielectric(&hit_info, dielectric_material, world->point_lights, world->point_light_count, false);
        }
        */
        default:
            RT_ASSERT(false && "Unhandled material type!");
            return rt_vec4_max_vec4(rt_world_compute_sky_color(world, ray),
                                    world->clear_color);
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
RT_API rt_status_t rt_framebuffer_create(uint32_t               width,
                                         uint32_t               height,
                                         rt_framebuffer_t*      framebuffer)
{
    RT_ASSERT(width         > 0);
    RT_ASSERT(height        > 0);
    RT_ASSERT(framebuffer   != NULL);

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
RT_API rt_status_t rt_framebuffer_resize(uint32_t               new_width,
                                         uint32_t               new_height,
                                         rt_framebuffer_t*      framebuffer)
{
    RT_ASSERT(new_width     > 0);
    RT_ASSERT(new_height    > 0);
    RT_ASSERT(framebuffer   != NULL);

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
    framebuffer->rgb_buffer[index] = rt_vec4_to_uint32_alpha(color, RT_FLOAT(1.0));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_framebuffer_free(rt_framebuffer_t*       framebuffer)
{
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
    rt_float_t  smoothing;
    rt_float_t  movement_speed;
    rt_float_t  rotation_speed;

    rt_float_t  pitch_delta;
    rt_float_t  yaw_delta;

    rt_idx_t    depth;

    rt_float_t  near;
    rt_float_t  far;

    rt_float_t  exposure;
    rt_float_t  exposure_new;
    rt_float_t  exposure_smoothing;
    bool        auto_exposure;
}
rt_fps_camera_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_fps_camera_t rt_fps_camera_create()
{
    rt_fps_camera_t camera = {

        .y_axis                 = { RT_FLOAT(0.0),
                                    RT_FLOAT(1.0),
                                    RT_FLOAT(0.0),
                                    RT_FLOAT(0.0), },

        .z_axis                 = { RT_FLOAT(0.0),
                                    RT_FLOAT(0.0),
                                    RT_FLOAT(1.0),
                                    RT_FLOAT(0.0), },

        .position               = { RT_FLOAT( 0.0),
                                    RT_FLOAT( 0.0),
                                    RT_FLOAT(-5.0),
                                    RT_FLOAT( 1.0), },

        .rotation               = {},
        .rotation_new           = {},

        .velocity_x             = {},
        .velocity_y             = {},
        .velocity_z             = {},

        .velocity_new_x         = {},
        .velocity_new_y         = {},
        .velocity_new_z         = {},

        .fovy                   = RT_PI * RT_FLOAT(0.5),

        .smoothing              = RT_FLOAT(10.0),
        .movement_speed         = RT_FLOAT(20.0),
        .rotation_speed         = RT_FLOAT(2.0),

        .pitch_delta            = RT_FLOAT(0.0),
        .yaw_delta              = RT_FLOAT(0.0),

        .depth                  = 3,

        .near                   = RT_FLOAT(0.3),
        .far                    = RT_FLOAT(1000.0),

        .exposure               = RT_FLOAT(1.0),
        .exposure_new           = RT_FLOAT(1.0),
        .exposure_smoothing     = RT_FLOAT(0.5),
        .auto_exposure          = true,
    };

    return camera;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_forward(rt_fps_camera_t*     camera,
                                       rt_float_t           delta_time)
{
    RT_ASSERT(camera != NULL);

    camera->velocity_new_z = rt_vec4_negate(rt_vec4_mul_scalar(camera->z_axis,
                                                               camera->movement_speed * delta_time));
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_forward(rt_fps_camera_t* camera)
{
    camera->velocity_new_z = (rt_vec4_t){};
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_backward(rt_fps_camera_t*    camera,
                                        rt_float_t          delta_time)
{
    camera->velocity_new_z = rt_vec4_mul_scalar(camera->z_axis,
                                                camera->movement_speed * delta_time);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_backward(rt_fps_camera_t* camera)
{
    camera->velocity_new_z = (rt_vec4_t){};
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_left(rt_fps_camera_t*    camera,
                                    rt_float_t          delta_time)
{
    rt_vec4_t strafe = rt_vec4_mul_scalar(rt_vec4_norm(rt_vec4_cross(camera->y_axis,
                                                                     camera->z_axis)),
                                          camera->movement_speed * delta_time);

    camera->velocity_new_x = rt_vec4_negate(strafe);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_left(rt_fps_camera_t* camera)
{
    camera->velocity_new_x = (rt_vec4_t){};
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_right(rt_fps_camera_t*   camera,
                                     rt_float_t         delta_time)
{
    rt_vec4_t strafe = rt_vec4_mul_scalar(rt_vec4_norm(rt_vec4_cross(camera->y_axis,
                                                                     camera->z_axis)),
                                          camera->movement_speed * delta_time);

    camera->velocity_new_x = strafe;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_right(rt_fps_camera_t* camera)
{
    camera->velocity_new_x = (rt_vec4_t){};
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_left(rt_fps_camera_t*      camera,
                                      rt_float_t            delta_time)
{
    camera->yaw_delta = camera->rotation_speed * delta_time;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_left(rt_fps_camera_t* camera)
{
    camera->yaw_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_right(rt_fps_camera_t*     camera,
                                       rt_float_t           delta_time)
{
    camera->yaw_delta = -camera->rotation_speed * delta_time;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_right(rt_fps_camera_t* camera)
{
    camera->yaw_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_down(rt_fps_camera_t*  camera,
                                      rt_float_t        delta_time)
{
    camera->pitch_delta = -camera->rotation_speed * delta_time;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_down(rt_fps_camera_t* camera)
{
    camera->pitch_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_rotate_up(rt_fps_camera_t*    camera,
                                    rt_float_t          delta_time)
{
    camera->pitch_delta = camera->rotation_speed * delta_time;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_rotating_up(rt_fps_camera_t* camera)
{
    camera->pitch_delta = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_up(rt_fps_camera_t*  camera,
                                  rt_float_t        delta_time)
{
    camera->velocity_new_y.y = camera->movement_speed * delta_time;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_up(rt_fps_camera_t* camera)
{
    camera->velocity_new_y.y = RT_FLOAT(0.0);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_move_down(rt_fps_camera_t*    camera,
                                    rt_float_t          delta_time)
{
    camera->velocity_new_y.y = -camera->movement_speed * delta_time;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_stop_moving_down(rt_fps_camera_t* camera)
{
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

    rt_float_t camera_h         = tan(camera->fovy * RT_FLOAT(0.5));

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

    if (camera->auto_exposure) {
        rt_float_t camera_exposure  = camera->exposure;
        camera_exposure             = camera_exposure + delta_time *
                                                        camera->exposure_smoothing *
                                                        (camera->exposure_new - camera_exposure);
        camera->exposure            = camera_exposure;
    }

    rt_float_t hdr_color_avg_accum = RT_FLOAT(0.0);

    for (uint32_t row = 0; row < rows; ++row) {

        rt_vec4_t ver = rt_vec4_mul_scalar(pixel_delta_v, (rt_float_t)row);

        for (uint32_t col = 0; col < cols; ++col) {

            rt_vec4_t hor           = rt_vec4_mul_scalar(pixel_delta_u,
                                                         (rt_float_t)col);

            rt_vec4_t pixel_center  = rt_vec4_add(pixel00_loc,
                                                  rt_vec4_add(hor, ver));

            rt_vec4_t ray_dir       = rt_vec4_norm(rt_vec4_sub(pixel_center,
                                                               camera->position));

            rt_ray_t r = {
                .dir = ray_dir,
                .org = camera->position,
            };

            rt_vec4_t hdr_color = rt_world_compute_color(world,
                                                         r,
                                                         camera->near,
                                                         camera->far,
                                                         camera->depth);

            hdr_color_avg_accum += (hdr_color.x + hdr_color.y + hdr_color.z);;

            rt_vec4_t tonemapped_color = rt_vec4_apply_1(rt_vec4_mul_scalar(hdr_color,
                                                                            -camera->exposure),
                                                         exp);

            tonemapped_color = rt_vec4_sub((rt_vec4_t){ RT_FLOAT(1.0),
                                                        RT_FLOAT(1.0),
                                                        RT_FLOAT(1.0),
                                                        RT_FLOAT(1.0)}, tonemapped_color);

            rt_vec4_t gamma_correct_color = rt_vec4_apply_2(tonemapped_color,
                                                            rt_apply_gamma_custom,
                                                            RT_GAMMA_INVERSE);

            rt_framebuffer_write(row,
                                 col,
                                 gamma_correct_color,
                                 framebuffer);
        }
    }

    hdr_color_avg_accum     = hdr_color_avg_accum / (rt_float_t)(cols * rows);
    camera->exposure_new    = RT_FLOAT(1.0) / (hdr_color_avg_accum + RT_FLOAT(1.0));
}

///////////////////////////////////////////////////////////////////////////
//////////////////// CAMERA MOVEMENT WITH NOTCURSES ///////////////////////
///////////////////////////////////////////////////////////////////////////
#ifdef RT_USE_NOTCURSES

#include <notcurses/notcurses.h>

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    uint32_t quit_key;
    uint32_t forward_key;
    uint32_t backward_key;
    uint32_t left_key;
    uint32_t right_key;
    uint32_t up_key;
    uint32_t down_key;
    uint32_t rotate_left_key;
    uint32_t rotate_right_key;
    uint32_t rotate_up_key;
    uint32_t rotate_down_key;
}
rt_fps_camera_notcurses_keybindings_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_fps_camera_notcurses_keybindings_t rt_fps_camera_notcurses_default_keybindings(void)
{
    rt_fps_camera_notcurses_keybindings_t keybindings = {

        .quit_key                = NCKEY_ESC,
        .forward_key             = 'w',
        .backward_key            = 's',
        .left_key                = 'a',
        .right_key               = 'd',
        .up_key                  = 'q',
        .down_key                = 'e',
        .rotate_left_key         = NCKEY_LEFT,
        .rotate_right_key        = NCKEY_RIGHT,
        .rotate_up_key           = NCKEY_UP,
        .rotate_down_key         = NCKEY_DOWN,

    };

    return keybindings;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_update_with_notcurses(rt_fps_camera_t*                                camera,
                                                rt_float_t                                      delta_time,
                                                const rt_fps_camera_notcurses_keybindings_t*    keybindings,
                                                struct notcurses*                               nc,
                                                bool*                                           is_running)
{
    RT_ASSERT(camera      != NULL);
    RT_ASSERT(keybindings != NULL);
    RT_ASSERT(nc          != NULL);
    RT_ASSERT(is_running  != NULL);

    struct ncinput input   = {};
    uint32_t input_id      = 0;

    while ((input_id = notcurses_get_nblock(nc, &input))) {

        if (input.id == keybindings->quit_key && input.evtype == NCTYPE_RELEASE) {
            *is_running = false;
            return;
        }

        if (input.id == keybindings->forward_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_move_forward(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_moving_forward(camera);
                    break;
                default:
                    break;
            }
        }
        else if (input.id == keybindings->backward_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_move_backward(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_moving_backward(camera);
                    break;
                default:
                    break;
            }
        }

        if (input.id == keybindings->left_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_move_left(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_moving_left(camera);
                    break;
                default:
                    break;
            }
        }
        else if (input.id == keybindings->right_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_move_right(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_moving_right(camera);
                    break;
                default:
                    break;
            }
        }

        if (input.id == keybindings->rotate_left_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_rotate_left(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_rotating_left(camera);
                    break;
                default:
                    break;
            }
        }
        else if (input.id == keybindings->rotate_right_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_rotate_right(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_rotating_right(camera);
                    break;
                default:
                    break;
            }
        }

        if (input.id == keybindings->rotate_down_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_rotate_down(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_rotating_down(camera);
                    break;
                default:
                    break;
            }
        }
        else if (input.id == keybindings->rotate_up_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_rotate_up(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_rotating_up(camera);
                    break;
                default:
                    break;
            }
        }

        if (input.id == keybindings->up_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_move_up(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_moving_up(camera);
                    break;
                default:
                    break;
            }
        }
        else if (input.id == keybindings->down_key) {
            switch (input.evtype) {
                case NCTYPE_PRESS:
                    rt_fps_camera_move_down(camera, delta_time);
                    break;
                case NCTYPE_RELEASE:
                    rt_fps_camera_stop_moving_down(camera);
                    break;
                default:
                    break;
            }
        }
    }
}

#endif

///////////////////////////////////////////////////////////////////////////
///////////////// CAMERA MOVEMENT WITH JOYSTICK (SDL3) ////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef RT_USE_SDL3

#include <SDL3/SDL.h>

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    int32_t quit_button;
    int32_t horizontal_movement_axis;
    int32_t vertical_movement_axis;
    int32_t horizontal_rotation_axis;
    int32_t vertical_rotation_axis;
    int32_t up_button;
    int32_t down_button;
}
rt_fps_camera_sdl3_joystick_keybindings_t;

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    SDL_Gamepad*    gamepad;
    int32_t         deadzone;
}
rt_sdl3_gamepad_info_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_fps_camera_sdl3_joystick_keybindings_t rt_fps_camera_sdl3_default_joystick_keybindings(void)
{
    rt_fps_camera_sdl3_joystick_keybindings_t keybindings = {
        
        .quit_button                = SDL_GAMEPAD_BUTTON_START,

        .horizontal_movement_axis   = SDL_GAMEPAD_AXIS_LEFTX,
        .vertical_movement_axis     = SDL_GAMEPAD_AXIS_LEFTY,

        .horizontal_rotation_axis   = SDL_GAMEPAD_AXIS_RIGHTX,
        .vertical_rotation_axis     = SDL_GAMEPAD_AXIS_RIGHTY,

        .up_button                  = SDL_GAMEPAD_BUTTON_LEFT_SHOULDER,
        .down_button                = SDL_GAMEPAD_BUTTON_RIGHT_SHOULDER,
    };

    return keybindings;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_sdl3_retrieve_gamepad(rt_sdl3_gamepad_info_t* gamepad_info)
{
    RT_ASSERT(gamepad_info != NULL);

    if (gamepad_info->gamepad) {
        return;
    }

    int32_t gamepad_cnt = 0;

    SDL_JoystickID* ids = SDL_GetGamepads(&gamepad_cnt);

    if (!gamepad_cnt) {
        SDL_free(ids);
        return;
    }

    for (int32_t i = 0; i < gamepad_cnt; ++i) {

        SDL_Gamepad* selected = SDL_OpenGamepad(ids[i]);
        if (!selected) {
            continue;
        }

        if (!gamepad_info->gamepad) {
            gamepad_info->gamepad = selected;
        }

        if (i > 0) {
            SDL_CloseGamepad(selected);
        }
    }

    SDL_free(ids);
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_fps_camera_update_with_sdl3_joystick(rt_fps_camera_t*                                 camera,
                                                    const rt_fps_camera_sdl3_joystick_keybindings_t* keybindings,
                                                    rt_sdl3_gamepad_info_t*                          gamepad_info,
                                                    rt_float_t                                       delta_time,
                                                    bool*                                            is_running)
{
    RT_ASSERT(camera        != NULL);
    RT_ASSERT(keybindings   != NULL);
    RT_ASSERT(gamepad_info  != NULL);
    RT_ASSERT(is_running    != NULL);

    rt_sdl3_retrieve_gamepad(gamepad_info);

    SDL_Gamepad* gamepad    = gamepad_info->gamepad;
    int32_t deadzone        = gamepad_info->deadzone;

    SDL_Event event         = {};

    while (SDL_PollEvent(&event)) {
        switch (event.type) {

            case SDL_EVENT_QUIT: {

                *is_running = false;

                break;
            }
            case SDL_EVENT_GAMEPAD_ADDED: {

                if (!gamepad) {
                    rt_sdl3_retrieve_gamepad(gamepad_info)  ;
                }

                break; 
            }
            case SDL_EVENT_GAMEPAD_REMOVED: {

                if (gamepad) {
                    SDL_CloseGamepad(gamepad);
                    gamepad_info->gamepad = NULL;
                }

                break;
            }
            case SDL_EVENT_GAMEPAD_BUTTON_DOWN: {

                if (event.gbutton.button == keybindings->quit_button) {
                    *is_running = false;
                }
                else if (event.gbutton.button == keybindings->up_button) {

                    camera->velocity_new_y = rt_vec4_mul_scalar(camera->y_axis,
                                                                camera->movement_speed * delta_time);

                }
                else if (event.gbutton.button == keybindings->down_button) {

                    camera->velocity_new_y = rt_vec4_mul_scalar(camera->y_axis,
                                                                -camera->movement_speed * delta_time);

                }

                break;
            }
            case SDL_EVENT_GAMEPAD_BUTTON_UP: {

                if (event.gbutton.button == keybindings->up_button) {
                    camera->velocity_new_y = (rt_vec4_t){};
                }
                else if (event.gbutton.button == keybindings->down_button) {
                    camera->velocity_new_y = (rt_vec4_t){};
                }

                break;
            }
            case SDL_EVENT_GAMEPAD_AXIS_MOTION: {

                if (!gamepad) {
                    break;
                }

                int32_t leftX = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_LEFTX);

                if (leftX < -deadzone || leftX > deadzone) {
                    camera->yaw_delta = -(rt_float_t)leftX / RT_FLOAT(32768.0)
                                                           * camera->rotation_speed
                                                           * delta_time;
                }
                else {
                    camera->yaw_delta = RT_FLOAT(0.0);
                }

                int32_t leftY = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_LEFTY);
                if (leftY < -deadzone || leftY > deadzone) {
                    camera->pitch_delta = -((rt_float_t)leftY / RT_FLOAT(32768.0)
                                                              * camera->rotation_speed
                                                              * delta_time);
                }
                else {
                    camera->pitch_delta = RT_FLOAT(0.0);
                }

                int32_t rightY = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_RIGHTY);

                if (rightY < -deadzone || rightY > deadzone) {
                    camera->velocity_new_z = rt_vec4_mul_scalar(camera->z_axis,
                                                                ((rt_float_t)rightY / RT_FLOAT(32768.0)
                                                                 * camera->movement_speed
                                                                 * delta_time));
                }
                else {
                    camera->velocity_new_z = (rt_vec4_t){};
                }

                int32_t rightX = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_RIGHTX);

                if (rightX < -deadzone || rightX > deadzone) {

                    rt_vec4_t strafe = rt_vec4_norm(rt_vec4_cross(camera->y_axis, camera->z_axis));

                    camera->velocity_new_x = rt_vec4_mul_scalar(strafe,
                                                                ((rt_float_t)rightX / RT_FLOAT(32768.0)
                                                                 * camera->movement_speed
                                                                 * delta_time));

                }
                else {
                    camera->velocity_new_x = (rt_vec4_t){};
                }
                break;
            }
            default:
                break;
        }
    }
}

#endif

///////////////////////////////////////////////////////////////////////////
/////////////////////////// NOTCURSES SURFACE /////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef RT_USE_NOTCURSES

#include <notcurses/notcurses.h>

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    struct ncplane*             render_surface;
    struct ncvisual_options     visual_options;

    uint32_t                    rows;
    uint32_t                    cols;

    ncblitter_e                 blitter;
}
rt_notcurses_surface_t;

///////////////////////////////////////////////////////////////////////////
RT_API rt_notcurses_surface_t rt_notcurses_surface_create(struct ncplane* ncplane)
{
    RT_ASSERT(ncplane != NULL);

    rt_notcurses_surface_t surface = {};

    surface.render_surface = ncplane;

    surface.rows = 0;
    surface.cols = 0;

    ncplane_dim_yx(ncplane, &surface.rows, &surface.cols);

    surface.rows *= 2;
    surface.cols *= 2;

    surface.visual_options.n         = ncplane;
    surface.visual_options.leny      = surface.rows;
    surface.visual_options.lenx      = surface.cols;
    surface.visual_options.blitter   = NCBLIT_2x2;
    surface.visual_options.scaling   = NCSCALE_NONE;
    surface.visual_options.flags     = NCVISUAL_OPTION_NODEGRADE;

    return surface;
}

///////////////////////////////////////////////////////////////////////////
RT_API bool rt_notcurses_surface_resize(rt_notcurses_surface_t*     surface,
                                        rt_framebuffer_t*           framebuffer)
{
    RT_ASSERT(surface      != NULL);
    RT_ASSERT(framebuffer  != NULL);

    uint32_t new_rows = 0;
    uint32_t new_cols = 0;

    ncplane_dim_yx(surface->render_surface, &new_rows, &new_cols);

    new_rows *= 2;
    new_cols *= 2;

    if (new_rows != surface->rows || new_cols != surface->cols) {

        RT_ASSERT(RT_STATUS_success == rt_framebuffer_resize(new_cols,
                                                             new_rows,
                                                             framebuffer));

        surface->visual_options.leny = new_rows;
        surface->visual_options.lenx = new_cols;

        surface->rows = new_rows;
        surface->cols = new_cols;

        return true;
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////
RT_API void rt_notcurses_surface_blit(const rt_notcurses_surface_t* surface,
                                      const rt_framebuffer_t*       framebuffer)
{
    RT_ASSERT(surface      != NULL);
    RT_ASSERT(framebuffer  != NULL);

    RT_ASSERT(-1 != ncblit_rgba(framebuffer->rgb_buffer,
                                (int32_t)(surface->cols * sizeof(uint32_t)),
                                &surface->visual_options));
}

#endif

///////////////////////////////////////////////////////////////////////////
/////////////////////// SDL INPUT (FOR JOYSTICKS) /////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef RT_USE_SDL3

#endif

#endif // RT_RAYTERM_H
