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

#include "rt_vec3.h"

///////////////////////////////////////////////////////////////////////////
typedef struct {
    rt_vec3_t org;
    rt_vec3_t dir;
} rt_ray_t;

///////////////////////////////////////////////////////////////////////////
static inline rt_ray_t rt_ray_create(rt_vec3_t org, rt_vec3_t dir)
{
    rt_ray_t ray = { org, dir };
    return ray;
}

///////////////////////////////////////////////////////////////////////////
static inline rt_vec3_t rt_ray_at(rt_ray_t ray, float t)
{
    rt_vec3_t dir_scaled = rt_vec3_mul_scalar(ray.dir, t);
    return rt_vec3_add(ray.org, dir_scaled);
}
