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

#include "rayterm.h"

#include <SDL3/SDL.h>
#include <notcurses/notcurses.h>

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <locale.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

///////////////////////////////////////////////////////////////////////////
static rt_vec3_t compute_color(rt_ray_t ray,
                               const rt_mesh_t* meshes,
                               uint32_t model_count,
                               const rt_point_light_t* point_light,
                               uint32_t depth)
{
    if (!depth) {
        goto compute_sky_color;
    }

    const float t_min = 1.0f;
    const float t_max = 1000.0f;

    rt_hit_info closest_hit_info = {};
    uint32_t closest_hint_index = {};

    bool has_hit = false;
    for (uint32_t i = 0; i < model_count; ++i) {
        rt_hit_info hit_info = {};
        if (rt_sphere_hit(meshes[i].sphere, ray, t_min, t_max, &hit_info)) {
            if (!hit_info.is_front_facing) {
                continue;
            }
            else if (!has_hit || hit_info.t < closest_hit_info.t) {
                has_hit = true;
                closest_hit_info = hit_info;
                closest_hint_index = i;
            }
        }
    }

    if (has_hit) {
        bool is_in_shadow = false;
        rt_ray_t new_ray = {
            .dir = rt_vec3_norm(rt_vec3_sub(point_light->position, closest_hit_info.position)),
            .org = closest_hit_info.position,
        };
        for (uint32_t i = 0; i < model_count; ++i) {
            rt_hit_info hit_info = {};
            if (i != closest_hint_index && rt_sphere_hit(meshes[i].sphere, new_ray, t_min, t_max, &hit_info)) {
                if (!is_in_shadow && meshes[i].material.type != RT_MATERIAL_TYPE_EMISSIVE) {
                    is_in_shadow = true;
                    break;
                }
            }
        }

        const rt_sphere_t* sphere = &meshes[closest_hint_index].sphere;
        const rt_material_t* material = &meshes[closest_hint_index].material;

        rt_vec3_t fragment_output = meshes[closest_hint_index].fragment_shader(sphere,
                                                                               material,
                                                                               point_light,
                                                                               &closest_hit_info);
        if (RT_MATERIAL_TYPE_METALLIC == material->type) {
            rt_vec3_t reflected_ray_dir = rt_vec3_reflect(ray.dir, closest_hit_info.normal);
            reflected_ray_dir = rt_vec3_norm(reflected_ray_dir);

            ray.dir = reflected_ray_dir;
            ray.org = closest_hit_info.position;

            rt_vec3_t reflected_color = compute_color(ray, meshes, model_count, point_light, depth - 1);
            return reflected_color;
        }
        // TODO: fix the broken dielectric code
        /*
        else if (RT_MATERIAL_TYPE_DIELECTRIC == material->type) {
            float refract_ind = 1.5f;
            if (!closest_hit_info.is_front_facing) {
                refract_ind = 1.0f / refract_ind;
            }

            rt_vec3_t raydir = rt_vec3_norm(ray.dir);
            rt_vec3_t refracted_ray_dir = rt_vec3_refract(raydir, closest_hit_info.normal, refract_ind);
            refracted_ray_dir = rt_vec3_norm(refracted_ray_dir);

            ray.dir = refracted_ray_dir;
            ray.org = closest_hit_info.position;

            rt_vec3_t refracted_color = compute_color(ray, meshes, model_count, point_light, depth - 1);
            return rt_vec3_mul(refracted_color, fragment_output);
        }
        */

        if (is_in_shadow) {
            return rt_vec3_mul_scalar(fragment_output, 0.1f);
        }

        return fragment_output;
    }

compute_sky_color:
    rt_vec3_t ray_dir_norm = rt_vec3_norm(ray.dir);
    rt_vec3_t white_color = { 1.0f, 0.2f, 0.2f };
    rt_vec3_t blue_color = { 0.0f, 0.0f, 1.0f };
    float k = (ray_dir_norm.y + 1.0f) * 0.5f;
    return rt_vec3_lerp(white_color, blue_color, k);
}

///////////////////////////////////////////////////////////////////////////
static SDL_Gamepad* retrieveGamepad()
{
    int32_t cnt = 0;
    SDL_JoystickID* ids = SDL_GetGamepads(&cnt);
    if (!cnt) {
        return NULL;
    }

    SDL_Gamepad* gamepad = NULL;
    for (int32_t i = 0; i < cnt; ++i) {
        SDL_Gamepad* selected = SDL_OpenGamepad(ids[i]);
        if (!gamepad) {
            gamepad = selected;
        }

        if (i > 0) {
            SDL_CloseGamepad(selected);
        }
    }
    SDL_free(ids);
    return gamepad;
}

///////////////////////////////////////////////////////////////////////////
int main([[maybe_unused]] int argc, char** argv)
{
    if (!SDL_Init(SDL_INIT_GAMEPAD)) {
        fprintf(stderr, "SDL_Init() failed: %s\n", SDL_GetError());
        return EXIT_FAILURE;
    }

    if (!SDL_HasGamepad()) {
        fprintf(stderr, "No gamepads are connected!\n");
        SDL_Quit();
        return EXIT_FAILURE;
    }

    SDL_Gamepad* gamepad = retrieveGamepad();

    if (!setlocale(LC_ALL, "")) {
        fprintf(stderr, "%s: setlocale(LC_ALL, \"\") failed\n", argv[0]);
        SDL_Quit();
        return EXIT_FAILURE;
    }

    struct notcurses_options opts = {};
    struct notcurses* nc = notcurses_core_init(&opts, stdout);
    if (!nc) {
        fprintf(stderr, "%s: notcurses_init() failed\n", argv[0]);
        SDL_Quit();
        return EXIT_FAILURE;
    }

    uint32_t rows = 0;
    uint32_t cols = 0;
    struct ncplane* nstd = notcurses_stddim_yx(nc, &rows, &cols);
    if (!nstd) {
        fprintf(stderr, "%s: notcurses_stdplane() failed\n", argv[0]);
        SDL_Quit();
        return EXIT_FAILURE;
    }

    rows *= 2;
    cols *= 2;

    size_t pixelbuffer_size = rows * cols * sizeof(uint32_t);
    uint32_t* rgb = malloc(pixelbuffer_size);
    if (!rgb) {
        fprintf(stderr, "%s: malloc() failed\n", argv[0]);
        SDL_Quit();
        return EXIT_FAILURE;
    }

    // Camera
    float aspect = (float)cols / (float)(rows * 2);

    rt_vec3_t camera_center = {
        .x = 0.0f,
        .y = 0.0f,
        .z = 0.0f,
    };

    rt_vec3_t camera_dir = {
        .x = 0.0f,
        .y = 0.0f,
        .z = 1.0f,
    };

    rt_vec3_t camera_up = {
        .x = 0.0f,
        .y = 1.0f,
        .z = 0.0f,
    };

    float camera_fovy = RT_PI * 0.5f;

    rt_mesh_t meshes[6];

    rt_vec3_t colors[] = {
        { 1.0f, 0.5f, 0.5f },
        { 0.5f, 1.0f, 0.5f },
        { 0.5f, 0.5f, 1.0f },
        { 0.25f, 0.5f, 0.9f }
    };

    for (uint32_t i = 0; i < 2; ++i) {
        for (uint32_t j = 0; j < 2; ++j) {
            const float sphere_radius = 3.0f;
            meshes[i * 2 + j].sphere = rt_sphere_create(rt_vec3_create((float)i * 10.0f - 5.0f,
                                                                      -2.0f,
                                                                      -(float)j * 10.0f - 5.0f),
                                                        sphere_radius);

            rt_vec3_t diffuse_color = colors[i * 2 + j];
            rt_vec3_t ambient_color = rt_vec3_mul_scalar(diffuse_color, 0.1f);

            enum rt_material_type material_type_choice = RT_MATERIAL_TYPE_LAMBERTIAN;
            if (i == j) {
                material_type_choice = RT_MATERIAL_TYPE_METALLIC;
            }

            meshes[i * 2 + j].material = rt_material_create(material_type_choice, ambient_color, diffuse_color);
            meshes[i * 2 + j].fragment_shader = diffuse_fragment_shader;
        }
    }

    // big sphere
    meshes[4].sphere = rt_sphere_create(rt_vec3_create(0.0f, -1005.5f, 0.0f), 1000.0f);
    meshes[4].material.ambient = rt_vec3_create(1.0f, 1.0f, 1.0f);
    meshes[4].material.diffuse = rt_vec3_create(1.0f, 1.0f, 1.0f);
    meshes[4].material.type = RT_MATERIAL_TYPE_LAMBERTIAN;
    meshes[4].fragment_shader = checkerboard_fragment_shader;

    rt_point_light_t point_light = {
        .position = rt_vec3_create(0.0f, 10.0f, -5.0f),
        .diffuse = rt_vec3_create(0.25f, 0.5f, 1.0f),
        .intensity = 2.0f,
    };

    // light sphere
    meshes[5].sphere = rt_sphere_create(point_light.position, 0.5f);
    meshes[5].material.ambient = rt_vec3_create(1.0f, 1.0f, 1.0f);
    meshes[5].material.diffuse = point_light.diffuse;
    meshes[5].material.type = RT_MATERIAL_TYPE_EMISSIVE;
    meshes[5].fragment_shader = emissive_fragment_shader;

    struct ncvisual_options vopts = {};
    vopts.n = nstd;
    vopts.leny = rows;
    vopts.lenx = cols;
    vopts.blitter = NCBLIT_2x2;
    vopts.scaling = NCSCALE_NONE;
    vopts.flags = NCVISUAL_OPTION_NODEGRADE;

    bool is_running = true;
    rt_vec3_t camera_velocity_forward = {};
    rt_vec3_t camera_velocity_strafe = {};
    rt_vec3_t camera_rotation = {};

    float last_time = 0.0f;
    float total_time = 0.0f;

    while (is_running) {

        float current_time = (float)SDL_GetTicks();
        float delta_time = (last_time > 0.0f) ? (current_time - last_time) : 0.0f;
        last_time = current_time;

        total_time += delta_time * 0.001f;

        struct ncinput input = {};
        uint32_t input_id = 0;

        while ((input_id = notcurses_get_nblock(nc, &input))) {

            if (input.id == 'q' && input.evtype == NCTYPE_RELEASE) {
                is_running = false;
            }
        }

        SDL_Event event = {};
        int32_t gamepad_deadzone = 1000;
        float camera_sensitivity = 0.1f;
        float camera_speed = 0.25f;
        float point_light_radius = 10.0f;

        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_EVENT_QUIT:
                    is_running = false;
                    break;
                case SDL_EVENT_GAMEPAD_BUTTON_DOWN:
                    if (SDL_GetGamepadButton(gamepad, SDL_GAMEPAD_BUTTON_RIGHT_STICK)) {
                        is_running = false;
                    }
                    break;
                case SDL_EVENT_GAMEPAD_BUTTON_UP:
                    break;
                case SDL_EVENT_GAMEPAD_AXIS_MOTION:
                    int leftX = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_LEFTX);
                    if (leftX < -gamepad_deadzone || leftX > gamepad_deadzone) {
                        camera_rotation.y = -(float)leftX / 32768.0f * camera_sensitivity;
                    }
                    else {
                        camera_rotation.y = 0.0f;
                    }

                    int32_t rightY = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_RIGHTY);
                    if (rightY < -gamepad_deadzone || rightY > gamepad_deadzone) {
                        camera_velocity_forward = rt_vec3_mul_scalar(camera_dir, ((float)rightY / 32768.0f));
                    }
                    else {
                        camera_velocity_forward = rt_vec3_create(0.0f, 0.0f, 0.0f);
                    }

                    int32_t rightX = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_RIGHTX);
                    if (rightX < -gamepad_deadzone || rightX > gamepad_deadzone) {
                        rt_vec3_t strafe = rt_vec3_norm(rt_vec3_cross(camera_up, camera_dir));
                        camera_velocity_strafe = rt_vec3_mul_scalar(strafe, ((float)rightX / 32768.0f));
                    }
                    else {
                        camera_velocity_strafe = rt_vec3_create(0.0f, 0.0f, 0.0f);
                    }
                    break;
                default:
                    break;
            }
        }

        camera_dir = rt_vec3_rotate_y(camera_dir, camera_rotation.y);

        rt_vec3_t camera_look_at = rt_vec3_sub(camera_center, camera_dir);
        float focal_length = rt_vec3_dist(camera_center, camera_look_at);
        float camera_h = tanf(camera_fovy * 0.5f);
        float viewport_height = 2.0f * camera_h * focal_length;
        float viewport_width = viewport_height * aspect;

        rt_vec3_t camera_w = rt_vec3_norm(rt_vec3_sub(camera_center, camera_look_at));
        rt_vec3_t camera_u = rt_vec3_norm(rt_vec3_cross(camera_up, camera_w));
        rt_vec3_t camera_v = rt_vec3_norm(rt_vec3_cross(camera_w, camera_u));

        rt_vec3_t viewport_u = rt_vec3_mul_scalar(camera_u, viewport_width);
        rt_vec3_t viewport_v = rt_vec3_mul_scalar(rt_vec3_negate(camera_v), viewport_height);

        rt_vec3_t pixel_delta_u = rt_vec3_div_scalar(viewport_u, (float)cols);
        rt_vec3_t pixel_delta_v = rt_vec3_div_scalar(viewport_v, (float)rows);
        rt_vec3_t pixel_delta_diag = rt_vec3_mul_scalar(rt_vec3_add(pixel_delta_u, pixel_delta_v), 0.5f);

        rt_vec3_t viewport_u_half = rt_vec3_mul_scalar(viewport_u, 0.5f);
        rt_vec3_t viewport_v_half = rt_vec3_mul_scalar(viewport_v, 0.5f);

        rt_vec3_t camera_total_velocity = rt_vec3_add(camera_velocity_forward, camera_velocity_strafe);
        float camera_total_velocity_len = rt_vec3_len(camera_total_velocity);
        if (camera_total_velocity_len > 1.0f) {
            camera_total_velocity = rt_vec3_div_scalar(camera_total_velocity, camera_total_velocity_len);
        }
        camera_center = rt_vec3_add(camera_center, rt_vec3_mul_scalar(camera_total_velocity, camera_speed * delta_time * 0.05f));

        rt_vec3_t viewport_upper_left = rt_vec3_sub(camera_center, rt_vec3_mul_scalar(camera_w, focal_length));
        viewport_upper_left = rt_vec3_sub(viewport_upper_left, viewport_u_half);
        viewport_upper_left = rt_vec3_sub(viewport_upper_left, viewport_v_half);

        rt_vec3_t pixel00_loc = rt_vec3_add(viewport_upper_left, pixel_delta_diag);

        point_light.position.x = point_light_radius * cosf(total_time);
        point_light.position.z = point_light_radius * sinf(total_time) - 5.0f;
        meshes[5].sphere.center = point_light.position;

        for (uint32_t row = 0; row < rows; ++row) {

            rt_vec3_t ver = rt_vec3_mul_scalar(pixel_delta_v, (float)row);

            for (uint32_t col = 0; col < cols; ++col) {

                rt_vec3_t hor = rt_vec3_mul_scalar(pixel_delta_u, (float)col);
                rt_vec3_t pixel_center = rt_vec3_add(pixel00_loc, rt_vec3_add(hor, ver));
                rt_vec3_t ray_dir = rt_vec3_sub(pixel_center, camera_center);

                rt_ray_t r = {
                    .dir = ray_dir,
                    .org = camera_center,
                };

                rt_vec3_t pixel_color = compute_color(r, meshes, RT_BUFFER_LEN(meshes), &point_light, 3);
                pixel_color = rt_vec3_apply_2(pixel_color, rt_apply_gamma_custom, 1.0f / 2.2f);
                rgb[row * cols + col] = rt_vec3_to_uint32_alpha(pixel_color, 1.0f);
            }
        }

        if (-1 == ncblit_rgba(rgb, (int32_t)(cols * sizeof(uint32_t)), &vopts)) {
            fprintf(stderr, "%s: ncblit_rgba() failed\n", argv[0]);
            SDL_Quit();
            return EXIT_FAILURE;
        }

        if (-1 == notcurses_render(nc)) {
            fprintf(stderr, "%s: notcurses_render() failed\n", argv[0]);
            SDL_Quit();
            return EXIT_FAILURE;
        }

        SDL_Delay(16);
    }
    free(rgb);
    SDL_Quit();
    return notcurses_stop(nc);
}
