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
#include <SDL3/SDL_keyboard.h>
#include <notcurses/notcurses.h>

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <locale.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

///////////////////////////////////////////////////////////////////////////
static SDL_Gamepad* retrieveGamepad()
{
    int32_t cnt = 0;
    SDL_JoystickID* ids = SDL_GetGamepads(&cnt);
    if (!cnt) {
        SDL_free(ids);
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

    rt_vec3_t camera_center     = rt_vec3_create(0.0f, 0.0f, 0.0f);
    rt_vec3_t camera_dir        = rt_vec3_create(0.0f, 0.0f, 1.0f);
    rt_vec3_t camera_up         = rt_vec3_create(0.0f, 1.0f, 0.0f);

    float camera_fovy = RT_PI * 0.5f;

    rt_world_t world = {};
    world.clear_color = rt_vec3_create(0.02f, 0.03f, 0.03f);
    world.face_cull_mode = RT_FACE_CULL_BACK;

    rt_vec3_t colors[] = {
        { 1.0f, 0.5f, 0.5f },
        { 0.5f, 1.0f, 0.5f },
        { 0.5f, 0.5f, 1.0f },
        { 0.25f, 0.5f, 0.9f }
    };

    const float sphere_radius = 3.0f;
    for (uint32_t i = 0; i < 2; ++i) {
        for (uint32_t j = 0; j < 2; ++j) {
            uint32_t sphere_index = rt_world_push_sphere(&world);
            RT_ASSERT_PUSH(sphere_index);

            rt_sphere_t* sphere = &world.spheres[sphere_index];
            sphere->center.x = (float)i * 10.0f - 5.0f;
            sphere->center.y = -2.0f;
            sphere->center.z = -(float)j * 10.0f - 5.0f;
            sphere->radius = sphere_radius;

            uint32_t material_index = rt_world_push_diffuse_material(&world);
            RT_ASSERT_PUSH(material_index);
            rt_sphere_set_diffuse_material(&world, sphere_index, material_index);

            rt_diffuse_material_t* material = &world.diffuse_materials[material_index];
            material->ambient = rt_vec3_mul_scalar(colors[i * 2 + j], 0.1f);
            material->diffuse = colors[i * 2 + j];
            material->receives_shadows = true;
        }
    }

    uint32_t big_sphere_index = rt_world_push_sphere(&world);
    RT_ASSERT_PUSH(big_sphere_index);

    rt_sphere_t* big_sphere = &world.spheres[big_sphere_index];
    big_sphere->center = rt_vec3_create(0.0f, -500.0f, 0.0f);
    big_sphere->radius = 495.0f;

    uint32_t big_sphere_material_index = rt_world_push_checkerboard_material(&world);
    RT_ASSERT_PUSH(big_sphere_material_index);
    rt_sphere_set_checkerboard_material(&world, big_sphere_index, big_sphere_material_index);

    rt_checkerboard_material_t* big_sphere_material = &world.checkerboard_materials[big_sphere_material_index];
    big_sphere_material->color_0 = rt_vec3_create(0.1f, 0.1f, 0.1f);
    big_sphere_material->color_1 = rt_vec3_create(0.9f, 0.9f, 0.9f);
    big_sphere_material->receives_shadows = true;
    big_sphere_material->shadow_factor = 0.2f;

    uint32_t point_light_index = rt_world_push_point_light(&world);
    RT_ASSERT_PUSH(point_light_index);

    rt_point_light_t* point_light = &world.point_lights[point_light_index];
    point_light->position = rt_vec3_create(0.0f, 10.0f, -5.0f);
    point_light->color = rt_vec3_create(0.25f, 0.5f, 1.0f);
    point_light->intensity = 2.0f;
    point_light->casts_shadows = true;

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

        point_light->position.x = point_light_radius * cosf(total_time);
        point_light->position.z = point_light_radius * sinf(total_time) - 5.0f;

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

                rt_vec3_t pixel_color = rt_world_compute_color(&world,
                                                               r,
                                                               (rt_vec2_t){0.1f, 1000.0f},
                                                               10);
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

    if (gamepad) {
        SDL_CloseGamepad(gamepad);
    }

    free(rgb);
    rt_world_free(&world);
    SDL_Quit();
    return notcurses_stop(nc);
}
