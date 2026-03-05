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

    rt_framebuffer_t framebuffer = {};
    RT_ASSERT(RT_STATUS_success == rt_framebuffer_create(cols, rows, &framebuffer));

    rt_world_t world            = {};

    rt_vec4_set(&world.clear_color,
                RT_FLOAT(0.02),
                RT_FLOAT(0.03),
                RT_FLOAT(0.03),
                RT_FLOAT(1.00));

    world.face_cull_mode = RT_FACE_cull_back;

    rt_idx_t plane_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_plane(&world,
                                                       &plane_index));

    rt_vec4_t plane_position = { RT_FLOAT(0.0),
                                 RT_FLOAT(-5.0),
                                 RT_FLOAT(0.0),
                                 RT_FLOAT(1.0) };

    rt_vec4_t plane_normal = { RT_FLOAT(0.0),
                               RT_FLOAT(1.0),
                               RT_FLOAT(0.0),
                               RT_FLOAT(0.0) };

    rt_world_set_plane_params(&world,
                              plane_index,
                              plane_position,
                              plane_normal);

    rt_idx_t plane_material_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_checkerboard_material(&world,
                                                                       &plane_material_index));

    rt_vec4_t mat_color_0 = { RT_FLOAT(0.1),
                              RT_FLOAT(0.1),
                              RT_FLOAT(0.1),
                              RT_FLOAT(1.0) };

    rt_vec4_t mat_color_1 = { RT_FLOAT(0.9),
                              RT_FLOAT(0.9),
                              RT_FLOAT(0.9),
                              RT_FLOAT(1.0) };

    rt_world_set_checkerboard_material_params(&world,
                                              plane_material_index,
                                              mat_color_0,
                                              mat_color_1,
                                              RT_FLOAT(0.2),
                                              true);

    rt_plane_link_checkerboard_material(&world,
                                        plane_index,
                                        plane_material_index);

    rt_idx_t point_light_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_point_light(&world,
                                                             &point_light_index));

    rt_vec4_t point_light_pos = { RT_FLOAT(0.0),
                                  RT_FLOAT(10.0),
                                  RT_FLOAT(-5.0),
                                  RT_FLOAT(1.0) };

    rt_vec4_t point_light_color = { RT_FLOAT(0.25),
                                    RT_FLOAT(0.5),
                                    RT_FLOAT(1.0),
                                    RT_FLOAT(1.0) };

    rt_world_set_point_light_params(&world,
                                    point_light_index,
                                    point_light_color,
                                    point_light_pos,
                                    RT_FLOAT(2.0),
                                    true);

    struct ncvisual_options vopts = {};
    vopts.n         = nstd;
    vopts.leny      = rows;
    vopts.lenx      = cols;
    vopts.blitter   = NCBLIT_2x2;
    vopts.scaling   = NCSCALE_NONE;
    vopts.flags     = NCVISUAL_OPTION_NODEGRADE;

    bool is_running = true;

    rt_vec4_t camera_z          = { RT_FLOAT(0.0),
                                    RT_FLOAT(0.0),
                                    RT_FLOAT(-1.0),
                                    RT_FLOAT(0.0) };

    rt_vec4_t camera_center     = { RT_FLOAT(0.0),
                                    RT_FLOAT(0.0),
                                    RT_FLOAT(0.0),
                                    RT_FLOAT(1.0) };

    rt_vec4_t camera_up         = { RT_FLOAT(0.0),
                                    RT_FLOAT(1.0),
                                    RT_FLOAT(0.0),
                                    RT_FLOAT(0.0) };

    rt_float_t camera_fovy      = RT_PI * RT_FLOAT(0.5);


    rt_vec4_t camera_velocity_forward   = {};
    rt_vec4_t camera_velocity_strafe    = {};
    rt_vec4_t camera_velocity_up        = {};
    rt_vec4_t camera_rotation           = {};

    rt_vec4_t camera_velocity_forward_new   = {};
    rt_vec4_t camera_velocity_strafe_new    = {};
    rt_vec4_t camera_velocity_up_new        = {};
    rt_vec4_t camera_rotation_new           = {};

    rt_float_t camera_sensitivity       = RT_FLOAT(0.1);
    rt_float_t camera_smoothing         = RT_FLOAT(0.01);
    rt_float_t camera_speed             = RT_FLOAT(0.5);
    rt_float_t camera_rotation_speed    = RT_FLOAT(0.05);
    rt_float_t camera_pitch_delta       = RT_FLOAT(0.0);
    rt_float_t camera_yaw_delta         = RT_FLOAT(0.0);

    rt_float_t point_light_radius       = RT_FLOAT(1.0);

    rt_float_t last_time     = RT_FLOAT(0.0);
    rt_float_t total_time    = RT_FLOAT(0.0);

    while (is_running) {

        uint32_t new_rows = 0;
        uint32_t new_cols = 0;

        ncplane_dim_yx(nstd, &new_rows, &new_cols);

        new_rows *= 2;
        new_cols *= 2;

        if (new_rows != vopts.leny || new_cols != vopts.lenx) {

            RT_ASSERT(RT_STATUS_success == rt_framebuffer_resize(new_cols,
                                                                 new_rows,
                                                                 &framebuffer));

            vopts.leny = new_rows;
            vopts.lenx = new_cols;

            rows = new_rows;
            cols = new_cols;

            continue;
        }

        // Camera
        rt_float_t aspect = (rt_float_t)cols / (rt_float_t)(rows * 2);

        rt_float_t current_time = (rt_float_t)SDL_GetTicks();
        rt_float_t delta_time = (last_time > RT_FLOAT(0.0)) ? (current_time - last_time)
                                                            : RT_FLOAT(0.0);
        last_time = current_time;

        total_time += delta_time * RT_FLOAT(0.001);

        struct ncinput input = {};
        uint32_t input_id = 0;

        rt_vec4_t camera_dir        = { RT_FLOAT(0.0),
                                        RT_FLOAT(0.0),
                                        RT_FLOAT(-1.0),
                                        RT_FLOAT(0.0) };


        while ((input_id = notcurses_get_nblock(nc, &input))) {

            if (input.id == NCKEY_ESC && input.evtype == NCTYPE_RELEASE) {
                is_running = false;
            }

            if (input.id == 'w') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_velocity_forward_new = rt_vec4_negate(rt_vec4_mul_scalar(camera_z,
                                                                                        camera_speed));
                        break;
                    case NCTYPE_RELEASE:
                        rt_vec4_zero(&camera_velocity_forward_new);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == 's') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_velocity_forward_new = rt_vec4_mul_scalar(camera_z,
                                                                         camera_speed);
                        break;
                    case NCTYPE_RELEASE:
                        rt_vec4_zero(&camera_velocity_forward_new);
                        break;
                    default:
                        break;
                }
            }

            rt_vec4_t strafe = rt_vec4_mul_scalar(rt_vec4_norm(rt_vec4_cross(camera_up,
                                                                             camera_z)), camera_speed);

            if (input.id == 'a') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_velocity_strafe_new = rt_vec4_negate(strafe);
                        break;
                    case NCTYPE_RELEASE:
                        rt_vec4_zero(&camera_velocity_strafe_new);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == 'd') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_velocity_strafe_new = strafe;
                        break;
                    case NCTYPE_RELEASE:
                        rt_vec4_zero(&camera_velocity_strafe_new);
                        break;
                    default:
                        break;
                }
            }

            if (input.id == NCKEY_LEFT) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_yaw_delta = camera_rotation_speed;
                        break;
                    case NCTYPE_RELEASE:
                        camera_yaw_delta = RT_FLOAT(0.0);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == NCKEY_RIGHT) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_yaw_delta = -camera_rotation_speed;
                        break;
                    case NCTYPE_RELEASE:
                        camera_yaw_delta = RT_FLOAT(0.0);
                        break;
                    default:
                        break;
                }
            }

            if (input.id == NCKEY_UP) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_pitch_delta = -camera_rotation_speed;
                        break;
                    case NCTYPE_RELEASE:
                        camera_pitch_delta = RT_FLOAT(0.0);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == NCKEY_DOWN) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_pitch_delta = camera_rotation_speed;
                        break;
                    case NCTYPE_RELEASE:
                        camera_pitch_delta = RT_FLOAT(0.0);
                        break;
                    default:
                        break;
                }
            }

            if (input.id == 'q') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_velocity_up_new.y = camera_speed;
                        break;
                    case NCTYPE_RELEASE:
                        camera_velocity_up_new.y = RT_FLOAT(0.0);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == 'e') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        camera_velocity_up_new.y = -camera_speed;
                        break;
                    case NCTYPE_RELEASE:
                        camera_velocity_up_new.y = RT_FLOAT(0.0);
                        break;
                    default:
                        break;
                }
            }
        }

        SDL_Event event                 = {};
        int32_t gamepad_deadzone        = 1000;

        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_EVENT_QUIT:
                    is_running = false;
                    break;
                case SDL_EVENT_GAMEPAD_AXIS_MOTION: {
                    int leftX = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_LEFTX);
                    if (leftX < -gamepad_deadzone || leftX > gamepad_deadzone) {
                        camera_rotation_new.y = -(rt_float_t)leftX / RT_FLOAT(32768.0) * camera_sensitivity;
                    }
                    else {
                        camera_rotation_new.y = RT_FLOAT(0.0);
                    }

                    int32_t rightY = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_RIGHTY);
                    if (rightY < -gamepad_deadzone || rightY > gamepad_deadzone) {
                        camera_velocity_forward_new = rt_vec4_mul_scalar(camera_dir, ((rt_float_t)rightY / RT_FLOAT(32768.0)));
                    }
                    else {
                        rt_vec4_zero(&camera_velocity_forward_new);
                    }

                    int32_t rightX = SDL_GetGamepadAxis(gamepad, SDL_GAMEPAD_AXIS_RIGHTX);
                    if (rightX < -gamepad_deadzone || rightX > gamepad_deadzone) {
                        rt_vec4_t strafe = rt_vec4_norm(rt_vec4_cross(camera_up, camera_dir));
                        camera_velocity_strafe_new = rt_vec4_mul_scalar(strafe, ((rt_float_t)rightX / RT_FLOAT(32768.0)));
                    }
                    else {
                        rt_vec4_zero(&camera_velocity_strafe_new);
                    }
                    break;
                }
                default:
                    break;
            }
        }

        camera_rotation_new.x       += camera_pitch_delta;
        camera_rotation_new.x       = rt_clamp(camera_rotation_new.x,
                                               RT_TO_RADIANS(RT_FLOAT(-89.0)),
                                               RT_TO_RADIANS(RT_FLOAT( 89.0)));

        camera_rotation_new.y       += camera_yaw_delta;

        camera_rotation             = rt_vec4_lerp(camera_rotation,
                                                   camera_rotation_new,
                                                   camera_smoothing * delta_time);

        camera_dir                  = rt_vec4_rotate_x(camera_dir, camera_rotation.x);
        camera_dir                  = rt_vec4_rotate_y(camera_dir, camera_rotation.y);

        camera_z                    = camera_dir;

        camera_velocity_forward     = rt_vec4_lerp(camera_velocity_forward,
                                                   camera_velocity_forward_new,
                                                   camera_smoothing * delta_time);

        camera_velocity_strafe      = rt_vec4_lerp(camera_velocity_strafe,
                                                   camera_velocity_strafe_new,
                                                   camera_smoothing * delta_time);

        camera_velocity_up          = rt_vec4_lerp(camera_velocity_up,
                                                   camera_velocity_up_new,
                                                   camera_smoothing * delta_time);

        rt_vec4_t camera_total_velocity = rt_vec4_add(camera_velocity_forward,
                                                      camera_velocity_strafe);

        camera_total_velocity       = rt_vec4_add(camera_total_velocity,
                                                  camera_velocity_up);

        camera_center               = rt_vec4_add(camera_center,
                                                  camera_total_velocity);

        rt_vec4_t camera_look_at    = rt_vec4_sub(camera_center,
                                                  camera_dir);

        rt_float_t focal_length     = rt_vec4_dist(camera_center,
                                                   camera_look_at);

        rt_float_t camera_h         = tanf(camera_fovy * RT_FLOAT(0.5));

        rt_float_t viewport_height  = RT_FLOAT(2.0) * camera_h * focal_length;

        rt_float_t viewport_width   = viewport_height * aspect;

        rt_vec4_t camera_w          = rt_vec4_norm(rt_vec4_sub(camera_center,
                                                               camera_look_at));

        rt_vec4_t camera_u          = rt_vec4_norm(rt_vec4_cross(camera_up,
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

        rt_vec4_t viewport_u_half   = rt_vec4_mul_scalar(viewport_u, RT_FLOAT(0.5));

        rt_vec4_t viewport_v_half   = rt_vec4_mul_scalar(viewport_v, RT_FLOAT(0.5));

        rt_vec4_t viewport_upper_left = rt_vec4_sub(camera_center,
                                                    rt_vec4_mul_scalar(camera_w,
                                                                       focal_length));

        viewport_upper_left         = rt_vec4_sub(viewport_upper_left, viewport_u_half);
        viewport_upper_left         = rt_vec4_sub(viewport_upper_left, viewport_v_half);

        rt_vec4_t pixel00_loc       = rt_vec4_add(viewport_upper_left, pixel_delta_diag);

        rt_point_light_t* point_light = &world.point_light_buffer[point_light_index];

        point_light->position.x     = point_light_radius * cosf(total_time);
        point_light->position.z     = point_light_radius * sinf(total_time);

        for (uint32_t row = 0; row < rows; ++row) {

            rt_vec4_t ver = rt_vec4_mul_scalar(pixel_delta_v, (rt_float_t)row);

            for (uint32_t col = 0; col < cols; ++col) {

                rt_vec4_t hor           = rt_vec4_mul_scalar(pixel_delta_u,
                                                             (rt_float_t)col);

                rt_vec4_t pixel_center  = rt_vec4_add(pixel00_loc,
                                                      rt_vec4_add(hor, ver));

                rt_vec4_t ray_dir       = rt_vec4_sub(pixel_center,
                                                      camera_center);

                rt_ray_t r = {
                    .dir = ray_dir,
                    .org = camera_center,
                };

                rt_vec4_t pixel_color = rt_world_compute_color(&world,
                                                               r,
                                                               RT_FLOAT(0.1),
                                                               RT_FLOAT(1000.0),
                                                               10);
                pixel_color = rt_vec4_apply_2(pixel_color,
                                              rt_apply_gamma_custom, RT_GAMMA_INVERSE);

                rt_framebuffer_write(row, col, pixel_color, &framebuffer);
            }
        }

        if (-1 == ncblit_rgba(framebuffer.rgb_buffer,
                              (int32_t)(cols * sizeof(uint32_t)),
                              &vopts)) {

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

    rt_framebuffer_free(&framebuffer);
    rt_world_free(&world);
    SDL_Quit();
    return notcurses_stop(nc);
}
