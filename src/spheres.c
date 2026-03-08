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

/*

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

*/

///////////////////////////////////////////////////////////////////////////
int main([[maybe_unused]] int argc, char** argv)
{
    if (!SDL_Init(SDL_INIT_GAMEPAD)) {
        fprintf(stderr, "SDL_Init() failed: %s\n", SDL_GetError());
        return EXIT_FAILURE;
    }

    /*
    SDL_Gamepad* gamepad = retrieveGamepad();
    */

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
                RT_FLOAT(0.02),
                RT_FLOAT(0.03),
                RT_FLOAT(1.00));

    world.face_cull_mode = RT_FACE_cull_back;

    rt_idx_t plane_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_plane(&world,
                                                       &plane_index));

    rt_plane_params_t plane_params = {

        .position = { RT_FLOAT( 0.0),
                      RT_FLOAT(-5.0),
                      RT_FLOAT( 0.0),
                      RT_FLOAT( 1.0), },

        .normal = { RT_FLOAT(0.0),
                    RT_FLOAT(1.0),
                    RT_FLOAT(0.0),
                    RT_FLOAT(0.0), },

        .side_length = RT_FLOAT(50.0),
    };

    rt_world_set_plane_params(&world, plane_index, &plane_params);

    rt_idx_t plane_material_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_checkerboard_material(&world,
                                                                       &plane_material_index));

    rt_checkerboard_material_t mat_params = {

        .color_0 = { RT_FLOAT(0.5),
                     RT_FLOAT(0.5),
                     RT_FLOAT(0.5),
                     RT_FLOAT(1.0), },

        .color_1 = { RT_FLOAT(0.05),
                     RT_FLOAT(0.05),
                     RT_FLOAT(0.05),
                     RT_FLOAT(1.0), },

        .shadow_factor      = RT_FLOAT(0.2),

        .receives_shadows   = true,
    };

    rt_world_set_checkerboard_material_params(&world,
                                              plane_material_index,
                                              &mat_params);

    rt_plane_link_checkerboard_material(&world,
                                        plane_index,
                                        plane_material_index);

    rt_idx_t sphere_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_sphere(&world, &sphere_index));

    rt_sphere_params_t sphere_params = {
        .center = { RT_FLOAT( 0.0),
                    RT_FLOAT( 3.0),
                    RT_FLOAT(-10.0),
                    RT_FLOAT( 1.0), },

        .radius = RT_FLOAT(3.0),
    };

    rt_world_set_sphere_params(&world,
                               sphere_index,
                               &sphere_params);

    rt_idx_t sphere_material_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_metallic_material(&world, &sphere_material_index));

    rt_metallic_material_t sphere_mat_params = {
        
        .ambient = { RT_FLOAT(0.02),
                     RT_FLOAT(0.04),
                     RT_FLOAT(0.08),
                     RT_FLOAT(1.0), },

        .specular = { RT_FLOAT(0.2),   
                      RT_FLOAT(0.4),
                      RT_FLOAT(0.8),
                      RT_FLOAT(1.0), },

        .receives_shadows = true,
    };


    rt_world_set_metallic_material_params(&world,
                                          sphere_material_index,
                                          &sphere_mat_params);

    rt_sphere_link_metallic_material(&world,
                                     sphere_index,
                                     sphere_material_index);

    rt_idx_t point_light_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_point_light(&world,
                                                             &point_light_index));

    rt_point_light_t point_light_params = {

        .color = { RT_FLOAT(0.25),
                   RT_FLOAT(0.5),
                   RT_FLOAT(1.0),
                   RT_FLOAT(1.0), },

        .position = { RT_FLOAT(0.0),
                      RT_FLOAT(10.0),
                      RT_FLOAT(-5.0),
                      RT_FLOAT(1.0), },

        .casts_shadows = true,
        .intensity = RT_FLOAT(3.0),
    };

    rt_world_set_point_light_params(&world,
                                    point_light_index,
                                    &point_light_params);

    rt_idx_t point_light_sphere_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_sphere(&world,
                                                        &point_light_sphere_index));

    rt_sphere_params_t point_light_sphere_params = {
        .center = point_light_params.position,
        .radius = RT_FLOAT(0.5),
    };

    rt_world_set_sphere_params(&world,
                               point_light_sphere_index,
                               &point_light_sphere_params);

    rt_emissive_material_t emissive_params = {
        .color = { RT_FLOAT(1.0),
                   RT_FLOAT(1.0),
                   RT_FLOAT(1.0),
                   RT_FLOAT(1.0), },
    };

    rt_idx_t point_light_emissive_material_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_emissive_material(&world,
                                                                   &point_light_emissive_material_index));

    rt_world_set_emissive_material_params(&world,
                                          point_light_emissive_material_index,
                                          &emissive_params);

    rt_sphere_link_emissive_material(&world,
                                     point_light_sphere_index,
                                     point_light_emissive_material_index);

    struct ncvisual_options vopts = {};
    vopts.n         = nstd;
    vopts.leny      = rows;
    vopts.lenx      = cols;
    vopts.blitter   = NCBLIT_2x2;
    vopts.scaling   = NCSCALE_NONE;
    vopts.flags     = NCVISUAL_OPTION_NODEGRADE;

    bool is_running = true;

    rt_fps_camera_t fps_camera = rt_fps_camera_create();

    rt_float_t point_light_radius       = RT_FLOAT(10.0);

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

        rt_float_t current_time = (rt_float_t)SDL_GetTicks();
        rt_float_t delta_time = (last_time > RT_FLOAT(0.0)) ? (current_time - last_time)
                                                            : RT_FLOAT(0.0);
        last_time = current_time;
        total_time += delta_time * RT_FLOAT(0.001);

        point_light_params.position.x = sin(total_time) * point_light_radius;
        point_light_params.position.z = cos(total_time) * point_light_radius;

        rt_world_set_point_light_params(&world,
                                        point_light_index,
                                        &point_light_params);

        point_light_sphere_params.center = point_light_params.position;

        rt_world_set_sphere_params(&world,
                                   point_light_sphere_index,
                                   &point_light_sphere_params);

        struct ncinput input = {};
        uint32_t input_id = 0;

        while ((input_id = notcurses_get_nblock(nc, &input))) {

            if (input.id == NCKEY_ESC && input.evtype == NCTYPE_RELEASE) {
                is_running = false;
            }

            if (input.id == 'w') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_move_forward(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_moving_forward(&fps_camera);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == 's') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_move_backward(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_moving_backward(&fps_camera);
                        break;
                    default:
                        break;
                }
            }

            if (input.id == 'a') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_move_left(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_moving_left(&fps_camera);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == 'd') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_move_right(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_moving_right(&fps_camera);
                        break;
                    default:
                        break;
                }
            }

            if (input.id == NCKEY_LEFT) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_rotate_left(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_rotating_left(&fps_camera);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == NCKEY_RIGHT) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_rotate_right(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_rotating_right(&fps_camera);
                        break;
                    default:
                        break;
                }
            }

            if (input.id == NCKEY_UP) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_rotate_down(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_rotating_down(&fps_camera);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == NCKEY_DOWN) {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_rotate_up(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_rotating_up(&fps_camera);
                        break;
                    default:
                        break;
                }
            }

            if (input.id == 'q') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_move_up(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_moving_up(&fps_camera);
                        break;
                    default:
                        break;
                }
            }
            else if (input.id == 'e') {
                switch (input.evtype) {
                    case NCTYPE_PRESS:
                        rt_fps_camera_move_down(&fps_camera);
                        break;
                    case NCTYPE_RELEASE:
                        rt_fps_camera_stop_moving_down(&fps_camera);
                        break;
                    default:
                        break;
                }
            }
        }

        /*
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
        */

        rt_fps_camera_render(&fps_camera,
                              &world,
                              &framebuffer,
                              delta_time);

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

    /*
    if (gamepad) {
        SDL_CloseGamepad(gamepad);
    }
    */

    rt_framebuffer_free(&framebuffer);
    rt_world_free(&world);
    SDL_Quit();
    return notcurses_stop(nc);
}
