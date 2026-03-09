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

#define RT_USE_NOTCURSES
#define RT_USE_SDL3
#include "rayterm.h"

#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <unistd.h>
#include <assert.h>

///////////////////////////////////////////////////////////////////////////
typedef struct
{
    struct notcurses*       nc;
    rt_framebuffer_t        framebuffer;
    rt_world_t              world;
    rt_sdl3_gamepad_info_t  gamepad_info;
}
app_state_t;

///////////////////////////////////////////////////////////////////////////
static void app_destroy(app_state_t* state)
{
    if (state) {

        rt_framebuffer_free(&state->framebuffer);
        rt_world_free(&state->world);

        if (state->nc) {
            notcurses_stop(state->nc);
        }

        if (state->gamepad_info.gamepad) {
            SDL_CloseGamepad(state->gamepad_info.gamepad);
        }
    }

    SDL_Quit();
}

///////////////////////////////////////////////////////////////////////////
static void error_callback(const char*      filename,
                           uint32_t         line,
                           const char*      function_name,
                           const char*      message,
                           void*            userParam)
{
    app_state_t* state = (app_state_t*)userParam;
    app_destroy(state);

    fprintf(stderr, "ERROR: %s:%" PRIu32 " in %s: %s\n", filename,
                                                         line,
                                                         function_name,
                                                         message);
    exit(EXIT_FAILURE);
}

///////////////////////////////////////////////////////////////////////////
int main(void)
{
    app_state_t app = {};
    app.gamepad_info.deadzone = 1000;

    rt_set_error_callback(error_callback, &app);

    RT_ASSERT(SDL_Init(SDL_INIT_GAMEPAD));

    RT_ASSERT(setlocale(LC_ALL, ""));

    struct notcurses_options opts = {};
    app.nc = notcurses_core_init(&opts, stdout);
    RT_ASSERT(app.nc != NULL);

    struct ncplane* std = notcurses_stdplane(app.nc);
    RT_ASSERT(std != NULL);

    rt_notcurses_surface_t notcurses_surface = rt_notcurses_surface_create(std);

    RT_ASSERT(RT_STATUS_success == rt_framebuffer_create(notcurses_surface.cols,
                                                         notcurses_surface.rows,
                                                         &app.framebuffer));

    app.world.clear_color = (rt_vec4_t){ RT_FLOAT(0.02),
                                         RT_FLOAT(0.03),
                                         RT_FLOAT(0.03),
                                         RT_FLOAT(1.00), };

    app.world.face_cull_mode = RT_FACE_cull_back;

    rt_idx_t plane_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_plane(&app.world,
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

    rt_world_set_plane_params(&app.world, plane_index, &plane_params);

    rt_idx_t plane_material_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_checkerboard_material(&app.world,
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

    rt_world_set_checkerboard_material_params(&app.world,
                                              plane_material_index,
                                              &mat_params);

    rt_plane_link_checkerboard_material(&app.world,
                                        plane_index,
                                        plane_material_index);

    rt_idx_t sphere_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_sphere(&app.world,
                                                        &sphere_index));

    rt_sphere_params_t sphere_params = {

        .center = { RT_FLOAT( 0.0),
                    RT_FLOAT( 3.0),
                    RT_FLOAT(-10.0),
                    RT_FLOAT( 1.0), },

        .radius = RT_FLOAT(3.0),

    };

    rt_world_set_sphere_params(&app.world,
                               sphere_index,
                               &sphere_params);

    rt_idx_t sphere_material_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_metallic_material(&app.world,
                                                                   &sphere_material_index));

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


    rt_world_set_metallic_material_params(&app.world,
                                          sphere_material_index,
                                          &sphere_mat_params);

    rt_sphere_link_metallic_material(&app.world,
                                     sphere_index,
                                     sphere_material_index);

    rt_idx_t point_light_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_point_light(&app.world,
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

    rt_world_set_point_light_params(&app.world,
                                    point_light_index,
                                    &point_light_params);

    rt_idx_t point_light_sphere_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_sphere(&app.world,
                                                        &point_light_sphere_index));

    rt_sphere_params_t point_light_sphere_params = {
        .center = point_light_params.position,
        .radius = RT_FLOAT(0.5),
    };

    rt_world_set_sphere_params(&app.world,
                               point_light_sphere_index,
                               &point_light_sphere_params);

    rt_emissive_material_t emissive_params = {
        .color = { RT_FLOAT(1.0),
                   RT_FLOAT(1.0),
                   RT_FLOAT(1.0),
                   RT_FLOAT(1.0), },
    };

    rt_idx_t point_light_emissive_material_index = 0;
    RT_ASSERT(RT_STATUS_success == rt_world_push_emissive_material(&app.world,
                                                                   &point_light_emissive_material_index));

    rt_world_set_emissive_material_params(&app.world,
                                          point_light_emissive_material_index,
                                          &emissive_params);

    rt_sphere_link_emissive_material(&app.world,
                                     point_light_sphere_index,
                                     point_light_emissive_material_index);

    bool is_running = true;

    rt_fps_camera_t fps_camera = rt_fps_camera_create();

    rt_float_t point_light_radius       = RT_FLOAT(10.0);

    // rt_fps_camera_notcurses_keybindings_t keybindings = rt_fps_camera_notcurses_default_keybindings();
    rt_fps_camera_sdl3_joystick_keybindings_t keybindings = rt_fps_camera_sdl3_default_joystick_keybindings();

    rt_timer_t timer = {};
    rt_float_t total_time = RT_FLOAT(0.0);

    while (is_running) {

        if (rt_notcurses_surface_resize(&notcurses_surface,
                                        &app.framebuffer)) {
            continue;
        }

        rt_float_t delta_time = rt_timer_update(&timer, NULL, NULL);
        total_time += delta_time;

        point_light_params.position.x       = sin(total_time) * point_light_radius;
        point_light_params.position.z       = cos(total_time) * point_light_radius;
        point_light_sphere_params.center    = point_light_params.position;

        rt_world_set_point_light_params(&app.world,
                                        point_light_index,
                                        &point_light_params);

        rt_world_set_sphere_params(&app.world,
                                   point_light_sphere_index,
                                   &point_light_sphere_params);

        // rt_fps_camera_update_with_notcurses(&fps_camera,
        //                                     delta_time,
        //                                     &keybindings,
        //                                     app.nc,
        //                                     &is_running);

        rt_fps_camera_update_with_sdl3_joystick(&fps_camera,
                                                 &keybindings,
                                                 &app.gamepad_info,
                                                 delta_time,
                                                 &is_running);

        rt_fps_camera_render(&fps_camera,
                             &app.world,
                             &app.framebuffer,
                             delta_time);

        rt_notcurses_surface_blit(&notcurses_surface,
                                  &app.framebuffer);

        RT_ASSERT(-1 != notcurses_render(app.nc));

        rt_timer_wait(&timer, RT_FLOAT(60.0));
    }

    app_destroy(&app);

    return EXIT_SUCCESS;
}
