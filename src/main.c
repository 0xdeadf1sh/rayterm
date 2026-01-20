#include "rt_vec3.h"
#include "rt_ray.h"

#include <notcurses/notcurses.h>

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <locale.h>
#include <unistd.h>
#include <time.h>

#define CAMERA_SPEED 0.1f
#define BUFFER_LEN(BUFFER) (sizeof(BUFFER) / sizeof((BUFFER)[0]))

static float compute_gamma(float x, float gamma)
{
    return powf(x, gamma);
}

static bool hit_sphere(rt_vec3_t center, float radius, rt_ray_t ray)
{
    rt_vec3_t oc = rt_vec3_sub(center, ray.org);
    float a = rt_vec3_dot(ray.dir, ray.dir);
    float b = -2.0f * rt_vec3_dot(ray.dir, oc);
    float c = rt_vec3_dot(oc, oc);
    c -= radius * radius;
    float discriminant = b * b - 4 * a * c;
    return discriminant >= 0.0f;
}

static rt_vec3_t compute_color(rt_ray_t ray)
{
    typedef struct {
        rt_vec3_t center;
        rt_vec3_t color;
        float radius;
    }
    sphere_t;

    sphere_t spheres[100];
    for (uint32_t i = 0; i < 10; ++i) {
        for (uint32_t j = 0; j < 10; ++j) {
            sphere_t sphere = {
                .center = rt_vec3_create((float)i * 2.0f - 10.0f, -2.0f, -(float)j * 2.0f - 2.0f),
                .color = rt_vec3_create((float)i / 9.0f, (float)j / 9.0f, 1.0f),
                .radius = 0.5f,
            };
            spheres[i * 10 + j] = sphere;
        }
    }

    for (uint32_t i = 0; i < BUFFER_LEN(spheres); ++i) {
        if (hit_sphere(spheres[i].center, spheres[i].radius, ray)) {
            return spheres[i].color;
        }
    }

    rt_vec3_t ray_dir_norm = rt_vec3_norm(ray.dir);
    rt_vec3_t white_color = { 1.0f, 1.0f, 1.0f };
    rt_vec3_t blue_color = { 0.0f, 0.0f, 1.0f };
    float k = (ray_dir_norm.y + 1.0f) * 0.5f;
    return rt_vec3_lerp(white_color, blue_color, k);
}

#define GAMMA       2.2f
#define GAMMA_INV   (1.0f / 2.2f)

int main([[maybe_unused]] int argc, char** argv)
{
    if (!setlocale(LC_ALL, "")) {
        fprintf(stderr, "%s: setlocale(LC_ALL, \"\") failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    struct notcurses_options opts = {};
    struct notcurses* nc = notcurses_core_init(&opts, stdout);
    if (!nc) {
        fprintf(stderr, "%s: notcurses_init() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    uint32_t rows = 0;
    uint32_t cols = 0;
    struct ncplane* nstd = notcurses_stddim_yx(nc, &rows, &cols);
    if (!nstd) {
        fprintf(stderr, "%s: notcurses_stdplane() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t pixelbuffer_size = rows * cols * sizeof(uint32_t);
    uint32_t* rgb = malloc(pixelbuffer_size);
    if (!rgb) {
        fprintf(stderr, "%s: malloc() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Camera
    float aspect = (float)cols / (float)(rows * 2);
    float focal_length = 1.0f;
    float viewport_height = 2.0f;
    float viewport_width = viewport_height * aspect;

    rt_vec3_t camera_center = {
        .x = 0.0f,
        .y = 0.0f,
        .z = 0.0f,
    };

    rt_vec3_t viewport_u = {
        .x = viewport_width,
        .y = 0.0f,
        .z = 0.0f,
    };

    rt_vec3_t viewport_v = {
        .x = 0.0f,
        .y = -viewport_height,
        .z = 0.0f,
    };

    rt_vec3_t pixel_delta_u = rt_vec3_div_scalar(viewport_u, (float)cols);
    rt_vec3_t pixel_delta_v = rt_vec3_div_scalar(viewport_v, (float)rows);
    rt_vec3_t pixel_delta_diag = rt_vec3_mul_scalar(rt_vec3_add(pixel_delta_u, pixel_delta_v), 0.5f);

    rt_vec3_t viewport_u_half = rt_vec3_mul_scalar(viewport_u, 0.5f);
    rt_vec3_t viewport_v_half = rt_vec3_mul_scalar(viewport_v, 0.5f);

    struct ncvisual_options vopts = {};
    vopts.n = nstd;
    vopts.leny = rows;
    vopts.lenx = cols;
    vopts.blitter = NCBLIT_1x1;
    vopts.scaling = NCSCALE_NONE;

    bool is_moving_forward = false;
    bool is_moving_backward = false;
    bool is_moving_left = false;
    bool is_moving_right = false;
    bool is_moving_up = false;
    bool is_moving_down = false;

    for (;;) {
        struct ncinput input = {};
        notcurses_get(nc, NULL, &input);

        if (input.id == 'z' && input.evtype == NCTYPE_RELEASE) {
            break;
        }

        if (input.id == 'w') {
            if (input.evtype == NCTYPE_PRESS) {
                is_moving_forward = true;
                is_moving_backward = false;
            }
            else if (input.evtype == NCTYPE_RELEASE) {
                is_moving_forward = false;
            }
        }
        else if (input.id == 's') {
            if (input.evtype == NCTYPE_PRESS) {
                is_moving_backward = true;
                is_moving_forward = false;
            }
            else if (input.evtype == NCTYPE_RELEASE) {
                is_moving_backward = false;
            }
        }

        if (input.id == 'a') {
            if (input.evtype == NCTYPE_PRESS) {
                is_moving_left = true;
                is_moving_right = false;
            }
            else if (input.evtype == NCTYPE_RELEASE) {
                is_moving_left = false;
            }
        }
        else if (input.id == 'd') {
            if (input.evtype == NCTYPE_PRESS) {
                is_moving_right = true;
                is_moving_left = false;
            }
            else if (input.evtype == NCTYPE_RELEASE) {
                is_moving_right = false;
            }
        }

        if (input.id == 'q') {
            if (input.evtype == NCTYPE_PRESS) {
                is_moving_up = true;
                is_moving_down = false;
            }
            else if (input.evtype == NCTYPE_RELEASE) {
                is_moving_up = false;
            }
        }
        else if (input.id == 'e') {
            if (input.evtype == NCTYPE_PRESS) {
                is_moving_down = true;
                is_moving_up = false;
            }
            else if (input.evtype == NCTYPE_RELEASE) {
                is_moving_down = false;
            }
        }

        if (is_moving_forward) {
            camera_center.z -= CAMERA_SPEED;
        }
        else if (is_moving_backward) {
            camera_center.z += CAMERA_SPEED;
        }

        if (is_moving_left) {
            camera_center.x -= CAMERA_SPEED;
        }
        else if (is_moving_right) {
            camera_center.x += CAMERA_SPEED;
        }

        if (is_moving_up) {
            camera_center.y += CAMERA_SPEED;
        }
        else if (is_moving_down) {
            camera_center.y -= CAMERA_SPEED;
        }

        rt_vec3_t viewport_upper_left = rt_vec3_sub(camera_center, rt_vec3_create(0.0f, 0.0f, focal_length));
        viewport_upper_left = rt_vec3_sub(viewport_upper_left, viewport_u_half);
        viewport_upper_left = rt_vec3_sub(viewport_upper_left, viewport_v_half);

        rt_vec3_t pixel00_loc = rt_vec3_add(viewport_upper_left, pixel_delta_diag);

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

                rt_vec3_t pixel_color = compute_color(r);
                pixel_color = rt_vec3_apply_2(pixel_color, compute_gamma, GAMMA_INV);
                rgb[row * cols + col] = rt_vec3_to_uint32_alpha(pixel_color, 1.0f);
            }
        }

        if (-1 == ncblit_rgba(rgb, (int32_t)(cols * sizeof(uint32_t)), &vopts)) {
            fprintf(stderr, "%s: ncblit_rgba() failed\n", argv[0]);
            return EXIT_FAILURE;
        }

        if (-1 == notcurses_render(nc)) {
            fprintf(stderr, "%s: notcurses_render() failed\n", argv[0]);
            return EXIT_FAILURE;
        }

        usleep(10000);
    }
    free(rgb);
    return notcurses_stop(nc);
}
