workspace "rayterm"

    configurations {
        "debug",
        "release",
        "tsan",
    }

    location "build"

project "rayterm"
    kind "ConsoleApp"
    language "C"

    targetdir "build/%{cfg.buildcfg}"

    files {
        "src/main.c",
    }

    includedirs {
        "include"
    }

    links {
        "notcurses-core",
    }

    defines {
        "_DEFAULT_SOURCE",
        "_XOPEN_SOURCE=600",
    }

    -- minimum: GCC 15
    filter { "system:linux", "action:gmake" }
        buildoptions {
            "`pkg-config --cflags notcurses-core`",
            "-std=c23",
            "-g",
            "-save-temps",
            "-march=native",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            "-Whardened",
            "-Werror",
            "-Wformat",
            "-Wformat=2",
            "-Wformat-security",
            "-Wconversion",
            "-Wpadded",
            "-Wsign-conversion",
            "-Wimplicit-fallthrough",
            "-Wfloat-equal",
            "-Wtrampolines -fzero-init-padding-bits=all",
            "-Wbidi-chars=any",
            "-Wundef",
            "-Wshadow",
            "-Wpointer-arith",
            "-Wcast-align",
            "-Wstrict-prototypes",
            "-Wmissing-prototypes",
            "-Wmissing-declarations",
            "-Wredundant-decls",
            "-Wold-style-definition",
            "-Wwrite-strings",
            "-Waggregate-return",
            "-Wcast-qual",
            "-Wswitch-default",
            "-Wunreachable-code",
            "-Winit-self",
            "-Wimplicit-function-declaration",
            "-Wdeprecated",
            "-U_FORTIFY_SOURCE",
            "-D_FORTIFY_SOURCE=3",
            "-fstrict-flex-arrays=3",
            "-fstack-clash-protection",
            "-fstack-protector-strong",
            "-Wl,-z,nodlopen",
            "-Wl,-z,noexecstack",
            "-Wl,-z,relro",
            "-Wl,-z,now",
            "-Wl,--as-needed",
            "-Wl,--no-copy-dt-needed-entries",
            "-fPIE -pie",
            "-fcf-protection=full", -- x86_64
            -- "-mbrach-protection=standard", -- aarch64
            "-fno-delete-null-pointer-checks",
            "-fno-strict-overflow",
            "-fno-strict-aliasing",
            "-ftrivial-auto-var-init=zero",
            "-fverbose-asm",
            "-masm=intel",
        }

        linkoptions {
            "`pkg-config --libs notcurses-core`",
            "-Wl,-z,nodlopen",
            "-Wl,-z,noexecstack",
            "-Wl,-z,relro",
            "-Wl,-z,now",
            "-Wl,--as-needed",
            "-Wl,--no-copy-dt-needed-entries",
            "-fPIE -pie",
        }

    filter "configurations:debug"
        symbols "On"

        defines {
            "_DEBUG",
        }

        buildoptions {
            "-O1",
            "-fno-omit-frame-pointer",
            "-fno-optimize-sibling-calls",
            "-fno-common",
            "-fsanitize=address,undefined",
        }

        linkoptions {
            "-fsanitize=address,undefined",
        }

    filter "configurations:release"
        symbols "On"

        defines {
            "NDEBUG",
        }

        buildoptions {
            "-O2",
            "-fsanitize-trap=undefined",
        }

        linkoptions {
            "-fsanitize-trap=undefined",
            "-flto",
        }

    filter "configurations:tsan"
        symbols "On"

        defines {
            "NDEBUG",
        }

        buildoptions {
            "-O2",
            "-fsanitize=thread",
        }

        linkoptions {
            "-fsanitize=thread",
        }
