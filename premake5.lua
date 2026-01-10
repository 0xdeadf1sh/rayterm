workspace "rayterm"
    configurations { "debug", "release" }
    location "build"

project "rayterm"
    kind "ConsoleApp"
    language "C"
    targetdir "build/%{cfg.buildcfg}"

    files {
        "src/main.c"
    }

    links {
        "notcurses-core"
    }

    filter "configurations:debug"
        defines { "_DEBUG" }
        symbols "On"

    filter "configurations:release"
        defines { "NDEBUG" }
        optimize "On"
