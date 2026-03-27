#!/bin/bash

watchexec -r -c -e c,h --wrap-process=none -- "make -C build config=release && ./build/release/$1"
