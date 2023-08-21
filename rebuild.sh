#!/bin/bash -xe

git submodule sync
git submodule update --init

case "$1" in
    "fmem")
        CFLAGS=' -g -Wall -Wpedantic -Wextra -DRUNMEMTEST '
        LDFLAGS=' -fsanitize=leak  '
        ;;
    "leak")
        CFLAGS=' -g -Wall -Wpedantic -Wextra '
        LDFLAGS=' -fsanitize=leak  '
        ;;
    "opt")
        CFLAGS=' -O3 -Wall -Wpedantic -Wextra '
        ;;
    "undefined")
        CFLAGS=' -g -Wall -Wpedantic -Wextra -fsanitize=undefined '
        LDFLAGS=' -fsanitize=undefined  '
        ;;
    *)
        CFLAGS=' -g -Wall -Wpedantic -Wextra '
        ;;
esac

./autogen.sh
./configure CFLAGS="${CFLAGS}" \
            LDFLAGS="${LDFLAGS}" \
             --enable-openmp \
             --localstatedir=/var \
             --prefix=/usr \
             --sysconfdir=/etc

make clean
make && sudo make install
