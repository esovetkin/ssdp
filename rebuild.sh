#!/bin/bash -xe

git submodule sync
git submodule update --init

case "$1" in
    "leak")
        CFLAGS=' -g -Wall -Wpedantic -Wextra '
        LDFLAGS=' -fsanitize=leak  '
        ;;
    "analyzer")
        CFLAGS=' -g -Wall -Wpedantic -Wextra -fanalyzer -Wno-analyzer-imprecise-fp-arithmetic '
        LDFLAGS=' -fsanitize=leak  '
        ;;
    "address")
        CFLAGS=' -g -Wall -Wpedantic -Wextra -fsanitize=address '
        LDFLAGS=' -fsanitize=address  '
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
make
sudo make install
