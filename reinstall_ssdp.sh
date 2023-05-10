sudo make uninstall
make clean
make distclean
 ./configure CFLAGS="-Wextra -Wpedantic -Wall -Og -g" CC=h5cc
 make
 sudo make install
