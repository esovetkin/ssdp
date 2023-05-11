sudo make uninstall
make clean
make distclean
 ./configure CFLAGS="-Wextra -Wpedantic -Wall -Og -g -I /usr/include/hdf5/serial" CC=h5cc
 make
 sudo make install
