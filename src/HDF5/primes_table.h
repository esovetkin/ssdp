#pragma once
#include <stddef.h>

extern unsigned const int PRIMES[];
extern const size_t N_PRIMES;
/*
    get the n-th prime number

    args:
        index: index of the prime e.g.
        prime: 2,3,5,7
        index: 0,1,2,3
    return:
        prime if successfull else 0

    This function reads from a table therefor it may fail.
*/
unsigned int get_nth_prime(unsigned int index);

/*
    get the next biggest prime number which comes after value
    
    args:
        value: a positive number
    return:
        next biggest prime if successfull else 0
    
    This function reads from a table therefor it may fail.

    example input-output pairs:
    input:  0   1   2   8   
    output: 2   2   2   11
*/
unsigned int get_next_biggest_prime(unsigned int value);