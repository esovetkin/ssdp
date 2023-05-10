#include "primes_table.h"
#include "_primes.h"

unsigned int get_nth_prime(unsigned int index){
    if(index >= N_PRIMES){
        return 0;
    }
    return PRIMES[index];
}

unsigned int get_next_biggest_prime(unsigned int value){
    for(unsigned int i = 0; i< N_PRIMES; i++){
        if(value<PRIMES[i]){
            return PRIMES[i];
        }
    }
    return 0;
}