#include <xxhash.h>
#include <string.h>

unsigned long long xxh3_hash_c(const char* key, int length) {
    return XXH3_64bits(key, length);
}