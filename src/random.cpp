#include "include/random.h"

// Seed random number generator from hardware device.
Engine engine{std::random_device{}()};
