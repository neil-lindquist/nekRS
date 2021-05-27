
#include "gslib.h"
#include "nrs.hpp"

#define D 2
#define WHEN_3D(a)
#include "interp_imp.h"
#undef WHEN_3D
#undef D

#define D 3
#define WHEN_3D(a) a
#include "interp_imp.h"
#undef WHEN_3D
#undef D
