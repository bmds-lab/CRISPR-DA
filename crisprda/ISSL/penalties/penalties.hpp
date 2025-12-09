#ifndef otScorePenalties
#define otScorePenalties
#include <unordered_map>
#include "phmap/phmap.h"

extern const phmap::parallel_flat_hash_map<uint64_t, double> precalculatedMITScores;
extern const double cfdPosPenalties[320];
extern const double cfdPamPenalties[16];

#endif // !otScorePenalties
