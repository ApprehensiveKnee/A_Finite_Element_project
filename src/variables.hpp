//===========================================================
//             HEADER FILE FOR SOME GLOBAL VARIABLES
//===========================================================

#ifndef VAR
#define VAR

extern inline constexpr double tol = 1.e-16;
// Degree of the FE spaces, along all possible dimensions
extern inline constexpr unsigned r = 3;
// Dimenstionality of the problem
extern inline constexpr unsigned short DIM = 2;
// Number of threads to be used in the parallel implementations
inline unsigned short THREADS = 4;
//  A flag for preformance comparison
inline bool PERFORMANCE_COMPARISON = false;

#endif