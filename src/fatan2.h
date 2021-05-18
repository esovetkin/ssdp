double fatan2(double y, double x);
double fatan(double x);
#ifdef FAST_ATAN2
	#define ATAN2(y,x) (fatan2((y),(x)))
	#define ATAN(x) (fatan(x))
#else
	#define ATAN2(y,x) (atan2((y),(x)))
	#define ATAN(x) (atan(x))
#endif	
