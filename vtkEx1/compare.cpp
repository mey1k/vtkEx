int CompareDoubleAbsoulte(double *x, double *y)
{
	double absTolerance = 0.00000000000000001;

	double diffX = x[0] - y[0];
	double diffY = x[1] - y[1];
	double diffZ = x[2] - y[2];

	if (diffX < -absTolerance || diffX > absTolerance || diffY < -absTolerance || diffY > absTolerance || diffZ < -absTolerance || diffZ > absTolerance)
		return 0;
	else
		return 1;
}