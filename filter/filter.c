#include <math.h>
#include "parameters.h"

/* Key:
 * del - change
 *
 * t - time
 *
 * r - roll
 * p - pitch
 * y - yaw
 *
 * g - gyroscope
 * a - accelerometer
 * m - magnetometer
*/

double r = 0;
double p = 0;
double y = 0;

double dr = 0;
double dp = 0;
double dy = 0;

/*
 * Update rate.  This is how often our state is updated with
 * gyro rate measurements.
 * eg .02 = 50 Hz rate
 */
static const double	dt	= dt_PARAM;

/*
 * Covariance matrix.  This is updated at every time step to
 * determine how well the sensors are tracking the actual state.
 */
static double P[3][3] = {	{ 1, 0, 0 },
							{ 0, 1, 0 },
							{ 0, 0, 1 },	};

/*
 * Our two states, the angle and the gyro bias.  As a byproduct of computing
 * the angle, we also have an unbiased angular rate available.   These are
 * read-only to the user of the module.
 */
double angle=0.0;
double q_bias=0.0;
double rate=0.0;

/*
 * R represents the measurement covariance noise in radians.
 * it is a 1x1 matrix that says that we expect 0.3 rad jitter
 * from the accelerometer.
 */
// !!! Adjusting R changes the speed at which the filter converges.
//     It is an indication of trust in the measurement
static const double	R_angle	= R_angle_PARAM;

/*
 * Q is a 3x3 matrix that represents the process covariance noise.
 * In this case, it indicates how much we trust the acceleromter
 * relative to the gyros.
 */
 // originally .001 and .003
static const double	Q_angle	= Q_angle_PARAM;
static const double	Q_gyro	= Q_gyro_PARAM;

/*
 * state_update is called every dt with a biased gyro measurement
 * by the user of the module.  It updates the current angle and
 * rate estimate.
 *
 * The pitch gyro measurement should be scaled into real units, but
 * does not need any bias removal.  The filter will track the bias.
 *
 * Our state vector is:
 *
 *	X = [ angle, gyro_bias ]
 *
 * It runs the state estimation forward via the state functions:
 *
 *	Xdot = [ angle_dot, gyro_bias_dot ]
 *
 *	angle_dot	= gyro - gyro_bias
 *	gyro_bias_dot	= 0
 *
 * And updates the covariance matrix via the function:
 *
 *	Pdot = A*P + P*A' + Q
 *
 * A is the Jacobian of Xdot with respect to the states:
 *
 *	A = [ d(angle_dot)/d(angle)     d(angle_dot)/d(gyro_bias) ]
 *	    [ d(gyro_bias_dot)/d(angle) d(gyro_bias_dot)/d(gyro_bias) ]
 *
 *	  = [ 0 -1 ]
 *	    [ 0  0 ]
 *
 * Due to the small CPU available on the microcontroller, we've
 * hand optimized the C code to only compute the terms that are
 * explicitly non-zero, as well as expanded out the matrix math
 * to be done in as few steps as possible.  This does make it harder
 * to read, debug and extend, but also allows us to do this with
 * very little CPU time.
 */
void state_update(const double q_m) /* Pitch gyro measurement */
{
	/* Unbias our gyro */
	const double q = q_m - q_bias;

	/*
	 * Compute the derivative of the covariance matrix
	 *
	 *	Pdot = A*P + P*A' + Q
	 *
	 * We've hand computed the expansion of A = [ 0 -1, 0 0 ] multiplied
	 * by P and P multiplied by A' = [ 0 0, -1, 0 ].  This is then added
	 * to the diagonal elements of Q, which are Q_angle and Q_gyro.
	 */
	const double Pdot[2 * 2] = {
		Q_angle - P[0][1] - P[1][0],	/* 0,0 */
		        - P[1][1],		/* 0,1 */
		        - P[1][1],		/* 1,0 */
		Q_gyro					/* 1,1 */
	};

	/* Store our unbias gyro estimate */
	rate = q;

	/*
	 * Update our angle estimate
	 * angle += angle_dot * dt
	 *       += (gyro - gyro_bias) * dt
	 *       += q * dt
	 */
	angle += q * dt;

	/* Update the covariance matrix */
	P[0][0] += Pdot[0] * dt;
	P[0][1] += Pdot[1] * dt;
	P[1][0] += Pdot[2] * dt;
	P[1][1] += Pdot[3] * dt;
}


/*
 * kalman_update is called by a user of the module when a new
 * accelerometer measurement is available.  ax_m and az_m do not
 * need to be scaled into actual units, but must be zeroed and have
 * the same scale.
 *
 * This does not need to be called every time step, but can be if
 * the accelerometer data are available at the same rate as the
 * rate gyro measurement.
 *
 * For a two-axis accelerometer mounted perpendicular to the rotation
 * axis, we can compute the angle for the full 360 degree rotation
 * with no linearization errors by using the arctangent of the two
 * readings.
 *
 * As commented in state_update, the math here is simplified to
 * make it possible to execute on a small microcontroller with no
 * floating point unit.  It will be hard to read the actual code and
 * see what is happening, which is why there is this extensive
 * comment block.
 *
 * The C matrix is a 1x2 (measurements x states) matrix that
 * is the Jacobian matrix of the measurement value with respect
 * to the states.  In this case, C is:
 *
 *	C = [ d(angle_m)/d(angle)  d(angle_m)/d(gyro_bias) ]
 *	  = [ 1 0 ]
 *
 * because the angle measurement directly corresponds to the angle
 * estimate and the angle measurement has no relation to the gyro
 * bias.
 */

double		angle_err;
double		C_0;
//double 		temp1, temp2;

/* angle_m: tilt determined from acceleration sensor */
void kalman_update(const double angle_m)
{
/* Compute our measured angle and the error in our estimate */
//const double		angle_m = atan2( -az_m, ax_m );
// !!!Changed since only 1 accel DOF
//	double angle_m;
//	angle_m = scale_accel((int)ax_m); // asin( ax_m / ACCEL_SCALE );
///////////////////////////////////////////////////////////////////

//#define GYRO_OFFSET 613		// value at 0 degrees / sec - unimportant, as the Kalman filter corrects this
//#define ACCEL_OFFSET 650	// value at 0 degrees

//#define GYRO_SCALE 0.8		// gyro is 2 mv/(deg/sec), ADC is 2.5mv/tick: 2/2.5 ; degrees/sec = adc * .8
//#define ACCEL_SCALE 135.0	// ADC bits per 90 degrees

//#define RAD_TO_DEG (180 / 3.1415926535897932384626433832795)

///////////////////////////////////////////////////////////////////
//	temp1 = ax_m - ACCEL_OFFSET;
//	if (temp1 < -ACCEL_SCALE) temp1 = -ACCEL_SCALE;
//	if (temp1 > ACCEL_SCALE)  temp1 =  ACCEL_SCALE;
//	temp2 = temp1 / ACCEL_SCALE;
//	angle_m = asin(temp2);
//	angle_m *= RAD_TO_DEG;
///////////////////////////////////////////////////////////////////

	angle_err = angle_m - angle;
	/*
	 * C_0 shows how the state measurement directly relates to
	 * the state estimate.
	 *
	 * The C_1 shows that the state measurement does not relate
	 * to the gyro bias estimate.  We don't actually use this, so
	 * we comment it out.
	 */
	C_0 = 1;
	/* const double		C_1 = 0; */

	/*
	 * PCt<2,1> = P<2,2> * C'<2,1>, which we use twice.  This makes
	 * it worthwhile to precompute and store the two values.
	 * Note that C[0,1] = C_1 is zero, so we do not compute that
	 * term.
	 */
	const double		PCt_0 = C_0 * P[0][0]; /* + C_1 * P[0][1] = 0 */
	const double		PCt_1 = C_0 * P[1][0]; /* + C_1 * P[1][1] = 0 */

	/*
	 * Compute the error estimate.  From the Kalman filter paper:
	 *
	 *	E = C P C' + R
	 *
	 * Dimensionally,
	 *
	 *	E<1,1> = C<1,2> P<2,2> C'<2,1> + R<1,1>
	 *
	 * Again, note that C_1 is zero, so we do not compute the term.
	 */
	const double		E = R_angle + C_0 * PCt_0;
	/*	+ C_1 * PCt_1 = 0 */

	/*
	 * Compute the Kalman filter gains.  From the Kalman paper:
	 *
	 *	K = P C' inv(E)
	 *
	 * Dimensionally:
	 *
	 *	K<2,1> = P<2,2> C'<2,1> inv(E)<1,1>
	 *
	 * Luckilly, E is <1,1>, so the inverse of E is just 1/E.
	 */
	const double		K_0 = PCt_0 / E;
	const double		K_1 = PCt_1 / E;

	/*
	 * Update covariance matrix.  Again, from the Kalman filter paper:
	 *
	 *	P = P - K C P
	 *
	 * Dimensionally:
	 *
	 *	P<2,2> -= K<2,1> C<1,2> P<2,2>
	 *
	 * We first compute t<1,2> = C P.  Note that:
	 *
	 *	t[0,0] = C[0,0] * P[0,0] + C[0,1] * P[1,0]
	 *
	 * But, since C_1 is zero, we have:
	 *
	 *	t[0,0] = C[0,0] * P[0,0] = PCt[0,0]
	 *
	 * This saves us a floating point multiply.
	 */
	const double		t_0 = PCt_0; /* C_0 * P[0][0] + C_1 * P[1][0] */
	const double		t_1 = C_0 * P[0][1]; /* + C_1 * P[1][1]  = 0 */

	P[0][0] -= K_0 * t_0;
	P[0][1] -= K_0 * t_1;
	P[1][0] -= K_1 * t_0;
	P[1][1] -= K_1 * t_1;

	/*
	 * Update our state estimate.  Again, from the Kalman paper:
	 *
	 *	X += K * err
	 *
	 * And, dimensionally,
	 *
	 *	X<2> = X<2> + K<2,1> * err<1,1>
	 *
	 * err is a measurement of the difference in the measured state
	 * and the estimate state.  In our case, it is just the difference
	 * between the two accelerometer measured angle and our estimated
	 * angle.
	 */
	angle	+= K_0 * angle_err;
	q_bias	+= K_1 * angle_err;
}
