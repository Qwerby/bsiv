package rubris.io;




/* 
 * A fast and efficient calculator for Black's implied volatility based on reduced iterations.
 * 
 * The original paper: http://jaeckel.16mb.com/LetsBeRational.pdf
 * 
 
 The Java code is a direct port of the c++ code http://jaeckel.16mb.com/LetsBeRational.7z
 
 Original reference from:
 
   https://github.com/vollib/vollib/blob/master/vollib/black_scholes_merton/implied_volatility.py

*************** Original copyright notice ********************

Copyright © 2013-2014 Peter Jäckel.

Permission to use, copy, modify, and distribute this software is freely granted,
provided that this notice is preserved.

WARRANTY DISCLAIMER
The Software is provided "as is" without warranty of any kind, either express or implied,
including without limitation any implied warranties of condition, uninterrupted use,
merchantability, fitness for a particular purpose, or non-infringement.

*/

public class BlackScholesJaekelIVCalculator{

	
	static final double ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887;
	static final double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;
	static final double SQRT_TWO_PI = 2.506628274631000502415765284811045253006986740610;
	static final double DBL_EPSILON = Math.ulp(1.0);
	static final double SQRT_DBL_EPSILON = Math.sqrt(DBL_EPSILON);
	static final double FOURTH_ROOT_DBL_EPSILON = Math.sqrt(SQRT_DBL_EPSILON);
	static final double EIGHTH_ROOT_DBL_EPSILON = Math.sqrt(FOURTH_ROOT_DBL_EPSILON);
	static final double SIXTEENTH_ROOT_DBL_EPSILON = Math.sqrt(EIGHTH_ROOT_DBL_EPSILON);
	static final double SQRT_DBL_MIN = Math.sqrt(Double.MIN_VALUE);
	static final double SQRT_DBL_MAX = Math.sqrt(Double.MAX_VALUE);

	static final double DBL_MAX = Double.MAX_VALUE;
	static final double DBL_MIN = Double.MIN_VALUE;

	static final double DENORMALIZATION_CUTOFF = 0;

	static final double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC = -DBL_MIN;
	static final double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM = DBL_MAX;

	static final double asymptotic_expansion_accuracy_threshold = -10;
	static final double small_t_expansion_of_normalized_black_threshold = 2 * SIXTEENTH_ROOT_DBL_EPSILON;

	static final double norm_cdf_asymptotic_expansion_first_threshold = -10.0;
	static final double norm_cdf_asymptotic_expansion_second_threshold = -1 / Math.sqrt(DBL_EPSILON);

	

	 
	/* Dividend yield and rate conversions from implied_volatility.py */
	
	/**
	 * Generates an Implied Volatility from the inputs defined.
	 * 
	 * 
	 * @param optionPrice - price of the option from the market
	 * @param underlyingPrice - price of the underlying
	 * @param strike - strike price of the option
	 * @param days - number of days until expiry
	 * @param riskFreeRate - interest rate (expressed where 1 == 100%). i.e 5% is 0.05.
	 * @param dividendYield - from the market (expressed where 1 == 100%). i.e 3% is 0.03.
	 * @param call - whether it is a call or a put
	 * @return the Implied Volatility of the option
	 */
	public double calculate(double optionPrice, double underlyingPrice, double strike, int days, double riskFreeRate, double dividendYield, boolean call ){
		double t = ((double)days)/365.0;
		double conversion_factor = Math.exp(-riskFreeRate*t);
		double adjustedPrice = optionPrice / conversion_factor;
		underlyingPrice = underlyingPrice * Math.exp((riskFreeRate-dividendYield)*t);
		double type = ((call)?1.0:-1.0);
		return implied_volatility_from_a_transformed_rational_guess(adjustedPrice, underlyingPrice, strike, t,type);
	}
	 
	
	
	
	/* Methods from LetsBeRational.cpp */
	protected double implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(double price, double F,
			double K, double T, double q /* q=±1 */, int N) {
		final double intrinsic = Math.abs(Math.max((q < 0 ? K - F : F - K), 0.0));
		if (price < intrinsic)
			return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC;
		final double max_price = (q < 0 ? K : F);
		if (price >= max_price)
			return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
		final double x = Math.log(F / K);
		// Map in-the-money to out-of-the-money
		if (q * x > 0) {
			price = Math.abs(Math.max(price - intrinsic, 0.0));
			q = -q;
		}
		return unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
				price / (Math.sqrt(F) * Math.sqrt(K)), x, q, N) / Math.sqrt(T);
	}

	public double implied_volatility_from_a_transformed_rational_guess(double price, double F, double K, double T,
			double q /* q=±1 */) {
		return implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price, F, K, T, q, 2);
	}

	protected double unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(double beta, double x, double q /* q=±1 */, int N){
		   // Subtract intrinsic.
		   if (q*x>0) {
		      beta = Math.abs(Math.max(beta-normalised_intrinsic(x, q),0.));
		      q = -q;
		   }
		   // Map puts to calls
		   if (q<0){
		      x = -x;
		      q = -q;
		   }
		   if (beta<=0) // For negative or zero prices we return 0.
		      return 0;
		   if (beta<DENORMALIZATION_CUTOFF) // For positive but denormalized (a.k.a. 'subnormal') prices, we return 0 since it would be impossible to converge to full machine accuracy anyway.
		      return 0;
		   final double b_max = Math.exp(0.5*x);
		   if (beta>=b_max)
		      return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
		   int iterations=0, direction_reversal_count = 0;
		   double f=-DBL_MAX, s=-DBL_MAX, ds=s, ds_previous=0, s_left=DBL_MIN, s_right=DBL_MAX;
		   // The temptation is great to use the optimised form b_c = exp(x/2)/2-exp(-x/2)·Phi(sqrt(-2·x)) but that would require implementing all of the above types of round-off and over/underflow handling for this expression, too.
		   final double s_c=Math.sqrt(Math.abs(2*x)), b_c = normalised_black_call(x,s_c), v_c = normalised_vega(x, s_c);
		   // Four branches.
		   if ( beta<b_c ) {
		      final double s_l = s_c - b_c/v_c, b_l = normalised_black_call(x,s_l);
		      if (beta<b_l){
		    	  // amended for input params in originl C++
		         double f_lower_map_l =f_lower_map(x,s_l);
		         double d_f_lower_map_l_d_beta =d_f_lower_map_d_beta(x,s_l);
		         double d2_f_lower_map_l_d_beta2 =d2_f_lower_map_d_beta2(x,s_l);
		        
		         
		       //  compute_f_lower_map_and_first_two_derivatives(x,s_l,f_lower_map_l,d_f_lower_map_l_d_beta,d2_f_lower_map_l_d_beta2);
		       
		         final double r_ll=convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0.,b_l,0.,f_lower_map_l,1.,d_f_lower_map_l_d_beta,d2_f_lower_map_l_d_beta2,true);
		         f = rational_cubic_interpolation(beta,0.,b_l,0.,f_lower_map_l,1.,d_f_lower_map_l_d_beta,r_ll);
		         if (!(f>0)) { // This can happen due to roundoff truncation for extreme values such as |x|>500.
		            // We switch to quadratic interpolation using f(0)≡0, f(b_l), and f'(0)≡1 to specify the quadratic.
		            final double t = beta/b_l;
		            f = (f_lower_map_l*t + b_l*(1-t)) * t;
		         }
		         s = inverse_f_lower_map(x,f);
		         s_right = s_l;
		         //
		         // In this branch, which comprises the lowest segment, the objective function is
		         //     g(s) = 1/ln(b(x,s)) - 1/ln(beta)
		         //          ≡ 1/ln(b(s)) - 1/ln(beta)
		         // This makes
		         //              g'               =   -b'/(b·ln(b)²)
		         //              newton = -g/g'   =   (ln(beta)-ln(b))·ln(b)/ln(beta)·b/b'
		         //              halley = g''/g'  =   b''/b'  -  b'/b·(1+2/ln(b))
		         //              hh3    = g'''/g' =   b'''/b' +  2(b'/b)²·(1+3/ln(b)·(1+1/ln(b)))  -  3(b''/b)·(1+2/ln(b))
		         //
		         // The Householder(3) iteration is
		         //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ]
		         //
		         for (; iterations<N && Math.abs(ds)>DBL_EPSILON*s; ++iterations){
		            if (ds*ds_previous<0)
		               ++direction_reversal_count;
		            if ( iterations>0 && ( 3==direction_reversal_count || !(s>s_left && s<s_right) ) ) {
		               // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
		               // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
		               s = 0.5*(s_left+s_right);
		               if (s_right-s_left<=DBL_EPSILON*s) break;
		               direction_reversal_count = 0;
		               ds = 0;
		            }
		            ds_previous=ds;
		            final double b = normalised_black_call(x,s), bp = normalised_vega(x, s);
		            if ( b>beta && s<s_right ) s_right=s; else if ( b<beta && s>s_left ) s_left=s; // Tighten the bracket if applicable.
		            if (b<=0||bp<=0) // Numerical underflow. Switch to binary nesting for this iteration.
		               ds = 0.5*(s_left+s_right)-s;
		            else {
		               final double ln_b=Math.log(b), ln_beta=Math.log(beta), bpob=bp/b, h=x/s, b_halley = h*h/s-s/4, newton = (ln_beta-ln_b)*ln_b/ln_beta/bpob, halley = b_halley-bpob*(1+2/ln_b);
		               final double b_hh3 = b_halley*b_halley-3* ((h/s) * (h/s))-0.25, hh3 = b_hh3+2* (bpob*bpob)*(1+3/ln_b*(1+1/ln_b))-3*b_halley*bpob*(1+2/ln_b);
		               ds = newton * householder_factor(newton,halley,hh3);
		            }
		            s += ds = Math.max(-0.5*s , ds );
		         }
		         return s;
		      } else {
		         final double v_l = normalised_vega(x, s_l), r_lm = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l,b_c,s_l,s_c,1/v_l,1/v_c,0.0,false);
		         s = rational_cubic_interpolation(beta,b_l,b_c,s_l,s_c,1/v_l,1/v_c,r_lm);
		         s_left = s_l;
		         s_right = s_c;
		      }
		   } else {
		      final double s_h = v_c>DBL_MIN ? s_c+(b_max-b_c)/v_c : s_c, b_h = normalised_black_call(x,s_h);
		      if(beta<=b_h){
		         final double v_h = normalised_vega(x, s_h), r_hm = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_c,b_h,s_c,s_h,1/v_c,1/v_h,0.0,false);
		         s = rational_cubic_interpolation(beta,b_c,b_h,s_c,s_h,1/v_c,1/v_h,r_hm);
		         s_left = s_c;
		         s_right = s_h;
		      } else {
		    	// changed from original c++ due to input parameters  
		         double f_upper_map_h = f_upper_map(s_h);
		         double d_f_upper_map_h_d_beta =d_f_upper_map_d_beta(x,s_h);
		         double d2_f_upper_map_h_d_beta2 =d2_f_upper_map_d_beta2(x,s_h);
		         
		     //    compute_f_upper_map_and_first_two_derivatives(x,s_h,f_upper_map_h,d_f_upper_map_h_d_beta,d2_f_upper_map_h_d_beta2);
		         
		         
		         if ( d2_f_upper_map_h_d_beta2>-SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2<SQRT_DBL_MAX ){
		            final double r_hh = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_h,b_max,f_upper_map_h,0.,d_f_upper_map_h_d_beta,-0.5,d2_f_upper_map_h_d_beta2,true);
		            f = rational_cubic_interpolation(beta,b_h,b_max,f_upper_map_h,0.,d_f_upper_map_h_d_beta,-0.5,r_hh);
		         }
		         if (f<=0) {
		            final double h=b_max-b_h, t=(beta-b_h)/h;
		            f = (f_upper_map_h*(1-t) + 0.5*h*t) * (1-t); // We switch to quadratic interpolation using f(b_h), f(b_max)≡0, and f'(b_max)≡-1/2 to specify the quadratic.
		         }
		         s = inverse_f_upper_map(f);
		         s_left = s_h;
		         if (beta>0.5*b_max) { // Else we better drop through and let the objective function be g(s) = b(x,s)-beta. 
		            //
		            // In this branch, which comprises the upper segment, the objective function is
		            //     g(s) = ln(b_max-beta)-ln(b_max-b(x,s))
		            //          ≡ ln((b_max-beta)/(b_max-b(s)))
		            // This makes
		            //              g'               =   b'/(b_max-b)
		            //              newton = -g/g'   =   ln((b_max-b)/(b_max-beta))·(b_max-b)/b'
		            //              halley = g''/g'  =   b''/b'  +  b'/(b_max-b)
		            //              hh3    = g'''/g' =   b'''/b' +  g'·(2g'+3b''/b')
		            // and the iteration is
		            //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ].
		            //
		            for (; iterations<N && Math.abs(ds)>DBL_EPSILON*s; ++iterations){
		               if (ds*ds_previous<0)
		                  ++direction_reversal_count;
		               if ( iterations>0 && ( 3==direction_reversal_count || !(s>s_left && s<s_right) ) ) {
		                  // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
		                  // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
		                  s = 0.5*(s_left+s_right);
		                  if (s_right-s_left<=DBL_EPSILON*s) break;
		                  direction_reversal_count = 0;
		                  ds = 0;
		               }
		               ds_previous=ds;
		               final double b = normalised_black_call(x,s), bp = normalised_vega(x, s);
		               if ( b>beta && s<s_right ) s_right=s; else if ( b<beta && s>s_left ) s_left=s; // Tighten the bracket if applicable.
		               if (b>=b_max||bp<=DBL_MIN) // Numerical underflow. Switch to binary nesting for this iteration.
		                  ds = 0.5*(s_left+s_right)-s;
		               else {
		                  final double b_max_minus_b = b_max-b, g = Math.log((b_max-beta)/b_max_minus_b), gp = bp/b_max_minus_b;
		                  final double b_halley = ( square(x/s) )/s-s/4, b_hh3 = b_halley*b_halley-3* ( square(x/(s*s)) )-0.25;
		                  final double newton = -g/gp, halley = b_halley+gp, hh3 = b_hh3+gp*(2*gp+3*b_halley);
		                  ds = newton * householder_factor(newton,halley,hh3);
		               }
		               s += ds = Math.max(-0.5*s , ds );
		            }
		            return s;
		         }
		      }
		   }
		   // In this branch, which comprises the two middle segments, the objective function is g(s) = b(x,s)-beta, or g(s) = b(s) - beta, for short.
		   // This makes
		   //              newton = -g/g'   =  -(b-beta)/b'
		   //              halley = g''/g'  =    b''/b'    =  x²/s³-s/4
		   //              hh3    = g'''/g' =    b'''/b'   =  halley² - 3·(x/s²)² - 1/4
		   // and the iteration is
		   //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ].
		   //
		   for (; iterations<N && Math.abs(ds)>DBL_EPSILON*s; ++iterations){
		      if (ds*ds_previous<0)
		         ++direction_reversal_count;
		      if ( iterations>0 && ( 3==direction_reversal_count || !(s>s_left && s<s_right) ) ) {
		         // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
		         // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
		         s = 0.5*(s_left+s_right);
		         if (s_right-s_left<=DBL_EPSILON*s) break;
		         direction_reversal_count = 0;
		         ds = 0;
		      }
		      ds_previous=ds;
		      final double b = normalised_black_call(x,s), bp = normalised_vega(x, s);
		      if ( b>beta && s<s_right ) s_right=s; else if ( b<beta && s>s_left ) s_left=s; // Tighten the bracket if applicable.
		      final double newton = (beta-b)/bp, halley = (square(x/s))/s-s/4, hh3 = halley*halley-3* square(x/(s*s) )-0.25;
		      s += ds = Math.max(-0.5*s , newton * householder_factor(newton,halley,hh3) );
		   }
		   return s;
		}
	
	
	
	protected boolean is_below_horizon(double x) {
		return Math.abs(x) < DENORMALIZATION_CUTOFF;
	}

	/* Normalisation functions */

	protected double normalised_intrinsic(double x, double q /* q=±1 */) {
		if (q * x <= 0)
			return 0;
		final double x2 = x * x;
		if (x2 < 98 * FOURTH_ROOT_DBL_EPSILON) {// The factor 98 is computed
												// from last coefficient:
												// √√92897280 = 98.1749
			return Math
					.abs(Math.max(
							(q < 0 ? -1 : 1) * x
									* (1 + x2 * ((1.0 / 24.0) + x2
											* ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2)))),
							0.0));
		}
		double b_max = Math.exp(0.5 * x), one_over_b_max = 1 / b_max;
		return Math.abs(Math.max((q < 0 ? -1 : 1) * (b_max - one_over_b_max), 0.));
	}

	// Calculation of
	//
	//	              b  =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
		//
	//	                    exp(-(h²+t²)/2)
	//	                 =  --------------- ·  [ Y(h+t) - Y(h-t) ]
	//	                        √(2π)
		// with
	//	           Y(z) := Φ(z)/φ(z)
	//
	// using an expansion of Y(h+t)-Y(h-t) for small t to twelvth order in t.
	// Theoretically accurate to (better than) precision  ε = 2.23E-16  when  h<=0  and  t < τ  with  τ := 2·ε^(1/16) ≈ 0.21.
	// The main bottleneck for precision is the coefficient a:=1+h·Y(h) when |h|>1 .
	protected double small_t_expansion_of_normalized_black_call(double h, double t){
	   // Y(h) := Φ(h)/φ(h) = √(π/2)·erfcx(-h/√2)
	   // a := 1+h·Y(h)  --- Note that due to h<0, and h·Y(h) -> -1 (from above) as h -> -∞, we also have that a>0 and a -> 0 as h -> -∞
	   // w := t² , h2 := h²
	   final double a = 1+h*(0.5*SQRT_TWO_PI)*erfcx_cody(-ONE_OVER_SQRT_TWO*h), w=t*t, h2=h*h;
	   final double expansion = 2*t*(a+w*((-1+3*a+a*h2)/6+w*((-7+15*a+h2*(-1+10*a+a*h2))/120+w*((-57+105*a+h2*(-18+105*a+h2*(-1+21*a+a*h2)))/5040+w*((-561+945*a+h2*(-285+1260*a+h2*(-33+378*a+h2*(-1+36*a+a*h2))))/362880+w*((-6555+10395*a+h2*(-4680+17325*a+h2*(-840+6930*a+h2*(-52+990*a+h2*(-1+55*a+a*h2)))))/39916800+((-89055+135135*a+h2*(-82845+270270*a+h2*(-20370+135135*a+h2*(-1926+25740*a+h2*(-75+2145*a+h2*(-1+78*a+a*h2))))))*w)/6227020800.0))))));
	   final double b = ONE_OVER_SQRT_TWO_PI*Math.exp((-0.5*(h*h+t*t)))*expansion;
	   return Math.abs(Math.max(b,0.0));
	}
	
	// Asymptotic expansion of
	//
	//	              b  =  Φ(h+t)·exp(x/2) - Φ(h-t)·exp(-x/2)
		// with
	//	              h  =  x/s   and   t  =  s/2
		// which makes
	//	              b  =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
		//
	//	                    exp(-(h²+t²)/2)
	//	                 =  ---------------  ·  [ Y(h+t) - Y(h-t) ]
	//	                        √(2π)
	// with
	//	           Y(z) := Φ(z)/φ(z)
		//
		// for large negative (t-|h|) by the aid of Abramowitz & Stegun (26.2.12) where Φ(z) = φ(z)/|z|·[1-1/z^2+...].
		// We define
	//	                     r
	//	         A(h,t) :=  --- · [ Y(h+t) - Y(h-t) ]
	//	                     t
	//
	// with r := (h+t)·(h-t) and give an expansion for A(h,t) in q:=(h/r)² expressed in terms of e:=(t/h)² .
	protected double asymptotic_expansion_of_normalized_black_call(double h, double t){
	   final double e=(t/h)*(t/h), r=((h+t)*(h-t)), q=(h/r)*(h/r);
	   // 17th order asymptotic expansion of A(h,t) in q, sufficient for Φ(h) [and thus y(h)] to have relative accuracy of 1.64E-16 for h <= η  with  η:=-10.
	   final double asymptotic_expansion_sum = (2.0+q*(-6.0E0-2.0*e+3.0*q*(1.0E1+e*(2.0E1+2.0*e)+5.0*q*(-1.4E1+e*(-7.0E1+e*(-4.2E1-2.0*e))+7.0*q*(1.8E1+e*(1.68E2+e*(2.52E2+e*(7.2E1+2.0*e)))+9.0*q*(-2.2E1+e*(-3.3E2+e*(-9.24E2+e*(-6.6E2+e*(-1.1E2-2.0*e))))+1.1E1*q*(2.6E1+e*(5.72E2+e*(2.574E3+e*(3.432E3+e*(1.43E3+e*(1.56E2+2.0*e)))))+1.3E1*q*(-3.0E1+e*(-9.1E2+e*(-6.006E3+e*(-1.287E4+e*(-1.001E4+e*(-2.73E3+e*(-2.1E2-2.0*e))))))+1.5E1*q*(3.4E1+e*(1.36E3+e*(1.2376E4+e*(3.8896E4+e*(4.862E4+e*(2.4752E4+e*(4.76E3+e*(2.72E2+2.0*e)))))))+1.7E1*q*(-3.8E1+e*(-1.938E3+e*(-2.3256E4+e*(-1.00776E5+e*(-1.84756E5+e*(-1.51164E5+e*(-5.4264E4+e*(-7.752E3+e*(-3.42E2-2.0*e))))))))+1.9E1*q*(4.2E1+e*(2.66E3+e*(4.0698E4+e*(2.3256E5+e*(5.8786E5+e*(7.05432E5+e*(4.0698E5+e*(1.08528E5+e*(1.197E4+e*(4.2E2+2.0*e)))))))))+2.1E1*q*(-4.6E1+e*(-3.542E3+e*(-6.7298E4+e*(-4.90314E5+e*(-1.63438E6+e*(-2.704156E6+e*(-2.288132E6+e*(-9.80628E5+e*(-2.01894E5+e*(-1.771E4+e*(-5.06E2-2.0*e))))))))))+2.3E1*q*(5.0E1+e*(4.6E3+e*(1.0626E5+e*(9.614E5+e*(4.08595E6+e*(8.9148E6+e*(1.04006E7+e*(6.53752E6+e*(2.16315E6+e*(3.542E5+e*(2.53E4+e*(6.0E2+2.0*e)))))))))))+2.5E1*q*(-5.4E1+e*(-5.85E3+e*(-1.6146E5+e*(-1.77606E6+e*(-9.37365E6+e*(-2.607579E7+e*(-4.01166E7+e*(-3.476772E7+e*(-1.687257E7+e*(-4.44015E6+e*(-5.9202E5+e*(-3.51E4+e*(-7.02E2-2.0*e))))))))))))+2.7E1*q*(5.8E1+e*(7.308E3+e*(2.3751E5+e*(3.12156E6+e*(2.003001E7+e*(6.919458E7+e*(1.3572783E8+e*(1.5511752E8+e*(1.0379187E8+e*(4.006002E7+e*(8.58429E6+e*(9.5004E5+e*(4.7502E4+e*(8.12E2+2.0*e)))))))))))))+2.9E1*q*(-6.2E1+e*(-8.99E3+e*(-3.39822E5+e*(-5.25915E6+e*(-4.032015E7+e*(-1.6934463E8+e*(-4.1250615E8+e*(-6.0108039E8+e*(-5.3036505E8+e*(-2.8224105E8+e*(-8.870433E7+e*(-1.577745E7+e*(-1.472562E6+e*(-6.293E4+e*(-9.3E2-2.0*e))))))))))))))+3.1E1*q*(6.6E1+e*(1.0912E4+e*(4.74672E5+e*(8.544096E6+e*(7.71342E7+e*(3.8707344E8+e*(1.14633288E9+e*(2.07431664E9+e*(2.33360622E9+e*(1.6376184E9+e*(7.0963464E8+e*(1.8512208E8+e*(2.7768312E7+e*(2.215136E6+e*(8.184E4+e*(1.056E3+2.0*e)))))))))))))))+3.3E1*(-7.0E1+e*(-1.309E4+e*(-6.49264E5+e*(-1.344904E7+e*(-1.4121492E8+e*(-8.344518E8+e*(-2.9526756E9+e*(-6.49588632E9+e*(-9.0751353E9+e*(-8.1198579E9+e*(-4.6399188E9+e*(-1.6689036E9+e*(-3.67158792E8+e*(-4.707164E7+e*(-3.24632E6+e*(-1.0472E5+e*(-1.19E3-2.0*e)))))))))))))))))*q)))))))))))))))));
	   final double b = ONE_OVER_SQRT_TWO_PI*Math.exp((-0.5*(h*h+t*t)))*(t/r)*asymptotic_expansion_sum;
	   return Math.abs(Math.max(b , 0.));
	}
	
	protected double normalised_intrinsic_call(double x) {
		return normalised_intrinsic(x, 1);
	}
	
	// b(x,s) = Φ(x/s+s/2)·exp(x/2) - Φ(x/s-s/2)·exp(-x/2)
	// = Φ(h+t)·exp(x/2) - Φ(h-t)·exp(-x/2)
	// with
	// h = x/s and t = s/2
	protected double normalized_black_call_using_norm_cdf(double x, double s){
		final double h = x/s, t = 0.5*s, b_max = Math.exp(0.5*x), b = norm_cdf(h + t) * b_max - norm_cdf(h - t) / b_max;
		return Math.abs(Math.max(b,0.0));
	}
	
	protected double normalised_black_call(double x, double s) {
		   if (x>0){
		      return normalised_intrinsic_call(x)+normalised_black_call(-x,s);
		   }
		   final double ax = Math.abs(x);
		   if (s<=ax*DENORMALIZATION_CUTOFF){
		      return normalised_intrinsic_call(x);
		   }
		   // Denote h := x/s and t := s/2.
		   // We evaluate the condition |h|>|η|, i.e., h<η  &&  t < τ+|h|-|η|  avoiding any divisions by s , where η = asymptotic_expansion_accuracy_threshold  and τ = small_t_expansion_of_normalized_black_threshold .
		   if ( x < s*asymptotic_expansion_accuracy_threshold  &&  0.5*s*s+x < s*(small_t_expansion_of_normalized_black_threshold+asymptotic_expansion_accuracy_threshold) )
		      // Region 1.
		      return asymptotic_expansion_of_normalized_black_call(x/s,0.5*s);
		   if ( 0.5*s < small_t_expansion_of_normalized_black_threshold ){
		      // Region 2.
		      return small_t_expansion_of_normalized_black_call(x/s,0.5*s);
		   }
		   // When b is more than, say, about 85% of b_max=exp(x/2), then b is dominated by the first of the two terms in the Black formula, and we retain more accuracy by not attempting to combine the two terms in any way.
		   // We evaluate the condition h+t>0.85  avoiding any divisions by s.
		   if ( x+0.5*s*s > s*0.85 ){
		      // Region 3.
		      return normalized_black_call_using_norm_cdf(x,s);
		   }
		   // Region 4.
		   return normalised_black_call_using_erfcx(x/s,0.5*s);
		}
	
	protected double normalised_black_call_using_erfcx(double h, double t) {
		   // Given h = x/s and t = s/2, the normalised Black function can be written as
		   //
		   //     b(x,s)  =  Φ(x/s+s/2)·exp(x/2)  -   Φ(x/s-s/2)·exp(-x/2)
		   //             =  Φ(h+t)·exp(h·t)      -   Φ(h-t)·exp(-h·t) .                     (*)
		   //
		   // It is mentioned in section 4 (and discussion of figures 2 and 3) of George Marsaglia's article "Evaluating the
		   // Normal Distribution" (available at http://www.jstatsoft.org/v11/a05/paper) that the error of any cumulative normal
		   // function Φ(z) is dominated by the hardware (or compiler implementation) accuracy of exp(-z²/2) which is not
		   // reliably more than 14 digits when z is large. The accuracy of Φ(z) typically starts coming down to 14 digits when
		   // z is around -8. For the (normalised) Black function, as above in (*), this means that we are subtracting two terms
		   // that are each products of terms with about 14 digits of accuracy. The net result, in each of the products, is even
		   // less accuracy, and then we are taking the difference of these terms, resulting in even less accuracy. When we are
		   // using the asymptotic expansion asymptotic_expansion_of_normalized_black_call() invoked in the second branch at the
		   // beginning of this function, we are using only *one* exponential instead of 4, and this improves accuracy. It
		   // actually improves it a bit more than you would expect from the above logic, namely, almost the full two missing
		   // digits (in 64 bit IEEE floating point).  Unfortunately, going higher order in the asymptotic expansion will not
		   // enable us to gain more accuracy (by extending the range in which we could use the expansion) since the asymptotic
		   // expansion, being a divergent series, can never gain 16 digits of accuracy for z=-8 or just below. The best you can
		   // get is about 15 digits (just), for about 35 terms in the series (26.2.12), which would result in an prohibitively
		   // long expression in function asymptotic expansion asymptotic_expansion_of_normalized_black_call(). In this last branch,
		   // here, we therefore take a different tack as follows.
		   //     The "scaled complementary error function" is defined as erfcx(z) = exp(z²)·erfc(z). Cody's implementation of this
		   // function as published in "Rational Chebyshev approximations for the error function", W. J. Cody, Math. Comp., 1969, pp.
		   // 631-638, uses rational functions that theoretically approximates erfcx(x) to at least 18 significant decimal digits,
		   // *without* the use of the exponential function when x>4, which translates to about z<-5.66 in Φ(z). To make use of it,
		   // we write
		   //             Φ(z) = exp(-z²/2)·erfcx(-z/√2)/2
		   //
		   // to transform the normalised black function to
		   //
		   //   b   =  ½ · exp(-½(h²+t²)) · [ erfcx(-(h+t)/√2) -  erfcx(-(h-t)/√2) ]
		   //
		   // which now involves only one exponential, instead of three, when |h|+|t| > 5.66 , and the difference inside the
		   // square bracket is between the evaluation of two rational functions, which, typically, according to Marsaglia,
		   // retains the full 16 digits of accuracy (or just a little less than that).
		   //
		   final double b = 0.5 * Math.exp(-0.5*(h*h+t*t)) * ( erfcx_cody(-ONE_OVER_SQRT_TWO*(h+t)) - erfcx_cody(-ONE_OVER_SQRT_TWO*(h-t)) );
		   return Math.abs(Math.max(b,0.0));
		}

	
	protected double normalised_vega(double x, double s) {
		final double ax = Math.abs(x);
		return (ax <= 0) ? ONE_OVER_SQRT_TWO_PI * Math.exp(-0.125 * s * s)
				: ((s <= 0 || s <= ax * SQRT_DBL_MIN) ? 0
						: ONE_OVER_SQRT_TWO_PI * Math.exp(-0.5 * (((x / s) * (x / s)) + ((0.5 * s) * (0.5 * s)))));
	}

	protected double householder_factor(double newton, double halley, double hh3) {
		return (1 + 0.5 * halley * newton) / (1 + newton * (halley + hh3 * newton / 6));
	}

	protected double inverse_f_lower_map(final double x, final double f) {
		return is_below_horizon(f) ? 0
				: Math.abs(x / (SQRT_THREE
						* inverse_norm_cdf(Math.pow(f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * Math.abs(x)), 1. / 3.))));
	}
	 

	
	/* Originally names space calls */ 
	static final double TWO_PI =                        6.283185307179586476925286766559005768394338798750;
	static final double SQRT_PI_OVER_TWO     =         1.253314137315500251207882642405522626503493370305 ; // sqrt(pi/2) to avoid misinterpretation.
	static final double SQRT_THREE          =          1.732050807568877293527446341505872366942805253810;
	static final double  SQRT_ONE_OVER_THREE    =       0.577350269189625764509148780501957455647601751270;
	static final double  TWO_PI_OVER_SQRT_TWENTY_SEVEN = 1.209199576156145233729385505094770488189377498728; // 2*pi/sqrt(27)
	static final double  PI_OVER_SIX       =            0.523598775598298873077107230546583814032861566563;
	
	protected double f_lower_map(final double x, final double s) {
		if (is_below_horizon(x)){
			return 0;
		}
		if (is_below_horizon(s)){
			return 0;
		}
		final double z = SQRT_ONE_OVER_THREE * Math.abs(x) / s, Phi = norm_cdf(-z);
		return TWO_PI_OVER_SQRT_TWENTY_SEVEN * Math.abs(x) * (Phi * Phi * Phi);
	}

	protected double d_f_lower_map_d_beta(final double x, final double s) {
		if (is_below_horizon(s)){
			return 1;
		}
		final double z = SQRT_ONE_OVER_THREE * Math.abs(x) / s, y = z * z, Phi = norm_cdf(-z);
		return TWO_PI * y * (Phi * Phi) * Math.exp(y + 0.125 * s * s);
	}

	protected double d2_f_lower_map_d_beta2(final double x, final double s) {
		final double ax = Math.abs(x), z = SQRT_ONE_OVER_THREE * ax / s, y = z * z, s2 = s * s, Phi = norm_cdf(-z),
				phi = norm_pdf(z);
		return PI_OVER_SIX * y / (s2 * s) * Phi
				* (8 * SQRT_THREE * s * ax + (3 * s2 * (s2 - 8) - 8 * x * x) * Phi / phi) * Math.exp(2 * y + 0.25 * s2);
	}
		
	protected double f_upper_map(final double s) {
		return norm_cdf(-0.5 * s);
	}

	protected double d_f_upper_map_d_beta(final double x, final double s) {
		return is_below_horizon(x) ? -0.5 : -0.5 * Math.exp(0.5 * square(x / s));
	}

	protected double d2_f_upper_map_d_beta2(final double x, final double s) {
		if (is_below_horizon(x)){
			return 0;
		}
		final double w = square(x / s);
		return SQRT_PI_OVER_TWO * Math.exp(w + 0.125 * s * s) * w / s;
	}

	
	protected double square(double s) {
		return s * s;
	}
	
	protected double inverse_f_upper_map(double f) {
		return -2. * inverse_norm_cdf(f);
	}
	

	
	
	
	
	/* ERCF functions from erf_cody.cpp*/
	  
	   /* S    REAL FUNCTION ERFCX(X) */
	   /*<       DOUBLE PRECISION FUNCTION DERFCX(X) >*/
	 protected double erfcx_cody(double x) {
	      /* ------------------------------------------------------------------ */
	      /* This subprogram computes approximate values for exp(x*x) * erfc(x). */
	      /*   (see comments heading CALERF). */
	      /*   Author/date: W. J. Cody, March 30, 1987 */
	      /* ------------------------------------------------------------------ */
	      /*<       INTEGER JINT >*/
	      /* S    REAL             X, RESULT */
	      /*<       DOUBLE PRECISION X, RESULT >*/
	      /* ------------------------------------------------------------------ */
	      /*<       JINT = 2 >*/
	      /*<       CALL CALERF(X,RESULT,JINT) >*/
	      return calerf(x, 2);
	      /* S    ERFCX = RESULT */
	      /*<       DERFCX = RESULT >*/
	      /*<       RETURN >*/
	      /* ---------- Last card of DERFCX ---------- */
	      /*<       END >*/
	  } /* derfcx_ */
	
	   
	/* S    REAL FUNCTION ERFC(X) */
	/*<       DOUBLE PRECISION FUNCTION DERFC(X) >*/
	 protected double erfc_cody(double x) {
	   /* -------------------------------------------------------------------- */
	   /* This subprogram computes approximate values for erfc(x). */
	   /*   (see comments heading CALERF). */
	   /*   Author/date: W. J. Cody, January 8, 1985 */
	   /* -------------------------------------------------------------------- */
	   /*<       INTEGER JINT >*/
	   /* S    REAL             X, RESULT */
	   /*<       DOUBLE PRECISION X, RESULT >*/
	   /* ------------------------------------------------------------------ */
	   /*<       JINT = 1 >*/
	   /*<       CALL CALERF(X,RESULT,JINT) >*/
	   return calerf(x, 1);
	   /* S    ERFC = RESULT */
	   /*<       DERFC = RESULT >*/
	   /*<       RETURN >*/
	   /* ---------- Last card of DERFC ---------- */
	   /*<       END >*/
	} /* derfc_ */
	
	protected double d_int(final double x){ return( (x>0) ? Math.floor(x) : -Math.floor(-x) ); }
	
	 final double[] a =  { 3.1611237438705656,113.864154151050156,377.485237685302021,3209.37758913846947,.185777706184603153 };
	   final double[] b = { 23.6012909523441209,244.024637934444173,1282.61652607737228,2844.23683343917062 };
	   final double[] c__ =  { .564188496988670089,8.88314979438837594,66.1191906371416295,298.635138197400131,881.95222124176909,1712.04761263407058,2051.07837782607147,1230.33935479799725,2.15311535474403846e-8 };
	   final double[] d__ = { 15.7449261107098347,117.693950891312499,537.181101862009858,1621.38957456669019,3290.79923573345963,4362.61909014324716,3439.36767414372164,1230.33935480374942 };
	   final double[] p =  { .305326634961232344,.360344899949804439,.125781726111229246,.0160837851487422766,6.58749161529837803e-4,.0163153871373020978 };
	   final double[] q = { 2.56852019228982242,1.87295284992346047,.527905102951428412,.0605183413124413191,.00233520497626869185 };

	   final double zero = 0.;
	   final double half = .5;
	   final double one = 1.;
	   final double two = 2.;
	   final double four = 4.;
	   final double sqrpi = 0.56418958354775628695;
	   final double thresh = .46875;
	   final double sixten = 16.;

	   /* ------------------------------------------------------------------ */
	   /* This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x) */
	   /*   for a real argument  x.  It contains three FUNCTION type */
	   /*   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX), */
	   /*   and one SUBROUTINE type subprogram, CALERF.  The calling */
	   /*   statements for the primary entries are: */
	   /*                   Y=ERF(X)     (or   Y=DERF(X)), */
	   /*                   Y=ERFC(X)    (or   Y=DERFC(X)), */
	   /*   and */
	   /*                   Y=ERFCX(X)   (or   Y=DERFCX(X)). */
	   /*   The routine  CALERF  is intended for internal packet use only, */
	   /*   all computations within the packet being concentrated in this */
	   /*   routine.  The function subprograms invoke  CALERF  with the */
	   /*   statement */
	   /*          CALL CALERF(ARG,RESULT,JINT) */
	   /*   where the parameter usage is as follows */
	   /*      Function                     Parameters for CALERF */
	   /*       call              ARG                  Result          JINT */
	   /*     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0 */
	   /*     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1 */
	   /*     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2 */
	   /*   The main computation evaluates near-minimax approximations */
	   /*   from "Rational Chebyshev approximations for the error function" */
	   /*   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This */
	   /*   transportable program uses rational functions that theoretically */
	   /*   approximate  erf(x)  and  erfc(x)  to at least 18 significant */
	   /*   decimal digits.  The accuracy achieved depends on the arithmetic */
	   /*   system, the compiler, the intrinsic functions, and proper */
	   /*   selection of the machine-dependent constants. */
	   /* ******************************************************************* */
	   /* ******************************************************************* */
	   /* Explanation of machine-dependent constants */
	   /*   XMIN   = the smallest positive floating-point number. */
	   /*   XINF   = the largest positive finite floating-point number. */
	   /*   XNEG   = the largest negative argument acceptable to ERFCX; */
	   /*            the negative of the solution to the equation */
	   /*            2*exp(x*x) = XINF. */
	   /*   XSMALL = argument below which erf(x) may be represented by */
	   /*            2*x/sqrt(pi)  and above which  x*x  will not underflow. */
	   /*            A conservative value is the largest machine number X */
	   /*            such that   1.0 + X = 1.0   to machine precision. */
	   /*   XBIG   = largest argument acceptable to ERFC;  solution to */
	   /*            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where */
	   /*            W(x) = exp(-x*x)/[x*sqrt(pi)]. */
	   /*   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to */
	   /*            machine precision.  A conservative value is */
	   /*            1/[2*sqrt(XSMALL)] */
	   /*   XMAX   = largest acceptable argument to ERFCX; the minimum */
	   /*            of XINF and 1/[sqrt(pi)*XMIN]. */
	   // The numbers below were preselected for IEEE .
	 final double xinf = 1.79e308;
	 final double xneg = -26.628;
	 final double xsmall = 1.11e-16;
	 final double xbig = 26.543;
	 final double xhuge = 6.71e7;
	 final double xmax = 2.53e307;
	 
	/*<       SUBROUTINE CALERF(ARG,RESULT,JINT) >*/
	protected double calerf(double x, final int jint) {

	  
	   double y=0.0, del=0.0, ysq=0.0, xden=0.0, xnum=0, result=0.0;

	 
	   /*   Approximate values for some important machines are: */
	   /*                          XMIN       XINF        XNEG     XSMALL */
	   /*  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15 */
	   /*  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15 */
	   /*  IEEE (IBM/XT, */
	   /*    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8 */
	   /*  IEEE (IBM/XT, */
	   /*    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16 */
	   /*  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17 */
	   /*  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18 */
	   /*  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17 */
	   /*  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16 */
	   /*                          XBIG       XHUGE       XMAX */
	   /*  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293 */
	   /*  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465 */
	   /*  IEEE (IBM/XT, */
	   /*    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37 */
	   /*  IEEE (IBM/XT, */
	   /*    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307 */
	   /*  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75 */
	   /*  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307 */
	   /*  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38 */
	   /*  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307 */
	   /* ******************************************************************* */
	   /* ******************************************************************* */
	   /* Error returns */
	   /*  The program returns  ERFC = 0      for  ARG .GE. XBIG; */
	   /*                       ERFCX = XINF  for  ARG .LT. XNEG; */
	   /*      and */
	   /*                       ERFCX = 0     for  ARG .GE. XMAX. */
	   /* Intrinsic functions required are: */
	   /*     ABS, AINT, EXP */
	   /*  Author: W. J. Cody */
	   /*          Mathematics and Computer Science Division */
	   /*          Argonne National Laboratory */
	   /*          Argonne, IL 60439 */
	   /*  Latest modification: March 19, 1990 */
	   /* ------------------------------------------------------------------ */
	   /*<       INTEGER I,JINT >*/
	   /* S    REAL */
	   /*<    >*/
	   /*<       DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5) >*/
	   /* ------------------------------------------------------------------ */
	   /*  Mathematical constants */
	   /* ------------------------------------------------------------------ */
	   /* S    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/, */
	   /* S   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/, */
	   /* S   2     SIXTEN/16.0E0/ */
	   /*<    >*/
	   /* ------------------------------------------------------------------ */
	   /*  Machine-dependent constants */
	   /* ------------------------------------------------------------------ */
	   /* S    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/, */
	   /* S   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/ */
	   /*<    >*/
	   /* ------------------------------------------------------------------ */
	   /*  Coefficients for approximation to  erf  in first interval */
	   /* ------------------------------------------------------------------ */
	   /* S    DATA A/3.16112374387056560E00,1.13864154151050156E02, */
	   /* S   1       3.77485237685302021E02,3.20937758913846947E03, */
	   /* S   2       1.85777706184603153E-1/ */
	   /* S    DATA B/2.36012909523441209E01,2.44024637934444173E02, */
	   /* S   1       1.28261652607737228E03,2.84423683343917062E03/ */
	   /*<    >*/
	   /*<    >*/
	   /* ------------------------------------------------------------------ */
	   /*  Coefficients for approximation to  erfc  in second interval */
	   /* ------------------------------------------------------------------ */
	   /* S    DATA C/5.64188496988670089E-1,8.88314979438837594E0, */
	   /* S   1       6.61191906371416295E01,2.98635138197400131E02, */
	   /* S   2       8.81952221241769090E02,1.71204761263407058E03, */
	   /* S   3       2.05107837782607147E03,1.23033935479799725E03, */
	   /* S   4       2.15311535474403846E-8/ */
	   /* S    DATA D/1.57449261107098347E01,1.17693950891312499E02, */
	   /* S   1       5.37181101862009858E02,1.62138957456669019E03, */
	   /* S   2       3.29079923573345963E03,4.36261909014324716E03, */
	   /* S   3       3.43936767414372164E03,1.23033935480374942E03/ */
	   /*<    >*/
	   /*<    >*/
	   /* ------------------------------------------------------------------ */
	   /*  Coefficients for approximation to  erfc  in third interval */
	   /* ------------------------------------------------------------------ */
	   /* S    DATA P/3.05326634961232344E-1,3.60344899949804439E-1, */
	   /* S   1       1.25781726111229246E-1,1.60837851487422766E-2, */
	   /* S   2       6.58749161529837803E-4,1.63153871373020978E-2/ */
	   /* S    DATA Q/2.56852019228982242E00,1.87295284992346047E00, */
	   /* S   1       5.27905102951428412E-1,6.05183413124413191E-2, */
	   /* S   2       2.33520497626869185E-3/ */
	   /*<    >*/
	   /*<    >*/
	   /* ------------------------------------------------------------------ */
	   /*<       X = ARG >*/
	   // x = *arg;
	   /*<       Y = ABS(X) >*/
	   y = Math.abs(x);
	   /*<       IF (Y .LE. THRESH) THEN >*/
	   if (y <= thresh) {
	      /* ------------------------------------------------------------------ */
	      /*  Evaluate  erf  for  |X| <= 0.46875 */
	      /* ------------------------------------------------------------------ */
	      /*<             YSQ = ZERO >*/
	      ysq = zero;
	      /*<             IF (Y .GT. XSMALL) YSQ = Y * Y >*/
	      if (y > xsmall) {
	         ysq = y * y;
	      }
	      /*<             XNUM = A(5)*YSQ >*/
	      xnum = a[4] * ysq;
	      /*<             XDEN = YSQ >*/
	      xden = ysq;
	      /*<             DO 20 I = 1, 3 >*/
	      for (int i__ = 1; i__ <= 3; ++i__) {
	         /*<                XNUM = (XNUM + A(I)) * YSQ >*/
	         xnum = (xnum + a[i__ - 1]) * ysq;
	         /*<                XDEN = (XDEN + B(I)) * YSQ >*/
	         xden = (xden + b[i__ - 1]) * ysq;
	         /*<    20       CONTINUE >*/
	         /* L20: */
	      }
	      /*<             RESULT = X * (XNUM + A(4)) / (XDEN + B(4)) >*/
	      result = x * (xnum + a[3]) / (xden + b[3]);
	      /*<             IF (JINT .NE. 0) RESULT = ONE - RESULT >*/
	      if (jint != 0) {
	         result = one - result;
	      }
	      /*<             IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT >*/
	      if (jint == 2) {
	         result = Math.exp(ysq) * result;
	      }
	      /*<             GO TO 800 >*/
	      return result;
	      /* ------------------------------------------------------------------ */
	      /*  Evaluate  erfc  for 0.46875 <= |X| <= 4.0 */
	      /* ------------------------------------------------------------------ */
	      /*<          ELSE IF (Y .LE. FOUR) THEN >*/
	   } else if (y <= four) {
	      /*<             XNUM = C(9)*Y >*/
	      xnum = c__[8] * y;
	      /*<             XDEN = Y >*/
	      xden = y;
	      /*<             DO 120 I = 1, 7 >*/
	      for (int i__ = 1; i__ <= 7; ++i__) {
	         /*<                XNUM = (XNUM + C(I)) * Y >*/
	         xnum = (xnum + c__[i__ - 1]) * y;
	         /*<                XDEN = (XDEN + D(I)) * Y >*/
	         xden = (xden + d__[i__ - 1]) * y;
	         /*<   120       CONTINUE >*/
	         /* L120: */
	      }
	      /*<             RESULT = (XNUM + C(8)) / (XDEN + D(8)) >*/
	      result = (xnum + c__[7]) / (xden + d__[7]);
	      /*<             IF (JINT .NE. 2) THEN >*/
	      if (jint != 2) {
	         /*<                YSQ = AINT(Y*SIXTEN)/SIXTEN >*/
	         double d__1 = y * sixten;
	         ysq = d_int(d__1) / sixten;
	         /*<                DEL = (Y-YSQ)*(Y+YSQ) >*/
	         del = (y - ysq) * (y + ysq);
	         /*<                RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT >*/
	         d__1 = Math.exp(-ysq * ysq) * Math.exp(-del);
	         result = d__1 * result;
	         /*<             END IF >*/
	      }
	      /* ------------------------------------------------------------------ */
	      /*  Evaluate  erfc  for |X| > 4.0 */
	      /* ------------------------------------------------------------------ */
	      /*<          ELSE >*/
	   } else {
	      /*<             RESULT = ZERO >*/
	      result = zero;
	      /*<             IF (Y .GE. XBIG) THEN >*/
	      if (y >= xbig) {
	         /*<                IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300 >*/
	         if (jint != 2 || y >= xmax) {
	        	 // REPLACED GOTO with function call
	        	 return negativeArgFix(result, jint, x, ysq, del, y);
	         }
	         /*<                IF (Y .GE. XHUGE) THEN >*/
	         if (y >= xhuge) {
	            /*<                   RESULT = SQRPI / Y >*/
	            result = sqrpi / y;
	            /*<                   GO TO 300 >*/
	            // REPLACED GOTO with function call
	            return negativeArgFix(result, jint, x, ysq, del, y);
	            /*<                END IF >*/
	         }
	         /*<             END IF >*/
	      }
	      /*<             YSQ = ONE / (Y * Y) >*/
	      ysq = one / (y * y);
	      /*<             XNUM = P(6)*YSQ >*/
	      xnum = p[5] * ysq;
	      /*<             XDEN = YSQ >*/
	      xden = ysq;
	      /*<             DO 240 I = 1, 4 >*/
	      for (int i__ = 1; i__ <= 4; ++i__) {
	         /*<                XNUM = (XNUM + P(I)) * YSQ >*/
	         xnum = (xnum + p[i__ - 1]) * ysq;
	         /*<                XDEN = (XDEN + Q(I)) * YSQ >*/
	         xden = (xden + q[i__ - 1]) * ysq;
	         /*<   240       CONTINUE >*/
	         /* L240: */
	      }
	      /*<             RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5)) >*/
	      result = ysq * (xnum + p[4]) / (xden + q[4]);
	      /*<             RESULT = (SQRPI -  RESULT) / Y >*/
	      result = (sqrpi - result) / y;
	      /*<             IF (JINT .NE. 2) THEN >*/
	      if (jint != 2) {
	         /*<                YSQ = AINT(Y*SIXTEN)/SIXTEN >*/
	         double d__1 = y * sixten;
	         ysq = d_int(d__1) / sixten;
	         /*<                DEL = (Y-YSQ)*(Y+YSQ) >*/
	         del = (y - ysq) * (y + ysq);
	         /*<                RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT >*/
	         d__1 = Math.exp(-ysq * ysq) * Math.exp(-del);
	         result = d__1 * result;
	         /*<             END IF >*/
	      }
	      /*<       END IF >*/
	   }
	   /* ------------------------------------------------------------------ */
	   /*  Fix up for negative argument, erf, etc. */
	   /* ------------------------------------------------------------------ */
	   /*<   300 IF (JINT .EQ. 0) THEN >*/
	   // REPLACED GOTO with function call
	return negativeArgFix(result, jint, x, ysq, del, y);
	   /* ---------- Last card of CALERF ---------- */
	   /*<       END >*/
	} /* calerf_ */
	
	
	protected double negativeArgFix(double result, final int jint,  double x, double ysq, double del, double y){
		  if (jint == 0) {
		      /*<             RESULT = (HALF - RESULT) + HALF >*/
		      result = (half - result) + half;
		      /*<             IF (X .LT. ZERO) RESULT = -RESULT >*/
		      if (x < zero) {
		         result = -(result);
		      }
		      /*<          ELSE IF (JINT .EQ. 1) THEN >*/
		   } else if (jint == 1) {
		      /*<             IF (X .LT. ZERO) RESULT = TWO - RESULT >*/
		      if (x < zero) {
		         result = two - result;
		      }
		      /*<          ELSE >*/
		   } else {
		      /*<             IF (X .LT. ZERO) THEN >*/
		      if (x < zero) {
		         /*<                IF (X .LT. XNEG) THEN >*/
		         if (x < xneg) {
		            /*<                      RESULT = XINF >*/
		            result = xinf;
		            /*<                   ELSE >*/
		         } else {
		            /*<                      YSQ = AINT(X*SIXTEN)/SIXTEN >*/
		            double d__1 = x * sixten;
		            ysq = d_int(d__1) / sixten;
		            /*<                      DEL = (X-YSQ)*(X+YSQ) >*/
		            del = (x - ysq) * (x + ysq);
		            /*<                      Y = EXP(YSQ*YSQ) * EXP(DEL) >*/
		            y = Math.exp(ysq * ysq) * Math.exp(del);
		            /*<                      RESULT = (Y+Y) - RESULT >*/
		            result = y + y - result;
		            /*<                END IF >*/
		         }
		         /*<             END IF >*/
		      }
		      /*<       END IF >*/
		   }
		  return result;
	}
	
	
	
	
	/* Normalised distribution from normaliseddistribution.cpp */
	
	protected double norm_cdf(double z) {
		if (z <= norm_cdf_asymptotic_expansion_first_threshold) {
			// Asymptotic expansion for very negative z following (26.2.12) on
			// page 408
			// in M. Abramowitz and A. Stegun, Pocketbook of Mathematical
			// Functions, ISBN 3-87144818-4.
			double sum = 1;
			if (z >= norm_cdf_asymptotic_expansion_second_threshold) {
				double zsqr = z * z, i = 1, g = 1, x, y, a = DBL_MAX, lasta;
				do {
					lasta = a;
					x = (4 * i - 3) / zsqr;
					y = x * ((4 * i - 1) / zsqr);
					a = g * (x - y);
					sum -= a;
					g *= y;
					++i;
					a = Math.abs(a);
				} while (lasta > a && a >= Math.abs(sum * DBL_EPSILON));
			}
			return -norm_pdf(z) * sum / z;
		}
		return 0.5 * erfc_cody(-z * ONE_OVER_SQRT_TWO);
	}
	
	
	protected double norm_pdf(double x) {
		return ONE_OVER_SQRT_TWO_PI * Math.exp(-.5 * x * x);
	}
	
	protected double inverse_norm_cdf(double u){
		   //
		   // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
		   //
		   // Produces the normal deviate Z corresponding to a given lower
		   // tail area of u; Z is accurate to about 1 part in 10**16.
		   // see http://lib.stat.cmu.edu/apstat/241
		   //
		   final double split1 = 0.425;
		   final double split2 = 5.0;
		   final double const1 = 0.180625;
		   final double const2 = 1.6;

		   // Coefficients for P close to 0.5
		   final double A0 = 3.3871328727963666080E0;
		   final double A1 = 1.3314166789178437745E+2;
		   final double A2 = 1.9715909503065514427E+3;
		   final double A3 = 1.3731693765509461125E+4;
		   final double A4 = 4.5921953931549871457E+4;
		   final double A5 = 6.7265770927008700853E+4;
		   final double A6 = 3.3430575583588128105E+4;
		   final double A7 = 2.5090809287301226727E+3;
		   final double B1 = 4.2313330701600911252E+1;
		   final double B2 = 6.8718700749205790830E+2;
		   final double B3 = 5.3941960214247511077E+3;
		   final double B4 = 2.1213794301586595867E+4;
		   final double B5 = 3.9307895800092710610E+4;
		   final double B6 = 2.8729085735721942674E+4;
		   final double B7 = 5.2264952788528545610E+3;
		   // Coefficients for P not close to 0, 0.5 or 1.
		   final double C0 = 1.42343711074968357734E0;
		   final double C1 = 4.63033784615654529590E0;
		   final double C2 = 5.76949722146069140550E0;
		   final double C3 = 3.64784832476320460504E0;
		   final double C4 = 1.27045825245236838258E0;
		   final double C5 = 2.41780725177450611770E-1;
		   final double C6 = 2.27238449892691845833E-2;
		   final double C7 = 7.74545014278341407640E-4;
		   final double D1 = 2.05319162663775882187E0;
		   final double D2 = 1.67638483018380384940E0;
		   final double D3 = 6.89767334985100004550E-1;
		   final double D4 = 1.48103976427480074590E-1;
		   final double D5 = 1.51986665636164571966E-2;
		   final double D6 = 5.47593808499534494600E-4;
		   final double D7 = 1.05075007164441684324E-9;
		   // Coefficients for P very close to 0 or 1
		   final double E0 = 6.65790464350110377720E0;
		   final double E1 = 5.46378491116411436990E0;
		   final double E2 = 1.78482653991729133580E0;
		   final double E3 = 2.96560571828504891230E-1;
		   final double E4 = 2.65321895265761230930E-2;
		   final double E5 = 1.24266094738807843860E-3;
		   final double E6 = 2.71155556874348757815E-5;
		   final double E7 = 2.01033439929228813265E-7;
		   final double F1 = 5.99832206555887937690E-1;
		   final double F2 = 1.36929880922735805310E-1;
		   final double F3 = 1.48753612908506148525E-2;
		   final double F4 = 7.86869131145613259100E-4;
		   final double F5 = 1.84631831751005468180E-5;
		   final double F6 = 1.42151175831644588870E-7;
		   final double F7 = 2.04426310338993978564E-15;

		   if (u<=0)
		      return Math.log(u);
		   if (u>=1)
		      return Math.log(1-u);

		   final double q = u-0.5;
		   if (Math.abs(q) <= split1)
		   {
		      final double r = const1 - q*q;
		      return q * (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0) /
		         (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + 1.0);
		   }
		   else
		   {
		      double r = q<0.0 ? u : 1.0-u;
		      r = Math.sqrt(-Math.log(r));
		      double ret;
		      if (r < split2)
		      {
		         r = r - const2;
		         ret = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0) /
		            (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + 1.0);
		      }
		      else
		      {
		         r = r - split2;
		         ret = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0) /
		            (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + 1.0);
		      }
		      return q<0.0 ? -ret : ret;
		   }
	}
	
	/* End of NORMALISED FUNCTIONS*/
	
	
	
	
	
	
/* Rational Cubic methods from rationalcubic.cpp*/
	
	final static double minimum_rational_cubic_control_parameter_value = -(1 - Math.sqrt(DBL_EPSILON));
	final static double maximum_rational_cubic_control_parameter_value = 2 / (DBL_EPSILON * DBL_EPSILON);

	protected boolean is_zero(double x) {
		return Math.abs(x) < DBL_MIN;
	}

	protected double rational_cubic_interpolation(double x, double x_l, double x_r, double y_l, double y_r, double d_l, double d_r, double r) {
	   final double h = (x_r - x_l);
	   if (Math.abs(h)<=0)
	      return 0.5 * (y_l + y_r);
	   // r should be greater than -1. We do not use  assert(r > -1)  here in order to allow values such as NaN to be propagated as they should.
	   double t = (x - x_l) / h;
	   if ( ! (r >= maximum_rational_cubic_control_parameter_value) ) {
	       t = (x - x_l) / h;
	       double omt = 1 - t, t2 = t * t, omt2 = omt * omt;
	      // Formula (2.4) divided by formula (2.5)
	      return (y_r * t2 * t + (r * y_r - h * d_r) * t2 * omt + (r * y_l + h * d_l) * t * omt2 + y_l * omt2 * omt) / (1 + (r - 3) * t * omt);
	   }
	   // Linear interpolation without over-or underflow.
	   return y_r * t + y_l * (1 - t);
	}

	protected double rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(double x_l, double x_r, double y_l, double y_r, double d_l, double d_r, double second_derivative_l) {
	   final double h = (x_r-x_l), numerator = 0.5*h*second_derivative_l+(d_r-d_l);
	   if (is_zero(numerator))
	      return 0;
	   final double denominator = (y_r-y_l)/h-d_l;
	   if (is_zero(denominator))
	      return numerator>0 ? maximum_rational_cubic_control_parameter_value : minimum_rational_cubic_control_parameter_value;
	   return numerator/denominator;
	}

	protected double rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(double x_l, double x_r, double y_l, double y_r, double d_l, double d_r, double second_derivative_r) {
	   final double h = (x_r-x_l), numerator = 0.5*h*second_derivative_r+(d_r-d_l);
	   if (is_zero(numerator))
	      return 0;
	   final double denominator = d_r-(y_r-y_l)/h;
	   if (is_zero(denominator))
	      return numerator>0 ? maximum_rational_cubic_control_parameter_value : minimum_rational_cubic_control_parameter_value;
	   return numerator/denominator;
	}

	protected double minimum_rational_cubic_control_parameter(double d_l, double d_r, double s, boolean preferShapePreservationOverSmoothness) {
	   final boolean monotonic = d_l * s >= 0 && d_r * s >= 0, convex = d_l <= s && s <= d_r, concave = d_l >= s && s >= d_r;
	   if (!monotonic && !convex && !concave) // If 3==r_non_shape_preserving_target, this means revert to standard cubic.
	      return minimum_rational_cubic_control_parameter_value;
	   final double d_r_m_d_l = d_r - d_l, d_r_m_s = d_r - s, s_m_d_l = s - d_l;
	   double r1 = -DBL_MAX, r2 = r1;
	   // If monotonicity on this interval is possible, set r1 to satisfy the monotonicity condition (3.8).
	   if (monotonic){
	      if (!is_zero(s)) // (3.8), avoiding division by zero.
	         r1 = (d_r + d_l) / s; // (3.8)
	      else if (preferShapePreservationOverSmoothness) // If division by zero would occur, and shape preservation is preferred, set value to enforce linear interpolation.
	         r1 =  maximum_rational_cubic_control_parameter_value;  // This value enforces linear interpolation.
	   }
	   if (convex || concave) {
	      if (!(is_zero(s_m_d_l) || is_zero(d_r_m_s))) // (3.18), avoiding division by zero.
	         r2 = Math.max(Math.abs(d_r_m_d_l / d_r_m_s), Math.abs(d_r_m_d_l / s_m_d_l));
	      else if (preferShapePreservationOverSmoothness)
	         r2 = maximum_rational_cubic_control_parameter_value; // This value enforces linear interpolation.
	   } else if (monotonic && preferShapePreservationOverSmoothness)
	      r2 = maximum_rational_cubic_control_parameter_value; // This enforces linear interpolation along segments that are inconsistent with the slopes on the boundaries, e.g., a perfectly horizontal segment that has negative slopes on either edge.
	   return Math.max(minimum_rational_cubic_control_parameter_value, Math.max(r1, r2));
	}

	protected double convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(double x_l, double x_r, double y_l, double y_r, double d_l, double d_r, double second_derivative_l, boolean preferShapePreservationOverSmoothness) {
	   final double r = rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l);
	   final double r_min = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r-y_l)/(x_r-x_l), preferShapePreservationOverSmoothness);
	   return Math.max(r,r_min);
	}

	protected double convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(double x_l, double x_r, double y_l, double y_r, double d_l, double d_r, double second_derivative_r, boolean preferShapePreservationOverSmoothness) {
	   final double r = rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r);
	   final double r_min = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r-y_l)/(x_r-x_l), preferShapePreservationOverSmoothness);
	   return Math.max(r,r_min);
	}
	
	/* END of Rational cubic */
}
