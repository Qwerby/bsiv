
# BlackScholesJaekelIVCalculator

A fast and efficient calculator for Black Schole's implied volatility based on reduced iterations.
 
The original paper from Peter Jaeckel: <http://jaeckel.16mb.com/LetsBeRational.pdf>
 
The Java code is a direct port of the c++ code <http://jaeckel.16mb.com/LetsBeRational.7z>
 
Original reference from:
 
<https://github.com/vollib/vollib/blob/master/vollib/black_scholes_merton/implied_volatility.py>


See the main file for the original copyright reference.


## Usage

Requires:

* optionPrice - price of the option from the market
* underlyingPrice - price of the underlying fromt he market
* strike - strike price of the option
* days - number of days until expiry
* riskFreeRate - interest rate (expressed where 1 == 100%). i.e 5% is 0.05.
* dividendYield - from the market (expressed where 1 == 100%). i.e 3% is 0.03.
* call - whether it is a call or a put


