package rubris.io;

import org.junit.Assert;
import org.junit.Test;



public class BlackScholesJaekelIVCalculatorTest {

	// validity check against 
	// http://www.fintools.com/resources/online-calculators/options-calcs/binomial/
	// http://www.option-price.com/implied-volatility.php
	
	@Test
	public void testSimpleCall() {
		
		double underlying =100;
		double rpice =1;
		int days = 30;
		double rate =0.05;
		double yield =.01;
		
		BlackScholesJaekelIVCalculator provider = new BlackScholesJaekelIVCalculator();
		double iv =provider.calculate(rpice, underlying, 100, days, rate, yield, true);
		Assert.assertEquals(0.0723,iv,0.0001);
	}
	
	@Test
	public void testHPricedCall() {
		
		double underlying =100;
		double rpice =10;
		int days = 30;
		double rate =0.05;
		double yield =.01;
		
		BlackScholesJaekelIVCalculator provider = new BlackScholesJaekelIVCalculator();
		double iv =provider.calculate(rpice, underlying, 100, days, rate, yield, true);
		Assert.assertEquals(0.8642,iv,0.0001);
	}
	
	@Test
	public void testSimplePut() {
		
		double underlying =107.48;
		double rpice =1.47;
		int days = 21;
		double rate =0.0023;
		double yield =.0215;
		
		BlackScholesJaekelIVCalculator provider = new BlackScholesJaekelIVCalculator();
		double iv =provider.calculate(rpice, underlying, 107, days, rate, yield, false);
		Assert.assertEquals(0.1604,iv,0.0001);
	}

	
	@Test
	public void testStrikesForPut() {
		
		double[] expected = {.2556,0.2269,.1973,.1662,.1329,.0953,.0436};
		double[] strikes={99.0,100.0,101.0,102.0,103.0,104.0,105.00};
		double underlying =104;
		double rprice =1.03;
		int days = 29;
		double rate =0.05;
		double yield =.03;
		
		BlackScholesJaekelIVCalculator provider = new BlackScholesJaekelIVCalculator();
		for(int i=0;i<strikes.length;i++){
			double iv =provider.calculate(rprice, underlying, strikes[i], days, rate, yield, false);
			Assert.assertEquals(expected[i],iv,0.0001);
		}
	}

	@Test
	public void testStrikesForCall() {
		
		double[] expected = {.2584,0.3297,.3873,.4375,.4831,.5252,.5647};
		double[] strikes={99.0,100.0,101.0,102.0,103.0,104.0,105.00};
		double underlying =104;
		double rprice =6.20;
		int days = 29;
		double rate =0.05;
		double yield =.03;
		
		BlackScholesJaekelIVCalculator provider = new BlackScholesJaekelIVCalculator();
		for(int i=0;i<strikes.length;i++){
			double iv =provider.calculate(rprice, underlying, strikes[i], days, rate, yield, true);
			Assert.assertEquals(expected[i],iv,0.0001);
		}
	}
	
	@Test
	public void testPerformance(){
		double underlying =104;
		double rprice =6.20;
		int days = 29;
		double rate =0.05;
		double yield =.03;
		int iterations =100000;
		double expected =.2584;
		double strike=99.0;
		
		BlackScholesJaekelIVCalculator provider = new BlackScholesJaekelIVCalculator();
		
		long time = System.currentTimeMillis();
		for(int i=0;i<iterations;i++){
			double iv =provider.calculate(rprice, underlying, strike, days, rate, yield, true);
			Assert.assertEquals(expected,iv,0.0001);
		}
		Assert.assertTrue((System.currentTimeMillis() -time) <500);
	}
	
}
