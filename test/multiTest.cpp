#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include <iomanip>
#include "../src/blackscholesconstprocess.hpp"
#include "../src/mceuropeanconstengine.hpp"

using namespace QuantLib;

void flatCurveSimulation(const Date & maturity, boost::shared_ptr<BlackScholesMertonProcess> process,boost::shared_ptr<StrikedTypePayoff> payoff){
        std::cout << "\nMaturity = "<< maturity << std::endl;
        // european exercise
        boost::shared_ptr<Exercise> europeanExercise(
                new EuropeanExercise(maturity));

        // options
        VanillaOption europeanOption(payoff, europeanExercise);
        // Black-Scholes for European plat
        
        europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
                    new AnalyticEuropeanEngine(process)));


        clock_t t1,t2;  
        Real res;
        
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "Black-Scholes : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;     
}

void flatCurveMontCarloSimulation(const Date & maturity, boost::shared_ptr<BlackScholesMertonProcess> bsmProcess,boost::shared_ptr<StrikedTypePayoff> payoff){
        // european exercise
        boost::shared_ptr<Exercise> europeanExercise(
        new EuropeanExercise(maturity));

        // options
        VanillaOption europeanOption(payoff, europeanExercise);
        // Black-Scholes for European plat
            // Monte Carlo Method: MC (crude)
	        
        Size timeSteps = 1;
        Size mcSeed = 42;
	
        boost::shared_ptr<PricingEngine> mcengine1;
        mcengine1 = MakeMCEuropeanConstEngine<PseudoRandom>(bsmProcess, false)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);
        europeanOption.setPricingEngine(mcengine1);
        
        clock_t t1,t2;  
        Real res;
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC (crude) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;
        
        
        boost::shared_ptr<PricingEngine> mcengine1c;
        mcengine1c = MakeMCEuropeanConstEngine<PseudoRandom>(bsmProcess, true)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);
        europeanOption.setPricingEngine(mcengine1c);
        
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC const(crude) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;
        
	
        
        // Monte Carlo Method: QMC (Sobol)
        Size nSamples = 32768;  // 2^15
	
        boost::shared_ptr<PricingEngine> mcengine2;
        mcengine2 = MakeMCEuropeanConstEngine<LowDiscrepancy>(bsmProcess, false)
            .withSteps(timeSteps)
            .withSamples(nSamples);
                 
        europeanOption.setPricingEngine(mcengine2);
        
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC (Sobol) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;
        
        boost::shared_ptr<PricingEngine> mcengine2c;
        mcengine2c = MakeMCEuropeanConstEngine<LowDiscrepancy>(bsmProcess, true)
            .withSteps(timeSteps)
            .withSamples(nSamples);
                 
        europeanOption.setPricingEngine(mcengine2c);
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC const(Sobol) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;

}

void forwardCurveSimulation(const Date &  maturity,
                              const boost::shared_ptr<StrikedTypePayoff> payoff,
                              boost::shared_ptr<BlackScholesMertonProcess> bsmProcess){
        std::cout << "\nMaturity = "<< maturity << std::endl;
                
                
        boost::shared_ptr<Exercise> europeanExercise(
        new EuropeanExercise(maturity));      
        VanillaOption europeanOption(payoff, europeanExercise);
        // Black-Scholes for European forward curve
        europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
                    new AnalyticEuropeanEngine(bsmProcess)));

        clock_t t1,t2;
        Real res;
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "Black-Scholes(forward curve) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;
                // Black-Scholes for European plat
            // Monte Carlo Method: MC (crude)
	        
        Size timeSteps = 1;
        Size mcSeed = 42;
	
        boost::shared_ptr<PricingEngine> mcengine1;
        mcengine1 = MakeMCEuropeanConstEngine<PseudoRandom>(bsmProcess, false)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);
        europeanOption.setPricingEngine(mcengine1);
          
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC (crude) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;
        
        
        boost::shared_ptr<PricingEngine> mcengine1c;
        mcengine1c = MakeMCEuropeanConstEngine<PseudoRandom>(bsmProcess, true)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);
        europeanOption.setPricingEngine(mcengine1c);
        
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC const(crude) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;
        
	
        
        // Monte Carlo Method: QMC (Sobol)
        Size nSamples = 32768;  // 2^15
	
        boost::shared_ptr<PricingEngine> mcengine2;
        mcengine2 = MakeMCEuropeanConstEngine<LowDiscrepancy>(bsmProcess, false)
            .withSteps(timeSteps)
            .withSamples(nSamples);
                 
        europeanOption.setPricingEngine(mcengine2);
        
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC (Sobol) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;
        
        boost::shared_ptr<PricingEngine> mcengine2c;
        mcengine2c = MakeMCEuropeanConstEngine<LowDiscrepancy>(bsmProcess, true)
            .withSteps(timeSteps)
            .withSamples(nSamples);
                 
        europeanOption.setPricingEngine(mcengine2c);
        t1 = clock();
        res = europeanOption.NPV();     
        t2 = clock();
        std::cout << "MC const(Sobol) : " << res << " (" << (float)(t2-t1)/(double(CLOCKS_PER_SEC)*1000) << "ms)"<<std::endl;

}





int main(int argc, char* argv[]){
    
    try{
        
        boost::timer timer;
        std::cout << std::endl;

        // set up dates
        Calendar calendar = TARGET();
        Date todaysDate(15, May, 1998);
        Date settlementDate(17, May, 1998);
        Settings::instance().evaluationDate() = todaysDate;

        // our option parameters
        Option::Type type(Option::Put);
        Real underlying = 36;
        Real strike = 40;
        Spread dividendYield = 0.00;
        Rate riskFreeRate = 0.06;
        Volatility volatility = 0.20;


        DayCounter dayCounter = Actual365Fixed();

        std::cout << "Option type = "  << type << std::endl;
        std::cout << "Underlying price = "        << underlying << std::endl;
        std::cout << "Strike = "                  << strike << std::endl;
        std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)
                  << std::endl;
        std::cout << "Dividend yield = " << io::rate(dividendYield)
                  << std::endl;
        std::cout << "Volatility = " << io::volatility(volatility)
                  << std::endl;
        std::cout << std::endl;
        std::string method;
        std::cout << std::endl ;

        // underlying handler
        Handle<Quote> underlyingH(
                boost::shared_ptr<Quote>(new SimpleQuote(underlying)));

        // bootstrap the yield/dividend/vol curves
        Handle<YieldTermStructure> flatTermStructure(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
        Handle<YieldTermStructure> flatDividendTS(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter)));
        Handle<BlackVolTermStructure> flatVolTS(
            boost::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(settlementDate, calendar, volatility,
                                     dayCounter)));
        // payoff
        boost::shared_ptr<StrikedTypePayoff> payoff(
                new PlainVanillaPayoff(type, strike));                
                                                                                                              
        // BlackScholes Merton Process platForward
        boost::shared_ptr<BlackScholesMertonProcess> flatbsmProcess(
                new BlackScholesMertonProcess(underlyingH, flatDividendTS, flatTermStructure, flatVolTS));
                          
        std::cout<< " \n\nflat curve tests :"<<std::endl;
        Date maturity(17, May, 2001);
        flatCurveSimulation(maturity, flatbsmProcess,payoff);
        flatCurveMontCarloSimulation(maturity, flatbsmProcess, payoff);
        
        maturity=Date(17, May, 2003);
        flatCurveSimulation(maturity, flatbsmProcess,payoff);
        flatCurveMontCarloSimulation(maturity, flatbsmProcess, payoff);
        
        maturity=Date(17, May, 2008);
        flatCurveSimulation(maturity, flatbsmProcess,payoff);
        flatCurveMontCarloSimulation(maturity, flatbsmProcess, payoff);
        
        //
        //
        //
        //fowardCurve part
        std::vector<Date> dates1(7);
        std::vector<Rate> rates(7);
        std::vector<Volatility> vols(6);
        std::vector<Date> dates2(6);

        dates1[0] = Date(17, May, 1998);    
        dates1[1] = Date(17, May, 1999);  
        dates1[2] = Date(17, May, 2001); 
        dates1[3] = Date(17, May, 2003);    
        dates1[4] = Date(17, May, 2005);  
        dates1[5] = Date(17, May, 2007);
        dates1[6] = Date(17, May, 2009);
        
        rates[0] = 0.06;
        rates[1] = 0.05;
        rates[2] = 0.04;
        rates[3] = 0.04;
        rates[4] = 0.06;
        rates[5] = 0.05;
        rates[6] = 0.04;
        
        
        dates2[0] = Date(17, May, 1999);    
        dates2[1] = Date(17, May, 2001);
        dates2[2] = Date(17, May, 2003);    
        dates2[3] = Date(17, May, 2005);
        dates2[4] = Date(17, May, 2007);    
        dates2[5] = Date(17, May, 2009);
        
        vols[0] = 0.20;
        vols[1] = 0.25;
        vols[2] = 0.30;
        vols[3] = 0.35;
        vols[4] = 0.40;
        vols[5] = 0.45;
                    // bootstrap the yield/dividend/vol forward curves 
        Handle<YieldTermStructure> fowardTermStructure(
            boost::shared_ptr<YieldTermStructure>(
                new ForwardCurve(dates1, rates, dayCounter)));
                
        Handle<YieldTermStructure> fowardDividendTS(
            boost::shared_ptr<YieldTermStructure>(
                new ForwardCurve(dates1, rates, dayCounter)));
                
        
        Handle<BlackVolTermStructure> fowardVolTS(
            boost::shared_ptr<BlackVolTermStructure>(
                new BlackVarianceCurve(todaysDate, dates2, vols,
                                     dayCounter)));
        // BlackScholes Merton Process forward curve
        boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
                new BlackScholesMertonProcess(underlyingH, fowardDividendTS, fowardTermStructure, fowardVolTS));
        
        std::cout<< "\n\nfoward curve tests : "<<std::endl;
        maturity=Date(17, May, 2001);
        forwardCurveSimulation(maturity, payoff,bsmProcess);
          
        maturity=Date(17, May, 2003);
        forwardCurveSimulation(maturity,payoff,bsmProcess);
        
        maturity=Date(17, May, 2008);
        forwardCurveSimulation(maturity, payoff,bsmProcess);
        
        
        
        
        
        // End test
        double seconds = timer.elapsed();
        Integer hours = int(seconds/3600);
        seconds -= hours * 3600;
        Integer minutes = int(seconds/60);
        seconds -= minutes * 60;
        std::cout << " \nRun completed in ";
        if (hours > 0)
            std::cout << hours << " h ";
        if (hours > 0 || minutes > 0)
            std::cout << minutes << " m ";
        std::cout << std::fixed << std::setprecision(0)
                  << seconds << " s\n" << std::endl;
        return 0;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }

}

